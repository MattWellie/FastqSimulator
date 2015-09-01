import os
import random
__author__ = 'mwelland'


class Sampler:
    """
        This class will simulate paired end reads across an input sequence. This uses a simulated fragment size of 140
        and 5-base intervals. A reverse complement is used to sample an every alternate sequence as reverse strand

        This class is highly dependent on dictionary in the format output by the bundled GBK and LRG parsers
    """

    def __init__(self, plain, modified, x, y, tile):
        """
        :param plain: The unaltered dictionary
        :param modified: The dictionary containing modified sequence
        :param x: The current x-coordinate for the read ID generator
        :param y: The current y-coordinate for the read ID generator
        :param tile: The current tile for the read ID generator
        """
        self.tile = tile
        self.x = x
        self.y = y
        self.plain_dict = plain
        self.mod_dict = modified
        self.gene_name = plain['genename']
        self.Read_1 = []
        self.Read_2 = []
        self.interval = 140  # Change for different fragment lengths
        self.alternate = 5  # read intervals
        self.reverse = 10  # should be double self.alternate
        self.read_length = 70  # Change for read sizes
        self.qual_string = ''.join(['I']*self.read_length)  # Temporary quality string
        self.out_fileR1 = ''
        self.out_fileR2 = ''
        self.R1_list = []
        self.R2_list = []

    def run(self):
        """
        Main control method.
        First the plain dictionary is processed, then the modified:
          -  For each transcript in the file
          -  Create a new pair of lists for reads 1 & 2
          -  Run the sampler which will generate reads and populate the lists
          -  Create filenames to write the output to
          -  Send the output to be written
        """
        print 'Gene: %s' % self.gene_name
        for transcript in self.plain_dict['transcripts']:
            self.R1_list = []
            self.R2_list = []
            self.run_plain(transcript)
            r1_filename = os.path.join('fastQs', self.gene_name + '0%sR1.fq' % transcript)
            r2_filename = os.path.join('fastQs', self.gene_name + '0%sR2.fq' % transcript)
            self.write_out(r1_filename, r2_filename)
            self.R1_list = []
            self.R2_list = []
            self.run_mod(transcript)
            r1_filename = os.path.join('fastQs', self.gene_name + '%sR1.fq' % transcript)
            r2_filename = os.path.join('fastQs', self.gene_name + '%sR2.fq' % transcript)
            self.write_out(r1_filename, r2_filename)
        return self.x, self.y, self.tile

    def run_plain(self, transcript):
        """
        Processes the dictionary containing an unaltered reference sequence
        Alternative transcripts may still occur, so each transcript is written separately
        Dictionary format:
        Transcripts > Transcript > List of exons
                                   Exons         > exon numbers > genomic_start
                                                                  genomic end
                                                                  padded sequence
                                                                  length
        """

        # Find the list of exon numbers using the selected transcript (method arg)
        exon_list = self.plain_dict['transcripts'][transcript]['list_of_exons']

        # Use only the dictionary section which applies for selected transcript (method arg)
        plain_dict = self.plain_dict['transcripts'][transcript]['exons']

        # For each exon
        for exon in exon_list:
            # Grab the exon region with padding added in list form
            exon_seq = plain_dict[exon]['padded seq']
            length = plain_dict[exon]['padded length']

            # Start from start of sequence
            offset = 0
            while offset <= length - (self.interval+1):  # To prevent index errors, flanking seq still leaves plenty of coverage

                # Using offset to move through the sequence, select a section of length *interval* (see __init__)
                sub_list = exon_seq[offset:offset + self.interval]

                # For every alternate sequence, work on the reverse complement instead
                if self.interval % self.reverse == 0:
                    sub_list = self.reverse_complement(sub_list)

                # Use extract method to populate list
                self.extract_to_list(sub_list)
                # Update offset
                offset += self.alternate

    def run_mod(self, transcript):
        """
        Processes the altered form of the dictionary, which includes a slightly different structure
        """

        # Use only the relevant transcript portion of the dictionary
        tran_dict = self.mod_dict[transcript]

        # Select the part which contains the exon numbers and coordinates
        exon_list = tran_dict['exon list']

        for exon in exon_list:
            exon_dict = tran_dict['exons'][exon]

            # A variable to count reads created per exon
            count = 0

            sequence = exon_dict['padded seq']
            length = exon_dict['padded length']
            offset = 0
            while offset < length - (self.interval+1):  # Stopping point to prevent index errors
                sub_list = sequence[offset:offset + self.interval]
                # Reverse every other 'fragment'
                if self.interval % self.reverse == 0:
                    sub_list = self.reverse_complement(sub_list)

                # Send for extraction to list
                self.extract_to_list(sub_list)
                # Increase the count of read pairs for this exon
                count += 1
                # Update the offset value
                offset += self.alternate


    def extract_to_list(self, sequence):
        """
        :param sequence: The 'fragment'
        For each 'fragment' this should take reads from either end and send to an output list
        Use a perfect Q40 quality String (reads do not need to simulate true qualities for this project)
        """

        read1 = ''.join(sequence[:self.read_length])

        # Produce a reverse complement of the end of the fragment
        read2 = ''.join(self.reverse_complement(sequence[self.interval-self.read_length:]))

        read_id = self.generate_seq_id()
        self.R1_list.append(read_id % 1)
        self.R1_list.append(read1)
        self.R1_list.append('+')
        self.R1_list.append(self.qual_string)
        self.R2_list.append(read_id % 2)
        self.R2_list.append(read2)
        self.R2_list.append('+')
        self.R2_list.append(self.qual_string)


    def write_out(self, r1_filename, r2_filename):
        # Open files and write contents of the Read files out
        with open(r1_filename, 'w') as outfile:
            for line in self.R1_list:
                print >> outfile, line
        with open(r2_filename, 'w') as outfile:
            for line in self.R2_list:
                print >> outfile, line


    def generate_seq_id(self):
        """
        :return: Creates a valid format Illumina FastQ header; generated number pairs are used for coordinates
        """
        instrument = 'MATTW_ART1'
        run_number = 00001
        flow_id = 'FLOW%d' % 01
        lane = 1
        tile = self.tile
        #  Don't repeat same numbers
        x_pos = self.x
        y_pos = self.y
        end = '%d:N:0:1'
        self.x += 1
        if self.x > 2500:
            self.y += 1
            self.x = 1
        if self.y == 2500:
            self.tile += 1
            self.x = 1
            self.y = 1

        return '@%s:%d:%s:%d:%d:%d:%d %s' % (instrument, run_number, flow_id, lane, tile, x_pos, y_pos, end)

    @staticmethod
    def reverse_complement(sequence):
        """
        :param sequence: DNA substring
        :return: reverse complement of the substring
        """
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        new_bases = []
        for base in sequence:
            new_bases.append(complement[base])
        return new_bases[::-1]
import os
import random
__author__ = 'mwelland'


class Sampler:
    """
        This class will simulate paired end reads across an input sequence. This uses a simulated fragment size of 220
        and 1-base intervals. A reverse complement will also be used to sample an every alternate sequence as a reverse
        strand with a 1-base offset.
    """

    def __init__(self, plain, modified):
        """
        :param plain: The unaltered dictionary
        :param modified: The dictionary containing modified sequence
        """
        self.plain_dict = plain
        self.mod_dict = modified
        self.gene_name = plain['genename']
        self.Read_1 = []
        self.Read_2 = []
        self.interval = 220  # Change for different fragment lengths
        self.alternate = 1  # read intervals
        self.reverse = 2  # should be double self.alternate
        self.read_length = 101  # Change for read sizes
        self.qual_string = ''.join(['I']*101)  # Temporary quality string
        self.out_fileR1 = ''
        self.out_fileR2 = ''
        self.R1_list = []
        self.R2_list = []

    def run(self):
        """
        Main control method.
        First the plain dictionary is processed, the the modified:
          -  For each transcript in the file
          -  Create a new pair of lists of reads 1 & 2
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
            #  print 'at regular print r1 length: %d' % len(self.R1_list)
            #  print 'at regular print r2 length: %d' % len(self.R2_list)
            #  this = raw_input()
            self.write_out(r1_filename, r2_filename)
        for transcript in self.mod_dict:
            self.R1_list = []
            self.R2_list = []
            self.run_mod(transcript)
            r1_filename = os.path.join('fastQs', self.gene_name + '%sR1.fq' % transcript)
            r2_filename = os.path.join('fastQs', self.gene_name + '%sR2.fq' % transcript)
            #  print 'at mod print r1 length: %d' % len(self.R1_list)
            #  print 'at mod print r2 length: %d' % len(self.R2_list)
            #  this = raw_input()
            self.write_out(r1_filename, r2_filename)

    def run_plain(self, transcript):
        """
        Processes the dictionary containing an unaltered reference sequence
        Alternative transcripts may still occur, so each transcript is written separately
        """
        #  Grab the full genetic sequence
        sequence = list(self.plain_dict['full genomic sequence'])
        #  Find the list of exon numbers using the selected transcript (method arg)
        exon_list = self.plain_dict['transcripts'][transcript]['list_of_exons']
        print exon_list
        this = raw_input()
        #  Use only the dictionary section which applies for selected transcript (method arg)
        plain_dict = self.plain_dict['transcripts'][transcript]['exons']

        #  Foreach exon
        for exon in exon_list:
            #  Gather start and stop coordinates from the dict
            start = plain_dict[exon]['genomic_start']
            end = plain_dict[exon]['genomic_end']

            #  Subselect part of the overall sequence with a large area of overlap
            exon_seq = sequence[start - 330:end + 330]  # Maybe change substring
            length = len(exon_seq)
            print 'Exon %d length: %d' % (exon, length)

            #  Start from start of sequence
            offset = 0
            while offset <= length - 205: #  To prevent index errors, flanking seq still leaves plenty of coverage

                #  Using offset to move through the sequence, select a section of length *interval* (see __init__)
                sub_list = exon_seq[offset:offset + self.interval]

                #  For every alternate sequence, work on the reverse complement instead
                if self.interval % self.reverse == 0:
                    sub_list = self.reverse_complement(sub_list)

                #  Use extract method to populate list
                self.extract_to_list(sub_list)
                #  Update offset
                offset += self.alternate
            #  print 'at exon %d end, r1 length: %d' % (exon, len(self.R1_list))
            #  print 'at exon %d end, r2 length: %d' % (exon, len(self.R2_list))
            #  this = raw_input()

    def run_mod(self, transcript):
        """
        Processes the altered form of the dictionary, which includes a slightly different structure
        """

        #  Use only the relevant transcript portion of the dictionary
        tran_dict = self.mod_dict[transcript]
        sequence = list(tran_dict['sequence'])
        #  Select the part which contains the exon numbers and coordinates
        exon_dict = tran_dict['exons']
        for exon in exon_dict:
            #  A variable to count reads created per exon
            count = 0
            #  Gather start and stop sequences
            start = exon_dict[exon]['start']
            end = exon_dict[exon]['end']
            exon_seq = sequence[start - 330:end + 330]
            length = len(exon_seq)
            print 'Exon %d length: %d' % (exon, length)
            offset = 0
            while offset < length - 205:  # Arbitrary stopping point
                sub_list = sequence[offset:offset + self.interval]
                if self.interval % self.reverse == 0:
                    sub_list = self.reverse_complement(sub_list)
                self.extract_to_list(sub_list)
                count += 1
                offset += self.alternate
            print 'Read pairs for exon %d = %d' % (exon, count)

    def extract_to_list(self, sequence):
        """
        :param sequence: The 'fragment'
        For each 'fragment' this should take reads from either end and send to an output list
        """
        read1 = ''.join(sequence[:self.read_length])
        read2 = ''.join(self.reverse_complement(sequence[self.interval-self.read_length:]))
        #read2 = ''.join(sequence[self.interval-self.read_length:])
        #  print 'read1: %s' % read1
        #  print 'read2: %s' % read2
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
        with open(r1_filename, 'w') as outfile:
            for line in self.R1_list:
                print >> outfile, line
        with open(r2_filename, 'w') as outfile:
            for line in self.R2_list:
                print >> outfile, line

    @staticmethod
    def generate_seq_id():
        """
        :return: Creates a valid format Illumina FastQ header
        At least I hope the format is valid
        """
        instrument = 'MATTW_ART1'
        run_number = 00001
        flow_id = 'FLOW%d' % 01
        lane = 1
        tile = 1
        x_pos = random.randint(0, 99999)
        y_pos = random.randint(0, 99999)
        end = '%d:N:0:1'

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
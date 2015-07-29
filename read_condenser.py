import os

__author__ = 'mwelland'


class Condenser:

    """ A class which will be used to take the outputs of several fastQ files with different contents and condense
        them into a format ready for alignment
    """

    def __init__(self, gene_list):
        self.gene_list = gene_list
        self.file_list = os.listdir('fastQs')
        self.gene_dictionary = {}

    def run(self):
        """
        Control method for condenser
        :return:
        """
        for gene in self.gene_list:
            name_length = len(gene)
            self.populate_dict(gene, name_length)
            r1_pairs = self.create_file_pairings(gene, name_length, 1)
            r2_pairs = self.create_file_pairings(gene, name_length, 2)
            self.combine_files(r1_pairs)
            self.combine_files(r2_pairs)
        self.erase_old_files()

    def populate_dict(self, gene, name_length):
        """
        The dictionary will contain the names of all files relating to the selected gene. This method uses the gene
        name and the transcript numbers to identify the unchanged sequence and the accompanying changed transcripts.
        This is done to allow each separate transcript to be paired with the unchanged sequence, so as to represent
        a homozygous variant as well as improving the read depth for each area.
        :param gene:
        :param name_length:
        :return:
        """
        self.gene_dictionary = {1: {'ref': '', 'transcripts': []}, 2: {'ref': '', 'transcripts': []}}
        read1s = []
        read2s = []
        fq_list = [name for name in self.file_list if name[:name_length] == gene]
        #  Separate read 1s from read 2s
        for filename in fq_list:
            first_part_of_name = filename.split('.')[0]
            if first_part_of_name[-1] == '1':
                read1s.append(filename)
            elif first_part_of_name[-1] == '2':
                read2s.append(filename)

        #  For each of Read 1 and 2, separate reference from altered
        for filename in read1s:
            transcript = filename[name_length:name_length+1]
            if transcript == '0':
                self.gene_dictionary[1]['ref'] = filename
            else:
                self.gene_dictionary[1]['transcripts'].append(filename)

        for filename in read2s:
            transcript = filename[name_length:name_length+1]
            if transcript == '0':
                self.gene_dictionary[2]['ref'] = filename
            else:
                self.gene_dictionary[2]['transcripts'].append(filename)

    def create_file_pairings(self, gene, name_length, read):
        """
        :param gene: gene name
        :param name_length: length of gene name
        :param read: 1 or 2
        :return: a list of 3-element tuples
                element 1: reference file name (unchanged seq)
                element 2: changed file name
                element 3: a name for the file once combined
        """
        file_pairs = []
        read_dict = self.gene_dictionary[read]
        for filename in read_dict['transcripts']:
            file_pairs.append([read_dict['ref'], filename, self.create_file_name(gene, name_length, filename, read)])
        return file_pairs

    @staticmethod
    def create_file_name(gene, name_length, filename, read):
        """
        Creates a new file name which combines the gene name, the transcript number and the read number
        """
        transcript = filename[name_length:name_length+1]
        filename = '%s_transcript%s_R%d.fq' % (gene, transcript, read)
        return filename

    @staticmethod
    def combine_files(read_pairs):
        """
        This will combine the created files in memory and write the output to a new file
        the filename for the new file will be the one created in create_file_name
        :param read_pairs:
        :return:
        """

        for triple_tuple in read_pairs:

            outfile_name = triple_tuple[2]

            with open(os.path.join('fastQs', outfile_name), 'w') as output_file:
                reference_file_name = triple_tuple[0]
                reference_file = open(os.path.join('fastQs', reference_file_name), "r")
                contents = reference_file.readlines()
                output_file.writelines(contents)
                alt_file = triple_tuple[1]
                contents = open(os.path.join('fastQs', alt_file), 'r').readlines()
                output_file.writelines(contents)

    def erase_old_files(self):
        """
        Deletes all the old files which have now been combined
        """
        for filename in self.file_list:
            os.remove(os.path.join('fastQs', filename))

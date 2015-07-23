import random
__author__ = 'mwelland'


class Modifier:

    """ This class creates single nucleotide changes in each exon of the input sequence
        The exons are located using the exon coordinates from the input file
        Random is used to choose a random location within the exon and a new base
        The base substitution is made and the position is recorded in a dict
    """

    def __init__(self, file_dict):
        self.dict = file_dict
        self.output_dict = {}

    def run_modifier(self):
        """
        Main control method for the class
        :return: the completed dictionary object
        """

        for transcript in self.dict['transcripts']:
            self.output_dict[transcript] = {'variants': {}, 'sequence': '', 'exons': {}}
            exon_list = self.dict['transcripts'][transcript]['list_of_exons']
            self.output_dict[transcript]['exon list'] = exon_list
            for exon in exon_list:
                self.output_dict[transcript]['exons'][exon] = {'start': '', 'end': ''}
            try:
                for exon in exon_list:
                    exon_seq = self.modify(exon, transcript)
                    self.output_dict[transcript]['exons'][exon]['seq'] = exon_seq
                    self.output_dict[transcript]['exons'][exon]['seq length'] = len(exon_seq)
            except IndexError:
                print 'The index was out of range, line 35 seq_mod'
        return self.output_dict

    def modify(self, exon_number, transcript):
        """
        Method to execute substitutions. This is called for each exon number.
        This creates a single substitution per exon.
        :param exon_number:
        :param transcript:
        """

        padding = self.dict['offset']
        exondict = self.dict['transcripts'][transcript]['exons'][exon_number]
        exon_seq = exondict['padded sequence']
        start = padding
        end = exondict['length'] + padding
        edit_coord = random.randint(start+1, end-1)
        # edit_coord = random.randint(start, end)
        base = str(random.sample(['A', 'C', 'G', 'T'], 1)[0])
        old_base = exon_seq[edit_coord]
        while base == old_base:
            base = str(random.sample(['A', 'C', 'G', 'T'], 1)[0])
        exon_seq[edit_coord] = base
        self.output_dict[transcript]['variants'][exon_number] = {'position': edit_coord, 'new base': base}
        return exon_seq
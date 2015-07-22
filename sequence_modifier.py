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
        self.full_seq = list(self.dict['full genomic sequence'])
        self.output_dict = {}

    def run_modifier(self):
        """
        Main control method for the class
        :return: the completed dictionary object
        """

        for transcript in self.dict['transcripts']:
            transcript_seq = self.full_seq
            self.output_dict[transcript] = {'variants': {}, 'sequence': '', 'exons': {}}
            exon_list = self.dict['transcripts'][transcript]['list_of_exons']
            for exon in exon_list:
                self.output_dict[transcript]['exons'][exon] = {'start': '', 'end': ''}
            try:
                for exon_number in exon_list:
                    transcript_seq = self.modify(exon_number, transcript, transcript_seq)
            except IndexError:
                print 'The index was out of range, line 34 seq_mod'
            self.output_dict[transcript]['sequence'] = ''.join(transcript_seq)
        return self.output_dict

    def modify(self, exon_number, transcript, transcript_seq):
        """
        Method to execute substitutions. This is called for each exon number.
        This creates a single substitution per exon.
        :param exon_number:
        :param transcript:
        """
        exondict = self.dict['transcripts'][transcript]['exons'][exon_number]
        start = exondict['genomic_start']
        self.output_dict[transcript]['exons'][exon_number]['start'] = start
        end = exondict['genomic_end']
        self.output_dict[transcript]['exons'][exon_number]['end'] = end
        edit_coord = random.randint(start+1, end-1)
        # edit_coord = random.randint(start, end)
        base = str(random.sample(['A', 'C', 'G', 'T'], 1)[0])
        old_base = transcript_seq[edit_coord]
        while base == old_base:
            base = str(random.sample(['A', 'C', 'G', 'T'], 1)[0])
        transcript_seq[edit_coord] = base
        self.output_dict[transcript]['variants'][exon_number] = {'position': edit_coord, 'new base': base}
        return transcript_seq
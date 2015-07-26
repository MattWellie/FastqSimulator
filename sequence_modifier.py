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
        self.genomic = file_dict['full sequence']
        self.padding = file_dict['offset']
        self.transcript_pos = 0

    def run_modifier(self):
        """
        Main control method for the class
        :return: the completed dictionary object
        """

        for transcript in self.dict['transcripts']:
            self.transcript_pos = 0

            self.output_dict[transcript] = {'variants': {}, 'exons': {}}
            exon_list = self.dict['transcripts'][transcript]['list_of_exons']

            #print 'exon list: %s' % str(exon_list)
            #this = raw_input()
            self.output_dict[transcript]['exon list'] = exon_list
            for exon in exon_list:

                try:
                    start = self.dict['transcripts'][transcript]['exons'][exon]['genomic_start']
                    end = self.dict['transcripts'][transcript]['exons'][exon]['genomic_end']

                    self.modify(exon, transcript, start, end)
                except IndexError:
                    print 'The index was out of range, line 35 seq_mod'
            for exon in exon_list:
                start = self.dict['transcripts'][transcript]['exons'][exon]['genomic_start']
                end = self.dict['transcripts'][transcript]['exons'][exon]['genomic_end']
                exon_seq = self.genomic[start - self.padding: end + self.padding]
                self.output_dict[transcript]['exons'][exon] = {'start': start, 'end': end,
                                                               'padded seq': exon_seq, 'padded length': len(exon_seq),
                                                               'length': end - start}
        return self.output_dict

    def modify(self, exon_number, transcript, start, end):
        """
        Method to execute substitutions. This is called for each exon number.
        This creates a single substitution per exon.
        :param exon_number:
        :param transcript:
        """

        edit_coord = random.randint(start+1, end-1)
        # edit_coord = random.randint(start, end)
        base = str(random.sample(['A', 'C', 'G', 'T'], 1)[0])
        old_base = self.genomic[edit_coord]
        while base == old_base:
            base = str(random.sample(['A', 'C', 'G', 'T'], 1)[0])
        self.genomic[edit_coord] = base

        ### Identify the position of the variant within the transcript (For HGVS)
        exon_length = end - start
        variant_pos = end - edit_coord
        hgvs_pos = self.transcript_pos + variant_pos
        hgvs = '%d%s>%s' %(hgvs_pos, old_base, base)
        self.transcript_pos += exon_length
        self.output_dict[transcript]['variants'][exon_number] = {'position': edit_coord, 'new base': base, 'hgvs': hgvs}

import random
__author__ = 'mwelland'


class Modifier:

    """
    This class creates single nucleotide changes in each exon of the input
    sequence. The exons are located using the exon coordinates from the input
    file. Random is used to choose a random location within the exon and a new
    base. The base substitution is made and the position is recorded in a dict
    """

    def __init__(self, file_dict, type):
        self.type = type
        self.dict = file_dict
        self.output_dict = {}
        self.genomic = file_dict['full sequence']
        self.padding = file_dict['offset']
        self.transcript_dict = {}
        self.transcript_pos = 0

    def run_modifier(self):
        """
        Main control method for the class
        :return: the completed dictionary object
        """

        for transcript in self.dict['transcripts']:
            self.transcript_pos = 0
            self.transcript_dict = self.dict['transcripts'][transcript]
            self.output_dict[transcript] = {'variants': {}, 'exons': {}}
            exon_list = self.transcript_dict['list_of_exons']
            self.output_dict[transcript]['exon list'] = exon_list
            for exon in exon_list:
                try:
                    start = self.transcript_dict['exons'][exon]['genomic_start']
                    end = self.transcript_dict['exons'][exon]['genomic_end']
                    b_or_a = self.transcript_dict['exons'][exon]['cds']
                    old_cds = self.transcript_dict['old_cds_offset']
                    self.modify(exon, transcript, start, end, b_or_a, old_cds)
                except IndexError:
                    print 'The index was out of range, line 43 seq_mod'
            for exon in exon_list:
                start = self.transcript_dict['exons'][exon]['genomic_start']
                end = self.transcript_dict['exons'][exon]['genomic_end']
                exon_seq = self.genomic[
                           start - self.padding:
                           end + self.padding
                           ]
                self.output_dict[transcript]['exons'][exon] = {
                    'start': start,
                    'end': end,
                    'padded seq': exon_seq,
                    'padded length': len(exon_seq),
                    'length': end - start
                }
        return self.output_dict

    def modify(self, exon_number, transcript, start, end, b_or_a, old_cds):
        """
        Method to execute substitutions. This is called for each exon number.
        This creates a single substitution per exon.
        :param exon_number:
        :param transcript:
        """
        if self.type == 'lrg':
            start -= 1
        old_exon = self.genomic[start:end]
        edit_coord = random.randint(start, end-1)
        base = str(random.sample(['A', 'C', 'G', 'T'], 1)[0])
        old_base = self.genomic[edit_coord]
        while base == old_base:
            base = str(random.sample(['A', 'C', 'G', 'T'], 1)[0])
        self.genomic[edit_coord] = base
        # Identify the position of the variant within the transcript (For HGVS)
        if b_or_a == 'before':
            cds_delay = old_cds
        else:
            cds_delay = self.transcript_dict['cds_offset']
        p_length = self.transcript_dict['protein_length'] + 3  # len(stop codon)
        exon_length = end - start
        """
        Put in something here to deal with if the exon is pre-coding
        b_or_a is 'before' if the exon is pre-coding sequence
        """
        # Edit coord + 1 to compensate for 0-base counting sequence index
        if b_or_a == 'before':
            hgvs_pos = (edit_coord + 1) - (old_cds - 1)
        else:
            variant_pos = (edit_coord + 1) - start
            hgvs_pos = (self.transcript_pos + variant_pos) - cds_delay

        """
        This is the correct hgvs nomenclature for annotation by annovar
        """
        if hgvs_pos <= 0:
            hgvs = 'c.%d%s>%s' % (hgvs_pos - 1, old_base, base)
        elif hgvs_pos > p_length:
            after_coding = hgvs_pos - p_length  # compensate for no 0s
            hgvs = 'c.*%d%s>%s' % (after_coding, old_base, base)
        else:
            hgvs = 'c.%s%d%s' % (old_base, hgvs_pos, base)
        self.transcript_pos += exon_length

        """
        # Un comment this section for positional debugging
        print 'start: %d, end: %d' % (start, end)
        print 'Length: %d' % (exon_length)
        print 'for gene: %s, file: %s cds: %d' % (self.dict['genename'],
         self.dict['filename'], cds_delay)
        print 'position %d in exon %d' % (variant_pos, exon_number)
        print 'changed %s to %s' % (old_base, base)
        print 'HGVS: %s' % hgvs
        print 'Exon: \n%s' % self.genomic[start:end]
        new_exon = self.genomic[start:end]
        i = 0
        while i <= len(new_exon):
            base1 = old_exon[i]
            base2 = new_exon[i]
            if base1 == base2:
                i += 1
            else:
                print base1
                print base2
                print 'position: %d' % (i + 1)
                j = (i + 1) - cds_delay
                if j <= 0:
                    j -= 1
                print 'Assumed HGVS: c.%s%d%s' % (base1, j, base2)
                break
        this = raw_input()
        """
        self.output_dict[transcript]['variants'][exon_number] = {
            'position': edit_coord,
            'new base': base,
            'hgvs': hgvs,
            'transcript': self.transcript_dict['NM_number']
        }
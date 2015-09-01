"""
This class is mostly reused from a previous project, and contains methods for
extracting contents from genomic reference files in GenBank file format.

Dict { full genomic sequence
       genename
       refseqname
       transcripts {  transcript { cds_offset
                                   exons      {  exon_number {   genomic_start
                                                                 genomic_stop
"""

import Bio
from Bio import SeqIO

__author__ = 'mwelland'
__version__ = 1.3
__version_date__ = '11/02/2015'


class GbkParser:


    def __init__(self, file_name, padding):
        """
        This class is created by instantiating with a file name and a padding
        value. These are used to locate the appropriate target file, and to
        select the amount of flanking sequence to be appended to exons.
        """
        '''
        :param file_name: the location/identity of the target input file
        :param padding: the required amount of intronic padding
        '''
        self.genomic = ''
        self.padding = padding
        self.exons = []
        self.cds = []
        self.mrna = []
        self.fileName = file_name

        # Read in the specified input file into a variable
        try:
            self.t_dict = {
                'transcripts': {},
                'input': SeqIO.to_dict(SeqIO.parse(file_name, 'genbank')),
                'offset': self.padding,
                'filename': file_name}
            self.t_dict['refseqname'] = \
                self.t_dict['input'].keys()[0]
            self.is_matt_awesome = True
        except IOError as fileNotPresent:
            print "The specified file cannot be located: "\
                  + fileNotPresent.filename
            exit()

    @property
    def get_version(self):
        """
        Quick function to grab version details for final printing
        :return:
        """
        return 'Version: {0}, Version Date: {1}'.format(str(__version__),
                                                        __version_date__)
          
    def get_mrna_exons(self):
        """
        This uses the list of exon start and stop positions to populate the
        exon positions in the dictionary
        """

        # Each transcript is handled separately
        for alternative in self.t_dict['Alt transcripts']:
            alt_dict = {'list_of_exons': [], 'exons': {}}
            selected_mrna = self.mrna[alternative-1]
            # The location of this feature is different
            # for NCBI and Ensembl genbank formats
            try:
                alt_dict['NM_number'] = \
                    selected_mrna.qualifiers['transcript_id'][0]
            except KeyError:
                alt_dict['NM_number'] = \
                    self.t_dict['genename']
                self.t_dict['refseqname'] = \
                    self.t_dict['genename']
                self.t_dict['genename'] = \
                    self.cds[0].qualifiers['gene'][0]
            exon = 1

            # The parser handles single- and multiple-exon genes differently
            if len(self.exons) == 1:
                alt_dict['exons'][exon] = {}
                alt_dict['list_of_exons'].append(exon)
                alt_dict['exons'][exon]['genomic_end'] = \
                    selected_mrna.location.end
                alt_dict['exons'][exon]['genomic_start'] = \
                    selected_mrna.location.start
            else:
                subfeatures = selected_mrna._get_sub_features()
                for coords in subfeatures:
                    alt_dict['list_of_exons'].append(exon)
                    start = coords.location.start
                    end = coords.location.end
                    exon_seq = list(self.genomic[
                                    start - self.padding:
                                    end + self.padding
                                    ])
                    length = len(exon_seq)
                    minidict = {
                        'genomic_start': start,
                        'genomic_end': end,
                        'padded seq': exon_seq,
                        'length': end - start,
                        'padded length': length
                    }
                    alt_dict['exons'][exon] = minidict
                    exon += 1
            self.t_dict['transcripts'][alternative] = alt_dict


    def get_protein(self):
        """
        This method takes the CDS tagged block from the GenBank features
        section and parses the contents to retrieve the protein sequence.
        This is added to the appropriate section of dictionary used to
        hold all required details of the file.
        """
        '''
        :param cds: a list containing the cds element(s) of the genbank features
        '''
        for alternative in self.t_dict['Alt transcripts']:
            selected_cds = self.cds[alternative-1]
            self.t_dict['transcripts'][alternative]['protein_length'] =\
                len(selected_cds.qualifiers['translation'][0])*3
            self.t_dict['transcripts'][alternative]['old_cds_offset'] =\
                selected_cds.location.start


    def find_cds_delay(self):
        """
        Method to find the actual start of the translated sequence - introduced
        to sort out non-coding exon problems. This involves taking the
        'transcript start location' variable and calculating how many exonic
        bases that is from the first exonic base in the file.
        The transcript start location indicates the genomic position of the
        first coding base relative to the start of the reference file, whilst
        the required value is the number of exonic bases from first exon to
        first codon.
        For each exon, if the coding sequence start point coord. is after the
        exon end coord., add the length of the exon to the offset value. If the
        coding sequence starts between the start and end of the exon, add the
        length of the start point from the start of the exon

        Sorry if this doesn't make any sense. I'll add a diagram to the Git repo
        """
        for t in self.t_dict['transcripts']:
            offset_total = 0
            offset = self.t_dict['transcripts'][t]['old_cds_offset']
            exon_list = self.t_dict['transcripts'][t]['list_of_exons']
            for e in exon_list:
                g_start = \
                    self.t_dict['transcripts'][t]['exons'][e]['genomic_start']
                g_stop =\
                    self.t_dict['transcripts'][t]['exons'][e]['genomic_end']
                if offset > g_stop:
                    offset_total = offset_total + (g_stop - g_start)
                    self.t_dict['transcripts'][t]['exons'][e]['cds'] = 'before'
                elif g_stop > offset > g_start:
                    self.t_dict['transcripts'][t]['cds_offset'] = \
                        offset_total + (offset - g_start)
                    self.t_dict['transcripts'][t]['exons'][e]['cds'] = 'after'
                elif offset < g_start:
                    self.t_dict['transcripts'][t]['exons'][e]['cds'] = 'after'


    def fill_and_find_features(self):
        """
        Go through the file and add a range of basic features to the dictionary
        :return:
        """
        dictionary = self.t_dict['input'][self.t_dict['refseqname']]
        self.genomic = dictionary.seq
        features = dictionary.features
        for feature in features:
            # Multiple exons are expected, not explicitly used
            if feature.type == 'exon':
                self.exons.append(feature)

        """
        This section works on the assumption that each exon in the file will
        use the appropriate gene name and that the only relevant CDS and mRNA
        sections will also contain the same accession
        """
        try:
            self.t_dict['genename'] = self.exons[0].qualifiers['gene'][0]
            for feature in features:
                if feature.type == 'CDS':
                    if feature.qualifiers['gene'][0] == self.t_dict['genename']:
                        self.cds.append(feature)
                elif feature.type == 'mRNA':
                    if feature.qualifiers['gene'][0] == self.t_dict['genename']:
                        self.mrna.append(feature)
        except KeyError:
            for feature in features:
                if feature.type == 'CDS':
                    self.cds.append(feature)
                elif feature.type == 'mRNA':
                    self.mrna.append(feature)
            note = self.mrna[0].qualifiers['note'][0]
            self.t_dict['genename'] = note.split('=')[1]
        assert len(self.cds) == len(self.mrna),\
            "There are a different number of CDS and mRNA"
        return features


    def run(self):
        """
        This is the main method of the GBK Parser. This method is called after
        class instantiation and handles the operation of all the other functions
        to complete the dictionary which will hold all of the sequence and exon
        details of the gene file being parsed
        """
        print 'BioPython version: ' + str(Bio.__version__)
        # initial sequence grabbing and populating dictionaries
        self.fill_and_find_features()
        self.t_dict['full sequence'] = list(self.genomic)
        self.t_dict['Alt transcripts'] = range(1, len(self.cds)+1)
        self.get_mrna_exons()
        self.get_protein()
        self.find_cds_delay()
        del self.t_dict['input']
        return self.t_dict
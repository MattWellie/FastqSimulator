import Bio
from Bio import SeqIO

__author__ = 'mwelland'
__version__ = 1.3
__version_date__ = '11/02/2015'


class GbkParser:
    """
    Notes:
        Isolated class to deal exclusively with GBK files
        Should return dictionary, not write full output

    Parses the input file to find all the useful values
    This will populate a dictionary to be returned at completion

            Dict { full genomic sequence
                   genename
                   refseqname
                   transcripts {  transcript {   exons {  exon_number {   genomic_start
                                                                          genomic_stop
    """

    def __init__(self, file_name, padding):

        """
        This class is created by instantiating with a file name and a padding value.
        These are used to locate the appropriate target file, and to select the
        amount of flanking sequence to be appended to exons.
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
            self.transcriptdict = {'transcripts': {}, 'input': SeqIO.to_dict(SeqIO.parse(file_name, 'genbank')),
                                   'offset': self.padding}
            self.transcriptdict['refseqname'] = self.transcriptdict['input'].keys()[0]
            self.is_matt_awesome = True
        except IOError as fileNotPresent:
            print "The specified file cannot be located: " + fileNotPresent.filename
            exit()

    @property
    def get_version(self):
        """
        Quick function to grab version details for final printing
        :return:
        """
        return 'Version: {0}, Version Date: {1}'.format(str(__version__), __version_date__)
          
    def get_mrna_exons(self):
        """ This uses the list of exon start and stop positions to populate 
            the exon positions in the dictionary""" 

        for alternative in self.transcriptdict['Alt transcripts']:
            self.transcriptdict['transcripts'][alternative] = {}
            self.transcriptdict['transcripts'][alternative]['list_of_exons'] = []
            self.transcriptdict['transcripts'][alternative]['exons'] = {}
            selected_mrna = self.mrna[alternative-1]
            try:
                self.transcriptdict['transcripts'][alternative]['NM_number'] = selected_mrna.qualifiers['transcript_id'][0]
            except KeyError:
                self.transcriptdict['transcripts'][alternative]['NM_number'] = self.transcriptdict['genename']
                self.transcriptdict['refseqname'] = self.transcriptdict['genename']
                self.transcriptdict['genename'] = self.cds[0].qualifiers['gene'][0]
            exon = 1
            subfeatures = selected_mrna._get_sub_features()
            
            for coords in subfeatures:
                self.transcriptdict['transcripts'][alternative]['exons'][exon] = {}
                self.transcriptdict['transcripts'][alternative]['list_of_exons'].append(exon)
                start = coords.location.start
                end = coords.location.end
                exon_seq = list(self.genomic[start - self.padding: end + self.padding])
                length = len(exon_seq)
                minidict = {'genomic_start': start, 'genomic_end': end,
                            'padded seq': exon_seq, 'length': end - start, 'padded length': length}
                self.transcriptdict['transcripts'][alternative]['exons'][exon] = minidict
                exon += 1

    def fill_and_find_features(self):
        dictionary = self.transcriptdict['input'][self.transcriptdict['refseqname']]
        self.genomic = dictionary.seq
        features = dictionary.features
        for feature in features:
            # Multiple exons are expected, not explicitly used
            if feature.type == 'exon':
                self.exons.append(feature)

        """ This section works on the assumption that each exon in the file will use the appropriate gene name
            and that the only relevant CDS and mRNA sections will also contain the same accession
        """
        try:
            self.transcriptdict['genename'] = self.exons[0].qualifiers['gene'][0]
            for feature in features:
                if feature.type == 'CDS':
                    if feature.qualifiers['gene'][0] == self.transcriptdict['genename']:
                        self.cds.append(feature)
                elif feature.type == 'mRNA':
                    if feature.qualifiers['gene'][0] == self.transcriptdict['genename']:
                        self.mrna.append(feature)
        except KeyError:
            for feature in features:
                if feature.type == 'CDS':
                    self.cds.append(feature)
                elif feature.type == 'mRNA':
                    self.mrna.append(feature)
            note = self.mrna[0].qualifiers['note'][0]
            self.transcriptdict['genename'] = note.split('=')[1]
        assert len(self.cds) == len(self.mrna), "There are a different number of CDS and mRNA"
        return features

    def run(self):
        """
        This is the main method of the GBK Parser. This method is called after class instantiation
        and handles the operation of all the other functions to complete the dictionary which will
        hold all of the sequence and exon details of the gene file being parsed
        """
        '''
        :return transcriptdict: This function fills and returns the dictionary, contents
                explained in Class docstring above
        '''
        print 'BioPython version: ' + str(Bio.__version__)
        # initial sequence grabbing and populating dictionaries
        self.fill_and_find_features()
        self.transcriptdict['full sequence'] = list(self.genomic)
        self.transcriptdict['Alt transcripts'] = range(1, len(self.cds)+1)
        self.get_mrna_exons()
        return self.transcriptdict

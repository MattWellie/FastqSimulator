from xml.etree.ElementTree import parse

__author__ = 'mwelland'
__version__ = 0.9
__version_date__ = '11/02/2015'


class LrgParser:

    """
    Class version: 0.9
    Modified Date: 21/07/2015
    Author : Matt Welland

    Minimal version of previously developed class
    Only requires exon numbers and coordinates, and full sequence

    Parses the input file to find all the useful values
    This will populate a dictionary to be returned at completion

            Dict { filename
                   genename
                   transcripts {  transcript {   exons {  exon_number {   genomic_start
                                                                          genomic_stop
                                                                          padded sequence
                                                                          length
                                                                          padded length
    """

    def __init__(self, file_name):
        self.fileName = file_name
        self.padding = 0
        self.sequence = ''
        # Read in the specified input file into a variable
        try:
            self.tree = parse(self.fileName)
            self.transcriptdict = {'transcripts': {},
                                   'root': self.tree.getroot()}
            self.transcriptdict['fixannot'] = self.transcriptdict['root'].find(
                'fixed_annotation')  # ensures only exons from the fixed annotation will be taken
            self.transcriptdict['updatable'] = self.transcriptdict['root'].find(
                'updatable_annotation')
            self.transcriptdict['genename'] = self.transcriptdict['root'].find(
                'updatable_annotation/annotation_set/lrg_locus').text
            self.transcriptdict['refseqname'] = self.transcriptdict['root'].find(
                'fixed_annotation/sequence_source').text
            if self.transcriptdict['root'].attrib['schema_version'] != '1.9':
                print 'This LRG file is not the correct version for this script'
                print 'This is designed for v.1.8'
                print 'This file is v.' + self.transcriptdict['root'].attrib['schema_version']
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

    # Grabs the sequence string from the <sequence/> tagged block
    def grab_element(self, path):
        """ Grabs specific element from the xml file from a provided path """
        try:
            for item in self.transcriptdict['root'].findall(path):
                return item.text
        except:
            print "No sequence was identified"
            print self.transcriptdict['filename']
            exit()

    def get_exon_coords(self):
        """ Traverses the LRG ETree to find all the useful values
            This should allow more robust use of the stored values, and enhances
            transparency of the methods put in place. Absolute references should
            also make the program more easily extensible
        """

        for items in self.transcriptdict['fixannot'].findall('transcript'):
            t_number = int(items.attrib['name'][1:])
            # print 'first t number = ' + str(t_number)
            self.transcriptdict['transcripts'][t_number] = {}  # First should be indicated with '1'; 'p1' can write on
            self.transcriptdict['transcripts'][t_number]["exons"] = {}
            self.transcriptdict['transcripts'][t_number]['list_of_exons'] = []
            # Gene sequence main coordinates are required to take introns
            # Transcript coordinates wanted for output  
            genomic_start = 0
            genomic_end = 0
            for exon in items.iter('exon'):
                exon_number = exon.attrib['label']
                if exon_number[-1] in ('a', 'b', 'c', 'd'):
                    # print exon_number
                    exon_number = exon_number[:-1]
                    print 'exon number: %d' % exon_number
                exon_number = int(exon_number)
                self.transcriptdict['transcripts'][t_number]['list_of_exons'].append(exon_number)
                self.transcriptdict['transcripts'][t_number]["exons"][exon_number] = {}
                for coordinates in exon:
                    if coordinates.attrib['coord_system'][-2] not in ['t', 'p']:
                        genomic_start = int(coordinates.attrib['start'])
                        genomic_end = int(coordinates.attrib['end'])
                assert genomic_start >= 0, "Exon index out of bounds"
                self.transcriptdict['transcripts'][t_number]["exons"][exon_number]['genomic_start'] = genomic_start
                self.transcriptdict['transcripts'][t_number]["exons"][exon_number]['genomic_end'] = genomic_end
                self.transcriptdict['offset'] = self.padding
                exon_seq = list(self.sequence[genomic_start - self.padding: genomic_end + self.padding])
                self.transcriptdict['transcripts'][t_number]["exons"][exon_number]['padded seq'] = exon_seq
                self.transcriptdict['transcripts'][t_number]["exons"][exon_number]['length'] = genomic_end-genomic_start
                self.transcriptdict['transcripts'][t_number]["exons"][exon_number]['padded length'] = len(exon_seq)


    def run(self, padding):
        self.padding = padding
        # Initial sequence grabbing and populating dictionaries
        self.sequence = self.grab_element('fixed_annotation/sequence')
        self.transcriptdict['full sequence'] = list(self.sequence)
        self.get_exon_coords()

        for transcript in self.transcriptdict['transcripts'].keys():
            self.transcriptdict['transcripts'][transcript]['list_of_exons'].sort(key=float)

        return self.transcriptdict
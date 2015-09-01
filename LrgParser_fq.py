"""
Class version: 0.9
Modified Date: 21/07/2015
Author : Matt Welland
Minimal version of previously developed class, similar to GbkParser
Only requires exon numbers and coordinates, and full sequence
Parses the input file to find all the useful values
This will populate a dictionary to be returned at completion

Dict { filename
       genename
       transcripts { transcript { exons { exon_number { genomic_start
                                                        genomic_stop
                                                        padded sequence
                                                        length
                                                        padded length
"""
__author__ = 'mwelland'
__version__ = 0.9
__version_date__ = '11/02/2015'

from xml.etree.ElementTree import parse


class LrgParser:


    def __init__(self, file_name):
        self.fileName = file_name
        self.padding = 0
        self.sequence = ''
        # Read in the specified input file into a variable
        try:
            self.tree = parse(self.fileName)
            self.t_dict = {'transcripts': {},
                                   'root': self.tree.getroot(),
                                   'filename': file_name}
            self.t_dict['fixannot'] = \
                self.t_dict['root'].find('fixed_annotation')
            self.t_dict['updatable'] = \
                self.t_dict['root'].find('updatable_annotation')
            self.t_dict['genename'] = \
                self.t_dict['root'].find(
                'updatable_annotation/annotation_set/lrg_locus').text
            self.t_dict['refseqname'] = \
                self.t_dict['root'].find(
                'fixed_annotation/sequence_source').text
            if self.t_dict['root'].attrib['schema_version'] != '1.9':
                print 'This LRG file is not the correct version for this script'
                print 'This is designed for v.1.8'
                print 'This file is v.' + \
                      self.t_dict['root'].attrib['schema_version']
            self.is_matt_awesome = True
        except IOError as fileNotPresent:
            print "The specified file cannot be located: " + \
                  fileNotPresent.filename
            exit()

    @property
    def get_version(self):
        """
        Quick function to grab version details for final printing
        :return:
        """
        return 'Version: {0}, Version Date: {1}'.format(str(__version__),
                                                        __version_date__)

    # Grabs the sequence string from the <sequence/> tagged block
    def grab_element(self, path):
        """ Grabs specific element from the xml file from a provided path """
        try:
            for item in self.t_dict['root'].findall(path):
                return item.text
        except:
            print "No sequence was identified"
            print self.t_dict['filename']
            exit()

    def get_exon_coords(self):
        """ Traverses the LRG ETree to find all the useful values
            This should allow more robust use of the stored values, and enhances
            transparency of the methods put in place. Absolute references should
            also make the program more easily extensible
        """

        for items in self.t_dict['fixannot'].findall('transcript'):
            t = int(items.attrib['name'][1:])
            # print 'first t number = ' + str(t_number)
            self.t_dict['transcripts'][t] = {}
            self.t_dict['transcripts'][t]["exons"] = {}
            self.t_dict['transcripts'][t]['list_of_exons'] = []
            # Gene sequence main coordinates are required to take introns
            # Transcript coordinates wanted for output  
            genomic_start = 0
            genomic_end = 0
            for e in items.iter('exon'):
                e_num = e.attrib['label']
                if e_num[-1] in ('a', 'b', 'c', 'd'):
                    # print exon_number
                    e_num = e_num[:-1]
                    print 'exon number: %d' % e_num
                e_num = int(e_num)
                self.t_dict['transcripts'][t]['list_of_exons'].append(e_num)
                self.t_dict['transcripts'][t]["exons"][e_num] = {}
                for coordinates in e:
                    if coordinates.attrib['coord_system'][-2] not in ['t', 'p']:
                        genomic_start = int(coordinates.attrib['start'])
                        genomic_end = int(coordinates.attrib['end'])
                assert genomic_start >= 0, "Exon index out of bounds"
                self.t_dict['transcripts'][t]["exons"][e_num]['genomic_start']=\
                    genomic_start
                self.t_dict['transcripts'][t]["exons"][e_num]['genomic_end'] =\
                    genomic_end
                self.t_dict['offset'] = self.padding
                exon_seq = list(self.sequence[
                                genomic_start - self.padding:
                                genomic_end + self.padding
                                ])
                self.t_dict['transcripts'][t]["exons"][e_num]['padded seq'] = \
                    exon_seq
                self.t_dict['transcripts'][t]["exons"][e_num]['length'] = \
                    genomic_end-genomic_start
                self.t_dict['transcripts'][t]["exons"][e_num]['padded length']=\
                    len(exon_seq)

    def find_cds_delay(self, t):
        """ Method to find the actual start of the translated sequence
            introduced to sort out non-coding exon problems """
        offset_total = 0
        offset = self.t_dict['transcripts'][t]['old_cds_offset']
        for e in self.t_dict['transcripts'][t]['list_of_exons']:
            g_start = self.t_dict['transcripts'][t]['exons'][e]['genomic_start']
            g_stop = self.t_dict['transcripts'][t]['exons'][e]['genomic_end']
            if offset > g_stop:
                self.t_dict['transcripts'][t]['exons'][e]['cds'] = 'before'
                offset_total = offset_total + (g_stop - g_start) + 1
            elif g_stop >= offset >= g_start:
                self.t_dict['transcripts'][t]['cds_offset'] = offset_total + \
                                                              (offset - g_start)
                self.t_dict['transcripts'][t]['exons'][e]['cds'] = 'after'
            elif offset < g_start:
                self.t_dict['transcripts'][t]['exons'][e]['cds'] = 'after'


    def get_nm(self):
        annotation_sets = self.t_dict['updatable'].findall('annotation_set')
        for annotation_set in annotation_sets:
            if annotation_set.attrib['type'] == 'ncbi':
                features = annotation_set.find('features')
                genes = features.findall('gene')
                for gene in genes:
                    transcripts = gene.findall('transcript')
                    for t_block in transcripts:
                        try:
                            t_no = t_block.attrib['fixed_id'][1:]
                            # print transcript_block.attrib['accession']
                            self.t_dict['transcripts'][int(t_no)]['NM_number']=\
                                t_block.attrib['accession']
                        except KeyError:
                            pass
                            # print 'found redundant transcript'

    def get_protein_exons(self):
        """ Collects full protein sequence for the appropriate transcript """
        for item in self.t_dict['fixannot'].findall('transcript'):
            p_number = int(item.attrib['name'][1:])
            coding_region = item.find('coding_region')
            coordinates = coding_region.find('coordinates')
            self.t_dict['transcripts'][p_number]['old_cds_offset'] = \
                int(coordinates.attrib['start'])
            translation = coding_region.find('translation')
            sequence = translation.find('sequence').text
            self.t_dict['transcripts'][p_number]['protein_length'] = \
                len(sequence)*3

    def run(self, padding):
        self.padding = padding
        # Initial sequence grabbing and populating dictionaries
        self.sequence = self.grab_element('fixed_annotation/sequence')
        self.t_dict['full sequence'] = list(self.sequence)
        self.get_exon_coords()
        self.get_nm()
        self.get_protein_exons()

        for t in self.t_dict['transcripts'].keys():
            self.find_cds_delay(t)
            self.t_dict['transcripts'][t]['list_of_exons'].sort(key=float)

        return self.t_dict

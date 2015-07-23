import os, cPickle
from subprocess import call
__author__ = 'mwelland'

class Aligner:

    def __init__(self, genelist):

        self.file_list = os.listdir('fastQs')
        self.read1_list= []
        self.read2_list= []
        self.read1name = 'all_refsR1.fastq'
        self.read2name = 'all_refsR2.fastq'
        

    def run(self):
        
        self.separate_files()
        self.add_file_contents_to_single_file()

        command_args = ['bwa', 'mem', '-t', 4, file_pair[0], file_pair[1]]
        output = '../SAMs/%s.sam' % gene
        call(command_args, stdout=output)
    
    @staticmethod        
    def add_file_contents_to_single_file():
        """ Takes each file from each read numbered list and adds all the contents
            to a single file
        """
        with open(self.read1name, 'w') as outfile:
            for fname in self.read1_list:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
        with open(self.read2name, 'w') as outfile:
            for fname in self.read2_list:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
    
    @staticmethod
    def separate_files():
        """ Read through the list of filenames and allocate each file to the
            appropriate list (read 1 or 2)
        """
        for filename in self.file_list:
            if filename.split('.')[-1] == 1:
                self.read1_list.append(filename)
            elif filename.split('.')[-1] == 2:
                self.read2_list.append(filename)
            else:
                print 'problem with filename %s' % filename

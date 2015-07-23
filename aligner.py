import os, cPickle
from subprocess import call
__author__ = 'mwelland'

class Aligner:

    def __init__(self, genelist):

        self.file_list = os.listdir('fastQs')

    def run(self):

        genelist_file = os.path.join('pickles', 'genelist.cPickle')
        with open(genelist_file, 'rU') as handle:
            genelist = cPickle.load(handle)

        for gene in genelist:
            file_pair = [name for name in self.file_list if name.split('_')[0] == gene]
            assert len(file_pair) == 2

            command_args = ['bwa', 'mem', '-t', 4, file_pair[0], file_pair[1]]
            output = '../SAMs/%s.sam' % gene
            call(command_args, stdout=output)
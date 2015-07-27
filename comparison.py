import cPickle
import os
import csv
__author__ = 'mwelland'


class VCF_Comparison:

    def __init__(self, run_number, pickle_dir, vcf_dir):
        self.run_number = run_number
        self.pickle_dir = pickle_dir
        self.vcf_dir = vcf_dir
        with open(os.path.join(self.pickle_dir, 'genelist.cPickle'), 'rU') as handle:
            self.genelist = cPickle.load(handle)

    def run(self):
        #for gene in self.genelist:
        print self.genelist[0]
        self.check_gene(self.genelist[0])

    def check_gene(self, gene):

        with open(os.path.join(self.pickle_dir, gene+'.cPickle'), 'rU') as handle:
            pickledict = cPickle.load(handle)


    def recursive_keys(self, dict):
        """
        Shitey method for printing out lists of keys
        :param dict:
        :return:
        """
        try:
            print str(dict.keys())
            for x in dict.keys():
                self.recursive_keys(dict[x])
                if x == 'hgvs':
                    print dict[x]
        except AttributeError:
            pass
            #print 'found a non-key'
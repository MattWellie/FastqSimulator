from excel_out_comparison import VCFComparison
import os

__author__ = 'mwelland'
run_number = 248
pickle_dict = os.path.join('pickles', '{}_variants.cPickle'.format(run_number))

vcomp = VCFComparison(run_number, pickle_dict, 'VCFs')
vcomp.run()

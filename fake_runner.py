from comparison import VCF_Comparison
__author__ = 'mwelland'

"""
Quick and dirty class for running the comparison part of the suite. This only needs the run number to be supplied as an
argument to get the process working
"""
vcomp = VCF_Comparison(707, 'pickles', 'VCFs')
vcomp.run()

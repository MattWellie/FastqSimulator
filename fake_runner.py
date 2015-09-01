"""
Quick and dirty class for running the comparison part of the suite. This only needs the run number to be supplied as an
argument to get the process working
"""

__author__ = 'mwelland'

from comparison import VCF_Comparison

vcomp = VCF_Comparison(707, 'pickles', 'VCFs')
vcomp.run()

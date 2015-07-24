import os
import cPickle
import random
from GbkParser_fq import GbkParser
from LrgParser_fq import LrgParser
from sequence_modifier import Modifier
from fq_sampler import Sampler
from read_condenser import Condenser
from aligner import Aligner

genelist = []
sam_directory = os.path.join(os.getcwd(), 'SAMs')
run_number = random.randint(1, 1000)
print 'This is Run number %d' % run_number
output_name = 'run_num_%d' % run_number
reference = '/DATA/references/hg19.fa'
#  these numbers will represent coordinates, and will be passed to the sampler class
#  X will be incremented up to a set value, then it will be reduced to 1 and Y will increase by 1
#  X will take on the value of Y + 1 to prevent coordinate clashes
#  This will create unique co-ordinates (hopefully)
#  Threshold may need to be increased for increased numbers of references
x_coord = 1
y_coord = 1

def check_file_type(file_name):
    """ This function takes the file name which has been selected
        as input. This will identify .xml and .gk/gbk files, and
        will print an error message and exit the application if
        a file is used which does not match either of these types
    """
    if file_name[-4:] == '.xml':
        return 'lrg'
    elif file_name[-3:] == '.gb':
        return 'gbk'
    elif file_name[-4:] == '.gbk':
        return 'gbk'
    else:
        print 'This program only works for GenBank and LRG files'
        print file_name
        exit()

#  This value will dictate how much sequence is taken from either side of the exons
padding = 150

input_files = os.listdir('input')
#  Clear the directory of old FastQs
#  These break the final condensing of files if left in the directory
file_list = os.listdir('fastQs')
for filename in file_list:
    os.remove(os.path.join('fastQs', filename))

for filename in input_files:
    try:
        file_type = check_file_type(filename)

        #  Use the PARSER project parsers to read required file contents
        #  Parsers have been modified (e.g. removing protein sequence)
        dictionary = {}
        if file_type == 'gbk':
            print 'Running parser'
            gbk_reader = GbkParser(os.path.join('input', filename), padding)
            dictionary = gbk_reader.run()
            parser_details = gbk_reader.get_version

        elif file_type == 'lrg':
            print 'Running parser'
            lrg_reader = LrgParser(os.path.join('input', filename))
            dictionary = lrg_reader.run(padding)
            parser_details = lrg_reader.get_version

        print 'FILE PARSED'

        #  Create a modifier instance and modify the dictionary
        #  This has been written to introduce one variant into each exon
        modifier = Modifier(dictionary)
        new_dict = modifier.run_modifier()
        #  Dump a copy of the changed dictionary using cPickle (troubleshooting)
        with open(os.path.join('pickles', dictionary['genename']+'_mod.cPickle'), 'wb') as handle:
            cPickle.dump(new_dict, handle)
        with open(os.path.join('pickles', dictionary['genename']+'_std.cPickle'), 'wb') as handle:
            cPickle.dump(dictionary, handle)
        print 'dict modified'

        #  Keep a record of the genes which have been processed
        #  This is used later on when combining files
        genelist.append(dictionary['genename'])

        #  Create a sampler instance to extract simulated reads
        sampler = Sampler(dictionary, new_dict,x_coord, y_coord)
        x_coord, y_coord = sampler.run()
    except TypeError:
        print 'problem with file %s' % filename

#  Dump a copy of the gene list
with open(os.path.join('pickles', 'genelist.cPickle'), 'wb') as handle:
    cPickle.dump(genelist, handle)

#  Create a condenser instance, and mix each variant transcript with the reference version
#  This is saved as a new .fq file
file_condenser = Condenser(genelist)
file_condenser.run()
aligner = Aligner(sam_directory, output_name, reference)
aligner.run()

print 'Run %s completed' % run_number
print 'successes:'
print genelist

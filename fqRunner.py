import os
import cPickle
import random
# Import classes from local files
from GbkParser_fq import GbkParser
from LrgParser_fq import LrgParser
from sequence_modifier import Modifier
from fq_sampler import Sampler
from read_condenser import Condenser
from aligner import Aligner
from subprocess import call
from comparison import VCF_Comparison

"""
This file controls the operation of the fastq generator program, and can be used to operate the full pipeline and
 results comparisons, based on specified arguments
"""

geneset = set()
sam_directory = os.path.join(os.getcwd(), 'SAMs')
# Randomly generate a new run number to keep each instance identifiable
# This is necessary as the variants created are unique to a run number
# The variant file and corresponding .fq will correspond to a given run number, included in the file names
run_number = random.randint(1, 1000)
print 'This is Run number %d' % run_number
output_name = 'run_num_%d' % run_number
vcf_name = 'run_%d.vcf' % run_number

# Location of required system resources
# These are not required to just use the fastq generator component of the program
reference = '/DATA/references/hg19.fa'
annovar = '/DATA/annovar/table_annovar.pl'
anno_db = '/DATA/annovar/humandb/'

#  Arguments in order are: reference, input, output
variant_call_string = 'samtools mpileup -g -f %s %s'
var_call_filter_string = 'bcftools call -vc %s'
anno_string = '%s --vcfinput %s %s --buildver hg19 --out %s --remove --protocol refGene --operation g --nastring .'

# An object to hold the names of all input files which fail
fail_list = []

# This value will dictate how much sequence is taken from either side of the exons
padding = 150

# These numbers will represent coordinates, and will be passed to the sampler class
# Thi will allow for reads to appear unique across a lareg number of input genes
# X will be incremented up to a set value, then it will be reduced to 1 and Y will increase by 1
# X will take on the value of Y + 1 to prevent coordinate clashes
# This will create unique co-ordinates (hopefully)
# Threshold may need to be increased for increased numbers of references
x_coord = 1
y_coord = 1
tile = 1


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

input_files = os.listdir('input')
# Clear the directory of old FastQs
# These break the final condensing of files if left in the directory
file_list = os.listdir('fastQs')
for filename in file_list:
    os.remove(os.path.join('fastQs', filename))
file_list = os.listdir('pickles')
for filename in file_list:
    os.remove(os.path.join('pickles', filename))
file_list = os.listdir('VCFs')
for filename in file_list:
    os.remove(os.path.join('VCFs', filename))

for filename in input_files:
    try:
        print 'File name: %s' % filename
        file_type = check_file_type(filename)

        # Use the PARSER project .gb and .lrg parsers to read required file contents
        # Parsers have been modified (e.g. removing protein sequence)
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

        print "FILE '%s' PARSED" % filename

        # Create a modifier instance and modify the dictionary
        # This has been written to introduce one variant into each exon
        modifier = Modifier(dictionary, file_type)
        new_dict = modifier.run_modifier()

        # Dump a copy of the changed dictionary using cPickle (troubleshooting/re-running)
        with open(os.path.join('pickles', '%s.cPickle' % dictionary['genename']), 'wb') as handle:
            cPickle.dump(new_dict, handle)

        # Keep a record of the genes which have been processed
        # This is used later on when combining files
        if dictionary['genename'] not in geneset:
            geneset.add(dictionary['genename'])

        # Create a sampler instance to extract simulated reads
        # This extracts reads for both the altered and unaltered dictionaries
        sampler = Sampler(dictionary, new_dict, x_coord, y_coord, tile)
        x_coord, y_coord, tile = sampler.run()

    # Allow for errors in parsing files without crashing program
    except TypeError:
        fail_list.append(filename)
        print 'problem with file %s' % filename
    except AttributeError:
        fail_list.append(filename)
        print 'problem with file %s' % filename

# Dump a copy of the gene list
with open(os.path.join('pickles', 'genelist.cPickle'), 'wb') as handle:
    cPickle.dump(geneset, handle)

# Create a condenser instance, and mix each variant transcript with the reference version
# At this point the output should contain one fq representing the unaltered gene, and one for each transcript
# The R1 and R2 files for each of these is merged, creating a single pair or R1 & R2 for each transcript
# The current method means that each .fq pair represents heterozygous variations, saved as a new .fq file
file_condenser = Condenser(geneset)
file_condenser.run()

# Creates an aligner instance and converts the multiple fq files into a single pair of files for conversion
aligner = Aligner(sam_directory, output_name, reference)
bam_filename = aligner.run()
bam_location = os.path.join('fastQs', bam_filename)

"""
This section is for the variant calling on the aligned files. From this point on, the process extends beyond the initial
goal of creating a fastq generator, and can be ignored/commented out depending on usage
Due to some aspect of the read generation, Platypus is unable to generate variant calls from the aligned input data.
The SAMtools mpileup feature, combined with the bcftools call function are used for the two-step cariant calling.
"""

temp_bcf = os.path.join('VCFs', 'temp.bcf')
vcf_location = os.path.join('VCFs', vcf_name)
variant_filled = variant_call_string % (reference, bam_location)
filled_filter = var_call_filter_string % temp_bcf

print variant_filled
with open(temp_bcf, 'w') as out_file:
    call(variant_filled.split(' '), stdout=out_file)

print 'filter command: %s' % filled_filter
with open(vcf_location, 'w') as vcf_out:
    call(filled_filter.split(' '), stdout=vcf_out)

# The next step is annotating the variant calls
anno_out = os.path.join('VCFs', '%d_anno.vcf' % run_number)
filled_anno = anno_string % (annovar, vcf_location, anno_db, anno_out)
call(filled_anno.split(' '))
os.remove(os.path.join('VCFs', 'temp.bcf'))

# And compare the VCF to the predicted:
vcf_comparison = VCF_Comparison(run_number, 'pickles', 'VCFs')
vcf_comparison.run()

print 'Run %s completed' % run_number
print 'successes:'
print geneset
if fail_list:
    print 'failures:'
    print fail_list
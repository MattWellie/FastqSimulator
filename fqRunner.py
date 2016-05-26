"""
The script used to run the FastQ simulation process. 

The program is designed to take input files in either GenBank or LRG formats. 
For each new run on the program, a number is randomly picked to be integrated into all file names for identification purposes.
From these input files, the process will use the transcript /exon details to isolate exonic regions of each transcript (+/- padding intron)
Using Sub-String methods, the program will 'walk' across the sequence segments in regions of a pre-defined length (hard-coded)
In each section (exon) a single base substitution will be made, with the intention of identifying which exons are 'seen' by the annotation software
For each inserted variant I do some magic to create the appropriate HGVS annotation, which is stored in a dictionary object
For each fragment, this will then simulate paired-end reads, writing the reads (complete with headers and quality strings) to R1 and R2 files

The program will then control the alignment, compression, indexing and variant calling of the files, with the files stored in appropriately names directories in the project root.
This currently uses hard coded references to the tools required (annovar, bwa..) during these steps
Annotation is minimal, only required to identify gene symbol and HGVS, no population data

Then the real magic happens

During the 'read' creation process variants were created and the HGVS annotations they were equivalent to were generated
The main process here is to compare the results from the VCF file with the 'expected' results 
On a per-gene basis (and per-transcript where appropriate), the VCF file and the expected results dictionary are parsed together
For every instance of matching HGVS being found, the VCF row is deleted and a match is recorded. 
At the end of parsing all the predictions, any remaining results in both the dict and VCF are classified as unmatched
An Excel report is written, with a summary page showing all mismatches, and an additional page to provide detail about each individual gene.
"""

import os
import cPickle
import random
from GbkParser_fq import GbkParser
from LrgParser_fq import LrgParser
from sequence_modifier import Modifier
from fq_sampler import Sampler
from read_condenser import Condenser
from aligner import Aligner
from subprocess import call
from excel_out_comparison import VCFComparison

geneset = set()
main_variant_dict = {}
sam_directory = os.path.join(os.getcwd(), 'SAMs')
run_number = random.randint(1, 1000)
print 'This is Run number %d' % run_number
output_name = 'run_num_%d' % run_number
vcf_name = 'run_%d.vcf' % run_number
reference = '/DATA/references/hg19.fa'
annovar = '/DATA/annovar/table_annovar.pl'
anno_db = '/DATA/annovar/humandb/'
#  Arguments in order are: reference, input, output
variant_call_string = 'samtools mpileup -g -f %s %s'
var_call_filter_string = 'bcftools call -vc %s'
anno_string = '%s --vcfinput %s %s --buildver hg19 --out %s --remove --protocol refGene --operation g --nastring .'
fail_list = []

#  these numbers will represent coordinates, and will be passed to the sampler class
#  X will be incremented up to a set value, then it will be reduced to 1 and Y will increase by 1
#  X will take on the value of Y + 1 to prevent coordinate clashes
#  This will create unique co-ordinates (hopefully)
#  Threshold may need to be increased for increased numbers of references
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

#  This value will dictate how much sequence is taken from either side of the exons
padding = 150

input_files = os.listdir('input')
#  Clear the directory of old FastQs
#  These break the final condensing of files if left in the directory
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
        modifier = Modifier(dictionary, file_type)
        new_dict = modifier.run_modifier()

        #  Keep a record of the genes which have been processed
        #  This is used later on when combining files
        if dictionary['genename'] not in geneset:
            geneset.add(dictionary['genename'])

        #  Create a sampler instance to extract simulated reads
        sampler = Sampler(dictionary, new_dict, x_coord, y_coord, tile)
        x_coord, y_coord, tile = sampler.run()

        main_variant_dict[dictionary['genename']] = {}
        for trans in new_dict:
            exon_range = new_dict[trans]['variants'].keys()
            # Extract the 'version free' transcript accession
            transcript = new_dict[trans]['variants'][exon_range[0]]['transcript'].split('.')[0]
            main_variant_dict[dictionary['genename']][transcript] = new_dict[trans]['variants']
        print 'dict modified'

    except TypeError:
        fail_list.append(filename)
        print 'problem with file %s' % filename
    except AttributeError:
        fail_list.append(filename)
        print 'problem with file %s' % filename

# Make sure everything in the gene list is accounted for in the variant dict
for gene in geneset:
    if gene not in main_variant_dict.keys():
        print 'Gene missing from variant dict: {}'.format(gene)
# Dump a copy of the main variant dictionary
variant_dict_pickle = os.path.join('pickles', '{}_variants.cPickle'.format(run_number))
with open(variant_dict_pickle, 'wb') as handle:
    cPickle.dump(main_variant_dict, handle)

# Count total number of inserted variants:
gene_count = 0
trans_count = 0
var_count = 0
for gene in main_variant_dict:
    gene_count += 1
    for transcript in main_variant_dict[gene]:
        trans_count += 1
        for exon in main_variant_dict[gene][transcript]:
            var_count += 1
print 'Genes: {0}\nTranscripts: {1}\nVariants: {2}'.format(gene_count, trans_count, var_count)

#  Create a condenser instance, and mix each variant transcript with the reference version
#  This is saved as a new .fq file
file_condenser = Condenser(geneset, run_number)
file_condenser.run()

#  Creates an aligner instance and converts the multiple fq files into a single pair of files for conversion
aligner = Aligner(sam_directory, output_name, reference, run_number)
bam_filename = aligner.run()
bam_location = os.path.join('fastQs', bam_filename)

"""
This section is for the variant calling on the aligned files.
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

#  The next step is annotating the variant calls
anno_out = os.path.join('VCFs', '%d_anno.vcf' % run_number)
filled_anno = anno_string % (annovar, vcf_location, anno_db, anno_out)
call(filled_anno.split(' '))
os.remove(os.path.join('VCFs', 'temp.bcf'))

# And compare the VCF to the predicted:
vcf_comparison = VCF_Comparison(run_number, variant_dict_pickle, 'VCFs')
vcf_comparison.run()

print 'Run %s completed' % run_number
if fail_list:
    print 'failures:'
    print fail_list

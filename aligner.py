"""
This class will transcribe the contents of every individual file's contents into a pair of files which will represent
reads 1 and 2 across all input reference transcripts. The output can be aligned in a single process to speed up the
creation of mapped output.
"""

__author__ = 'mwelland'

import os
from subprocess import call


class Aligner:

    def __init__(self, sam_directory, out_name, reference):
        self.file_list = os.listdir('fastQs')
        self.read1_list = []
        self.read2_list = []
        self.read1name = os.path.join('fastQs', 'all_refsR1.fq')
        self.read2name = os.path.join('fastQs', 'all_refsR2.fq')
        self.sam_directory = sam_directory
        self.out_name = out_name
        self.reference = reference

    def run(self):
        self.separate_files()
        self.add_file_contents_to_single_file()
        output_name = self.run_alignment()
        return output_name

    def add_file_contents_to_single_file(self):
        """
        Takes each file from each read numbered list and adds all the contents to a single file
        Delete the file after it has been written
        """
        print 'Combining file contents'
        with open(self.read1name, 'w') as outfile:
            for fname in self.read1_list:
                with open(os.path.join('fastQs', fname)) as infile:
                    for line in infile:
                        outfile.write(line)
                os.remove(os.path.join('fastQs', fname))

        with open(self.read2name, 'w') as outfile:
            for fname in self.read2_list:
                with open(os.path.join('fastQs', fname)) as infile:
                    for line in infile:
                        outfile.write(line)
                os.remove(os.path.join('fastQs', fname))

        # Each file should now be written into either R1 or R2
        # This method may conflict with the reference-transcript pairings

    def separate_files(self):
        """
        Read through the list of filenames and allocate each file to the
        appropriate list (read 1 or 2)
        """
        print 'Sorting fastQs'
        for filename in self.file_list:
            if filename.split('.')[0][-1] == '1':
                self.read1_list.append(filename)
            elif filename.split('.')[0][-1] == '2':
                self.read2_list.append(filename)
            else:
                print filename.split('.')[0][-1]
                print 'problem with filename %s' % filename
                this = raw_input()

    def run_alignment(self):
        """ Processes the command line alignment using BWA MEM
            Uses SAMtools to convert the file to BAM
        """
        print 'Starting the alignment'
        sam_name = os.path.join(self.sam_directory, self.out_name+'.sam')
        bam_name = os.path.join(self.sam_directory, self.out_name+'.bam')
        sorted_file = os.path.join(self.sam_directory, 'sorted'+self.out_name)

        with open(sam_name, 'w') as out_file:
            call(['bwa', 'mem', '-t', '4', self.reference, self.read1name, self.read2name], stdout=out_file)
        with open(bam_name, 'w') as bam_out:
            call(['samtools', 'view', '-b', sam_name], stdout=bam_out)
        call(['samtools', 'sort', bam_name, sorted_file])
        call(['samtools', 'index', sorted_file+'.bam'])
        return os.path.join(self.sam_directory, sorted_file+'.bam')

    def clear_out_sams(self):
        file_list = [x for x in os.listdir(self.sam_directory)
                     if x[:6] != 'sorted']
        for x in file_list:
            os.remove(os.path.join(self.sam_directory, x))
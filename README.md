# FastqSimulator
A program for taking reference sequences and creating fake fastq reads

## Process

- Using the parser modules created for the reference sequence generator project, this program takes input files in genbank (representing a genomic region, NG_ accession) or LRG file formats. It extracts coordinates and sequences from these files and stores the contents in a dictionary. 
- For each transcript found in the input file, a version of the original dictionary is modified using the sequence_modifier class. This class uses the full genomic sequence and the exon coordinates to create artificial variants using the following algorithm:
- - ddddd
- - d

 

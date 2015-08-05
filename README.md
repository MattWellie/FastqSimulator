# FastqSimulator
A program for taking reference sequences and creating fake fastq reads

## Process

- Using the parser modules created for the reference sequence generator project, this program takes input files in genbank (representing a genomic region, NG_ accession) or LRG file formats. It extracts coordinates and sequences from these files and stores the contents in a dictionary. 
- For each transcript found in the input file, a version of the original dictionary is modified using the sequence_modifier class. This class uses the full genomic sequence and the exon coordinates to create artificial variants using the following algorithm:
  - Find *start* and *end* coordinates and generate a random number between the two positions as the coordinate to change
  - Use a random selecter to pick one of A, C, G, or T (reselect if the new base is the same)
  - substitute the base in the genomic sequence with the new one
  - Use a combination of the old base, the new base, the coordinate position relative to the start of the coding sequence and the length of the protein sequence to  predict HGVS nomenclature of the variant annotation 
- Once the sequence has been modified, the exon sequence (along with a region of padding either side to allow for decent read coverage) is extracted, along with details about the exon and the variant. These are put in a new dictionary and returned.



# Known issues

- Annovar requires a separate database for Ensembl gene IDs. Use of this database (ensGene) is not currently in the code
- for non-coding exons, the variants may be mismatched by 1 base. This may be due to differences between LRG and GB coordinates, or may be a different issue. Hard coding a change to the offset will fix some and break others... investigate more.

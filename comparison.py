import cPickle
import os
import re
__author__ = 'mwelland'


class VCF_Comparison:

    def __init__(self, run_number, pickle_dir, vcf_dir):
        self.vcf_file = os.path.join(vcf_dir, '%d_anno.vcf.hg19_multianno.vcf' % run_number)
        self.pickle_dir = pickle_dir
        self.pickles = os.listdir(os.path.join(self.pickle_dir))
        self.vcf = {}
        self.tempvcf = os.path.join(vcf_dir, 'tempout.vcf')
        with open(os.path.join(self.pickle_dir, 'genelist.cPickle'), 'rU') as handle:
            self.geneset = cPickle.load(handle)
        self.titles = 'CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,OTHER\n'
        self.matches = 0
        self.rows = 0
        self.variants = 0

    def run(self):
        self.squish_vcf()
        self.open_vcf()
        for gene in self.geneset:
            self.check_gene(gene)
        if self.matches != self.rows:
            print 'Total variants counted: %s' % self.variants
            print 'Total variants in VCF: %d' % self.rows
            print 'Total matches: %d' % self.matches
        else:
            print 'All variants found'
        os.remove(self.tempvcf)

    def squish_vcf(self):
        """
        This mini method just writes out only the non-header information from the original vcf into a new file
        The file is written to a new output to make sure that it can be read again if required
        The output is written in CSV format so that the csv.DictWriter method can be used
        """
        with open(self.vcf_file, 'rU') as input_vcf:
            with open(self.tempvcf, 'wb') as output_vcf:
                for line in input_vcf:
                    if line[0] == '#':
                        pass
                    else:
                        output_vcf.write(line)
                        self.rows += 1

    def open_vcf(self):
        """
        Add all contents from the VCF into a dictionary object which can be sorted through by gene
        Use regex to capture the gene name, create a dictionary index which is the gene name (if not already an index)
        Add the row to a list in the dictionary
        Might be best to treat the whole 'INFO' block as a single string, as different variants are annotated in
        different columns, depending on whether they are 5'UTR, 3'UTR or exonic...
        Ready to begin matching against pickled contents
        """
        with open(self.tempvcf) as csvfile:

            for row in csvfile:
                search_string = row.split('\t')[7]
                match = re.search(';Gene\.refGene=(?P<gene_name>,?.*?);', search_string)
                if match:
                    gene = match.group('gene_name')
                    if gene in self.vcf:
                        self.vcf[gene].append(search_string)
                    else:
                        self.vcf[gene] = [search_string]
            for gene in self.vcf:
                print 'Gene located: %s, %d rows' % (gene, len(self.vcf[gene]))

    def check_gene(self, gene):
        gene_pickles = [x for x in self.pickles if x.split('.')[0] == gene]
        try:
            gene_vcf = self.vcf[gene]
            for gene_file in gene_pickles:
                with open(os.path.join(self.pickle_dir, gene_file), 'rU') as handle:
                    pickledict = cPickle.load(handle)
                for transcript in pickledict:
                    vars = pickledict[transcript]['variants']
                    exons = vars.keys()
                    for exon in exons:
                        self.variants += 1
                        transcript = vars[exon]['transcript'].split('.')[0]
                        hgvs = vars[exon]['hgvs']
                        found = False
                        if hgvs[2] == '-' or hgvs[2] == '*':
                            variant = '%s:%s:%s' % (gene, transcript, hgvs)
                            for row in gene_vcf:
                                match = re.search('(%s:)?%s:%s' % (gene, transcript, hgvs), row)
                                if match:
                                    found = True
                                    self.matches += 1
                        else:
                            variant = '%s:%s:exon%d:%s' % (gene, transcript, exon, hgvs)
                            for row in gene_vcf:
                                match = re.search('(%s:)?%s:(exon%d)?:%s' % (gene, transcript, exon, hgvs), row)
                                if match:
                                    found = True
                                    self.matches += 1
                        if not found:
                            print 'This variant was not found: %s' % variant
        except KeyError:
            print '%s not found as a key' % gene
            this = raw_input()

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
        self.missing = {}
        self.tempvcf = os.path.join(vcf_dir, 'tempout.vcf')
        with open(os.path.join(self.pickle_dir, 'genelist.cPickle'), 'rU') as handle:
            self.geneset = cPickle.load(handle)
        self.titles = 'CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,OTHER\n'
        self.matches = 0
        self.variants = 0
        self.unmatched_predictions = 0

    def run(self):
        self.squish_vcf()
        self.open_vcf()
        for gene in self.geneset:
            self.missing[gene] = []
            self.check_gene(gene)

        os.remove(self.tempvcf)
        print 'Remaining entries:'
        self.add_missing_variants()
        perfect_genes = []
        for gene in self.missing:
            if self.missing[gene]:
                print gene
                for row in self.missing[gene]:
                    print row
            else:
                perfect_genes.append(gene)
        if self.matches != self.variants:
            print 'Total variants counted: %s' % self.variants
            print 'Total matches: %d' % self.matches
            print 'Predicted and not found: %d' % self.unmatched_predictions
            print 'Perfect genes: %s' % str(perfect_genes)
        else:
            print 'All variants found'

    def add_missing_variants(self):
        for gene in self.vcf:
            for row in self.vcf[gene]:
                # Search for specific string in complete row
                # Example row:
                # GeneDetail.refGene=NM_002506:c.-6897A>G
                m = re.search('GeneDetail.refGene=(?P<HGVS>NM.*?);', row)
                if m:
                    try:
                        self.missing[gene].append('In VCF, not expected: %s: %s' % (gene, m.group('HGVS')))
                    except KeyError:
                        try:
                            self.missing[gene].append('In VCF, not expected: %s: %s' % (gene.split(',')[0], m.group('HGVS')))
                        except KeyError:
                            try:
                                self.missing[gene].append('In VCF, not expected: %s: %s' % (gene.split(',')[1], m.group('HGVS')))
                            except:
                                print 'Theres an error with a gene ID'

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
                else:
                    print "couldn't match the variant in %s" % row

    def check_gene(self, gene):
        gene_pickles = [x for x in self.pickles if x.split('.')[0] == gene]
        try:
            rows_to_delete = []
            gene_vcf = self.vcf[gene]
            rows = range(len(gene_vcf))
            for gene_file in gene_pickles:
                with open(os.path.join(self.pickle_dir, gene_file), 'rU') as handle:
                    pickledict = cPickle.load(handle)
                for trans in pickledict:
                    vars = pickledict[trans]['variants']
                    exons = vars.keys()
                    for exon in exons:
                        self.variants += 1
                        transcript = vars[exon]['transcript'].split('.')[0]
                        hgvs = vars[exon]['hgvs']
                        found = False
                        if hgvs[2] == '-' or hgvs[2] == '*':
                            variant = '{0}:{1}:{2}'.format(gene, transcript, hgvs)
                            for row in rows:
                                match = re.search('({0}:)?{1}:{2}'.format(gene, transcript, hgvs), gene_vcf[row])
                                if match:
                                    rows_to_delete.append(row)
                                    found = True
                                    self.matches += 1
                        else:
                            variant = '{0}:{1}:exon{2}:{3}'.format(gene, transcript, exon, hgvs)
                            for row in rows:
                                # match = re.search('(%s:)?%s:exon.{1,3}?:%s' % (gene, transcript, hgvs), row)
                                match = re.search('({0}:)?{1}:.*?:{2}'.format(gene, transcript, hgvs), gene_vcf[row])
                                # match = re.search('(%s:)?%s:exon%d?:%s' % (gene, transcript, exon, hgvs), row)
                                if match:
                                    rows_to_delete.append(row)
                                    found = True
                                    self.matches += 1
                        if not found:
                            self.missing[gene].append('Predicted, not found in VCF: %s' % variant)
                            self.unmatched_predictions += 1
            # Delete any rows which have been matched against
            # This is done in reverse, high indexes first
            # From low to high means the list shrinks and high indexes are invalid
            # print 'Delete list: %s' % sorted(rows_to_delete, reverse=True)
            for row in sorted(rows_to_delete, reverse=True):
                try:
                    del gene_vcf[row]
                except IndexError:
                    print 'problem with this list: %s' % str(sorted(rows_to_delete, reverse=True))
                    print 'Index: %d' % row
                    print 'vcf length: %d' % len(gene_vcf)
                    this = raw_input()
            #self.unmatched_predictions += len(gene_vcf)
            self.vcf[gene] = gene_vcf
        except KeyError:
            print 'Gene %s not found as a key' % gene
            this = raw_input()

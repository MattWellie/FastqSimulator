import cPickle
import os
import re
from openpyxl import load_workbook
import xlsxwriter

__author__ = 'mwelland'


class VCFComparison:

    def __init__(self, run_number, variant_dict, vcf_dir):
        self.vcf_file = os.path.join(vcf_dir, '{}_anno.vcf.hg19_multianno.vcf'.format(run_number))
        self.run_number = run_number
        with open(variant_dict, 'r') as handle:
            self.variant_dict = cPickle.load(handle)
        self.genes = self.variant_dict.keys()
        self.vcf = {}
        self.results = {}
        self.tempvcf = os.path.join(vcf_dir, 'tempout.vcf')
        self.matches = 0
        self.variants = 0
        self.unmatched_predictions = 0
        self.excel_dir = 'Excels'
        self.mystery_genes = set()
        self.transcripts = 0
        self.perfect_genes = 0

    def run(self):
        self.squish_vcf()
        self.open_vcf()
        for gene in self.genes:
            self.results[gene] = {'found':     set(),
                                  'not_found': {'in_vcf': set(),
                                                'in_fq':  set()}}
            self.check_gene(gene)

        os.remove(self.tempvcf)
        self.add_missing_variants()
        # Print all the output stuff to an excel document
        for gene in self.results:
            perfect = 0
            for section in self.results[gene]['not_found']:
                perfect += len(self.results[gene]['not_found'][section])
            if perfect == 0:
                self.perfect_genes += 1
        self.write_excel()
        if self.matches != self.variants:
            print 'Total variants counted: {}'.format(self.variants)
            print 'Total matches: {}'.format(self.matches)
            print 'Predicted and not found:{}'.format(self.unmatched_predictions)
            print 'Perfect genes: {}'.format(self.perfect_genes)
        else:
            print 'All variants found'

    def write_excel(self):
        # This method will take the results from the process and output to an excel
        # There will be a summary page to condense the main details of the comparison
        # Each gene will have a further page to describe results in detail
        excel_out_name = os.path.join(self.excel_dir, 'run_{}_results.xlsx'.format(self.run_number))
        workbook = xlsxwriter.Workbook(excel_out_name)
        format_bold = workbook.add_format({'bold': True})
        format_matched = workbook.add_format({'bg_color': '#ADFF2F'})
        format_missing_db = workbook.add_format({'bg_color': '#F4A460'})
        format_missing_excel = workbook.add_format({'bg_color': '#F08080'})
        format_hyperlink = workbook.add_format({'font_color': '#0000FF'})

        worksheet = workbook.add_worksheet('Summary')
        worksheet.set_column(0, 0, 20)
        worksheet.set_column(2, 2, 17)
        row = 0
        col = 0
        worksheet.write(row, col, 'Summary Page', format_bold); row =+ 2
        worksheet.write(row, 0, 'Genes featured:', format_bold)
        worksheet.write(row, 1, '{}'.format(len(self.genes)), format_bold)
        row += 1
        worksheet.write(row, 0, 'Transcripts featured:', format_bold)
        worksheet.write(row, 1, '{}'.format(self.transcripts), format_bold)
        row += 1
        worksheet.write(row, 0, 'Perfect genes:', format_bold)
        worksheet.write(row, 1, '{}'.format(self.perfect_genes), format_bold)
        row += 1
        worksheet.write(row, 0, 'Total Variants:', format_bold)
        worksheet.write(row, 1, '{}'.format(self.variants), format_bold)
        row += 1
        worksheet.write(row, 0, 'Variants Matched:', format_bold)
        worksheet.write(row, 1, '{}'.format(self.matches), format_bold)
        row += 1
        worksheet.write(row, 0, 'Variants not Matched:', format_bold)
        worksheet.write(row, 1, '{}'.format(self.unmatched_predictions), format_bold)
        row += 2
        worksheet.write(row, 0, 'Mismatches by Gene:', format_missing_excel)
        row += 1
        worksheet.write(row, 0, 'Gene', format_bold)
        worksheet.write(row, 1, 'FastQ Prediction', format_bold)
        worksheet.write(row, 2, 'VCF Prediction', format_bold)
        highest_row = row + 1
        for gene in self.results:
            worksheet.write(highest_row, 0, gene, format_bold); row += 1
            fq_row = highest_row
            vcf_row = highest_row
            if self.results[gene]['not_found']['in_fq']:
                for result in self.results[gene]['not_found']['in_fq']:
                    worksheet.write(fq_row, 1, result); fq_row += 1
                fq_row += 1
            if self.results[gene]['not_found']['in_vcf']:
                for result in self.results[gene]['not_found']['in_vcf']:
                    worksheet.write(vcf_row, 2, result); vcf_row += 1
                vcf_row += 1
            if vcf_row > fq_row:
                highest_row = vcf_row
            else:
                highest_row = fq_row
        worksheet.set_column(1, 1, 45)
        worksheet.set_column(2, 2, 100)


        for gene in self.results:
            matches = len(self.results[gene]['found'])
            mismatches = 0
            for section in self.results[gene]['not_found']:
                mismatches += len(self.results[gene]['not_found'][section])
            total = mismatches + matches

            worksheet = workbook.add_worksheet(gene)
            row = 0
            col = 0
            worksheet.write(row, col, gene, format_bold); row =+ 2
            worksheet.write(row, col, 'Total Variants:', format_bold); col += 1
            worksheet.write(row, col, '{}'.format(total), format_bold); row += 1; col -= 1
            worksheet.write(row, col, 'Matched:', format_bold); col += 1
            worksheet.write(row, col, '{}'.format(matches), format_bold); row += 1; col -= 1
            worksheet.write(row, col, 'Not Matched:', format_bold); col += 1
            worksheet.write(row, col, '{}'.format(mismatches), format_bold); row += 1

            row += 2

            worksheet.write(row, col, 'Variants Matched:', format_matched)
            row += 1
            for variant in self.results[gene]['found']:
                worksheet.write(row, col, variant, format_matched)
                row += 1
            row += 2

            if self.results[gene]['not_found']['in_vcf'] or self.results[gene]['not_found']['in_fq']:
                worksheet.write(row, col, 'Unmatched Variants:', format_missing_excel)
                row += 1
                if self.results[gene]['not_found']['in_fq']:
                    worksheet.write(row, col, 'Predicted:', format_missing_excel); row += 1
                    for variant in self.results[gene]['not_found']['in_fq']:
                        worksheet.write(row, col, variant, format_missing_excel); row += 1
                    row += 2
                else:
                    worksheet.write(row, col, 'No Predicted Variants:', format_missing_db); row += 2
                if self.results[gene]['not_found']['in_vcf']:
                    worksheet.write(row, col, 'Unexpected:', format_missing_excel); row += 1
                    for variant in self.results[gene]['not_found']['in_vcf']:
                        worksheet.write(row, col, variant, format_missing_excel); row += 1
                    row += 2
                else:
                    worksheet.write(row, col, 'No Unexpected Variants:', format_missing_db); row += 2

            else:
                worksheet.write(row, col, 'No Unmatched Variants:', format_missing_db)
            worksheet.set_column(0, 0, 15)
            worksheet.set_column(1, 1, 40)

        workbook.close()

    def add_missing_variants(self):
        for gene in self.vcf:
            for row in self.vcf[gene]:
                # Search for specific string in complete row
                # RegEx required as columns are not always in order
                # GeneDetail.refGene=NM_002506:c.-6897A>G
                active_match = 'Unknown Variant'
                matched = False
                if 'GeneDetail.refGene=.;' in row:
                    if 'AAChange.refGene=.;' not in row:
                        a = re.search('AAChange.refGene=.*?:(?P<HGVS>NM_.*?);', row)
                        b = re.search('AAChange.refGene=(?P<HGVS>NM.*?);', row)
                        if a:
                            active_match = a
                            matched = True
                        elif b:
                            active_match = b
                            matched = True
                else:
                    a = re.search('GeneDetail.refGene=.*?:(?P<HGVS>NM_.*?);', row)
                    b = re.search('GeneDetail.refGene=(?P<HGVS>NM_.*?);', row)
                    if a:
                        active_match = a
                        matched = True
                    elif b:
                        active_match = b
                        matched = True
                if matched:
                    filtered_list = self.filter_matches(active_match.group('HGVS'), gene)
                    self.results[gene]['not_found']['in_vcf'].add(', '.join(filtered_list))
                else:
                    if gene in self.results:
                        self.results[gene]['not_found']['in_vcf'].add('Variant unknown')
                    else:
                        self.mystery_genes.add(gene)

    def filter_matches(self, string, gene):
        output_list = []
        for element in string.split(';'):
            nm_number = element.split(':')
            if nm_number in self.variant_dict[gene].keys():
                output_list.append(element)
        if not output_list:
            output_list.append(string.split(';')[0])
        return output_list


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
        rows_to_delete = []
        gene_vcf = self.vcf[gene]
        rows = range(len(gene_vcf))
        for transcript in self.variant_dict[gene]:
            self.transcripts += 1
            variants = self.variant_dict[gene][transcript]
            exons = variants.keys()
            for exon in exons:
                self.variants += 1
                hgvs = variants[exon]['hgvs']
                found = False
                # If the variant is 3' UTR or 5' UTR, e.g. c.-69A>C:
                if hgvs[2] == '-' or hgvs[2] == '*':
                    # Use the exact sequence predicted to write the gene
                    variant = '{0}:{1}:{2}'.format(gene, transcript, hgvs)
                    for row in rows:
                        match = re.search('({0}:)?{1}:{2}'.format(gene, transcript, hgvs), gene_vcf[row])
                        if match:
                            rows_to_delete.append(row)
                            found = True
                            self.matches += 1
                            self.results[gene]['found'].add(variant)
                            break
                else:
                    # Use the exact sequence predicted to write the gene
                    variant = '{0}:{1}:exon{2}:{3}'.format(gene, transcript, exon, hgvs)
                    for row in rows:
                        match = re.search('({0}:)?{1}:.*?:{2}'.format(gene, transcript, hgvs), gene_vcf[row])
                        if match:
                            rows_to_delete.append(row)
                            found = True
                            self.matches += 1
                            self.results[gene]['found'].add(variant)
                            break
                """
                This section will allow matches to be made which are less specific.
                This may not be useful if exon numbers are required, but exon numbers
                may change between systems more easily than variant nomenclature.
                Matching only on nomenclature should be fine for this project.
                """
                if not found:
                    for row in rows:
                        if hgvs in gene_vcf[row]:
                            rows_to_delete.append(row)
                            found = True
                            self.matches += 1
                            self.results[gene]['found'].add(variant)
                            break
                if not found:
                    self.results[gene]['not_found']['in_fq'].add(variant)
                    self.unmatched_predictions += 1
        # Delete any rows which have been matched against
        # This is done in reverse, high indexes first
        # From low to high would mean the list shrinks and high indexes are invalid
        for row in sorted(rows_to_delete, reverse=True):
            del gene_vcf[row]
        self.vcf[gene] = gene_vcf


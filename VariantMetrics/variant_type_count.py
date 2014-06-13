#!/usr/bin/env python

#Written to get number of SNP, Indels and Mixed (SNP&Indel) sites

import argparse
import gzip
import sys

def main(args):
    fin = gzip.open(args.vcf,'rt') if args.vcf.endswith('.gz') else open(args.vcf)
    fout = open(args.out,'wt') if args.out != "" else sys.stdout

    #Probably add a function to check VCF header, etc.
    header = None
    first_sample_index = 9;

    var_site_types = ['SNP_BI','SNP_MULTI','INDEL_BI','INDEL_MULTI','MIXED']
    site_counts = dict(zip(var_site_types,[0,0,0,0,0]))
    total_var_sites = 0

    for line in fin:
        line = line.strip()

        if not line.startswith('#'):

            vcf_fields = line.split('\t')

            if vcf_fields[header['FILTER']] != 'PASS':
                continue
                
            alt_alleles = vcf_fields[header['ALT']].split(',')
            total_var_sites += 1

            if len(alt_alleles) == 1:
                if len(alt_alleles[0]) == 1 and len(vcf_fields[header['REF']]) == 1:
                    site_counts['SNP_BI'] += 1
                else:
                    site_counts['INDEL_BI'] += 1
            else:
                snp_flag = 0
                indel_flag = 0

                for a in alt_alleles:
                    if len(a) == len(vcf_fields[header['REF']]):
                        snp_flag = 1
                    else:
                        indel_flag = 1

                if snp_flag == 1 and indel_flag == 1:
                    site_counts['MIXED'] += 1
                elif snp_flag == 1:
                    site_counts['SNP_MULTI'] += 1
                elif indel_flag == 1:
                    site_counts['INDEL_MULTI'] += 1



        elif line.startswith('#CHROM'):
            line = line.lstrip('#')

            header = line.split('\t')
            header = dict(zip(header, range(len(header))))

    print("VarType\tCount\tPercent")
    
    for k in var_site_types:
        print("{0}\t{1}\t{2:.3f}".format(k,site_counts[k],100*site_counts[k]/total_var_sites))

    fin.close()

    if fout != None:
        fout.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', '--input', '-i', help='Input VCF file (.vcf) or gzipped VCF file (.vcf.gz)', required=True)
    parser.add_argument('--out', '-o', help='Output tsv file', default="")

    args = parser.parse_args()
    main(args)

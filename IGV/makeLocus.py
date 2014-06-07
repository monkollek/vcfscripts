#!/usr/bin/env python

#Written to produce input locus file for IGV_plotter

import argparse
import gzip
import sys

def main(args):
    fin = gzip.open(args.vcf) if args.vcf.endswith('.gz') else open(args.vcf)
    fout = open(args.out,'wt') if args.out != "" else sys.stdout

    #Probably add a function to check VCF header, etc.
    header = None
    first_sample_index = 9;

    for line in fin:
        line = line.strip()

        if not line.startswith('#'):

            vcf_fields = line.split('\t')

            gt_field_names = vcf_fields[header['FORMAT']].split(':')

            chr_pos = vcf_fields[header['CHROM']]+"_"+vcf_fields[header['POS']]        
            alleles = vcf_fields[header['ALT']].split(',')
            alleles.insert(0,vcf_fields[header['REF']])

            locus_line = ""    

            for sample_index in range(first_sample_index,len(vcf_fields)):
                #For some reason sample != './.' does not work here
                if not vcf_fields[sample_index].startswith('./.'):                
                    sample = dict(zip(gt_field_names,vcf_fields[sample_index].split(':')))

                    sample_alleles = sample['GT'].split('/')

                    if int(sample['GQ'])  >= args.GQ and int(sample['DP'])  >= args.DP and sample['GT'] != '0/0':
                        
                        sample_str = "{0}_{1}_{2}".format(sample_names[sample_index],alleles[int(sample_alleles[0])],alleles[int(sample_alleles[1])])
                        file_name = "{0}_{1}".format(chr_pos,sample_str)

                        print("{0}\t{1}+{2}\t{3}".format(file_name,chr_pos,args.locus_size,sample_names[sample_index]),file=fout)

            #info_field = dict([(x.split('=', 1)) for x in re.split(';(?=\w)', fields[header['INFO']]) if x.find('=') > -1])

        elif line.startswith('#CHROM'):
            line = line.lstrip('#')

            header = line.split('\t')
            sample_names = header
            header = dict(zip(header, range(len(header))))


    fin.close()

    if fout != None:
        fout.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', '--input', '-i', help='Input VCF file (.vcf) or gzipped VCF file (.vcf.gz)', required=True)
    parser.add_argument('--out', '-o', help='Output locus file for IGV', default="")
    parser.add_argument('--GQ', default=30, type=int, help='Minimum Sample Genotype Quality')
    parser.add_argument('--DP', default=10, type=int, help='Minimum Sample Depth')
    parser.add_argument('--locus_size','-lz',default=50, type=int, help='Locus Window Size')

    args = parser.parse_args()
    main(args)

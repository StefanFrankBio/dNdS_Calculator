import argparse
from Bio import SeqIO
import json
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--consensus" )
    parser.add_argument("-v", "--variants")
    parser.add_argument("-o", "--output")
    return parser.parse_args()

args = parse_args()
reference = list(SeqIO.parse(args.consensus, "fasta"))
variants = list(SeqIO.parse(args.variants, "fasta"))

with open('site_ratios.json', 'r') as infile:
    codon_data = json.load(infile)

with open('substitutions.json', 'r') as infile:
    sub_dict = json.load(infile)

ref_codons = re.findall(r'...', str(reference[0].seq))
N = sum(codon_data[codon][0] for codon in ref_codons)
S = sum(codon_data[codon][1] for codon in ref_codons)

with open(args.output, 'w') as outfile:
    for variant in variants:
        var_codons = re.findall(r'...', str(variant.seq))
        NS = 0
        SS = 0
        for ref, var in zip(ref_codons, var_codons):
            subs = sub_dict[ref][var]
            NS += subs[0]
            SS += subs[1]
        try:
            print(variant.id, NS, SS, (NS/N)/(SS/S), file=outfile, sep='\t')
        except ZeroDivisionError:
            print(variant.id, NS, SS, file=outfile, sep='\t')

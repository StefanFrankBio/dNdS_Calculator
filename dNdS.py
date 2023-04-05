#!/usr/bin/env python3
import argparse
import json
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
from collections import Counter

def site_counts(codon_table, outfile):
    with open(codon_table, 'r') as infile:
        codon_dict = json.load(infile)
    keys = codon_dict.keys()
    alphabet = {char for string in keys for char in string}
    site_counts_dict = {}
    for codon in keys:
        subs = single_substitutions(codon, alphabet)
        synonymous_sites = sum([codon_dict[sub] == codon_dict[codon] for sub in subs])
        nonsynonymous_sites = len(subs) - synonymous_sites
        site_counts_dict[codon] = (nonsynonymous_sites, synonymous_sites)
    with open(outfile, 'w') as outfile:
        json.dump(site_counts_dict, outfile, indent=4)

def single_substitutions(codon, alphabet):
    substitutions = []
    for i, base in enumerate(codon):
        for new_base in alphabet:
            if base != new_base:
                new_codon = codon[:i] + new_base + codon[i+1:]
                substitutions.append(new_codon)
    return substitutions

def sub_counts(codon_table, outfile):
    with open(codon_table, 'r') as infile:
        codon_dict = json.load(infile)
    keys = codon_dict.keys()
    codon_pairs = list(itertools.product(keys, repeat=2))
    sub_counts_dict = {}
    for r, v in codon_pairs:
        subs = [r[:i] + v[i] + r[i+1:] for i in range(3) if r[i] != v[i]]
        NS = sum([codon_dict[sub] != codon_dict[r] for sub in subs])
        SS = len(subs) - NS
        sub_counts_dict.setdefault(r, {})[v] = (NS, SS)
        sub_counts_dict[r]['---'] = (0,0)
    with open(outfile, 'w') as outfile:
        json.dump(sub_counts_dict, outfile, indent=4)

def transpose(infile, prefix):
    msa = list(SeqIO.parse(infile, "fasta"))
    codons_by_seq = []
    for record in msa:
        codons_by_seq.append(re.findall(r'...', str(record.seq)))
    codons_by_site = list(map(list, zip(*codons_by_seq)))
    transpose_dict = {}
    for i, site in enumerate(codons_by_site):
        count = Counter(site)
        del count['---']
        sorted_counts = sorted(count.items(), key=lambda x: (-x[1], x[0]))
        transpose_dict[i] = sorted_counts
    with open(f'{prefix}.json', 'w') as outfile:
        json.dump(transpose_dict, outfile, indent=4)

def consensus(infile, prefix):
    with open(infile, 'r') as infile:
        transpose_dict = json.load(infile)
    consensus = Seq(''.join([val[0][0] for val in transpose_dict.values()]))
    record = SeqRecord(consensus, id='prefix', description='')
    SeqIO.write(record, f'{prefix}.cons', 'fasta')

def per_sequence(infile, prefix, site_counts, sub_counts):
    print(infile)
    transpose(infile, prefix)
    consensus(f'{prefix}.json', prefix)

    reference = next(SeqIO.parse(f'{prefix}.cons', "fasta"))
    variants = list(SeqIO.parse(infile, "fasta"))

    with open(site_counts, 'r') as infile:
        site_counts_dict = json.load(infile)

    with open(sub_counts, 'r') as infile:
        sub_counts_dict = json.load(infile)

    ref_codons = re.findall(r'...', str(reference.seq))
    N = sum(site_counts_dict[codon][0] for codon in ref_codons)
    S = sum(site_counts_dict[codon][1] for codon in ref_codons)

    with open(f'{prefix}.dNdS', 'w') as outfile:
        print('site', 'N', 'S', 'NS', 'SS', 'dNdS', file=outfile, sep='\t')
        for variant in variants:
            var_codons = re.findall(r'...', str(variant.seq))
            NS = 0
            SS = 0
            for ref, var in zip(ref_codons, var_codons):
                subs = sub_counts_dict[ref][var]
                NS += subs[0]
                SS += subs[1]
            try:
                print(variant.id, N, S, NS, SS, (NS/N)/(SS/S), file=outfile, sep='\t')
            except ZeroDivisionError:
                print(variant.id, N, S, NS, SS, file=outfile, sep='\t')

def per_site(infile, prefix, site_counts, sub_counts):
    transpose(infile, prefix)

    with open(f'{prefix}.json', 'r') as infile:
        transpose_dict = json.load(infile)

    with open(site_counts, 'r') as infile:
        site_counts_dict = json.load(infile)

    with open(sub_counts, 'r') as infile:
        sub_counts_dict = json.load(infile)
    with open(f'{prefix}.dNdS', 'w') as outfile:
        print('site', 'N', 'S', 'NS', 'SS', 'dNdS', file=outfile, sep='\t')
        for key, value in transpose_dict.items():
            ref_codon = value[0][0]
            N = site_counts_dict[ref_codon][0]
            S = site_counts_dict[ref_codon][1]
            NS = 0
            SS = 0
            for v in value:
                NS += sub_counts_dict[ref_codon][v[0]][0] * v[1]
                SS += sub_counts_dict[ref_codon][v[0]][1] * v[1]
            if SS != 0:
                print(key, N, S, NS, SS, (NS/N)/(SS/S), file=outfile, sep='\t')
            elif NS != 0:
                print(key, N, S, NS, SS, file=outfile, sep='\t')
            else:
                pass

def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')
    
    site_counts_parser = subparsers.add_parser('site_counts')
    site_counts_parser.add_argument('-c', '--codon_table', default='resources/default_codon_table.json')
    site_counts_parser.add_argument('-o', '--outfile')
    
    sub_counts_parser = subparsers.add_parser('sub_counts')
    sub_counts_parser.add_argument('-c', '--codon_table', default='resources/default_codon_table.json')
    sub_counts_parser.add_argument('-o', '--outfile')

    transpose_parser = subparsers.add_parser('transpose')
    transpose_parser.add_argument('-i', '--infile')
    transpose_parser.add_argument('-p', '--prefix')

    consensus_parser = subparsers.add_parser('consensus')
    consensus_parser.add_argument('-i', '--infile')
    consensus_parser.add_argument('-p', '--prefix')

    per_sequence_parser = subparsers.add_parser('per_sequence')
    per_sequence_parser.add_argument('-i', '--infile')
    per_sequence_parser.add_argument('-p', '--prefix')
    per_sequence_parser.add_argument('-s', '--site-counts', default='resources/default_site_counts.json')
    per_sequence_parser.add_argument('-u', '--sub-counts', default='resources/default_sub_counts.json')

    per_site_parser = subparsers.add_parser('per_site')
    per_site_parser.add_argument('-i', '--infile')
    per_site_parser.add_argument('-p', '--prefix')
    per_site_parser.add_argument('-s', '--site-counts', default='resources/default_site_counts.json')
    per_site_parser.add_argument('-u', '--sub-counts', default='resources/default_sub_counts.json')

    return parser.parse_args()

def main():
    args = parse_args()
    if args.command == 'site_counts':
        site_counts(args.codon_table, args.outfile)
    elif args.command == 'sub_counts':
        sub_counts(args.codon_table, args.outfile)
    elif args.command == 'transpose':
        transpose(args.infile, args.prefix)
    elif args.command == 'consensus':
        consensus(args.infile, args.prefix)
    elif args.command == 'per_sequence':
        per_sequence(args.infile, args.prefix, args.site_counts, args.sub_counts)
    elif args.command == 'per_site':
        per_site(args.infile, args.prefix, args.site_counts, args.sub_counts)

if __name__ == '__main__':
    main()
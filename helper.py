import itertools
import json

nucleotides = 'ACGT'
# Length: 4 ** 3 = 64
codons = [''.join(i) for i in itertools.product(nucleotides, repeat=3)]
stop_codons = ['TAA', 'TAG', 'TGA']
# Length: 4 ** 3 - 3 = 61
codons = [codon for codon in codons if codon not in stop_codons]
# Length: 61 ** 2
codon_pairs = list(itertools.product(codons, repeat=2))

codon_to_aa = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

sub_dict = {}
for ref_codon, var_codon in codon_pairs:
    substitutions = [ref_codon[:i] + var_codon[i] + ref_codon[i+1:] for i in range(len(ref_codon)) if ref_codon[i] != var_codon[i]]
    # REMOVE STOP CODON INTERMEDIATED FROM SUBSTITUTIONS???
    ref_aa = codon_to_aa[ref_codon]
    NS_count = sum([codon_to_aa[sub] != ref_aa for sub in substitutions])
    SS_count = len(substitutions) - NS_count
    sub_dict.setdefault(ref_codon, {})[var_codon] = (NS_count, SS_count)

with open('substitutions.json', 'w') as outfile:
    json.dump(sub_dict, outfile, indent=4)

def get_substitutions(codon):
    substitutions = []
    for i, base in enumerate(codon):
        for new_base in 'ACGT':
            if base != new_base:
                new_codon = codon[:i] + new_base + codon[i+1:]
                substitutions.append(new_codon)
    return substitutions

sites_dict = {}
for codon in codons:
    single_subs = get_substitutions(codon)
    N_sites = sum([codon_to_aa[sub] != codon_to_aa[codon] for sub in single_subs])
    S_sites = len(single_subs) - N_sites
    sites_dict[codon] = (N_sites, S_sites)

with open('site_ratios.json', 'w') as outfile:
    json.dump(sites_dict, outfile, indent=4)



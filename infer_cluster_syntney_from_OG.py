#!/usr/bin/env python
# -*- encoding:utf-8 -*-
# Purpose: ：Input: Orthogroups, gene ID mamping, species ID mapping
#            and species tree
#            Output: Orthogroups for each species
# Author: 
# Note:
# Last updated on: 2019-11-03

import sys
import os
import re
from difflib import SequenceMatcher


if len(sys.argv) != 5:
    print("Usage: %s og_group.txt species_map gene_map species_tree" % sys.argv[0])
    exit(1)


og_group = sys.argv[1]
map_species = sys.argv[2]
map_genes = sys.argv[3]
species_tree = sys.argv[4]


# load species mapping file
d_genes = {}
with open(map_genes) as r:
    for l in r:
        items = l.rstrip().split()
        d_genes[items[1]] = items[2].split('.')[0].replace('_', '')
        
targets = "OPB37942 OPB37943 OPB37944 OPB37945 OPB37946 OPB37947 OPB37948 OPB37949 OPB37950 OPB37951"
l_tags = []
for item in targets.split():
    for key, value in d_genes.items():
        if item == value:
            l_tags.append(key)

# print(l_tags)

# load_species_tree
d_species_tree = {}
with open(map_species) as r:
    for l in r:
        species_code, species_file = l.rstrip().split(': ')
        species_code = species_code[:5]
        species_file = species_file.split('/')[-1][:-12]
        species_file = species_file.replace(' ', '_')

        d_species_tree[species_code] = species_file


# Update te species tree
with open(species_tree) as r, open(species_tree + '.nwk', 'w') as w:
    for l in r:
        tree = l.rstrip()
        for code, name in d_species_tree.items():
            if tree.count(code) > 1:
                sys.stderr.write('%s contain more than once. Should be: %s \n' % (code, name))
            else:
                tree = tree.replace(code, name)
    
    w.write(tree + '\n')


l_og = []
# Extract gene clusters
with open(og_group) as r, open('gene_clusters.txt', 'w') as w:
    for l in r:
        exist = any(tag in l for tag in l_tags)
        l_genes = []
        if exist:
            items = l.rstrip().split()
            l_genes.append(items[0])
            for item in items[1:]:
                l_genes.append(d_species_tree[item[:5]] + '_' + d_genes[item])
            # 
            w.write(' '.join(l_genes) + '\n')
            l_og.append(' '.join(l_genes))

# print(l_og[0])
with open("gene_clusters.new.txt", 'w') as w:
    for gene in targets.split():
        l_spes = []
        for og in l_og:
            if gene in og:
                for spe_gene in og.split()[1:]:
                    species_name = '_'.join(spe_gene.split('_')[:-1])
                    if species_name not in l_spes:
                        l_spes.append(species_name)
                continue
        w.write(" ".join(l_spes) + '\n')

# Looping for gene in each gene cluster
# l_og = l_og.sort(key=len, reverse=True)
for spe in d_species_tree.values():
    l_og_genes = []
    l_new_og_genes = []
    for og in l_og: # og is a string: OG0000005: Acremonium_chrysogenum_ATCC_11550_KFH40504 Acremonium_chrysogenum_ATCC_11550_KFH45
        t_idx = -1
        for i, tag in enumerate(targets.split()):
            if tag in og:
                t_idx = i
                break
        if spe in og:
            #for _ in l_tags: # tags are target cluster
            cont = ''
            
            for spe_gen in og.split()[1:]:
                if spe in spe_gen:
                    cont += str(i) + '-' + spe_gen.split('_')[-1] + ':'
            l_og_genes.append(cont.rstrip(':'))
        else:
            l_og_genes.append(str(i) + '-X')

    # Sorting...
    
    l_og_genes.sort()
    l_alone_genes = []
    for i, cont in enumerate(l_og_genes):
        for ii, cont_gene in enumerate(cont.split(':')): # Marking the gene having no similarity with all the other genes
            if len(cont_gene) < 4:
                continue
            alone = True
            concat = cont_gene
            for cont2 in l_og_genes[i+1:]:
                for cont2_gene in cont2.split(':'):
                    if len(cont2_gene) < 4:
                        continue
                    if cont_gene[2:-2] == cont2_gene[2:-2] and \
                       SequenceMatcher(None, cont_gene[2:], cont2_gene[2:]).ratio() >= 0.75 and \
                       abs(int(cont_gene[-2:]) - int(cont2_gene[-2:])) < 20:
                        alone = False
                        concat += ';' + cont2_gene
        if concat.count(';') > 0:
            l_new_og_genes.append(concat)
        if alone:
            l_alone_genes.append(cont_gene)
    # print(l_alone_genes)

    for gene in l_alone_genes:
        for cont in l_og_genes:
            cont.replace(gene, '')
    
    # print(spe)
    # print(l_og_genes)
    l_new_og_genes.sort(key=len, reverse=True)
    print(spe + ': ' + str(l_new_og_genes))



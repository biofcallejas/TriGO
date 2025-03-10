#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from __future__ import division
from collections import defaultdict
import sys
import collections
import os
from pronto import Ontology
import argparse
import operator
from tabulate import tabulate

program = 'enrichmentAnalysisGO.py'

parser = argparse.ArgumentParser(prog=program, formatter_class=argparse.RawTextHelpFormatter, 
    description="\n\n\tCalculate GO enrichment analysis using hypergeometric distribution\nVersion 2.0 for python 3 or higher\n\n")
requiredNamed = parser.add_argument_group('Required arguments')
requiredNamed.add_argument("-i", "--input", dest='input', required=True, type=str, help='Input database, 2 columns; GO terms and list of geneIDs associated (no headers)')
requiredNamed.add_argument("-g", "--geneset", dest='geneset', required=True, type=str, help='List of geneIDs to calculate enrichment, one geneID per line (no headers)')
requiredNamed.add_argument("-o", "--output", dest='output', required=True, type=str, help='Prefix: Output text file, results are also printed on screeen')
args = parser.parse_args()

go_values = []

##################################################################
#STEP1: Unique GOterms are identified
genes_in_db = 0
with  open (args.input, 'r') as infile:
	for line in infile:
		line = line.rstrip('\n').rstrip('\r').strip()
		if line.startswith('New_ID'):
			continue
		col = line.split('\t')
		genes_in_db += 1
		go_lists = col[1].split(',')
		go_values =  go_values + go_lists

uniq_go_vals = set(go_values)
print ('\n' + str('*' * 60) + '\nLoading GO database: ' + str(args.input))
print ('Total GO terms: ' + (str(len(go_values))))
print ('Unique GO terms: ' + (str(len(uniq_go_vals))))
print ('Genes in DB: ' + str(genes_in_db))

##################################################################
#STEP2: Using unique values, geneIDs associated are stored to the dictionary: dict_go_terms

dict_go_terms = {}

with  open (args.input, 'r') as infile:
	for line in infile:
		line = line.rstrip('\n').rstrip('\r').strip()
		if line.startswith('New_ID'):
			continue
		col = line.split('\t')
		for val in uniq_go_vals:
			if val in col[1].split(','):
				try:
					dict_go_terms[val].append(col[0])
				except KeyError:
					dict_go_terms[val] = [col[0]]
ge_ids_list = []

print ('Loading geneIDs: ' + str(args.geneset))
with  open (args.geneset, 'r') as in2_file:
	for line in in2_file:
		line = line.rstrip('\n').rstrip('\r').strip()
		ge_ids_list.append(line.strip())
		
##################################################################
#STEP3: Function, definition
import numpy as np
import pandas as pd
from scipy.stats import rankdata

def log_factorial(n):
    return np.log(np.arange(1,n+1)).sum()
def log_binomial(n,k):
    return log_factorial(n) - (log_factorial(k) + log_factorial(n-k))

def GOEA(target_genes,GENE_SETS,df_key='GO',goterms=None,fdr_thresh=0.25,p_thresh=1e-3): 
    
    if isinstance(GENE_SETS,pd.DataFrame):
        genes = np.array(list(GENE_SETS.index))
        agt = np.array(list(GENE_SETS[df_key].values))
        idx = np.argsort(agt)
        genes = genes[idx]
        agt = agt[idx]
        bounds = np.where(agt[:-1]!=agt[1:])[0]+1
        bounds = np.append(np.append(0,bounds),agt.size)
        bounds_left=bounds[:-1]
        bounds_right=bounds[1:]
        genes_lists = [genes[bounds_left[i]:bounds_right[i]] for i in range(bounds_left.size)]
        GENE_SETS = dict(zip(np.unique(agt),genes_lists))
    all_genes = np.unique(np.concatenate(list(GENE_SETS.values())))
    all_genes = np.array(all_genes)
    
    if goterms is None:
        goterms = np.unique(list(GENE_SETS.keys()))
    else:
        goterms = goterms[np.in1d(goterms,np.unique(list(GENE_SETS.keys())))]
    
    _,ix = np.unique(target_genes,return_index=True)
    target_genes=target_genes[np.sort(ix)]
    target_genes = target_genes[np.in1d(target_genes,all_genes)]
    
    N = all_genes.size

    probs=[]
    probs_genes=[]
    counter=0
    for goterm in goterms:
        if counter%1000==0:
            print(counter)
        counter+=1
        
        gene_set = np.array(GENE_SETS[goterm])
        
        B = gene_set.size
        
        gene_set_in_target = gene_set[np.in1d(gene_set,target_genes)]
        b = gene_set_in_target.size        
        if b != 0:
            n = target_genes.size
            num_iter = min(n,B)
            rng = np.arange(b,num_iter+1)
            probs.append(sum([np.exp(log_binomial(n,i)+log_binomial(N-n,B-i) - log_binomial(N,B)) for i in rng]))
        else:
            probs.append(1.0)
        
        probs_genes.append(gene_set_in_target)
        
    probs = np.array(probs)    
    probs_genes = np.array(probs_genes)
    
    fdr_q_probs = probs.size*probs / rankdata(probs,method='ordinal')
    
    filt = np.logical_and(fdr_q_probs<fdr_thresh,probs<p_thresh)
    enriched_goterms = goterms[filt]
    p_values = probs[filt]
    fdr_q_probs = fdr_q_probs[filt]    
    probs_genes=probs_genes[filt]
    
    gns = []
    for i in probs_genes:
        gns.append(';'.join(i))
    gns = np.array(gns)
    enriched_goterms = pd.DataFrame(data=fdr_q_probs,index=enriched_goterms,columns=['fdr_q_value'])
    enriched_goterms['p_value'] = p_values
    enriched_goterms['genes'] = gns
    
    enriched_goterms = enriched_goterms.sort_values('p_value')   
    return enriched_goterms

##################################################################
#STEP4: Functional enrichment, analysis
print ('\nCalculating functional enrichment.....pvalue and fdr')
#from light_goea import GOEA
target_genes = np.array(ge_ids_list)
GENE_SETS = dict_go_terms

#optional ,fdr_thresh=4.0,p_thresh=0.9
enriched_goterms = GOEA(target_genes,GENE_SETS)
enriched_goterms.to_csv(args.output + '_ENRICHMENT_RESULTS_temp.csv', encoding='utf-8')

go_result_1 = {}

##################################################################
#STEP5: Functional enrichment calculation
print ('\nCalculating functional enrichment..... Count..')
with  open (args.output + '_ENRICHMENT_RESULTS_temp.csv', 'r') as in_4_f:
	for line in in_4_f:
		line = line.rstrip('\n').rstrip('\r').strip()
		if line.startswith(','):
			continue
		col = line.split(',')
		go_result_1[col[0]] = len(col[3].split(';'))

enrich_res_go = {}

##genes in test sample
len_sample = len(ge_ids_list)
#total genes en la base de datos = genes_in_db (es un int)

with open('enrichment_values.txt', 'w') as out_temp_enr:
	out_temp_enr.write('GO\tCount\tExpected\tEnrichment\n')
	for k, v in go_result_1.items():
		if k in dict_go_terms.keys():
			expected_val = len(dict_go_terms[k])
			#enrich_res_go[k] = v/expected_val
			enrich_res_go[k] = v / (len_sample*(expected_val/genes_in_db))
			out_temp_enr.write(k + '\t' + str(v) + '\t' + str(expected_val) + '\t' + str(v / (len_sample*(expected_val/genes_in_db))) + '\n')
		else:
			print ('Fatal Error!!... GO term: ' + str(k) + ' is not recognized!, line 249')

with  open (args.output + '_ENRICHMENT_RESULTS_temp.csv', 'r') as in5_f, open(args.output + '_ENRICHMENT_RESULTS.tsv', 'w') as out_final:
	out_final.write('>GO\tEnrichment\tfdr_q_value\tp_value\tgenes\n')
	for line in in5_f:
		line = line.rstrip('\n').rstrip('\r').strip()
		if line.startswith(','):
			continue
		col = line.split(',')
		if col[0] in enrich_res_go.keys():
			enr_final_val = enrich_res_go[col[0]]
			out_final.write(col[0] + '\t' + str(enr_final_val) + '\t' + str(col[1]) + '\t' + str(col[2]) + '\t' + str(col[3]) + '\n')
		else:
			print ('Fatal Error!!... GO term: ' + str(col[0]) + ' is not recognized!, line 262')
os.remove(args.output + '_ENRICHMENT_RESULTS_temp.csv')

##################################################################
#STEP6: Ontology description
print ('\nLoading GO descriptions.....')
print(' ')
fil_go = './go.obo'
go = Ontology(fil_go)
fin_dct = {}

def extract_go_level(go_obo_file):
	# Dictionary to store GO terms with their corresponding ontology level
	go_levels = {}
	with open(go_obo_file, "r") as obo_file:
		current_id = None
		current_namespace = None

		for line in obo_file:
			line = line.strip()

			# Detect a new [Term] entry
			if line == "[Term]":
				current_id = None
				current_namespace = None
			elif line.startswith("id: GO:"):
				current_id = line.split(": ")[1]
			elif line.startswith("namespace:"):
				# Get the ontology level (BP, MF, CC)
				namespace = line.split(": ")[1]
				if namespace == "biological_process":
					current_namespace = "BP"
				elif namespace == "molecular_function":
					current_namespace = "MF"
				elif namespace == "cellular_component":
					current_namespace = "CC"

			# If both the ID and namespace have been found, store them
			if current_id and current_namespace:
				go_levels[current_id] = current_namespace
				current_id = None
				current_namespace = None
	return go_levels

go_lvls_dict = extract_go_level("go.obo")

print ('\nIdentifying GO descriptions.....')
with open(args.output + '_ENRICHMENT_RESULTS.tsv', 'r') as inf_6, open(args.output + '_ENRICHMENT_RESULTS_final.tsv', 'w') as out_final_2:
	out_final_2.write('GO\tGO_level\tDescription\tEnrichment\tfdr_q_value\tp_value\tgenes\n')
	for line in inf_6:
		line = line.rstrip('\n').rstrip('\r').strip()
		if line.startswith('>GO'):
			continue
		col = line.split('\t')
		try:
			fin_dct[col[0]] =[str(go[str(col[0])].name), str(go_lvls_dict[col[0]]) , str(col[1]), str(col[2]), str(col[3])]
			out_final_2.write(col[0] + '\t' + str(go_lvls_dict[col[0]]) + '\t' + str(go[str(col[0])].name) + '\t' + str(col[1]) + '\t' + str(col[2]) + '\t' + str(col[3]) + '\t' + str(col[4]) + '\n')
		except KeyError:
			fin_dct[col[0]] = ['NA', 'NA' ,str(enr_final_val), str(col[1]), str(col[2]), str(col[3])]
			print ('Warning!!... GO term does not have description: ' + str(col[0]) + ' is not recognized!, line 262')
			out_final_2.write(col[0] + '\tNA\tNA\t' + str(enr_final_val) + '\t' + str(col[1]) + '\t' + str(col[2]) + '\t' + str(col[3]) + '\n' + str(col[4]) + '\n')
os.remove(args.output + '_ENRICHMENT_RESULTS.tsv')

print ('\n\n\t\t-----------------------GO enrichment results-----------------------')
print ((tabulate([(k,) + tuple(v) for k,v in fin_dct.items()], headers=['GO','Description','GO_level','Enrichment','fdr_q_value','p_value'], tablefmt="fancy_grid")))
print ('***List of genes associated to GOids is available at the final result file')
print ('\nResul files: ' + args.output + '_ENRICHMENT_RESULTS_final.tsv \nenrichment_values.txt\n' + str('*'*60)+ '\nJob completed!...')
print ('\n\n')


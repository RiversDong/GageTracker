#!/usr/bin/env python3

import os
from gtfparse import read_gtf
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline


def seq_aln(protein, out):
    blastp_res = os.path.join(out, "blast_res")
    fwd_out = os.path.join(out, "fwd-results.tab")
    fwd_blastp = NcbiblastpCommandline(query=protein, subject=protein,out=fwd_out, outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",max_target_seqs=2)
    fwd_stdout, fwd_stderr = fwd_blastp()
    fwd_results = pd.read_csv(fwd_out, sep="\t", header=None)
    headers = ["query", "subject", "identity", "coverage",\
            "qlength", "slength", "alength",\
            "bitscore", "E-value"]
    fwd_results.columns = headers

    fwd_results = fwd_results[fwd_results["query"]!=fwd_results["subject"]]

    # the normalised score because bit socore is accociated with the sequence length 
    fwd_results['norm_bitscore'] = fwd_results.bitscore/fwd_results.qlength
    fwd_results['qcov'] = fwd_results.alength/fwd_results.qlength
    fwd_results['scov'] = fwd_results.alength/fwd_results.slength
    fwd_results['qcov'] = fwd_results['qcov'].clip(upper=1)
    fwd_results['scov'] = fwd_results['scov'].clip(upper=1)

    fwd_results.to_csv(blastp_res, sep="\t", index=False)
    return fwd_results

def origination(gtf, rbbh, age, age_list, young_list, mechanism_path):
    '''
    infer the origination mechanism
    '''
    hit_rbbh = list(rbbh["subject"]);
    df = read_gtf(gtf); exonInfo = df[df["feature"] == "exon"]
    gene_exon_count = exonInfo.groupby("gene").agg({'feature':'count'})
    age_info = pd.read_csv(age, sep="\t")
    genes = age_info["Gene"]
    mechanism_type = []
    parent_gene = []; gene2parent = {}

    for i in genes:
        query_branch=age_info[age_info["Gene"]==i]["Branch"].iloc[0]
        if i not in gene2parent:
            gene2parent[i]=[]
        '''
        only consider the genes in young branch
        '''
        if query_branch in young_list:
            if i in hit_rbbh:
                target = rbbh[rbbh["subject"]==i]["query"].iloc[0]
                target_branch=age_info[age_info["Gene"]==target]["Branch"].iloc[0]
                query_branch_index = age_list.index(query_branch)
                target_branch_index = age_list.index(target_branch)

                ## the age of query gene yonger than the hit genes
                # index of young gene in the behind of the old genes
                if query_branch_index > target_branch_index:
                    parent = target
                    query_exon_count = gene_exon_count.loc[i]["feature"]
                    target_exon_count = gene_exon_count.loc[target]["feature"]
                    identity = rbbh[rbbh["subject"]==i]["identity"].iloc[0]
                    coverage = rbbh[rbbh["subject"]==i]["coverage"].iloc[0]
                    evalue = rbbh[rbbh["subject"]==i]['E-value'].iloc[0]
                    gene2parent[i].append(parent)
                    # old gene has multi-exons, and young gene has a single exon-->retro
                    if evalue<=1e-10:
                        if query_exon_count == 1 and target_exon_count > 1:
                            if identity>=50 and coverage>=70: 
                                mechanism_type.append("Retro")
                            else:
                                mechanism_type.append("Retro-like")
                        else:
                            if identity>=50 and coverage>=70:
                                mechanism_type.append("Duplication")
                            else:
                                mechanism_type.append("Duplication-like")
                    elif evalue>0.05:
                        mechanism_type.append("NA")
                        gene2parent[i].append("NA")
                    else:
                        mechanism_type.append("NA")
                        gene2parent[i].append("NA")
                else:
                    mechanism_type.append("NA")
                    gene2parent[i].append("NA")
            else:
                mechanism_type.append("De novo")
                gene2parent[i].append("NA")
        else:
            '''
            if not in young gene list, do not consider the origination of such genes
            such gene is the old gene
            '''
            mechanism_type.append("NA")
            gene2parent[i].append("NA")
    for i in genes:
        parent_gene.append(gene2parent[i][0])

    age_info["Parent ID"]=parent_gene
    age_info["mechanism_type"]=mechanism_type
    age_info.to_csv(os.path.join(mechanism_path, "mechanism"), index=False,sep="\t")
    
def main_mechanism(age_table, age_list, young_list, protein_seq, out, annotation):
    mechanism_path = os.path.join(out, "mechanism")
    blastp_res = seq_aln(protein_seq, mechanism_path)
    origination(annotation, blastp_res, age_table, age_list, young_list, mechanism_path)

if __name__ == "__main__":
    main_mechanism("/home/chuand/new_gene/virilis/dating/dvirilis.age", \
            ["B-2", "B-1", "B0","B1","B2","B3","B4","B5","B6","B7","B8"],\
            ["B1","B2","B3","B4","B5","B6","B7","B8"],\
            "/home/chuand/new_gene/virilis/dating/mechanism/protein.faa",\
            "/home/chuand/new_gene/virilis/dating",\
            "/home/chuand/new_gene/virilis/dvirilis.gtf")



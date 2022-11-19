#!/usr/bin/env python3

import os
from gtfparse import read_gtf
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline

def splitGeneByAge(ageTable, OldBr, proteinSeq, out):
    '''
    Split sequences to young and old
    '''
    mechanism=os.path.join(out, "mechanism")
    if not os.path.exists(mechanism):
        os.mkdir(mechanism)
    old_faa = os.path.join(mechanism, "old.fa")
    young_faa = os.path.join(mechanism, "young.fa")

    age_f = pd.read_csv(ageTable, sep="\t");old_gene = list(age_f[age_f["Branch"].isin(OldBr)]["Gene"])
    young_gene = list(age_f[~age_f["Branch"].isin(OldBr)]["Gene"])
    records = SeqIO.parse(proteinSeq, "fasta")
    gene2proteinseq = {}; tmpseqLen=0; tmpseq=""
    for i in records:
        iid=str(i.id).split("|"); iprotein = iid[0]; igene = iid[1]
        iseq = str(i.seq)
        if igene not in gene2proteinseq:
            gene2proteinseq[igene] = [iseq]
        else:
            gene2proteinseq[igene].append(iseq)
    gene2longseq = {}
    for i in gene2proteinseq:
        seq_list = gene2proteinseq[i]
        if len(seq_list) == 1:
            gene2longseq[i] = seq_list[0]
        else:
            tmpseq = ""; tmplen=0
            for j in seq_list:
                if len(j)>tmplen:
                    tmpseq=j; tmplen=len(j)
            gene2longseq[i] = tmpseq

    young_recs=[]; old_recs=[]
    for i in gene2longseq:
        if i in old_gene:
            rec = SeqRecord(Seq(gene2longseq[i]), id=i, description="")
            old_recs.append(rec)
        if i in young_gene:
            rec = SeqRecord(Seq(gene2longseq[i]), id=i, description="")
            young_recs.append(rec)
    SeqIO.write(old_recs, old_faa, "fasta")
    SeqIO.write(young_recs, young_faa, "fasta")
    return old_faa, young_faa, mechanism

def RBH(young_fa, old_fa, out):
    '''
    Get the RBH result from BLASTp result
    '''
    rbh_result = os.path.join(out, "rbh")
    fwd_out = os.path.join(out, "young-fwd-results.tab")
    rev_out = os.path.join(out, "young-rev-results.tab")
    fwd_blastp = NcbiblastpCommandline(query=young_fa, subject=old_fa,\
                                        out=fwd_out,\
                    outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",\
                    max_target_seqs=1)
    rev_blastp = NcbiblastpCommandline(query=old_fa, subject=young_fa,\
                                        out=rev_out,\
                    outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",\
                    max_target_seqs=1)
    #fwd_stdout, fwd_stderr = fwd_blastp()
    #rev_stdout, rev_stderr = rev_blastp()

    fwd_results = pd.read_csv(fwd_out, sep="\t", header=None)
    rev_results = pd.read_csv(rev_out, sep="\t", header=None)
    headers = ["query", "subject", "identity", "coverage",\
            "qlength", "slength", "alength",\
            "bitscore", "E-value"]
    fwd_results.columns = headers
    rev_results.columns = headers

    fwd_results['norm_bitscore'] = fwd_results.bitscore/fwd_results.qlength
    rev_results['norm_bitscore'] = rev_results.bitscore/rev_results.qlength

    fwd_results['qcov'] = fwd_results.alength/fwd_results.qlength
    rev_results['qcov'] = rev_results.alength/rev_results.qlength
    fwd_results['scov'] = fwd_results.alength/fwd_results.slength
    rev_results['scov'] = rev_results.alength/rev_results.slength

    fwd_results['qcov'] = fwd_results['qcov'].clip(upper=1)
    rev_results['qcov'] = rev_results['qcov'].clip(upper=1)
    fwd_results['scov'] = fwd_results['scov'].clip(upper=1)
    rev_results['scov'] = rev_results['scov'].clip(upper=1)

    rbbh = pd.merge(fwd_results, rev_results[['query', 'subject']],\
                left_on='subject', right_on='query',\
                how='outer')
    rbbh = rbbh.loc[rbbh.query_x == rbbh.subject_y]
    rbbh = rbbh.groupby(['query_x', 'subject_x']).max()
    rbbh.to_csv(rbh_result, sep="\t", index=False)
    return rbbh

def origination(gtf, rbbh, age, old_br, mechanism_path):
    '''
    infer the origination mechanism
    '''
    young_rbbh = list(rbbh["subject_y"]);
    df = read_gtf(gtf); exonInfo = df[df["feature"] == "exon"]
    if "gene" in exonInfo.columns:
        gene_exon_count = exonInfo.groupby("gene").agg({'feature':'count'})
    else:
        gene_exon_count = exonInfo.groupby("gene_id").agg({'feature':'count'})
    age_info = pd.read_csv(age, sep="\t")
    genes = age_info["Gene"]
    mechanism_type = []
    parent_gene = []
    for i in genes:
        # iloc is used to get value from series
        branch=age_info[age_info["Gene"]==i]["Branch"].iloc[0]
        if branch in old_br:
            mechanism_type.append("NA")
            parent_gene.append("NA")
        else:
            if i in young_rbbh:
                parent = rbbh[rbbh["subject_y"]==i]["query_y"].iloc[0]
                parent_gene.append(parent)
                young_exon_count = gene_exon_count.loc[i]["feature"]
                old_exon_count = gene_exon_count.loc[parent]["feature"]
                identity = rbbh[rbbh["subject_y"]==i]["identity"].iloc[0]
                coverage = rbbh[rbbh["subject_y"]==i]["coverage"].iloc[0]
                if young_exon_count == 1 and old_exon_count>1:
                    if identity>=50 and coverage>=70: 
                        mechanism_type.append("RNA-based")
                    else:
                        mechanism_type.append("RNA-based like")
                else:
                    if identity>=50 and coverage>=70:
                        mechanism_type.append("DNA-based")
                    else:
                        mechanism_type.append("DNA-based like")
            else:
                mechanism_type.append("De novo")
                parent_gene.append("NA")
    age_info["Parent ID"]=parent_gene
    age_info["mechanism_type"]=mechanism_type
    age_info.to_csv(os.path.join(mechanism_path, "mechanism"), index=False,sep="\t")
    
def main_mechanism(age_table, old_br, protein_seq, out, annotation):
    # age_table: /home/chuand/new_gene/virilis/dating/dvirilis.age
    # old_br: ["B-2", "B-1", "B0"]
    # protein_seq: /home/chuand/new_gene/virilis/dvirilis.faa
    # out: /home/chuand/new_gene/virilis/dating
    # annotation: /home/chuand/new_gene/virilis/dvirilis.gtf
    old_faa,young_faa,mechanism_path = splitGeneByAge(age_table,old_br,protein_seq, out)
    rbbh = RBH(young_faa, old_faa, mechanism_path)
    origination(annotation, rbbh, age_table, old_br, mechanism_path)

if __name__ == "__main__":
    #main_mechanism("/home/chuand/new_gene/virilis/dating/dvirilis.age", ["B-2", "B-1", "B0"], "/home/chuand/new_gene/virilis/dvirilis.faa", "/home/chuand/new_gene/virilis/dating","/home/chuand/new_gene/virilis/dvirilis.gtf")
    main_mechanism("/home/chuand/new_gene/dm_new/dm_previous",\
            [-2,-1,0],\
            "/home/chuand/new_gene/dm_new/dm97_new.fa",\
            "/home/chuand/new_gene/dm_new/dating/","/home/chuand/new_gene/dm_new/dm97.gtf")






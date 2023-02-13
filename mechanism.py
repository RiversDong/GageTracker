#!/usr/bin/env python3

'''
This script is used to infer the origination mechanism
Prepare you all protein sequences (in fasta format) coded by all genes
'''
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def splitGeneByAge(ageTable, OldBr, proteinSeq, out):
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

def alignment(old_faa, young_faa, mechanism):
    result = os.path.join(mechanism, "alignment")
    dbname = os.path.join(mechanism, os.path.basename(old_faa))
    mkdb = "makeblastdb -in {} -dbtype prot -out {} &> /dev/null ".format(old_faa, dbname)
    os.system(mkdb)
    align = "blastp -query {} -db {} -out {}  \
            -outfmt \"6 qseqid sseqid pident length qlen slen evalue bitscore\"".format(young_faa,dbname,result)
    os.system(align)
    return result

def structure(youngTOold):
    f = open(youngTOold).read().split("\n")[0:-1]
    for i in f:
        info = i.split("\t"); young1 = info[0]; old1 = info[1]
        ident1 = float(info[2]);evalue=float(info[6])
        hsp_len=int(info[3]); qlen = int(info[4]); slen=int(info[5])
        longer_seq = qlen if qlen>slen else slen; coverage = hsp_len/longer_seq


if __name__ == "__main__":
    #old_faa, young_faa, mechanism = splitGeneByAge("/home/chuand/new_gene/virilis/dating/dvirilis.age",\
    #        ["B-2", "B-1", "B0"], "/home/chuand/new_gene/virilis/dvirilis.faa",\
    #        "/home/chuand/new_gene/virilis/dating")
    #result = alignment(old_faa, young_faa, mechanism)



















    


    


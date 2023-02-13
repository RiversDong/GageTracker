#!/usr/bin/env python3

import os
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def pathBuilt(outpath):
    masking = os.path.join(outpath, "masking")
    if not os.path.exists(masking):
        os.mkdir(masking)
    last = os.path.join(outpath, "last")
    if not os.path.exists(last):
        os.mkdir(last)
    lastdb = os.path.join(outpath, "lastdb")
    if not os.path.exists(lastdb):
        os.mkdir(lastdb)
    chain = os.path.join(outpath, "chain")
    if not os.path.exists(chain):
        os.mkdir(chain)
    twoBit = os.path.join(outpath, "2bit")
    if not os.path.exists(twoBit):
        os.mkdir(twoBit)
    size = os.path.join(outpath, "size")
    if not os.path.exists(size):
        os.mkdir(size)
    cleanChain = os.path.join(outpath, "cleanChain")
    if not os.path.exists(cleanChain):
        os.mkdir(cleanChain)
    rBest = os.path.join(outpath, "rBest")
    if not os.path.exists(rBest):
        os.mkdir(rBest)
    maf = os.path.join(outpath, "axt2maf")
    if not os.path.exists(maf):
        os.mkdir(maf)
    block = os.path.join(outpath, "block")
    if not os.path.exists(block):
        os.mkdir(block)
    psl = os.path.join(outpath, "psl")
    if not os.path.exists(psl):
        os.mkdir(psl)
    newfa = os.path.join(outpath, "nfasta")
    if not os.path.exists(newfa):
        os.mkdir(newfa)
    mechanism = os.path.join(outpath, "mechanism")
    if not os.path.exists(mechanism):
        os.mkdir(mechanism)
    return masking, last, lastdb, psl, twoBit, size, chain, cleanChain, rBest, maf, block, newfa

def get_focal_size(target, size_path):
    target_base = os.path.basename(target)
    isize=os.path.join(size_path, target_base+".chrom.sizes")
    cmd = "faSize {0} -detailed -tab | sort -k1 > {1}".format(target, isize)
    os.system(cmd)
    return isize

#def addPrefix(refDic, newFaPath, mask="false"):
def addPrefix(refDic, newFaPath, process_num=1, mask="false"):
    print("Step I: Add prefix for each ID of focal species and outgroup species...")
    new_reference = {}
    for i in refDic:
        new_reference[i]=[]
        ispecies = refDic[i]
        for j in ispecies:
            new_reference[i].append(os.path.join(newFaPath, os.path.basename(j)))
    tmp_ref_fasta = list(refDic.values()); ref_fasta=[]
    for i in tmp_ref_fasta:
        ref_fasta.extend(i)
    for i in ref_fasta:
        ibase = os.path.basename(i)
        new_fa = os.path.join(newFaPath, ibase)

        # --- following is the newly added for maked the reference genome ---
        masked_ibase = os.path.join(newFaPath, ibase+".msk")
        cmd_tantan = "tantan {0} > {1}".format(i, masked_ibase)
        os.system(cmd_tantan)
        # ---end ---

        # After mask reference genome by tantan, we should read the masked genome
        # And add "ref_" prefix to the masked reference genomes
        irecords = SeqIO.parse(masked_ibase, "fasta"); new_records = []
        for j in irecords:
            jid = "ref_"+str(j.id) 
            #jseq = str(j.seq).upper()
            jseq = str(j.seq)
            jrec = SeqRecord(Seq(jseq), id = jid, description="")
            new_records.append(jrec)
        SeqIO.write(new_records, new_fa, "fasta")
    if mask == "true":
        '''
        mask the genomes by multiprocessing
        '''
        new_reference_mask={}
        cmd1_list = []; cmd2_list = []
        for i in new_reference:
            refs = new_reference[i]
            new_reference_mask[i]=[]
            for j in refs:
                window_tmp = j+".tmp"
                window_mask=j+".w"
                cmd1 = "windowmasker -mk_counts -in {0} -infmt fasta -out {1} -sformat obinary &> /dev/null".format(j, window_tmp); cmd1_list.append(cmd1)
                cmd2 = "windowmasker -ustat {0} -in {1} -out {2} -outfmt fasta".format(window_tmp,j,window_mask); cmd2_list.append(cmd2)
                new_reference_mask[i].append(window_mask)
        pool = ThreadPool(process_num); pool.map(os.system, cmd1_list)
        pool.map(os.system, cmd2_list);
        # 以下代码将.w后缀去掉
        for i in new_reference:
            ifiles = new_reference[i]
            for j in ifiles:
                wfile = j + ".w"
                os.remove(j);
                os.rename(wfile, j)
        print(new_reference)
        return new_reference
    else:
        return new_reference
    
def repeatMasking(tfa, masking):
    '''
    repeat mask the focal species
    '''
    print("Step II: Mask the repeat sequences of focal species...")
    ibase = os.path.basename(tfa)
    itmp = os.path.join(masking, ibase+".tmp")
    cmd1 = "windowmasker -mk_counts -in {0} -infmt fasta -out {1} -sformat obinary &> /dev/null".format(tfa, itmp)
    os.system(cmd1)
    
    imask = os.path.join(masking, ibase+".tmpmask")
    cmd2 = "windowmasker -ustat {0} -in {1} -out {2} -outfmt fasta".format(itmp, tfa, imask) 
    os.system(cmd2)
    
    iintervals = os.path.join(masking, ibase+".intervals")
    cmd3 = "windowmasker -ustat {0} -in {1} -out {2} -outfmt interval".format(itmp, tfa, iintervals) 
    os.system(cmd3)
    return imask, iintervals
    
def repFunc(positions, seq):
    for j in positions:
        start = j[0]-1; end = j[1]
        exon_seq = seq[start:end]
        replace_str = exon_seq.upper()
        seq = seq[:start] + replace_str + seq[end:]
    return seq
    
def mask_gene_ratio(intervals, exonbed, pre_mask, chrsize):
    '''
    pre_mask is the masked genomes
    '''
    print("Step III: Change the masked exon to unmask state...")
    "get the mask bed file that storing the mask region of each chromosome"
    f = open(intervals).read().split("\n")[0:-1]
    intervals_bed = os.path.join(os.path.split(pre_mask)[0], "intervals.bed")
    out_intervals = open(intervals_bed, "w")
    out_overlap = intervals_bed+".wo"
    chr2pos = {}; chr_list = []
    for i in f:
        if ">" in i:
            chromosome = i.split()[0].replace(">",""); chr_list.append(chromosome)
            chr2pos[chromosome] = []
        else:
            start, end = i.split(" - ")
            start = str(int(start)); end = str(int(end)+2)
            chr2pos[chromosome].append([start, end]);
    chr_list.sort()
    for i in chr_list:
        ipositions = chr2pos[i]
        for j in ipositions:
            out_intervals.write(i+"\t" + "\t".join(j)+"\n")
    out_intervals.close()
    cmd = "bedtools intersect -a {0} -b {1} -wo > {2}".format(intervals_bed, exonbed, out_overlap)
    os.system(cmd)
    f = open(out_overlap).read().split("\n")[0:-1]
    overlapped_exon = []
    for i in f:
        exon = i.split("\t")[3:6]; exon = "\t".join(exon); overlapped_exon.append(exon)
    overlapped_exon_set = set(overlapped_exon)
    chr2exonpos = {}
    for i in overlapped_exon_set:
        info = i.split("\t"); chromosome = info[0]; exonstart,exonend=int(info[1]),int(info[2])
        if chromosome not in chr2exonpos:
            chr2exonpos[chromosome] = [[exonstart, exonend]]
        else:
            chr2exonpos[chromosome].append([exonstart, exonend])
    #records = SeqIO.parse(pre_mask, "fasta")
    #rec_list = []
    #for i in records:
    #    iid = str(i.id)
    #    sequence = str(i.seq);
    #    if iid in chr2exonpos:
    #        positions = chr2exonpos[iid]
    #         new_sequence = repFunc(positions, sequence);
    #        rec = SeqRecord(Seq(new_sequence), id = iid, description=""); rec_list.append(rec)
    #    else:
    #        new_sequence = sequence
    #        rec = SeqRecord(Seq(new_sequence), id = iid, description=""); rec_list.append(rec)
    new_mask = pre_mask.replace(".tmpmask", ".mask")
    #SeqIO.write(rec_list, new_mask, "fasta")
    #return new_mask
    '''
    replace the aboved python script by the perl script written by Fang
    '''
    perl_script = os.path.join(os.path.dirname(__file__), "gene.exon.upper.pl")
    cmd ="perl " + perl_script + " {0} {1} {2} > {3}".format(exonbed, pre_mask, chrsize, new_mask)
    os.system(cmd)
    return new_mask

def genomeAlignment(focalmask, reference, dbpath, last, threadNum=5, lg="false"):
    '''
    genome alignment by last program
    '''
    focus_basename = os.path.basename(focalmask)
    dbname = os.path.join(dbpath, focus_basename + ".db")
    focus_mask = focalmask
    reference_tmp = reference.values(); reference_mask = []
    for i in reference_tmp:
        reference_mask.extend(i)
    print("Step IV: Whole genome alignments...")
    print("  IV.1 construct the db index using {0}".format(focus_basename))
    if lg=="true":
        '''
        if the genome is large genome, perform the alignment by lastdb5 and lastal5
        '''
        makedb = "lastdb5  -P18 -c -R11 -uMAM8 {0} {1}".format(dbname ,focus_mask)
        os.system(makedb)
        po = Pool(threadNum)
        print("  IV.2 whole genome alignment between reference and focus in {} processes".format(threadNum))
        for i in reference_mask:
            ibase = os.path.basename(i);iout = os.path.join(last, ibase+".maf")
            cmd = "lastal5 -P 6 -C2 -u0 -m50 -p HOXD70 {0} {1} > {2}".format(dbname,i, iout)
            po.apply_async(os.system,(cmd, ))
        po.close(); po.join()
    else:
        makedb = "lastdb -P18 -c -R01 -uMAM8 {0} {1}".format(dbname ,focus_mask)
        os.system(makedb)
        po = Pool(threadNum)
        print("  IV.2 whole genome alignment between reference and focus in {} processes".format(threadNum))
        for i in reference_mask:
            ibase = os.path.basename(i);iout = os.path.join(last, ibase+".maf")
            cmd = "lastal -P 6 -C2 -u0 -m50 -p HOXD70 {0} {1} > {2}".format(dbname,i, iout)
            po.apply_async(os.system,(cmd, ))
        po.close(); po.join()

    
if __name__ == "__main__":
    chr_size = "/home/chuand/new_gene/virilis/dating/size/dvirilis.fasta.chrom.sizes"
    new_mask = mask_gene_ratio("/home/chuand/new_gene/virilis/dating/masking/dvirilis.fasta.intervals", "/home/chuand/new_gene/virilis/dating/exon.bed", "/home/chuand/new_gene/virilis/dating/masking/dvirilis.fasta.tmpmask",chr_size)
    print(new_mask)

    #refDic={"br1": ["/home/chuand/new_gene/human/data/fasta/x_tropicalis.fasta"], "br2":["/home/chuand/new_gene/human/data/fasta/chicken.fasta","/home/chuand/new_gene/human/data/fasta/zebra_finch.fasta","/home/chuand/new_gene/human/data/fasta/lizard.fasta"]}
    #new_reference = addPrefix(refDic, "/home/chuand/new_gene/human/dating/nfasta",mask="true")
    #new_reference = addPrefix(refDic, "/home/chuand/new_gene/human/dating/nfasta",mask="false")





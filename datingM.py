#!/usr/bin/env python3

from gtfparse import read_gtf
import os
import re
import pandas as pd
from pandas import Series
import polars as pl
import subprocess

def gene_mask_ratio(exon_bed, mask_interval):
    exon_f = open(exon_bed).read().split("\n")[0:-1]
    gene2exonregion = {}; gene2chr = {}; 
    for i in exon_f:
        info = i.split(); gene=info[4]; chr = info[0]
        exon_start = int(info[1]); exon_end = int(info[2])
        gene2chr[gene] = chr
        if gene not in gene2exonregion:
            # because python range does't contain the last position
            # to get the actial position of an exon, I add 1 manually
            gene2exonregion[gene] = list(range(exon_start, exon_end+1))
        else:
            gene2exonregion[gene].extend(list(range(exon_start, exon_end+1)))
    f = open(mask_interval).read().split("\n")[0:-1]
    chr2maskregion = {}
    for i in f:
        if ">" in i:
            chromosome = i.split()[0].replace(">","")
            chr2maskregion[chromosome] = []
        else:
            start, end = i.split(" - ")
            start = int(start)+1; end = int(end)+1
            # because python range does't contain the last position
            # to get the actual position of an interval, I add 1 manually
            chr2maskregion[chromosome].extend(list(range(start, end+1)))
    for i in chr2maskregion:
        chr2maskregion[i] = set(chr2maskregion[i])

    gene2ratio = {}
    for i in gene2exonregion:
        igene_exonset = set(gene2exonregion[i]); ichr = gene2chr[i]
        imaskset = chr2maskregion[ichr]
        ratio = round(len(igene_exonset.intersection(imaskset))/len(igene_exonset), 3)
        gene2ratio[i] = ratio
    return gene2ratio

def group_by_element(lst):
    index = []; result = []
    for i, _ in enumerate(lst):
        if i < len(lst) - 1 and lst[i + 1] != lst[i]:
            index.append(i + 1)
    result.append(lst[:index[0]])
    for i, item in enumerate(index):
        if i < len(index) - 1:
            result.append(lst[index[i]:index[i + 1]])
    result.append(lst[item:])
    return result

def codeName(tid, codefor):
    tid = list(tid); tid = group_by_element(tid); exonId = []; chrId=[]
    if codefor == "exon":
        for i in tid:
            eindex = 0
            for j in i:
                exonId.append("E_"+str(eindex));eindex = eindex+1
        return exonId
    elif codefor == "chr":
        for i in tid:
            gindex = 0
            for j in i:
                chrId.append(j+"_"+str(gindex)); gindex+=1
        return chrId


def getInfo(gtf, exonBed, geneBed):
    # 使用 gtfparse 解析 GTF 文件
    df_pandas = read_gtf(gtf)

    # 提取所需的列
    required_columns = ["seqname", "start", "end", "transcript_id", "gene_id", "feature"]
    df_pandas = df_pandas[required_columns]

    # 转换为 Polars DataFrame
    df = pl.DataFrame(df_pandas)

    # 处理 exon 信息
    exonInfo = df.filter(df["feature"] == "exon")
    exonInfo = exonInfo.with_columns(
        (exonInfo["end"] - exonInfo["start"]).alias("exonLen")
    )

    # 计算每个 transcript 的外显子长度总和
    t2len = exonInfo.groupby("transcript_id").agg(pl.sum("exonLen").alias("tlen"))

    # 添加 exonId
    tids = exonInfo["transcript_id"].to_list()
    exonId = codeName(tids, "exon")  # 使用现有的 codeName 函数
    exonInfo = exonInfo.with_columns(pl.Series(exonId).alias("exonId"))

    # 合并 transcript 长度信息
    exonInfo = exonInfo.join(t2len, on="transcript_id", how="left")

    # 保存 exon 信息到输出文件
    out_columns = ["seqname", "start", "end", "transcript_id", "gene_id", "exonId", "tlen"]
    exonInfo.select(out_columns).write_csv(exonBed, separator="\t", has_header=False)

    # 处理 gene 信息
    geneInfo = df.filter(df["feature"] == "gene")
    geneInfo.select(["seqname", "start", "end", "gene_id"]).write_csv(
        geneBed, separator="\t", has_header=False
    )


    
def blockGene(axt2maf, block, alignment_file, tarchrsize, gene_bed):
    '''
    deal with block gene functions
    '''
    data = pd.read_csv(tarchrsize, sep="\t", header=None); tarchrs = list(data[0]); mafs = os.listdir(axt2maf)
    pd_gap_species = pd.DataFrame()
    f = open(gene_bed).read().split("\n"); genes = [i.split("\t")[3] for i in f if i!=""]
    pd_gap_species["gene"] = genes
    gapGeneOut = os.path.join(block, "gene.gap")

    for i in mafs:
        inmaf = os.path.join(axt2maf, i)
        ifasta_f = i.replace(".maf", "")
        iblock = os.path.join(block, ifasta_f)
        cmd = ["maf2synteny", "-o", iblock, "-b", "100", "-s", alignment_file, inmaf]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print("Error executing maf2synteny command:")
            print(result.stderr)
            continue  # 或者根据需要处理错误

        # calculate the gap gene list for reference ifasta_f
        coord = os.path.join(iblock, "100"); coord = os.path.join(coord, "blocks_coords.txt")
        f = open(coord).read().split("\n"); f = f[0:-1]
        id2des = {} # sequence id to chromosome
        for j in f:
            iinfo = j.split()
            leninfo=len(iinfo)
            if leninfo!=0 and re.match('^[0-9]',iinfo[0]) and re.match('^[0-9]',iinfo[1]):
                if "ref_" not in iinfo[2]:
                    id2des[iinfo[0]] = iinfo[2]
        chrs=[]; starts=[]; ends=[]
        for j in f:
            iinfo = j.split(); leninfo=len(iinfo); iid = iinfo[0]
            if leninfo!=0 and leninfo >= 5 and iid in id2des:
                strand = iinfo[1]; length = iinfo[4]
                ichr = id2des[iid]
                if ichr in tarchrs: 
                    if strand == "+":
                        start,end=iinfo[2],iinfo[3]
                    else:
                        start,end=iinfo[3],iinfo[2]
                    chrs.append(id2des[iid]); starts.append(int(start)), ends.append(int(end))
        blockouttmp = os.path.join(block, "blocktmp")
        blockout = os.path.join(block, ifasta_f + ".block")
        blockrange={"chrs":chrs, "starts":starts,"ends":ends}; blockrange=pd.DataFrame(blockrange)
        blockrange.to_csv(blockouttmp, sep="\t", index=False, header=0)
        sortblock = "cat {0} |cut -f1-3 |sort -k1,1 -k2,2n > {1}".format(blockouttmp, blockout)
        os.system(sortblock)
        gap_list_file = os.path.join(block, "tmp.gaplist")
        cmd_gapgene = "bedtools complement -i {}".format(blockout)
        cmd_gapgene = cmd_gapgene + " -g {}".format(tarchrsize)
        cmd_gapgene = cmd_gapgene + " | awk '{if($3-$2>20000){print}}'"
        cmd_gapgene = cmd_gapgene + " | bedtools intersect -a {} -b stdin -f 1 -wa -wb | cut -f4 | sort -u > {}".format(gene_bed, gap_list_file)
        os.system(cmd_gapgene)
        
        # form the present and absent list for all genes and all species
        geneInGap = open(gap_list_file).read().split("\n"); geneInGap.remove(""); gap_table = []
        for gene in genes:
            if gene in geneInGap:
                gap_table.append(1)
            else:
                gap_table.append(0)
        pd_gap_species[ifasta_f] = gap_table
    pd_gap_species.to_csv(gapGeneOut, sep="\t", index=False)

def exonIntersectRbest(rBest, exon_bed, block):
    files = os.listdir(rBest);filespath = [os.path.join(rBest, i) for i in files]
    for i in filespath:
        f = open(i).read().split("\n")
        ioutfile = os.path.join(os.path.dirname(exon_bed),"rBest.tmp")
        iout = open(ioutfile ,"w")
        for j in f:
            info = j.split()
            if len(info)==9:
                outinfo = [info[1], info[2], info[3]]
                iout.write("\t".join(outinfo)+"\n")
        iout.close()
        interfile = os.path.join(os.path.dirname(exon_bed),"interExonRbest.tmp")
        cmdinter = "bedtools intersect -a {} -b {} -wo > {}".format(exon_bed, ioutfile, interfile)
        os.system(cmdinter)
        wo2out = os.path.join(block, os.path.basename(i).replace(".rBest.chain.net.axt","")+".wo2")
        cmd_wo_2 = "less %s |awk 'BEGIN{OFS=\"\t\"}{s[$4] += $11;g[$4]=$5;t[$4]=$7} END {for (i in s) {print g[i],i,t[i],s[i],s[i]/t[i]}}' |sort -k1,1 -k2,2 > %s" %(interfile, wo2out)
        os.system(cmd_wo_2)
        overlappedGene = os.path.join(block, os.path.basename(i).replace(".rBest.chain.net.axt","")+".rbhGene")
        cmd_wo2_gene = "awk '{if($5>=0.3){print $1}}' %s |sort -u > %s" %(wo2out,overlappedGene)
        os.system(cmd_wo2_gene)

def presentOrabsent(x):
    value = 1 if x>=1 else 0
    return value
    
def speciesGap2BranchGap(y):
    if y.sum() == len(y):
        value = 1
    else:
        value = 0
    return value
    
def getIndex(inlist):
    return [index for (index, value) in enumerate(inlist) if value == 1]
    
#def getIndex(inlist):
#    pos = [i for i in range(len(inlist)) if inlist[i]==1]
#    return pos

def dating(block, gene_bed, reference, branch, outpath, agefile, voting, gene2ratio):
    br2ref = {}
    for i in reference:
        ilist = reference[i]
        newList = [os.path.basename(j) for j in ilist]
        br2ref[i]=newList
    f = open(gene_bed).read().split("\n")
    genes = [i.split("\t")[3] for i in f if i!=""]
    chrs = [i.split("\t")[0] for i in f if i!=""]
    start_positoin = [i.split("\t")[1] for i in f if i!=""]
    end_positoin = [i.split("\t")[2] for i in f if i!=""]

    files = os.listdir(block)
    homoFile = [os.path.join(block, i) for i in files if ".rbhGene" in i]
    gene2present = {}
    gene2present["gene"]=genes
    for i in homoFile:
        ibase = os.path.basename(i).replace(".rbhGene","")
        homoLists = open(i).read().split("\n")
        homoLists.remove("")
        table = [1 if genei in homoLists else 0 for genei in genes]
        gene2present[ibase]=table
    pd_gene2present = pd.DataFrame(gene2present)
    brHomoPd = pd.DataFrame();  brHomoPd["gene"]=genes
    brBreakPd = pd.DataFrame(); brBreakPd["gene"]=genes
    speciesBreakPd = pd.read_csv(os.path.join(block,"gene.gap"),sep="\t")
    refbranch = branch[1:]
    for i in refbranch:
        iref = br2ref[i]; ihomo = pd_gene2present.loc[:,iref]
        isum = ihomo.apply(lambda x: presentOrabsent(x.sum()), axis=1)
        brHomoPd[i]=isum
        ibreak = speciesBreakPd.loc[:,iref]
        breaksum = ibreak.apply(lambda y:speciesGap2BranchGap(y), axis=1)
        brBreakPd[i]=breaksum
    #--test
    #brBreakPd.to_csv("del", sep="\t", index=False)
    homo_table_path = os.path.join(outpath, "homo.table")
    pd_gene2present.to_csv(homo_table_path, sep="\t", index=False)
    #--
    geneindex = range(0, len(genes))
    refBrNum = len(branch)-1
    age_out = open(os.path.join(outpath, agefile), "w")
    tmp_title1 = ["In_"+br for br in branch]
    tmp_title2 = ["Confidence","Branch", "Chromosome", "Start", "End", "Gene", "GeneMaskRatio"]
    out_title = "\t".join(tmp_title2) + "\t" + "\t".join(tmp_title1) + "\t" + "MC"
    age_out.write(out_title+"\n")
    for i in geneindex:
        iinfo = list(brHomoPd.loc[i,:]); genei=iinfo[0]
        tablei = iinfo[1:] # not include the gene ID in the focal species
        tmp_label = ""
        if sum(tablei) == 0:
            age = branch[0]
        else:
            homoindex = getIndex(tablei) # first homo index in references
            condiction_checking=[0] # -- newly added for checking if satisfying >= voting --
            for j in homoindex:
                jsublist = tablei[:j+1]

                if (jsublist.count(1)+1)/(len(jsublist)+1) >= voting:
                    # Missed out on handling of a situation, I have solved it by adding a check list "condiction_checking"
                    age = branch[j+1]
                    condiction_checking.append(1)
                        
            #-- newly added for handle the 'for loop' for the condiction that hasn't gene meant the internal 'if condiction'--
            
            if sum(condiction_checking) == 0:
                age = branch[0]
                tmp_label = "NT"
            # --
                    
        istart = str(start_positoin[i]); iend = str(end_positoin[i])
        ichr = chrs[i]
        tmp_table = [str(kk) for kk in tablei]
        gaplist = list(brBreakPd.loc[i,:]); gaplist = gaplist[1:]
        confidence = "CON"

        # remove genes in gap positions
        if age != branch[-1]:
            gene_age_index = branch.index(age)
            if gene_age_index == 0:
                gap_num = sum(gaplist[0:2])
                nb_num = 2
            elif gene_age_index == len(branch)-1:
                gap_num = sum(gaplist[gene_age_index])
                nb_num = 1
            else:
                gap_num = sum(gaplist[gene_age_index-1:gene_age_index+1])
                nb_num = 2
            if gap_num == nb_num:
                confidence = "NCON"
        if genei in gene2ratio:
            mask_ratio = str(gene2ratio[genei])
        else:
            mask_ratio = "NA"
        resultInfo = [confidence, age, ichr,istart,iend, genei,mask_ratio,"1"+"\t"+"\t".join(tmp_table), tmp_label]
        age_out.write("\t".join(resultInfo)+"\n")
    age_out.close()

if __name__ == "__main__":
    #getInfo("/home/chuand/new_gene/bin/gene_age_pridiction/example/dvirilis.gtf","/home/chuand/new_gene/tmp/exon.bed","/home/chuand/new_gene/tmp/gene.bed")
    #blockGene("/home/chuand/new_gene/tmp/axt2maf", "/home/chuand/new_gene/tmp/block", "/home/chuand/new_gene/bin/gene_age_pridiction/bin/maf2synteny.param.file", "/home/chuand/new_gene/tmp/size/dvirilis.fasta.chrom.sizes", "/home/chuand/new_gene/tmp/gene.bed")
    #exonIntersectRbest("/home/chuand/new_gene/tmp/rBest", "/home/chuand/new_gene/tmp/exon.bed","/home/chuand/new_gene/tmp/block")
    ctr = open("/home/chuand/new_gene/geneAger/ctl").read()
    exec(ctr)
    dating("/home/chuand/new_gene/tmp/block", "/home/chuand/new_gene/tmp/gene.bed",reference, branch,"/home/chuand/new_gene/tmp")
    

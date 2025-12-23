#!/usr/bin/env python3

import os
from multiprocessing import Pool
import subprocess


def maf2axt(last, psl):
    mafs = os.listdir(last)
    mafs_num = len(mafs); threadNum = mafs_num if mafs_num <= 5 else 5
    po = Pool(mafs_num)
    print("Step IV: transform the maf file to axt files...")
    for inmaf in mafs:
        outpsl = os.path.join(psl, inmaf+".psl")
        inmaf = os.path.join(last, inmaf)
        cmd = "maf-convert axt {0} > {1}".format(inmaf, outpsl)
        po.apply_async(os.system, (cmd,))
    po.close();po.join()

def faToTwoBit(focal_mask, reference, bit, size):
    allfa = [focal_mask]; reffa = reference.values()
    for i in reffa:
        allfa.extend(i)
    for i in allfa:
        ibase = os.path.basename(i)
        ibase = ibase.replace(".mask","")
        out2bit = os.path.join(bit, ibase+".2bit")
        infasta = i
        cmd = "faToTwoBit {0} {1}".format(infasta, out2bit)
        os.system(cmd)
        isize=os.path.join(size, ibase.replace(".mask","")+".chrom.sizes")
        cmd = "faSize {0} -detailed -tab | sort -k1 > {1}".format(infasta, isize)
        os.system(cmd)


def run_command(cmd):
    try:
        result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result.stdout.decode())
    except subprocess.CalledProcessError as e:
        print("Error executing command:", e.cmd)
        print("Return code:", e.returncode)
        print("Output:", e.stdout.decode())
        print("Error:", e.stderr.decode())

def axtChain(psl, chain, twoBit, size, cleanChain, output, target=""):
    target = target+".2bit"
    targetTwobit = os.path.join(twoBit, target)
    tmp_psld = os.listdir(psl)
    psl_num = len(tmp_psld); threadNum = psl_num if psl_num <= 5 else 5
    po = Pool(threadNum)
    print("Step V: chaining process...")
    for inpsl in tmp_psld:
        outchain = os.path.join(chain, inpsl.replace(".maf.psl","")+".chain")
        queryTwobit = os.path.join(twoBit, inpsl.replace(".maf.psl","")+".2bit")
        inpsl_path = os.path.join(psl, inpsl)
        cmd = f"axtChain -linearGap=loose {inpsl_path} {targetTwobit} {queryTwobit} {outchain} &>/dev/null"
        print(cmd)
        po.apply_async(run_command, (cmd,))
    po.close(); po.join()
    
    print("Step VI: chain clean...")
    chains = os.listdir(chain)
    tSizes = os.path.join(size, target.replace(".2bit","")+".chrom.sizes")
    t2bit = os.path.join(twoBit, target.replace(".2bit","")+".2bit")
    for i in chains:
        inchain = os.path.join(chain,i)
        outchain = os.path.join(cleanChain, i.replace(".chain","")+".chain.clean")
        qSize = os.path.join(size, i.replace(".chain","")+".chrom.sizes")
        q2bit = os.path.join(twoBit, i.replace(".chain","")+".2bit")
        removedSuspects = os.path.join(output+".removedSuspects.bed")
        cmd = f"chainCleaner {inchain} -tSizes={tSizes} -qSizes={qSize} {t2bit} {q2bit} {outchain} {removedSuspects} -linearGap=loose"
        run_command(cmd)

# Example usage
# axtChain(psl, chain, twoBit, size, cleanChain, output, target)

#def netting(cleanChain, size, target, twoBit, rBest):
def netting(cleanChain, size, target, twoBit, rBest, outpath):
    chain_files = os.listdir(cleanChain)
    clearnChains = [os.path.join(cleanChain, i) for i in chain_files if ".chain.clean" in i]
    tSizes = os.path.join(size, target+".chrom.sizes")
    t2bit = os.path.join(twoBit, target+".2bit")
    print("Step VII: chain netting...")
    for i in clearnChains:
        ibase = os.path.basename(i)
        qSize = os.path.join(size,ibase.replace(".chain.clean","")+".chrom.sizes") 
        q2bit = os.path.join(twoBit,ibase.replace(".chain.clean","")+".2bit") 
        ipath = os.path.dirname(i)
        tBestChain = os.path.join(ipath, ibase.replace(".chain.clean","")+".tBest.chain")
        rBestNet = os.path.join(rBest, ibase.replace(".chain.clean","")+".rBest.chain.net.axt")
        cmd_tbest = "chainStitchId {0} stdout\
                   | chainSwap stdin stdout \
                   | chainSort stdin {1}".format(i, tBestChain)
        os.system(cmd_tbest)
        
        focal_ref_rBest_chain = os.path.join(outpath, "focal_ref.rBest.chain")
        focal_ref_rBest_chain_net = os.path.join(outpath, "focal_ref.rBest.chain.net")
        focal_ref_rBest_chain_net_axt = os.path.join(outpath, "focal_ref.rBest.chain.net.axt")

        cmd_fr_rbest_chain = "chainPreNet {0} {1} {2} stdout\
                | chainNet -minSpace=1 -minScore=0 stdin {1} {2} stdout /dev/null \
                | netSyntenic stdin stdout \
                | netChainSubset stdin {0} stdout \
                | chainStitchId stdin stdout \
                | chainSwap stdin stdout \
                | chainSort stdin {3}".format(tBestChain, qSize, tSizes, focal_ref_rBest_chain)
        os.system(cmd_fr_rbest_chain)
        cmd_fr_rbest_chain_net = "chainPreNet {0} ".format(focal_ref_rBest_chain) + tSizes + " " + qSize
        cmd_fr_rbest_chain_net += " stdout \
                | chainNet -minSpace=1 -minScore=0 stdin " + tSizes + " " + qSize + " stdout /dev/null \
                | netSyntenic stdin " + focal_ref_rBest_chain_net
        os.system(cmd_fr_rbest_chain_net)

        cmd_netToAxt = "netToAxt " + focal_ref_rBest_chain_net + " " + focal_ref_rBest_chain +" {0} {1} {2} >/dev/null".format(t2bit, q2bit, focal_ref_rBest_chain_net_axt)
        os.system(cmd_netToAxt)
        cmd_axtsort = "axtSort " + focal_ref_rBest_chain_net_axt + " " + rBestNet
        os.system(cmd_axtsort)

def axtToMaf(rBest, size, tfa, maf):
    fsize = os.path.join(size, tfa + ".chrom.sizes")
    rfiles = os.listdir(rBest)
    threadNum = len(rfiles) if len(rfiles) <= 5 else 5
    print("Step VIII: transform the *.net.axt to *.maf...")
    po = Pool(threadNum)
    for i in rfiles:
        ibase = os.path.basename(i)
        qsize = os.path.join(size, ibase.replace(".rBest.chain.net.axt","")+".chrom.sizes")
        iin = os.path.join(rBest, i)
        iout = os.path.join(maf, ibase.replace(".rBest.chain.net.axt","")+".maf")
        cmd = "axtToMaf {0} {1} {2} {3}".format(iin, fsize, qsize, iout)
        po.apply_async(os.system,(cmd,))
    po.close(); po.join()

if __name__ == "__main__":
    #maf2psl("/home/chuand/new_gene/tmp/last", "/home/chuand/new_gene/tmp/psl")
    faToTwoBit("/home/chuand/new_gene/tmp/masking", "/home/chuand/new_gene/tmp/2bit", "/home/chuand/new_gene/tmp/size")
    #clearnChains, tSizes, t2bit = axtChain("/home/chuand/new_gene/tmp/psl", "/home/chuand/new_gene/tmp/chain", "/home/chuand/new_gene/tmp/twoBit", "/home/chuand/new_gene/tmp/size", "/home/chuand/new_gene/tmp/cleanChain",target="dvirilis.fasta")
    #netting(clearnChains, "/home/chuand/new_gene/tmp/size", tSizes, t2bit, "/home/chuand/new_gene/tmp/twoBit", "/home/chuand/new_gene/tmp/rBest")
    #axtToMaf("/home/chuand/new_gene/tmp/rBest", "/home/chuand/new_gene/tmp/size", "/home/chuand/new_gene/tmp/size/dvirilis.fasta.chrom.sizes", "/home/chuand/new_gene/tmp/axt2maf")



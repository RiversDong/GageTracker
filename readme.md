# 1. About GageTracker
GageTracker is a python package for dating gene age by micro- and macro-syteny with high speed and accuracy. It's can:
-  Data gene age according to phylogeny
- Run each command step-by-step, which facilitates to add of new outgroup species without having to repeat the comparison of the previously aligned species, thus can save more time
- GageTracker can handle the comparison of large genomes by masking outgroup species and utilizing lastdb5 alignment
- Easily add new reference genome RBH alignment without performing additional alignment that have done previously


<img src="https://github.com/RiversDong/GageTracker/assets/45725241/b3fff6d1-4f31-4235-968b-bd25e2f88827" width="80%" height="80%">


# 2. Installation

## 2.1. Dependencies
All the dependencies (listed in the following table) should be pre-installed. The users need to add all the corresponding executable programs to the environmental path before running GageTracker for dating gene age.
| Software | Links |
| --- | --- |
| python | Python >=3.8 |
| last | https://gitlab.com/mcfrith/last |
| tantan | https://gitlab.com/mcfrith/tantan |
| bedtools | https://github.com/arq5x/bedtools2/releases |
| maf2synteny | https://github.com/fenderglass/maf2synteny |
| ~~windowmasker~~ | ~~https://github.com/goeckslab/WindowMasker~~ |
| BLAST toolkit | ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ |
| bedtools (v2.29.2) | https://github.com/arq5x/bedtools2/releases |

WindowMasker is a subprogram of BLAST+, so it does not need to be installed separately.
And the UCSC tools are also needed, please downloaded these tools from http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/.
```
axtChain, axtSort, axtToMaf, axtToPsl, chainCleaner,
chainFilter, chainNet, chainPreNet, chainSort, chainStitchId,     
chainSwap, chainToPsl, faSize, faToTwoBit, mafToPsl,
netChainSubset, netSyntenic, netToAxt, pslToChain
```
**Please note that all dependencies must be added to the environment path before using GageTracker!!**

## 2.2. Download the latest released GageTracker version from github:
```
  git clone https://github.com/RiversDong/GageTracker.git
  cd GageTracker && chmod 755 ./*
  export PATH=$PATH:/path/to/GageTracker
```
/path/to/GageTracker is the path of unziped GageTracker file that downloaded by git clone https://github.com/RiversDong/GageTracker.git

## 2.3. Install package gtfparse, pandas, and biopython by typing the following command lines
```
pip install gtfparse==1.2.1
pip install pandas==1.4.3
pip install biopython==1.79
```

# 3. Usage
## 3.1. Prepare the users defined control file
There are 9 necessary parameters in the control file (ctl), which are used to designate the input annotation of target species, output path, the main branches and the reference genome list of each branch. To illustrate how to prepare the control file, we used the gene dating task of Drosophila melanogaster (D.melanogaster) as a case. In this example, there  is a total including 7 branches (from branch 0 to branch 6), and 12 Drosophila (Figure 1 from Zhang Y E et al. Genome research, 2010, 20(11): 1526-1533). The detailed parameters are listed in table 1, which totally contains three columns, the first is the required parameters, the second is an example to show how to set this parameter, and the final one is an explanation of the corresponding parameter.

![Tree](https://user-images.githubusercontent.com/45725241/202661978-5b76599b-e118-4f93-ba72-735686bfae6e.png "Phylogeny tree used for dating gene age of D.melanogaster. D.melanogaster is our target species and others are outgroup species (reference species). The tree showing here was cited from Zhang et al ")

The following table gives the detailed explanations for each of the parameters in user’s control file, in which the first column lists the necessary parameters, the second column lists the example for what value should be given and the final one gives its corresponding description.
<html>

<head>
<meta http-equiv=Content-Type content="text/html; charset=gb2312">
<meta name=Generator content="Microsoft Word 15 (filtered)">
</head>
<body lang=ZH-CN style='text-justify-trim:punctuation'>
<div class=WordSection1 style='layout-grid:15.6pt'>
<table class=MsoTable15Plain1 border=1 cellspacing=0 cellpadding=0 width=595
 style='width:446.3pt;border-collapse:collapse;border:none'>
 <tr>
  <td width=85 valign=top style='width:63.55pt;border:solid #BFBFBF 1.0pt;
  padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><b><span lang=EN-US
  style='font-family:"Times New Roman",serif'>Parameters</span></b></p>
  </td>
  <td width=236 valign=top style='width:177.2pt;border:solid #BFBFBF 1.0pt;
  border-left:none;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><b><span lang=EN-US
  style='font-family:"Times New Roman",serif'>Example</span></b></p>
  </td>
  <td width=400 valign=top style='width:205.55pt;border:solid #BFBFBF 1.0pt;
  border-left:none;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><b><span lang=EN-US
  style='font-family:"Times New Roman",serif'>Explanations</span></b></p>
  </td>
 </tr>
 <tr>
  <td width=85 valign=top style='width:63.55pt;border:solid #BFBFBF 1.0pt;
  border-top:none;background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>branch</span></p>
  </td>
  <td width=236 valign=top style='width:177.2pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal align=left style='margin-bottom:4.65pt;text-align:left'><span
  lang=EN-US style='font-size:9.0pt;font-family:"Times New Roman",serif'>branch
  = [&quot;br6&quot;, &quot;br5&quot;, &quot;br4&quot;, &quot;br3&quot;,
  &quot;br2&quot;, &quot;br1&quot;, &quot;br0&quot;]</span></p>
  </td>
  <td width=274 valign=top style='width:205.55pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>The parameter </span><span
  style='font-size:9.0pt'>“</span><span lang=EN-US style='font-size:9.0pt;
  font-family:"Times New Roman",serif'>branch</span><span style='font-size:
  9.0pt'>”</span><span lang=EN-US style='font-size:9.0pt;font-family:"Times New Roman",serif'>
  is a python list type, which is used for designating the branch name. <span
  style='color:red'>Please label the branch name from old to young. In this example
  br6 represents the youngest branch and br0 represents the oldest branch.</span></span></p>
  </td>
 </tr>
 <tr>
  <td width=85 valign=top style='width:63.55pt;border:solid #BFBFBF 1.0pt;
  border-top:none;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>old_br</span></p>
  </td>
  <td width=236 valign=top style='width:177.2pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal align=left style='margin-bottom:4.65pt;text-align:left'><span
  lang=EN-US style='font-size:9.0pt;font-family:"Times New Roman",serif'>old_br
  = [&quot;br0&quot;]</span></p>
  </td>
  <td width=274 valign=top style='width:205.55pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><a name="OLE_LINK3"><span
  lang=EN-US style='font-size:9.0pt;font-family:"Times New Roman",serif'>The
  parameter is a</span></a><span lang=EN-US style='font-size:9.0pt;font-family:
  "Times New Roman",serif'> python list type, which is used for designating the
  old branch.</span></p>
  </td>
 </tr>
 <tr>
  <td width=85 valign=top style='width:63.55pt;border:solid #BFBFBF 1.0pt;
  border-top:none;background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>protein_seq</span></p>
  </td>
  <td width=236 valign=top style='width:177.2pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal align=left style='margin-bottom:4.65pt;text-align:left'><span
  lang=EN-US style='font-size:9.0pt;font-family:"Times New Roman",serif'>protein_seq
  = &quot;/path/to/fasta/dmel_protein.fasta&quot;</span></p>
  </td>
  <td width=274 valign=top style='width:205.55pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><a name="OLE_LINK5"></a><a
  name="OLE_LINK4"><span lang=EN-US style='font-size:9.0pt;font-family:"Times New Roman",serif'> (Target protein sequence) The
  parameter is a</span></a><span lang=EN-US style='font-size:9.0pt;font-family:
  "Times New Roman",serif'> python string type. This parameter specifies the
  path of the longest protein sequence of each gene (stored in fasta format). Please note that the sequence
  names should be the same as their corresponding genes and the sequences are
  represented by the longest protein sequences. For example, gene g encodes two
  proteins: p1 and p2, among which p2 is the longest one, then the sequence
  should be organized as:</span></p>
  <p class=MsoNoSpacing><span lang=EN-US style='font-size:9.0pt;font-family:
  "Times New Roman",serif'>&gt;g</span></p>
  <p class=MsoNoSpacing><span lang=EN-US style='font-size:9.0pt;font-family:
  "Times New Roman",serif'>The protein sequence of p2</span></p>
  </td>
 </tr>
 <tr>
  <td width=85 valign=top style='width:63.55pt;border:solid #BFBFBF 1.0pt;
  border-top:none;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>outpath</span></p>
  </td>
  <td width=236 valign=top style='width:177.2pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal align=left style='margin-bottom:4.65pt;text-align:left'><span
  lang=EN-US style='font-size:9.0pt;font-family:"Times New Roman",serif'>outpath
  = &quot;/path/to/dating</span></p>
  </td>
  <td width=274 valign=top style='width:205.55pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>(output path) The parameter is
  a python string type. This parameter specifies the output of temporary files
  and gene age files</span></p>
  </td>
 </tr>
 <tr>
  <td width=85 valign=top style='width:63.55pt;border:solid #BFBFBF 1.0pt;
  border-top:none;background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>target</span></p>
  </td>
  <td width=236 valign=top style='width:177.2pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal align=left style='margin-bottom:4.65pt;text-align:left'><span
  lang=EN-US style='font-size:9.0pt;font-family:"Times New Roman",serif'>target
  = &quot;/path/to/fasta/dmel.fasta&quot;</span></p>
  </td>
  <td width=274 valign=top style='width:205.55pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>(Target genome) The parameter is
  a python string type. This parameter specifies the genome of our target
  species that we want to data gene age.</span></p>
  </td>
 </tr>
 <tr>
  <td width=85 valign=top style='width:63.55pt;border:solid #BFBFBF 1.0pt;
  border-top:none;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><span lang=EN-US
  style='font-family:"Times New Roman",serif'>annotation</span></p>
  </td>
  <td width=236 valign=top style='width:177.2pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal align=left style='margin-bottom:4.65pt;text-align:left'><span
  lang=EN-US style='font-size:9.0pt;font-family:"Times New Roman",serif'>annotation
  = &quot;/path/to/fasta/dmel.gtf&quot;</span></p>
  </td>
  <td width=274 valign=top style='width:205.55pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>(Target annotation) The parameter is
  a python string type. This parameter specifies the gene annotation of our
  target species (in gtf format). We have tested the annotation from Ensembl and
  NCBI. All of these annotations work well.  Note: the GTF file should contain these features: gene, exon, and CDS in the third column.</span></p>
  </td>
 </tr>
 <tr>
  <td width=85 valign=top style='width:63.55pt;border:solid #BFBFBF 1.0pt;
  border-top:none;background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>reference</span></p>
  </td>
  <td width=236 valign=top style='width:177.2pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal align=left style='text-align:left'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>reference[&quot;br0&quot;]
  = [&quot;/path/to/fasta/dvir.fasta&quot;,
  &quot;/path/to/fasta/dmoj.fasta&quot;, &quot;/path/to/fasta/dgri.fasta&quot;]</span></p>
  <p class=MsoNormal align=left style='text-align:left'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>reference[&quot;br1&quot;]
  = [&quot;/path/to/fasta/dwil.fasta&quot;]</span></p>
  <p class=MsoNormal align=left style='text-align:left'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>reference[&quot;br2&quot;]
  = [&quot;/path/to/fasta/ dper.fasta&quot;,
  &quot;/path/to/fasta/dpse.fasta&quot;]</span></p>
  <p class=MsoNormal align=left style='text-align:left'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>reference[&quot;br3&quot;]
  = [&quot;/path/to/fasta/dana.fasta&quot;]</span></p>
  <p class=MsoNormal align=left style='text-align:left'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>reference[&quot;br4&quot;]
  = [&quot;/path/to/fasta/dere.fasta&quot;, &quot;/path/to/fasta/dyak.fasta&quot;]</span></p>
  <p class=MsoNormal align=left style='margin-bottom:4.65pt;text-align:left'><span
  lang=EN-US style='font-size:9.0pt;font-family:"Times New Roman",serif'>reference[&quot;br5&quot;]
  = [&quot;/path/to/fasta/dsec.fasta&quot;,
  &quot;/path/to/fasta/dsim.fasta&quot;]</span></p>
  </td>
  <td width=274 valign=top style='width:205.55pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal align=left style='text-align:left'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>(Out-group species genomes) The parameter is
  a python dictionary type. This parameter designates the branch and its
  corresponding genome path. For example, branch br4 contains two species (Figure
  1) dere and dyak, which are stored in /path/to/fasta/. Therefore, the pair of
  key and values between br4 and its corresponding species should be written as
  reference[&quot;br4&quot;] = [&quot;/path/to/fasta/dere.fasta&quot;,
  &quot;/path/to/fasta/dyak.fasta&quot;]. <span style='color:red'>Please note that reference parameter should only contain the genome of outgroup species excluding the genome of our target species. </span></span></p>
  </td>
 </tr>
 <tr>
  <td width=85 valign=top style='width:63.55pt;border:solid #BFBFBF 1.0pt;
  border-top:none;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>age</span></p>
  </td>
  <td width=236 valign=top style='width:177.2pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal align=left style='text-align:left'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>age = &quot;dmel.age&quot;</span></p>
  </td>
  <td width=274 valign=top style='width:205.55pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal align=left style='text-align:left'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>The parameter is
  a python string type. This parameter specifies the file for storing gene age,
  which will be stored in the path designated by output parameter. In this
  example, the final output is /path/to/outpath/dmel.age.</span></p>
  </td>
 </tr>
 <tr>
  <td width=85 valign=top style='width:63.55pt;border:solid #BFBFBF 1.0pt;
  border-top:none;background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>voting</span></p>
  </td>
  <td width=236 valign=top style='width:177.2pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal align=left style='text-align:left'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>voting = 0.5</span></p>
  </td>
  <td width=274 valign=top style='width:205.55pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal align=left style='text-align:left'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>The parameter is
  a python int type, which is used to determine whether a gene is present or
  absent in outgroup species(default is 0.5).</span></p>
  </td>
 </tr>
</table>
<p class=MsoNormal><span lang=EN-US>&nbsp;</span></p>
</div>
</body>
</html>

## 3.2. Running the program after preparing the user-defined control file
### The basic command
```
GageTracker example.ctl [options]
    Options:
        -h, -help   show all options and their default settings, and exit
        -V         -version show version information, and exit
        -mos        mask the outgroup species using windowmasker (default: NOT mask the outgroup species)
        -lg         align the large genome using lastal5 (default: lastal)
        -p          running the program in multi-processes
        -ao or -step1   only performs alignments and get the rbh alignment according to the outgroup species. 
        -rbh or -step2  calculate the rBH regions by cosindering the genome alignment.
        -da or -step3   only performs gene dating according to the ctl files.
        -m          infer the originating mechanism for young genes (based on the BLASTp alignments)
```
### Several examples for running the command line
```
# Dating gene age without masking outgroups genomes (for genome size <1 Gb)
GageTracker dm.ctl

# Dating gene age and infer originating mechanism with 12 processes
GageTracker dm.ctl -m -p 12

# Dating gene age without masking outgroups genomes (for genome size >1 Gb)
GageTracker zebrafish.ctl -lg -p 2 

# Dating gene age with masking outgroups genomes (for large genome size ~3 Gb)
GageTracker human.ctl -lg -p 2 -mos 

# Only run the genome alignment with 5 processes
GageTracker dm.ctl -ao -p 5

# Only run the genome alignment with 5 processes and mask the outgroup species
GageTracker dm.ctl -ao -mos -p 5

# Only run the genome alignment with 5 processes and mask the outgroup species, and also handling with aligning with large genomes
GageTracker dm.ctl -ao -mos -lg -p 5

# Get the RBH alignments based on the genome alignments with 5 processes
GageTracker dm.ctl -rbh -p 5
or
GageTracker dm.ctl -step2 -p 5

# Get the gene age based on the results from the previous two steps (genome alignment and RBH results) with 5 processes
GageTracker dm.ctl -da -p 5
or
GageTracker dm.ctl -step3 -p 5
```
# 4. TIPs
## 4.1. TIP1: add new whole genome alignment as reference
We have provided a toolkit (addaln), which allow users add a new reference genome without performing additional alignments that have done. Just typing the following, a new RBH alignment will be done
```
addaln -add /home/chuand/new_gene/virilis/fasta/dbusckii.fasta -br B1 -ct dvirilis.ctl
-add: specify the genome path that user wants to add
-br: specify the branch that the newly add genome belong to
-ct: tell the path of control.
```

After this procedure is done, add the following to user’s previous control file
```
reference[“B1”] = [/home/chuand/new_gene/virilis/fasta/dbusckii.fasta]
```
Then, perform the age dating procedure (the third stage) by typing the following, Noting that dvirilis.ctl contains the newly added reference species and the updated branch list.
```
GageTracker dvirilis.ctl -da -p 5
or
GageTracker dvirilis.ctl -step3 -p 5
```
## 4.2. Tip2: get the gene age of different gene type
We provide a tool, gage_diff.py (get gene age of different annotation type), to filter gene age according to user’s needed. Just type the gene type according to the prompt message from the tool. Bellowing is an example
```
gage_diff.py test.ctl

You can select the following gene type:
1       tRNA
2       rRNA
3       transcribed_pseudogene
4       pseudogene
5       snRNA
6       lncRNA
7       guide_RNA
8       snoRNA
9       protein_coding

Select the gene type listed above: lncRNA

Your output is stored in the outpath that specifies in ctl.

```
The above example will give the gene age list of lncRNA, which stored in the path specified by ctl file. In the same way, users can also choose the age of protein-coding genes from the annotation results. Noting that the results is filtered by repeat sequences ratio among the total exon length of a gene. 

# 5. Output
The output contains four key columns (Confidence, Branch, Chromosome and GeneMaskRatio) and several supplementary columns. In “Confidence”, CON means the alignment is not detected in sequencing gaps and NCON means the alignment is detected in sequencing gaps, therefore such genes marked by NCON should be considered as unreliable, which means that such genes are deemed as young genes not because it can not be found in out group species, but because of the sequencing quality.

# 6. Some error messages and solutions
* ERROR NetFilterNonNested.perl (comes with the chaincleaner source code) is not a binary in sPATH, Either install it or provide the nets as input   
Please download NetFilterNonNested.perl from https://github.com/hillerlab/GenomeAlignmentTools/blob/master/src/NetFilterNonNested.perl, then move this script to your PATH and grant it executable permissions.

* gene.bed file is empty    
Check the GTF file to ensure it includes the "gene" feature. If the file lacks the "gene" feature, GageTracker won't be able to extract position information, as it relies on the presence of the "gene" keyword in the feature columns.

#  References
1. Zhang YE, Vibranovski MD, Krinsky BH et al. Age-dependent chromosomal distribution of male-biased genes in Drosophila, Genome Res 2010;20:1526-1533.
2. Shao Y, Chen C, Shen H, et al. GenTree, an integrated resource for analyzing the evolution and function of primate-specific coding genes[J]. Genome research, 2019, 29(4): 682-696.
3. Fang CC#, Dong C# et al. GageTracker: a tool for dating gene age by micro- and macro-synteny with high speed and accuracy (unpublished article)
4. Dong C, Zhang L, Xia S, et al. New gene evolution with subcellular expression patterns detected in PacBio-sequenced genomes of Drosophila genus[J]. bioRxiv, 2022: 2022.11. 30.518489.

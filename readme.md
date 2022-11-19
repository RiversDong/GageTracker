# Description
GageTracker is a python package for dating gene age by micro- and macro- collinearity in a high speed and accuracy. It's can:
-  Data gene age according to phyostratigraphy
-  Estimate the origination mechanism of new (young) genes
- Run each command step-by-step, which facilitate to add new outgroup species without having to repeat the comparison of the previous aligned species, thus cab save more time
- GageTracker can handle the comparison of large genomes by masking outgroup species and utilizing lastdb5 alignment
- Easily add new reference genome RBH alignment without performing additional alignment that have done previously

<img src="https://user-images.githubusercontent.com/45725241/202850768-d9fdc0a2-a9e6-4c13-a5a6-a296152e76d1.png" width="80%" align="middle">
# Dependencies
All the dependencies (listed in the following table) should be pre-installed. The users need to add all the corresponding executable programs to environmental path before running GageTracker for dating gene age.
| Software | Links |
| --- | --- |
| python | Python >=3.8 |
| last | https://gitlab.com/mcfrith/last |
| tantan | https://gitlab.com/mcfrith/tantan |
| bedtools | https://github.com/arq5x/bedtools2/releases |
| maf2synteny | https://github.com/fenderglass/maf2synteny |
| windowmasker | https://github.com/goeckslab/WindowMasker |
| BLAST toolkit | ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ |
# Installation
## Download the latest released GageTracker version from github:
```
  git clone https://github.com/xxxx/GageTracker_XX
  tar xzf GageTracker_XX.tar.gz
  cd GageTracker_XX && chmod 755 ./*
  export PATH=$PATH:/path/to/GageTracker
```
## Install package gtfparse, pandas and biopython by typing the following command lines
```
pip install gtfparse
pip install pandas
pip install biopython
```

#  Usage
## Prepare the users defined control file
There are 9 necessary parameters in control file (ctl), which are used to designate the input annotation of target species, output path, the main branches and the reference genome list of each branch. To illustrate how to preparation the control file, we used the gene dating task of Drosophila melanogaster (D.melanogaster) as a case. In this example, there are total including 7 branches (from branch 0 to branch 6), and 12 Drosophila (Figure 1 from Zhang Y E et al. Genome research, 2010, 20(11): 1526-1533). The detailed parameters are listed in table 1, which totally contains three columns, the first is the required parameters, the second is an example to show how to set this parameter, and the final one is an explanation about the corresponding parameter.

![Tree](https://user-images.githubusercontent.com/45725241/202661978-5b76599b-e118-4f93-ba72-735686bfae6e.png "Tree for dating D.melanogaster. D.melanogaster is our target species and others are outgroup species (reference species). The tree showing here was cited from Zhang et al ")


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
  <td width=274 valign=top style='width:205.55pt;border:solid #BFBFBF 1.0pt;
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
  style='color:red'>Please write the branch name from young to old. In this example
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
  name="OLE_LINK4"><span lang=EN-US style='font-size:9.0pt;font-family:"Times New Roman",serif'>The
  parameter is a</span></a><span lang=EN-US style='font-size:9.0pt;font-family:
  "Times New Roman",serif'> python string type. This parameter specifies the
  path of protein sequence (stored in fasta format). Please note that the sequence
  names should be the same with their corresponding genes and the sequences are
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
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>The parameter is
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
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>The parameter is
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
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>The parameter is
  a python string type. This parameter specifies the gene annotation of our
  target species (in gtf format). We have tested the annotation from Ensembl and
  NCBI. All of these annotations work well.</span></p>
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
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>The parameter is
  a python dictionary type. This parameter designates the branch and its
  corresponding genome path. For example, branch br4 contains two species (Figure
  1) dere and dyak, which are stored in /path/to/fasta/. Therefore, the pair of
  key and values between br4 and its corresponding species should be written as
  reference[&quot;br4&quot;] = [&quot;/path/to/fasta/dere.fasta&quot;,
  &quot;/path/to/fasta/dyak.fasta&quot;]. <span style='color:red'>Please note
  that reference parameter should only contain the genome of outgroup species
  with excluding the genome of our target species </span></span></p>
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
  example, the final output is /home/chuand/new_gene/test_1.1/ dmel.age.</span></p>
  </td>
 </tr>
 <tr>
  <td width=85 valign=top style='width:63.55pt;border:solid #BFBFBF 1.0pt;
  border-top:none;background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal style='margin-bottom:4.65pt'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>votings</span></p>
  </td>
  <td width=236 valign=top style='width:177.2pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal align=left style='text-align:left'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>voting = 3</span></p>
  </td>
  <td width=274 valign=top style='width:205.55pt;border-top:none;border-left:
  none;border-bottom:solid #BFBFBF 1.0pt;border-right:solid #BFBFBF 1.0pt;
  background:#F2F2F2;padding:0cm 5.4pt 0cm 5.4pt'>
  <p class=MsoNormal align=left style='text-align:left'><span lang=EN-US
  style='font-size:9.0pt;font-family:"Times New Roman",serif'>The parameter is
  a python int type, which is used to determine whether a gene is present or
  absent in outgroup species.</span></p>
  </td>
 </tr>
</table>
<p class=MsoNormal><span lang=EN-US>&nbsp;</span></p>
</div>
</body>
</html>

##  Running the program after preparing the user-defined control file
### The basic command
```
GageTracker example.ctl [options]
    Options:
        -h, -help   show all options and their default settings, and exit
        -V         -version show version information, and exit
        -mos        mask the outgroup species using windowmasker (default: NOT mask the outgroup species)
        -lg         align the large genome using lastal5 (default: lastal)
        -p          running the program in multi-processes
        -m          infer the originating mechanism for young genes
```
### Here is some examples
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

# Get the gene age based on the results from the previous two steps (genome alignment and RBH results) with 5 processes
GageTracker dm.ctl -da -p 5 
```
#  TIPs
## TIP1: add new whole genome alignment as reference
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
```
## Tip2: get the gene age of different gene type
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

Your output is stored in /home/chuand/new_gene/test/test.age.lncRNA

```
# Output
The output contains four keys columns (Confidence, Branch, Chromosome and GeneMaskRatio) and several supplementary columns. In “Confidence”, CON means the alignment is not detected in sequencing gaps and NCON means the alignment is detected in sequencing gaps, therefore such genes marked by NCON should be considered as unreliable, which means that such genes are deemed as young genes not because it can not be found in out group species, but because of the sequencing quality.


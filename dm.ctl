# Write the user-defined branch name from young branch to the old branch

branch = ["br6", "br5", "br4", "br3", "br2", "br1", "br0"]

old_br = ["br0"]

protein_seq = "/home/chuand/new_gene/dm/dm.faa"

outpath = "/home/chuand/new_gene/fangdata/dating3"

target = "/home/chuand/new_gene/fangdata/fasta/dm6_32.fasta"

annotation = "/home/chuand/new_gene/fangdata/dm6_32_108.gtf"


reference["br0"] = ["/home/chuand/new_gene/fangdata/fasta/dgri.fasta", "/home/chuand/new_gene/fangdata/fasta/dmoj.fasta", "/home/chuand/new_gene/fangdata/fasta/dvir.fasta"]

reference["br1"] = ["/home/chuand/new_gene/fangdata/fasta/dwil.fasta"]

reference["br2"] = ["/home/chuand/new_gene/fangdata/fasta/dper.fasta","/home/chuand/new_gene/fangdata/fasta/dpse.fasta"]

reference["br3"] = ["/home/chuand/new_gene/fangdata/fasta/dana.fasta"]

reference["br4"] = ["/home/chuand/new_gene/fangdata/fasta/dere.fasta", "/home/chuand/new_gene/fangdata/fasta/dyak.fasta"]

reference["br5"] = ["/home/chuand/new_gene/fangdata/fasta/dsec.fasta", "/home/chuand/new_gene/fangdata/fasta/dsim.fasta"]

age = "GageTracker_r1_0.5.age"

voting = 0.5

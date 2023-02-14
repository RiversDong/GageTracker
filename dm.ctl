# Write the user-defined branch name from young branch to the old branch

branch = ["br6", "br5", "br4", "br3", "br2", "br1", "br0"]

old_br = ["br0"]

protein_seq = "/home/chuand/new_gene/dm/dm.faa"

outpath = "/home/user/new_gene/dm_data/dating3"

target = "/home/user/new_gene/dm_data/fasta/dm6_32.fasta"

annotation = "/home/user/new_gene/dm_data/dm6_32_108.gtf"


reference["br0"] = ["/home/user/new_gene/dm_data/fasta/dgri.fasta", "/home/user/new_gene/dm_data/fasta/dmoj.fasta", "/home/user/new_gene/dm_data/fasta/dvir.fasta"]

reference["br1"] = ["/home/user/new_gene/dm_data/fasta/dwil.fasta"]

reference["br2"] = ["/home/user/new_gene/dm_data/fasta/dper.fasta","/home/user/new_gene/dm_data/fasta/dpse.fasta"]

reference["br3"] = ["/home/user/new_gene/dm_data/fasta/dana.fasta"]

reference["br4"] = ["/home/user/new_gene/dm_data/fasta/dere.fasta", "/home/user/new_gene/dm_data/fasta/dyak.fasta"]

reference["br5"] = ["/home/user/new_gene/dm_data/fasta/dsec.fasta", "/home/user/new_gene/dm_data/fasta/dsim.fasta"]

age = "GageTracker_r1_0.5.age"

voting = 0.5

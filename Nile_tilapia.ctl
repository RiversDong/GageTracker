# Write the user-defined branch name from young branch to the old branch

branch = ["br7", "br6", "br5", "br4", "br3", "br2", "br1", "br0"]

old_br = ["br0"]

protein_seq = "/newstorage/fangcc/genome/01.Nile_tilapia/genome/Onil.pro.fa"

outpath = "/newstorage/fangcc/genome/01.Nile_tilapia/dating"

target = "/newstorage/fangcc/genome/01.Nile_tilapia/genome/Onil.fasta"

annotation = "/newstorage/fangcc/genome/01.Nile_tilapia/genome/Onil.gtf"


reference["br0"] = ["/newstorage/fangcc/genome/01.Nile_tilapia/genome/Drer.fasta", "/newstorage/fangcc/genome/01.Nile_tilapia/genome/Spil.fasta", "/newstorage/fangcc/genome/01.Nile_tilapia/genome/Smerid.fasta", "/newstorage/fangcc/genome/01.Nile_tilapia/genome/Bbel.fasta", "/newstorage/fangcc/genome/01.Nile_tilapia/genome/Lcamt.fasta", "/newstorage/fangcc/genome/01.Nile_tilapia/genome/Aspat.fasta", "/newstorage/fangcc/genome/01.Nile_tilapia/genome/Rtyp.fasta"]

reference["br1"] = ["/newstorage/fangcc/genome/01.Nile_tilapia/genome/Csem.fasta"]

reference["br2"] = ["/newstorage/fangcc/genome/01.Nile_tilapia/genome/Olatip.fasta"]

reference["br3"] = ["/newstorage/fangcc/genome/01.Nile_tilapia/genome/Esur.fasta"]

reference["br4"] = ["/newstorage/fangcc/genome/01.Nile_tilapia/genome/Cgun.fasta"]

reference["br5"] = ["/newstorage/fangcc/genome/01.Nile_tilapia/genome/Tpoll.fasta"]

reference["br6"] = ["/newstorage/fangcc/genome/01.Nile_tilapia/genome/Oaur.fasta"]

age = "GageTracker.age"

voting = 0.5

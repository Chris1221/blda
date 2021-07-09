# Pseudobulk the single cell experiment that Ravza made.
#
# I have four clusters that have distinct cell IDs. So all I should 
# need to do is take the mapping, find where the cell ID maps to (in
# terms of the major one) then write it out to a connection there with
# a new cell barcode tag. 

# %%
from pdb import set_trace
import pandas as pd
import yaml
import pysam as ps
import os 
import os.path
import bulk_lda as bl
# %%

data = pd.read_csv("../data/test2.csv")

mapping = {}

for i, row in data.iterrows():
    name = row.cells[12:]
    if name not in mapping:
        mapping[name] = []
    
    mapping[name].append(name)

with open("../data/sc_mapping.yaml", 'w') as out:
    yaml.dump(mapping, out)

# %%
#bl.pb.split_bam("../data/BuenMerged.more150.bam", prefix = "../data/sc")

#set_trace()

#for cluster in mapping.keys():
    #os.system(f"samtools index ../data/{cluster}.bam")
#    os.system(f"macs2 callpeak -n {cluster} --format BAMPE -g hs -q 0.001 --treatment ../data/sc/{cluster}.bam -B --outdir ../data/sc/peaks/{cluster}/")
# %%

blah = {f"../data/sc/peaks/{k}/{k}_peaks.narrowPeak": f"../data/sc/{k}.bam" for k in mapping.keys()}

#prefix = "blah"

bl.cm.make_count_matrix(blah, merged_bed = "../data/sc/sc_norm.bed", output = "sc_norm", dummy = False, norm = "none")
bl.cm.make_count_matrix(blah, merged_bed = "../data/sc/sc_dummy.bed", output = "sc_dummy", dummy = True, norm = "none")
# %%
#bl.ct.run_cistopic(prefix)

#bl.ct.optimise_cistopic_parameters(prefix, output = "blah.yaml")

#bl.ct.run_cistopic(mtx_prefix = "blah", alpha = 50, beta = 0.1, dir = "sc", write_out=True, cores = 2)
# Pseudobulk the single cell experiment that Ravza made.
#
# I have four clusters that have distinct cell IDs. So all I should 
# need to do is take the mapping, find where the cell ID maps to (in
# terms of the major one) then write it out to a connection there with
# a new cell barcode tag. 

# %%
import pandas as pd
import yaml
import pysam as ps
# %%

data = pd.read_csv("../data/test2.csv")

mapping = {}
for i, row in data.iterrows():
    if row.Clusters not in mapping:
        mapping[row.Clusters] = []
    
    mapping[row.Clusters].append(row.cells[12:])

with open("../data/mapping.yaml", 'w') as out:
    yaml.dump(mapping, out)
# %%
import bulk_lda as bl

bl.pb.pb_bam("../data/BuenMerged.more150.bam", mapping, output_prefix="../data")
# %% 
# Do some simple peak calling on the bam files
# Construct the calls
for cluster in mapping.keys():
    os.system(f"macs2 callpeak -n {cluster} --format BAMPE -g hs -q 0.001 --treatment ../data/{cluster}.bam -B --outdir ../data/peaks/{cluster}/")
# %%

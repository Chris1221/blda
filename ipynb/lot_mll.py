# %% 
import bulk_lda as bl
from glob import glob
from pathlib import Path

# First need to construct the dictionary of peak calls to BAM files.
# Does it matter that they all have more than 3 columns? 

cov = glob("../data/mll/*.bigWig")
data = {}
for c in cov: 
    base = Path(c).stem
    peak = f"../data/mll/peaks/{base}/{base}.bed"
    data[peak] = c



#bl.cm.make_count_matrix(data,merged_bed="../data/mll/merged_bed.bed", output= "../data/mll/mll_dummy", format = "dummy")
#bl.cm.make_count_matrix(data, merged_bed="../data/mll/merged_bed.bed", output= "../data/mll/mll", format = "bigwig")


#bl.ct.optimise_cistopic_parameters("../data/mll/mll_dummy", output = "mll_dummy.yaml")
bl.ct.optimise_cistopic_parameters("../data/mll/mll", output = "mll.yaml")


#bl.ct.run_cistopic(mtx_prefix = "../data/mll/mll_dummy", alpha = 10, beta = 0.1, dir = "mll_dummy", write_out=True, cores = 2)
#bl.ct.run_cistopic(mtx_prefix = "../data/mll/mll", alpha = 10, beta = 0.1, dir = "mll", write_out=True, cores = 1)
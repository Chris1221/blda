# %%
import pandas as pd
import yaml
import pysam as ps
import os 
import os.path
import bulk_lda as bl

# %%
blah = {
    '../data/peaks/C1/C1_peaks.narrowPeak': "../data/C1.bam",
    '../data/peaks/C2/C2_peaks.narrowPeak': "../data/C2.bam",
    '../data/peaks/C3/C3_peaks.narrowPeak': "../data/C3.bam",
    '../data/peaks/C4/C4_peaks.narrowPeak': "../data/C4.bam",
}

bl.cm.make_count_matrix(blah, output = "mtx/rpkm_new", dummy = False)
bl.cm.make_count_matrix(blah, output = "mtx/dummy_new", dummy = True)

bl.ct.optimise_cistopic_parameters("mtx/rpkm_new", output = "output/rpkm_new.yaml", n_iter = 100)
bl.ct.optimise_cistopic_parameters("mtx/dummy_new", output = "output/dummy_new.yaml", n_iter = 100)
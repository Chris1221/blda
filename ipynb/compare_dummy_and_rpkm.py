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

bl.cm.make_count_matrix(blah, output = "mtx/rpkm", dummy = False)
bl.cm.make_count_matrix(blah, output = "mtx/dummy", dummy = True)

bl.ct.optimise_cistopic_parameters("mtx/rpkm", output = "output/rpkm.yaml", n_iter = 20)
bl.ct.optimise_cistopic_parameters("mtx/dummy", output = "output/dummy.yaml", n_iter = 20)
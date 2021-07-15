# %%
import pandas as pd
from plotnine import *

# %% 
mtx = pd.read_csv("mtx/rpkm.mtx", sep = " ", skiprows = 3, header = None)
mtx.columns = ["peak", "cell", "rpkm"]
# %%
(
    ggplot(mtx, aes(x = "rpkm")) + 
    geom_histogram(bins = 30, position = "dodge") + 
    xlim(0, 500) + 
    #facet_wrap("cell") +
    ggtitle("Normalised Read Count (RPKM)") + 
    xlab("RPKM") + 
    ylab("Count") + 
    theme_bw()
)
# %%
ns = []
for peak in mtx['peak'].unique():
    n = mtx.loc[mtx['peak'] == peak].shape[0]
    if n == 9:
        print(mtx.loc[mtx['peak'] == peak])
    ns.append(n)
# %%
import collections
import numpy as np
ns = np.array(ns)
collections.Counter(ns)
# %%
# Find a sub-collection with all the 2s
twos = pd.DataFrame()
for peak in mtx['peak'].unique():
    n = mtx.loc[mtx['peak'] == peak].shape[0]
    if n == 2:
        two = mtx.loc[mtx['peak'] == peak]
        two['cell'] = [0, 1]
        twos = pd.concat([twos, two])

# %%

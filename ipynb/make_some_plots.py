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
    xlim(0, 1000) + 
    facet_wrap("cell") +
    ggtitle("Normalised Read Count (RPKM)") + 
    xlab("RPKM") + 
    ylab("Count") + 
    theme_bw()
)
# %%

# %%
import pandas as pd
from plotnine import *

# %% 
mtx = pd.read_csv("mtx/rpkm.mtx", sep = " ", skiprows = 3, header = None)
mtx.columns = ["peak", "cell", "rpkm"]
# %%
(
    ggplot(mtx, aes(x = "rpkm", fill = "cell")) + 
    geom_histogram(bins = 30, position = "dodge") + 
    xlim(0, 1000) + 
    ggtitle("Normalised Read Count (RPKM)") + 
    xlab("RPKM") + 
    ylab("Count") + 
    theme_bw()
)
# %%

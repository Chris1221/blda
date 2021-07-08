import pandas as pd
import matplotlib
import seaborn as sns
from pdb import set_trace

def create_topic_heatmap(dir: str, output: str):
    matplotlib.use('Agg')
    cells = pd.read_csv(f"{dir}/cell-topic.tsv", index_col = 0, sep = " ")
    cells = cells.div(cells.sum(axis=0), axis = 1)
    g = sns.clustermap(cells)
    g.savefig(output) 

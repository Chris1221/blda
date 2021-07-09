import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from pdb import set_trace
import os

def create_topic_heatmap(dir: str, output: str, norm: str = "row"):
    matplotlib.use('Agg')
    cells = pd.read_csv(f"{dir}/cell-topic.tsv", index_col = 0, sep = " ")
    if norm == "col":
        cells = cells.div(cells.sum(axis=0), axis = 1)
    g = sns.clustermap(cells,cmap = "Blues")
    g.savefig(output) 

def region_scatterplot(dir: str, output_dir: str):
    df = pd.read_csv(f"{dir}/region_all.tsv", sep = "\t", index_col = 0)
    other_cols = [i for i in df.columns if not i.startswith("tSNE")]

    mkdir(output_dir)

    for i in other_cols:
        fig = sns.scatterplot(data = df, x = "tSNE1", y = "tSNE2", hue = i).get_figure()
        fig.savefig(f"{output_dir}/{i}.png")
        plt.clf()

def mkdir(dir: str):
    if not os.path.exists(dir):
        os.makedirs(dir)
import rpy2
import rpy2.robjects as robjects
from bayes_opt import BayesianOptimization
import os
import yaml

from .plots import create_topic_heatmap, region_scatterplot
from .constants import MTX_SUFFIX, CELLS_SUFFIX, REGIONS_SUFFIX

def run_cistopic(mtx_prefix: str = "notblah", alpha: float = 50, beta: float = 0.1, cores = 4, write_out = False, dir = None) -> float:
    """Runs cisTopic on a count matrix. Returns log likelihood of best model.

    Args:
        mtx_prefix (str): Path to the matrix without the file extension suffix.

    Returns:
        float: Log-likelihood of best fitting model.
    """
    print(mtx_prefix)

    mtx = mtx_prefix + MTX_SUFFIX
    cells = mtx_prefix + CELLS_SUFFIX
    regions = mtx_prefix + REGIONS_SUFFIX


    cistopic_script = (
        "library(cisTopic)\n"
        "library(Matrix)\n"
        "library(dplyr)\n"
        "print(R.home())\n"
        "find.package('cisTopic')\n"
        "options(device=pdf)\n"
        f"countMatrix <- readMM('{mtx}')\n"
        f"colnames(countMatrix) <- read.csv('{cells}',header=F)$V1\n"
        f"regions = read.csv('{regions}',header=F)$V1\n"
        "if(any(!grepl(x = regions, pattern = 'chr'))){"
        "   regions = paste0('chr', regions)"
        "}\n"
        "rownames(countMatrix) <- regions\n"
        "cisTopicObject <- createcisTopicObject(countMatrix,keepCountsMatrix=FALSE)\n"

        # parmeters
        #"topic <- c(2,4,8,10,12,14,16,18,20,22,24,26,28,30,40,50)\n"
        "topic <- c(1,2,3,4,5,6,7,8)\n"
        f"alpha = {alpha}\n"
        f"beta = {beta}\n"
        "alphaByTopic = TRUE\n"
        "iterations = 3000\n"

        "cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=topic, seed=1, alpha=alpha, alphaByTopic=TRUE, beta=beta,"
        f"nCores={cores}, iterations=3000, addModels=TRUE)\n"
        "cisTopicObject <- selectModel(cisTopicObject, type='derivative',plot=FALSE)\n"

        #"numTopics <- dim(cisTopicObject@selected.model$document_expects)[1]\n"
        "return(cisTopicObject@log.lik %>% dplyr::filter(second_derivative == max(second_derivative)))"
    )

    b = robjects.r(cistopic_script)
    print("Ran first R")

    if write_out:
        
        if not os.path.exists(f"{dir}"):
            os.mkdir(dir)

        if not os.path.exists(f"{dir}/region-topic_BW"):
            os.mkdir(f"{dir}/region-topic_BW")
        
        if not os.path.exists(f"{dir}/region-topic_BW"):
             os.mkdir(f"{dir}/selectedregion-topic_Bed")

        script = (
            "source('../r/cistopic_fns.R')\n"
            "export_dr_cisTopic(cisTopicObject,\n"
            " thrP=0.99,\n"
            " genome='hg19',\n"
            " request_interval='10',\n"
            " top='3',\n"
            " var.y='name',\n"
            " order.by='Binom_Adjp_BH',\n"
            f" out_path_cell_topic='{dir}/cell-topic.tsv',\n"
            f" out_path_region_topic='{dir}/region-topic.tsv',\n"
            f" out_path_BW='{dir}/region-topic_BW',\n"
            f" out_path_Bed='{dir}/selectedregion-topic_Bed',\n"
            f" out_path_geneId = '{dir}/gene.tsv',\n"
            f" test = FALSE)"
        )
        robjects.r(script) 
        print("Ran second R")

        perp_cells = 1
        perp_region = 10
        seed = 1

        clustering = (
            "source('../r/clustering_helper.R')\n"
            f"seed = {seed}\n"
            f"cell_topic <- as.data.frame(read.table('{dir}/cell-topic.tsv'))\n"
            f"region_topic <- as.data.frame(read.table('{dir}/region-topic.tsv'))\n"
            #f"cell_tsne <- run_tsne(cell_topic,seed=seed, check_duplicates = F, perplexity = {perp_cells})\n"
            #"cell_umap <- run_umap(cell_topic,seed=seed)\n"
            f"region_tsne <- run_tsne(region_topic,seed=seed,perplexity={perp_region},check_duplicates=F)\n"
            #f"cell_density <- density_clustering(input=cell_topic,seed=seed)\n"
            #"cell_louvain <- louvain_clustering(input=cell_topic,seed=seed,k=15)\n"
            #"write_result(cbind(cell_tsne,cell_umap,cell_density,cell_louvain), '{dir}/cell_all.csv')\n"
            f"write_result(cbind(region_topic, region_tsne), '{dir}/region_all.tsv')"
        )

        robjects.r(clustering)
        print("Ran third R")

        create_topic_heatmap(dir, f"{dir}/topic_heatmap.png")
        region_scatterplot(dir, f"{dir}/region_plots")

    return (-b[1][0])


def optimise_cistopic_parameters(mtx_prefix: str, output: str, n_iter: int = 10):
    os.environ['OPENBLAS_NUM_THREADS'] = '1'

    # Hacky but set the default based on the input parameter
    defaults = [i for i in run_cistopic.__defaults__]
    defaults[0] = mtx_prefix
    run_cistopic.__defaults__ = tuple(defaults)
 
    p_dict = {
        'alpha': [10, 500],
        'beta': [0.0001, 0.9999]
    }

    optimizer = BayesianOptimization(
        f=run_cistopic,
        pbounds=p_dict,
        random_state=1,
    )

    optimizer.maximize(
        init_points=5,
        n_iter=n_iter,
    )

    data = {
        "max": optimizer.max['params']
    }

    # Dump the parameters.
    with open(output, 'w') as o:
        yaml.dump(data, o)

    # And dump the results
    run_cistopic(**optimizer.max['params'], dir = output.replace(".yaml", ""), write_out=True, cores = 2)
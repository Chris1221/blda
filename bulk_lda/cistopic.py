import rpy2
import rpy2.robjects as robjects

from .constants import MTX_SUFFIX, CELLS_SUFFIX, REGIONS_SUFFIX

def run_cistopic(mtx_prefix: str) -> float:
    """Runs cisTopic on a count matrix. Returns log likelihood of best model.

    Args:
        mtx_prefix (str): Path to the matrix without the file extension suffix.

    Returns:
        float: Log-likelihood of best fitting model.
    """

    mtx = mtx_prefix + MTX_SUFFIX
    cells = mtx_prefix + CELLS_SUFFIX
    regions = mtx_prefix + REGIONS_SUFFIX


    cistopic_script = (
        "library(cisTopic)\n"
        "library(Matrix)\n"
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
        "topic <- c(2,4,8,10,12,14,16,18,20,22,24,26,28,30,40,50)\n"
        "alpha = 50\n"
        "beta = 0.1\n"
        "alphaByTopic = TRUE\n"
        "iterations = 3000\n"

        "cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=topic, seed=1, alpha=alpha, alphaByTopic=alphaByTopic, beta=beta,"
        "nCores=length(topic), iterations=iterations, addModels=TRUE)\n"
        "cisTopicObject <- selectModel(cisTopicObject, type='derivative',plot=FALSE)\n"

        #"numTopics <- dim(cisTopicObject@selected.model$document_expects)[1]\n"
        "return(cisTopicObject@log.lik)"
    )

    b = robjects.r(cistopic_script)
    print(b)
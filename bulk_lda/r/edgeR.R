# Have to do `conda activate edgeR` for this function

library(tidyr)
library(edgeR)
library(dplyr)
library(glue)

diffexp <- function(prefix, path_to_prefix) {

    # Shitty workaround
    genes = read.table(paste0(path_to_prefix, prefix, ".regions"))
    genes$V2 = genes$V1
    write.table(genes, paste0(path_to_prefix, prefix, ".regions2"), col.names = F, row.names = F, quote = F, sep = "\t")
    mtx = read10X(
       mtx = paste0(prefix, ".mtx"),
       genes = paste0(prefix, ".regions2"),
       barcodes = paste0(prefix, ".cells"),
       path = path_to_prefix
    ) 

    # Figure out the groupings from the names
    k562 = grep("K562", mtx$samples$Barcode)
    h1 = grep("H1ESC", mtx$samples$Barcode)
    gm = grep("GM12878", mtx$samples$Barcode)

    group = c()
    for(cell in mtx$samples$Barcode){
        if(grepl("K562", cell)){
            group <- c(group, 1)
        } else if(grepl("H1ESC", cell)){
            group <- c(group,2)
        } else if(grepl("GM12878", cell)){
            group <- c(group, 3)
        }  
    }

    design = model.matrix(~group)

    mtx = estimateDisp(mtx, design)
    fit = glmQLFit(mtx, design)

    # Get the top 100 differentially accessible genes
    tt = topTags(glmQLFTest(fit), n = 100) %>%
        as.data.frame %>%
        mutate(coords = rownames(.)) %>%
        separate(coords, into = c("chr", "start", "end"), sep = ":|-") %>%
        select(chr, start, end) %>%
        write.table(glue("{path_to_prefix}{prefix}.tt.txt"), sep = "\t", col.names = F, row.names = F, quote = F)

    #write.table(tt, paste0(path_to_prefix, prefix, ".tt.txt"), sep = "\t", col.names = T, row.names = F)

    #qlf.2vs1 <- glmQLFTest(fit, coef=3)
    #topTags(qlf.2vs1, n = 100)



    #goana(qlf.2vs1, species = "Hs")

    #library(biomaRt)
    #ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    #ensembl <- useEnsembl(biomart = "genes")
    
}

get_motifs = function(
        bed_file,
        featherFilePath = "~/well/03_FEATHER/hg19-regions-9species.all_regions.mc9nr.feather",
    ){

    regionsList <- rtracklayer::import.bed(bed_file)
    regionSets <- list(GATA1_peaks=regionsList)
    
    
    data("motifAnnotations_hgnc")
    motifAnnotation <- motifAnnotations_hgnc
    data(dbRegionsLoc_hg19)
    dbRegionsLoc <- dbRegionsLoc_hg19

    regionSets_db <- lapply(regionSets, function(x) convertToTargetRegions(queryRegions=x, targetRegions=dbRegionsLoc))

    allRegionsToImport <- unique(unlist(regionSets_db)); length(allRegionsToImport)
    motifRankings <- importRankings(featherFilePath, columns=allRegionsToImport)

    motifEnrichmentTable <- cisTarget(regionSets_db, motifRankings, aucMaxRank=0.005*getNumColsInDB(motifRankings))

    resultsToShow <- addLogo(motifEnrichmentTable)
}

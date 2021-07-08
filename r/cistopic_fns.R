export_dr_cisTopic <- function(input_cisTopicObject, thrP, genome, request_interval, top, var.y, order.by, out_path_cell_topic, out_path_region_topic, out_path_BW, out_path_Bed, out_path_geneId, test){
    suppressWarnings(library(cisTopic))
    # derive and export cell-topic distribution for each cell
    df_cell_topic <- as.data.frame(t(modelMatSelection(input_cisTopicObject, 'cell', 'Probability')))
    write.table(df_cell_topic, row.names = TRUE, col.names = TRUE, sep= " ",file = out_path_cell_topic)

    # derive and export normalized topic-region distribution for each genomic region
    cisTopicObject <- getRegionsScores(input_cisTopicObject, method='NormTop', scale=TRUE)
    colnames_region_data <- colnames(cisTopicObject@region.data)
    ScoresTopic_ids <- colnames_region_data[grepl("Scores_Topic",colnames_region_data)]
    df_region_topic <- as.data.frame(cisTopicObject@region.data[,ScoresTopic_ids])
    write.table(df_region_topic,row.names = TRUE, col.names = TRUE, sep= "\t",file = out_path_region_topic)

    # export the entrez id for each region
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    library(org.Hs.eg.db)
    cisTopicObject <- annotateRegions(cisTopicObject, txdb=txdb, annoDb='org.Hs.eg.db')
    write(cisTopicObject@region.data$geneId, file = out_path_geneId)

    # export BigWig files for observing the scored regions in the genome. (information on the length of the chromosomes has to be provided.)
    #library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    #txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    if(test){
	dir.create(out_path_BW)
        dir.create(out_path_Bed)
    } else{
	getBigwigFiles(cisTopicObject, path=out_path_BW, seqlengths=seqlengths(txdb))
	#select the most contributing regions to represent each topic
	cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=thrP, plot=FALSE)
    
    # export the binarized topics to .bed files for external analysis
    getBedFiles(cisTopicObject, path=out_path_Bed)
    }
}
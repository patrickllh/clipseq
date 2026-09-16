
library(DEWSeq)
library(IHW)
library(tidyverse)
library(data.table)
# the files get big so i run in the NAIST server
# run in server like this:
# nohup Rscript script.R > output.log 2>&1 &

setwd(getwd())

# I made the input files using HTseq according to DEWseq and HTseq's documentations

countFile <- "matrix.gz"
annotationFile <- "mapping.gz"

countData <- fread(countFile, sep = "\t") # i did head(countData) to know what to put in colData factors later
annotationData <- fread(annotationFile, sep = "\t")

# create dataset and metadata
colData <- data.frame(row.names = colnames(countData)[-1],
                      type = factor(
                        c("IN","IN","IP","IP"),
                        levels = c("IN","IP"))) # remember...the first element here is going to be the reference

ddw <- DESeqDataSetFromSlidingWindows(countData  = countData,
                                      colData    = colData,
                                      annotObj   = annotationData,
                                      tidy       = TRUE,
                                      design     = ~type)

# filter out low count windows
keep <- rowSums(counts(ddw)) >= 10 # i'm thinking of making it just 5
ddw <- ddw[keep,]

# calculate preliminary size factor for later calculations
ddw <- estimateSizeFactors(ddw)
# sizeFactors(ddw)

# filter an mRNA-only list then calculate its size factor then put back to dataset (this is accd to the docs...)
ddw_mRNAs <- ddw[ rowData(ddw)[,"gene_type"] == "protein_coding", ]
ddw_mRNAs <- estimateSizeFactors(ddw_mRNAs)
sizeFactors(ddw) <- sizeFactors(ddw_mRNAs)

# From the docs: "In general, when there is asymmetry in the data, like enrichment in IP over control, 
# it is a good pratice to call differentially expressed features with DESeq2 
# and then exclude them for normalisation"
# so I try to do that according to their code in the docs
# i am not sure yet though how to call differentially expressed "features"...
# how is it different from differentially expressed windows?
# anyway for now...

# make copy first, after getting the preliminary size factor avbove
ddw_tmp <- ddw
# then calculate some necessary statistical scores
ddw_tmp <- estimateDispersions(ddw_tmp, fitType = "local", quiet = TRUE) 
ddw_tmp <- nbinomWaldTest(ddw_tmp)
ddw_tmp_results <- results(ddw_tmp,contrast = c("type", "IP", "IN"),
                           tidy = TRUE, filterFun =  ihw)
# get those "differentially expressed features"
tmp_significant_windows <- ddw_tmp_results %>%
  dplyr::filter(padj < 0.05) %>% 
  dplyr::pull(row)

# remove those "differentially expressed features"
ddw_mRNAs <- ddw_mRNAs[ !rownames(ddw_mRNAs) %in% tmp_significant_windows, ]

# data is know filtered
# calculate final size factors
ddw_mRNAs <- estimateSizeFactors(ddw_mRNAs)
sizeFactors(ddw) <- sizeFactors(ddw_mRNAs)

# this part assesses the most appropriate dispersion model
# either "local" or "parametric"
# we generate some graphs then assess visually
# although at the end they have automatic way to decide
decide_fit <- TRUE

parametric_ddw  <- estimateDispersions(ddw, fitType="parametric")
if(decide_fit){
  local_ddw  <- estimateDispersions(ddw, fitType="local")
}

# plot parametric fit just for manual viewing
png("dispersion_plot.png", width = 800, height = 600)
plotDispEsts(parametric_ddw, main="Parametric fit")
dev.off()

# plot local fit just for manual viewing
if(decide_fit){
  # Save as PNG
  png("dispersion_plot_local.png", width = 800, height = 600)
  plotDispEsts(local_ddw, main="Local fit")
  dev.off()}

# calculate some stats for the dispersion.
# here the docs pay attention to the median value (how different is the dispersion from the fit)
parametricResid <- na.omit(with(mcols(parametric_ddw),abs(log(dispGeneEst)-log(dispFit))))
if(decide_fit){
  localResid <- na.omit(with(mcols(local_ddw),abs(log(dispGeneEst)-log(dispFit))))
  residDf <- data.frame(residuals=c(parametricResid,localResid),
                        fitType=c(rep("parametric",length(parametricResid)),
                                  rep("local",length(localResid))))
  summary(residDf)
}

# plot overlaid histogram to compare counts with high dispersion between the two fits
# again just for manual assessment
if(decide_fit){
  png("residuals_histogram.png", width = 1000, height = 800)
  
  print(
    ggplot(residDf, aes(x = residuals, fill = fitType)) + 
      scale_fill_manual(values = c("darkred", "darkblue")) + 
      geom_histogram(alpha = 0.5, position='identity', bins = 100) + 
      theme_bw()
  )
  
  dev.off()
}

# this is the important part. it will output which it decides to use!!!
if(decide_fit){
  summary(localResid)
  if (median(localResid) <= median(parametricResid)){
    cat("chosen fitType: local")
    ddw <- local_ddw
  }else{
    cat("chosen fitType: parametric")
    ddw <- parametric_ddw
  }
  rm(local_ddw,parametric_ddw,residDf,parametricResid,localResid)
}else{
  ddw <- parametric_ddw
  rm(parametric_ddw)
}

# based on the output above, we will do the DEWseq with the fitType either parametric or local
# in this run i found it to be parametric so...
ddw <- estimateDispersions(ddw, fitType = "parametric", quiet = TRUE)
ddw <- nbinomWaldTest(ddw)

resultWindows <- resultsDEWSeq(ddw,
                               contrast = c("type", "IP", "IN"),
                               tidy = TRUE) %>% as_tibble
resultWindows

resultWindows[,"p_adj_IHW"] <- adj_pvalues(ihw(pSlidingWindows ~ baseMean, 
                                               data = resultWindows,
                                               alpha = 0.05,
                                               nfolds = 10))

resultWindows <- resultWindows %>% 
  mutate(significant = resultWindows$p_adj_IHW < 0.01)

sum(resultWindows$significant)

resultWindows %>%
  filter(significant) %>%
  arrange(desc(log2FoldChange)) %>%
  .[["gene_name"]] %>%
  unique %>%
  head(200)

write.csv(resultWindows, "DEWSeq_results.csv", row.names = FALSE)

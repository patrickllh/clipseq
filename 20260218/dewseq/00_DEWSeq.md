# DEWSeq

Set-up dependencies for HTSeq and DEWSeq. Good luck!!!
I kept a record here at the end of this markdown as of 20260227.

### Workflow

0. cp `gencode.v38.annotation.gff3.gz` from `skipper/annotations/`. I moved the bams I want to analyze out of `skipper/20260218/output/bams/dedup/genome/`. Note that each skipper run may generate different bams because it chooses 1 random alignment per multimapped read, during the STAR step.

1. Make flattened annotation.

```bash
htseq-clip annotation -g gencode.v38.annotation.gff3.gz --splitExons --unsorted -o annotation.gz
```

2. Make sliding windows, default size of 50nt.

```bash
htseq-clip createSlidingWindows -i annotation.gz -o slidingwindows.gz
```

3. Make mapped annotation file.
### Note that this will be the annotation file input for DEWSeq!!!

```bash
htseq-clip mapToId -a slidingwindows.gz -o mapping.gz
```

4. Get the 5' alignments. These represent crosslink sites.

```bash
htseq-clip extract -e 1 -s s --ignore -c 52 -i SNRNP70_OSM_IP_1.genome.Aligned.sort.dedup.bam -o SNRNP70_OSM_IP_1.crosslinks.gz
htseq-clip extract -e 1 -s s --ignore -c 52 -i SNRNP70_OSM_IP_2.genome.Aligned.sort.dedup.bam -o SNRNP70_OSM_IP_2.crosslinks.gz
htseq-clip extract -e 1 -s s --ignore -c 52 -i SNRNP70_OSM_IN_1.genome.Aligned.sort.dedup.bam -o SNRNP70_OSM_IN_1.crosslinks.gz
htseq-clip extract -e 1 -s s --ignore -c 52 -i SNRNP70_OSM_IN_2.genome.Aligned.sort.dedup.bam -o SNRNP70_OSM_IN_2.crosslinks.gz
```

5. Count the crosslink sites according to the sliding window coordinates. This takes a while...

```bash
mkdir -p inputFolder
htseq-clip count -i SNRNP70_OSM_IP_1.crosslinks.gz -a slidingwindows.gz -o inputFolder/SNRNP70_OSM_IP_1_counts.gz
htseq-clip count -i SNRNP70_OSM_IP_2.crosslinks.gz -a slidingwindows.gz -o inputFolder/SNRNP70_OSM_IP_2_counts.gz
htseq-clip count -i SNRNP70_OSM_IN_1.crosslinks.gz -a slidingwindows.gz -o inputFolder/SNRNP70_OSM_IN_1_counts.gz
htseq-clip count -i SNRNP70_OSM_IN_2.crosslinks.gz -a slidingwindows.gz -o inputFolder/SNRNP70_OSM_IN_2_counts.gz
```

6. Reformat the counts into an R-friendly matrix.

```bash
htseq-clip createMatrix -i inputFolder -b SNRNP70 -o matrix.gz
```

7. DEWSeq analysis starts here. Note that I prepared this so we can run in the server.

Install the necessary packages. The versions were noted at the end of this markdown.
```R
library(DEWSeq)
library(IHW)
library(tidyverse)
library(data.table)
```

Load the input.
```R
countFile <- "matrix.gz"
annotationFile <- "mapping.gz"

countData <- fread(countFile, sep = "\t") # i did head(countData) to know what to put in colData factors later
annotationData <- fread(annotationFile, sep = "\t")
```

Create DESeq-like object and its metadata. Similar to DESeq2 workflows.
It is important to modify the metadata `type` according to the actual data.
```R
colData <- data.frame(row.names = colnames(countData)[-1],
                      type = factor(
                        c("IN","IN","IP","IP"),
                        levels = c("IN","IP"))) # remember...the first element here is going to be the reference

ddw <- DESeqDataSetFromSlidingWindows(countData  = countData,
                                      colData    = colData,
                                      annotObj   = annotationData,
                                      tidy       = TRUE,
                                      design     = ~type)
```

Filter low counts.
```R
keep <- rowSums(counts(ddw)) >= 10 # i'm thinking of making it just 5
ddw <- ddw[keep,]
```

Calculate preliminary size factor for later calculations
```R
ddw <- estimateSizeFactors(ddw)
```

Filter an mRNA-only list then calculate its size factor then put back to dataset (this is accd to the docs...)
```R
ddw_mRNAs <- ddw[ rowData(ddw)[,"gene_type"] == "protein_coding", ]
ddw_mRNAs <- estimateSizeFactors(ddw_mRNAs)
sizeFactors(ddw) <- sizeFactors(ddw_mRNAs)
```

From the docs: "In general, when there is asymmetry in the data, like enrichment in IP over control, it is a good pratice to call differentially expressed features with DESeq2 and then exclude them for normalisation" so I followed that for now.

```R
# make copy first, after getting the preliminary size factor above
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
```

What seems to happen is that we take our dataset with a size factor calculated only from mRNAs, then using that we calculate for dispersions and a Wald test to create a list of "differentially expressed features". We then remove these mRNAs in the mRNA dataset. From this filtered mRNA dataset, we calculate again a size factor, which we then put into the original, non-filtered dataset containing all RNA biotypes.

```R
# data is know filtered
# calculate final size factors
ddw_mRNAs <- estimateSizeFactors(ddw_mRNAs)
sizeFactors(ddw) <- sizeFactors(ddw_mRNAs)
```

We then assess the best statistical model fitting for our data, either a local or parametric fit.
We essentially make graphs to assess manually for ourselves...but later a dispersian median from either fit is generated.
The fit with the lower median is then used for the final, actual DEWSeq analysis.

```R
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
```

We know see which has better fit, based on the dispersion medians between parametric and local fit.
Pay attention to the output of R, and then use that later for DEWSeq analysis!
```R
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
```
We now do DEWSeq. For this run, parametric was recommended so I put `fitType = "parametric"`:
```R
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
```

How many significant windows did we get?
What genes did we get as targets?
```R
sum(resultWindows$significant)
```
resultWindows %>%
  filter(significant) %>%
  arrange(desc(log2FoldChange)) %>%
  .[["gene_name"]] %>%
  unique %>%
  head(200)
```

Save the results.
```R
write.csv(resultWindows, "DEWSeq_results.csv", row.names = FALSE)
```

# Dependencies here!
I also saved the .yml files in this folder. For R packages, they're in a .csv

### For Skipper
```bash
name: skipper
channels:
  - defaults
  - conda-forge
  - bioconda
  - https://repo.anaconda.com/pkgs/main
  - https://repo.anaconda.com/pkgs/r
dependencies:
  - _libgcc_mutex=0.1=conda_forge
  - _openmp_mutex=4.5=2_gnu
  - _r-mutex=1.0.1=anacondar_1
  - aioeasywebdav=2.4.0=pyha770c72_0
  - aiohappyeyeballs=2.6.1=pyhd8ed1ab_0
  - aiohttp=3.12.15=py311h3778330_0
  - aiosignal=1.4.0=pyhd8ed1ab_0
  - alsa-lib=1.2.9=hd590300_0
  - amply=0.1.6=pyhd8ed1ab_1
  - annotated-types=0.7.0=pyhd8ed1ab_1
  - appdirs=1.4.4=pyhd8ed1ab_1
  - attrs=25.3.0=pyh71513ae_0
  - bcrypt=4.3.0=py311hdae7d1d_1
  - bedtools=2.31.0=hf5e1c6e_3
  - binutils_impl_linux-64=2.44=h4bf12b8_1
  - bioconductor-biobase=2.54.0=r41hc0cfd56_2
  - bioconductor-biocgenerics=0.40.0=r41hdfd78af_0
  - bioconductor-biocio=1.4.0=r41hdfd78af_0
  - bioconductor-biocparallel=1.28.3=r41hc247a5b_1
  - bioconductor-biostrings=2.62.0=r41hc0cfd56_2
  - bioconductor-delayedarray=0.20.0=r41hc0cfd56_2
  - bioconductor-fgsea=1.20.0=r41hc247a5b_2
  - bioconductor-genomeinfodb=1.30.1=r41hdfd78af_0
  - bioconductor-genomeinfodbdata=1.2.7=r41hdfd78af_2
  - bioconductor-genomicalignments=1.30.0=r41hc0cfd56_2
  - bioconductor-genomicranges=1.46.1=r41hc0cfd56_1
  - bioconductor-iranges=2.28.0=r41hc0cfd56_2
  - bioconductor-matrixgenerics=1.6.0=r41hdfd78af_0
  - bioconductor-rhtslib=1.26.0=r41hc0cfd56_2
  - bioconductor-rsamtools=2.10.0=r41hc247a5b_2
  - bioconductor-rtracklayer=1.54.0=r41hd029910_0
  - bioconductor-s4vectors=0.32.4=r41hc0cfd56_0
  - bioconductor-summarizedexperiment=1.24.0=r41hdfd78af_0
  - bioconductor-xvector=0.34.0=r41hc0cfd56_2
  - bioconductor-zlibbioc=1.40.0=r41hc0cfd56_2
  - boto3=1.39.16=pyhd8ed1ab_0
  - botocore=1.39.17=pyge310_1234567_0
  - brotli-python=1.1.0=py311hfdbb021_3
  - bwidget=1.10.1=ha770c72_1
  - bzip2=1.0.8=h4bc722e_7
  - c-ares=1.34.5=hb9d3cd8_0
  - ca-certificates=2025.9.9=h06a4308_0
  - cachetools=5.5.2=pyhd8ed1ab_0
  - cairo=1.16.0=hbbf8b49_1016
  - certifi=2025.10.5=py311h06a4308_0
  - cffi=1.17.1=py311hf29c0ef_0
  - charset-normalizer=3.4.2=pyhd8ed1ab_0
  - click=8.2.1=pyh707e725_0
  - coin-or-cbc=2.10.12=h8b142ea_1
  - coin-or-cgl=0.60.7=h516709c_0
  - coin-or-clp=1.17.8=h1ee7a9c_0
  - coin-or-osi=0.108.10=haf5fa05_0
  - coin-or-utils=2.11.11=hee58242_0
  - coincbc=2.10.12=2_metapackage
  - colorama=0.4.6=pyhd8ed1ab_1
  - coloredlogs=15.0.1=pyhd8ed1ab_4
  - configargparse=1.7.1=pyhe01879c_0
  - connection_pool=0.0.3=pyhd3deb0d_0
  - cryptography=45.0.5=py311hafd3f86_0
  - curl=8.8.0=he654da7_1
  - datrie=0.8.2=py311h9ecbd09_8
  - defusedxml=0.7.1=pyhd8ed1ab_0
  - docutils=0.22=pyhd8ed1ab_0
  - dpath=2.2.0=pyha770c72_0
  - dropbox=12.0.2=pyhd8ed1ab_1
  - eido=0.2.4=pyhd8ed1ab_0
  - exceptiongroup=1.3.0=pyhd8ed1ab_0
  - expat=2.7.1=hecca717_0
  - fastp=0.23.4=hadf994f_3
  - fastqc=0.12.1=hdfd78af_0
  - filechunkio=1.8=py_2
  - font-ttf-dejavu-sans-mono=2.37=hab24e00_0
  - font-ttf-inconsolata=3.000=h77eed37_0
  - font-ttf-source-code-pro=2.038=h77eed37_0
  - font-ttf-ubuntu=0.83=h77eed37_3
  - fontconfig=2.14.2=h14ed4e7_0
  - fonts-conda-ecosystem=1=0
  - fonts-conda-forge=1=0
  - freetype=2.12.1=h267a509_2
  - fribidi=1.0.10=h36c2ea0_0
  - frozenlist=1.7.0=py311h52bc045_0
  - ftputil=5.1.0=pyhd8ed1ab_0
  - gcc_impl_linux-64=15.1.0=h4393ad2_4
  - gettext=0.25.1=h3f43e3d_1
  - gettext-tools=0.25.1=h3f43e3d_1
  - gfortran_impl_linux-64=15.1.0=h3b9cdf2_4
  - giflib=5.2.2=hd590300_0
  - gitdb=4.0.12=pyhd8ed1ab_0
  - gitpython=3.1.45=pyhff2d567_0
  - google-api-core=2.25.1=pyhd8ed1ab_0
  - google-api-python-client=2.177.0=pyhff2d567_0
  - google-auth=2.40.3=pyhd8ed1ab_0
  - google-auth-httplib2=0.2.0=pyhd8ed1ab_1
  - google-cloud-core=2.4.3=py311h06a4308_0
  - google-cloud-storage=3.2.0=pyhd8ed1ab_0
  - google-crc32c=1.7.1=py311h0973507_0
  - google-resumable-media=2.7.2=pyhd8ed1ab_2
  - googleapis-common-protos=1.70.0=pyhd8ed1ab_0
  - graphite2=1.3.14=h5888daf_0
  - gsl=2.7=he838d99_0
  - gxx_impl_linux-64=15.1.0=h6a1bac1_4
  - h2=4.2.0=pyhd8ed1ab_0
  - harfbuzz=7.3.0=hdb3a94d_0
  - homer=4.11=pl5262h4ac6f70_9
  - hpack=4.1.0=pyhd8ed1ab_0
  - htslib=1.20=h81da01d_0
  - httplib2=0.22.0=pyhd8ed1ab_1
  - humanfriendly=10.0=pyh707e725_8
  - hyperframe=6.1.0=pyhd8ed1ab_0
  - icu=72.1=hcb278e6_0
  - idna=3.10=pyhd8ed1ab_1
  - importlib-metadata=8.7.0=pyhe01879c_1
  - iniconfig=2.0.0=pyhd8ed1ab_1
  - isa-l=2.31.1=hb9d3cd8_1
  - jinja2=3.1.6=pyhd8ed1ab_0
  - jmespath=1.0.1=pyhd8ed1ab_1
  - jsonschema=4.25.0=pyhe01879c_0
  - jsonschema-specifications=2025.4.1=pyh29332c3_0
  - jupyter_core=5.8.1=pyh31011fe_0
  - kernel-headers_linux-64=5.14.0=he073ed8_2
  - keyutils=1.6.1=h166bdaf_0
  - krb5=1.21.3=h659f571_0
  - lcms2=2.15=h7f713cb_2
  - ld_impl_linux-64=2.44=h1423503_1
  - lerc=4.0.0=h0aef613_1
  - libasprintf=0.25.1=h3f43e3d_1
  - libasprintf-devel=0.25.1=h3f43e3d_1
  - libblas=3.9.0=32_h59b9bed_openblas
  - libcblas=3.9.0=32_he106b2a_openblas
  - libcrc32c=1.1.2=h9c3ff4c_0
  - libcups=2.3.3=h4637d8d_4
  - libcurl=8.8.0=hca28451_1
  - libdeflate=1.19=hd590300_0
  - libedit=3.1.20250104=pl5321h7949ede_0
  - libev=4.33=hd590300_2
  - libexpat=2.7.1=hecca717_0
  - libffi=3.4.6=h2dba641_1
  - libgcc=15.1.0=h767d61c_4
  - libgcc-devel_linux-64=15.1.0=h4c094af_104
  - libgcc-ng=15.1.0=h69a702a_4
  - libgettextpo=0.25.1=h3f43e3d_1
  - libgettextpo-devel=0.25.1=h3f43e3d_1
  - libgfortran=15.1.0=h69a702a_4
  - libgfortran-ng=15.1.0=h69a702a_4
  - libgfortran5=15.1.0=hcea5267_4
  - libglib=2.78.1=hebfc3b9_0
  - libgomp=15.1.0=h767d61c_4
  - libiconv=1.18=h4ce23a2_1
  - libidn2=2.3.8=ha4ef2c3_0
  - libjpeg-turbo=2.1.5.1=hd590300_1
  - liblapack=3.9.0=32_h7ac8fdf_openblas
  - liblapacke=3.9.0=32_he2f377e_openblas
  - liblzma=5.8.1=hb9d3cd8_2
  - liblzma-devel=5.8.1=hb9d3cd8_2
  - libnghttp2=1.58.0=h47da74e_1
  - libnsl=2.0.1=hb9d3cd8_1
  - libopenblas=0.3.30=pthreads_h94d23a6_1
  - libpng=1.6.43=h2797004_0
  - libsanitizer=15.1.0=h97b714f_4
  - libsodium=1.0.20=h4ab18f5_0
  - libsqlite=3.46.0=hde9e2c9_0
  - libssh2=1.11.0=h0841786_0
  - libstdcxx=15.1.0=h8f9b012_4
  - libstdcxx-devel_linux-64=15.1.0=h4c094af_104
  - libstdcxx-ng=15.1.0=h4852527_4
  - libtiff=4.6.0=h29866fb_1
  - libunistring=0.9.10=h7f98852_0
  - libuuid=2.38.1=h0b41bf4_0
  - libwebp-base=1.6.0=hd42ef1d_0
  - libxcb=1.15=h0b41bf4_0
  - libxcrypt=4.4.36=hd590300_1
  - libxml2=2.11.5=h0d562d8_0
  - libzlib=1.2.13=h4ab18f5_6
  - logmuse=0.2.8=pyhd8ed1ab_1
  - make=4.4.1=hb9d3cd8_2
  - markdown-it-py=3.0.0=pyhd8ed1ab_1
  - markupsafe=3.0.2=py311h2dc5d0c_1
  - mdurl=0.1.2=pyhd8ed1ab_1
  - multidict=6.6.3=py311h2dc5d0c_0
  - mysql-connector-c=6.1.11=h659d440_1008
  - nbformat=5.10.4=pyhd8ed1ab_1
  - ncurses=6.5=h2d0b736_3
  - numpy=2.3.2=py311h2e04523_0
  - oauth2client=4.1.3=pyhd8ed1ab_1
  - openjdk=20.0.0=h8e330f5_0
  - openssl=3.5.1=h7b32b05_0
  - packaging=25.0=pyh29332c3_1
  - pandas=2.3.1=py311hed34c8f_0
  - pandoc=2.19.2=h32600fe_2
  - pango=1.50.14=heaa33ce_1
  - paramiko=3.5.1=pyhd8ed1ab_0
  - pcre2=10.40=hc3806b6_0
  - pephubclient=0.4.4=pyhd8ed1ab_1
  - peppy=0.40.7=pyhd8ed1ab_2
  - perl=5.32.1=7_hd590300_perl5
  - pip=25.2=pyh8b19718_0
  - pixman=0.46.4=h537e5f6_0
  - plac=1.4.5=pyhd8ed1ab_0
  - platformdirs=4.3.8=pyhe01879c_0
  - pluggy=1.6.0=pyhd8ed1ab_0
  - ply=3.11=pyhd8ed1ab_3
  - prettytable=3.16.0=pyhd8ed1ab_0
  - propcache=0.3.1=py311h2dc5d0c_0
  - proto-plus=1.26.1=pyhd8ed1ab_0
  - protobuf=5.28.3=py311hfdbb021_0
  - psutil=7.0.0=py311h9ecbd09_0
  - pthread-stubs=0.4=hb9d3cd8_1002
  - pulp=2.7.0=py311h38be061_1
  - pyasn1=0.6.1=pyhd8ed1ab_2
  - pyasn1-modules=0.4.2=pyhd8ed1ab_0
  - pycparser=2.22=pyh29332c3_1
  - pydantic=2.11.7=pyh3cfb1c2_0
  - pydantic-core=2.33.2=py311hdae7d1d_0
  - pygments=2.19.2=pyhd8ed1ab_0
  - pynacl=1.5.0=py311h9ecbd09_4
  - pyopenssl=25.1.0=pyhd8ed1ab_0
  - pyparsing=3.2.3=pyhe01879c_2
  - pysftp=0.2.9=py_1
  - pysocks=1.7.1=pyha55dd90_7
  - pytest=8.4.1=pyhd8ed1ab_0
  - python=3.11.4=hab00c5b_0_cpython
  - python-dateutil=2.9.0.post0=pyhe01879c_2
  - python-fastjsonschema=2.21.1=pyhd8ed1ab_0
  - python-irodsclient=2.1.0=pyhd8ed1ab_0
  - python-tzdata=2025.2=pyhd8ed1ab_0
  - python_abi=3.11=8_cp311
  - pytz=2025.2=pyhd8ed1ab_0
  - pyu2f=0.1.5=pyhd8ed1ab_1
  - pyyaml=6.0.2=py311h2dc5d0c_2
  - r-askpass=1.1=r41h06615bd_3
  - r-assertthat=0.2.1=r41hc72bb7e_3
  - r-backports=1.4.1=r41h06615bd_1
  - r-base=4.1.3=h63daf7b_11
  - r-base64enc=0.1_3=r41h06615bd_1005
  - r-bh=1.81.0_1=r41hc72bb7e_0
  - r-bit=4.0.5=r41h06615bd_0
  - r-bit64=4.0.5=r41h06615bd_1
  - r-bitops=1.0_7=r41h06615bd_1
  - r-blob=1.2.4=r41hc72bb7e_0
  - r-broom=1.0.5=r41hc72bb7e_0
  - r-bslib=0.5.0=r41hc72bb7e_0
  - r-cachem=1.0.8=r41h57805ef_0
  - r-callr=3.7.3=r41hc72bb7e_0
  - r-cellranger=1.1.0=r41hc72bb7e_1005
  - r-cli=3.6.1=r41h38f115c_0
  - r-clipr=0.8.0=r41hc72bb7e_1
  - r-colorspace=2.1_0=r41h133d619_0
  - r-cowplot=1.1.1=r41hc72bb7e_1
  - r-cpp11=0.5.2=r41hc72bb7e_0
  - r-crayon=1.5.2=r41hc72bb7e_1
  - r-curl=4.3.3=r41hf9611b0_2
  - r-data.table=1.14.8=r41h133d619_0
  - r-dbi=1.1.3=r41hc72bb7e_1
  - r-dbplyr=2.3.2=r41hc72bb7e_0
  - r-digest=0.6.31=r41h38f115c_0
  - r-dplyr=1.1.2=r41ha503ecb_0
  - r-dtplyr=1.3.1=r41hc72bb7e_0
  - r-ellipsis=0.3.2=r41h06615bd_1
  - r-evaluate=0.21=r41hc72bb7e_0
  - r-fansi=1.0.4=r41h133d619_0
  - r-farver=2.1.1=r41h7525677_1
  - r-fastmap=1.1.1=r41h38f115c_0
  - r-fastmatch=1.1_3=r41h06615bd_1
  - r-fontawesome=0.5.1=r41hc72bb7e_0
  - r-forcats=1.0.0=r41hc72bb7e_0
  - r-formatr=1.14=r41hc72bb7e_0
  - r-fs=1.6.2=r41ha503ecb_0
  - r-futile.logger=1.4.3=r41hc72bb7e_1004
  - r-futile.options=1.0.1=r41hc72bb7e_1003
  - r-gargle=1.4.0=r41h785f33e_0
  - r-generics=0.1.3=r41hc72bb7e_1
  - r-ggdendro=0.1.23=r41hc72bb7e_1
  - r-ggplot2=3.4.2=r41hc72bb7e_0
  - r-ggrepel=0.9.3=r41h38f115c_0
  - r-ggupset=0.3.0=r41hc72bb7e_2
  - r-glue=1.6.2=r41h06615bd_1
  - r-googledrive=2.1.0=r41hc72bb7e_0
  - r-googlesheets4=1.1.0=r41h785f33e_0
  - r-gridextra=2.3=r41hc72bb7e_1004
  - r-gtable=0.3.3=r41hc72bb7e_0
  - r-haven=2.5.2=r41h38f115c_0
  - r-highr=0.10=r41hc72bb7e_0
  - r-hms=1.1.3=r41hc72bb7e_0
  - r-htmltools=0.5.5=r41h38f115c_0
  - r-httr=1.4.6=r41hc72bb7e_0
  - r-ids=1.0.1=r41hc72bb7e_2
  - r-isoband=0.2.7=r41h38f115c_1
  - r-jquerylib=0.1.4=r41hc72bb7e_1
  - r-jsonlite=1.8.5=r41h57805ef_0
  - r-knitr=1.43=r41hc72bb7e_0
  - r-labeling=0.4.2=r41hc72bb7e_2
  - r-lambda.r=1.2.4=r41hc72bb7e_2
  - r-lattice=0.21_8=r41h133d619_0
  - r-lifecycle=1.0.3=r41hc72bb7e_1
  - r-lubridate=1.9.2=r41h133d619_1
  - r-magrittr=2.0.3=r41h06615bd_1
  - r-mass=7.3_58.3=r41h133d619_0
  - r-matrix=1.5_4.1=r41h316c678_0
  - r-matrixstats=1.0.0=r41h57805ef_0
  - r-memoise=2.0.1=r41hc72bb7e_1
  - r-mgcv=1.8_42=r41he1ae0d6_0
  - r-mime=0.12=r41h06615bd_1
  - r-modelr=0.1.11=r41hc72bb7e_0
  - r-munsell=0.5.0=r41hc72bb7e_1005
  - r-nlme=3.1_162=r41hac0b197_0
  - r-openssl=2.0.6=r41habfbb5e_0
  - r-pillar=1.9.0=r41hc72bb7e_0
  - r-pkgconfig=2.0.3=r41hc72bb7e_2
  - r-prettyunits=1.1.1=r41hc72bb7e_2
  - r-processx=3.8.1=r41h133d619_0
  - r-progress=1.2.2=r41hc72bb7e_3
  - r-ps=1.7.5=r41h133d619_0
  - r-purrr=1.0.1=r41h133d619_0
  - r-r6=2.5.1=r41hc72bb7e_1
  - r-rappdirs=0.3.3=r41h06615bd_1
  - r-rcolorbrewer=1.1_3=r41h785f33e_1
  - r-rcpp=1.0.10=r41h38f115c_0
  - r-rcurl=1.98_1.12=r41hf9611b0_1
  - r-readr=2.1.4=r41h38f115c_0
  - r-readxl=1.4.2=r41h81ef4d7_0
  - r-rematch=1.0.1=r41hc72bb7e_1005
  - r-rematch2=2.1.2=r41hc72bb7e_2
  - r-reprex=2.0.2=r41hc72bb7e_1
  - r-restfulr=0.0.15=r41h73dbb54_0
  - r-rjson=0.2.21=r41h7525677_2
  - r-rlang=1.1.1=r41ha503ecb_0
  - r-rmarkdown=2.22=r41hc72bb7e_0
  - r-rstudioapi=0.14=r41hc72bb7e_1
  - r-rtsne=0.16=r41h37cf8d7_1
  - r-rvest=1.0.3=r41hc72bb7e_1
  - r-sass=0.4.6=r41ha503ecb_0
  - r-scales=1.2.1=r41hc72bb7e_1
  - r-selectr=0.4_2=r41hc72bb7e_2
  - r-snow=0.4_4=r41hc72bb7e_1
  - r-stringi=1.7.12=r41hc0c3e09_1
  - r-stringr=1.5.0=r41h785f33e_0
  - r-sys=3.4.2=r41h57805ef_0
  - r-tibble=3.2.1=r41h133d619_1
  - r-tidyr=1.3.0=r41h38f115c_0
  - r-tidyselect=1.2.0=r41hc72bb7e_0
  - r-tidyverse=1.3.2=r41hc72bb7e_1
  - r-timechange=0.2.0=r41h38f115c_0
  - r-tinytex=0.45=r41hc72bb7e_0
  - r-tzdb=0.4.0=r41ha503ecb_0
  - r-utf8=1.2.3=r41h133d619_0
  - r-uuid=1.1_0=r41h06615bd_1
  - r-vctrs=0.6.2=r41ha503ecb_0
  - r-vgam=1.1_8=r41hac0b197_0
  - r-viridis=0.6.3=r41hc72bb7e_0
  - r-viridislite=0.4.1=r41hc72bb7e_1
  - r-vroom=1.6.3=r41ha503ecb_0
  - r-withr=2.5.0=r41hc72bb7e_1
  - r-xfun=0.39=r41ha503ecb_0
  - r-xml=3.99_0.14=r41hc38eee6_1
  - r-xml2=1.3.4=r41h1ad5fc0_1
  - r-yaml=2.3.7=r41h133d619_0
  - readline=8.2=h8c095d6_2
  - referencing=0.36.2=pyh29332c3_0
  - requests=2.32.4=pyhd8ed1ab_0
  - reretry=0.11.8=pyhd8ed1ab_1
  - rich=14.1.0=pyhe01879c_0
  - rpds-py=0.26.0=py311hdae7d1d_0
  - rsa=4.9.1=pyhd8ed1ab_0
  - s3transfer=0.13.1=pyhd8ed1ab_0
  - samtools=1.17=hd87286a_2
  - sed=4.9=h6688a6e_0
  - seqkit=2.10.1=he881be0_0
  - setuptools=80.9.0=pyhff2d567_0
  - setuptools-scm=8.3.1=pyhd8ed1ab_0
  - shellingham=1.5.4=pyhd8ed1ab_1
  - six=1.17.0=pyhe01879c_1
  - skewer=0.2.2=hc9558a2_3
  - slacker=0.14.0=pyhd8ed1ab_1
  - smart_open=7.3.0.post1=pyhe01879c_0
  - smmap=5.0.2=pyhd8ed1ab_0
  - snakemake=7.32.3=hdfd78af_1
  - snakemake-minimal=7.32.3=pyhdfd78af_1
  - star=2.7.10b=h6b7c446_1
  - stone=3.3.2=pyhd8ed1ab_0
  - stopit=1.1.2=pyhd8ed1ab_1
  - sysroot_linux-64=2.34=h087de78_2
  - tabulate=0.9.0=pyhd8ed1ab_2
  - throttler=1.2.2=pyhd8ed1ab_0
  - tk=8.6.13=noxft_h4845f30_101
  - tktable=2.10=h8d826fa_7
  - tomli=2.2.1=pyhe01879c_2
  - toposort=1.10=pyhd8ed1ab_1
  - traitlets=5.14.3=pyhd8ed1ab_1
  - typer=0.16.0=pyh167b9f4_0
  - typer-slim=0.16.0=pyhe01879c_0
  - typer-slim-standard=0.16.0=hf964461_0
  - typing-extensions=4.14.1=h4440ef1_0
  - typing-inspection=0.4.1=pyhd8ed1ab_0
  - typing_extensions=4.14.1=pyhe01879c_0
  - tzdata=2025b=h78e105d_0
  - ubiquerg=0.8.0=pyhd8ed1ab_0
  - ucsc-bedgraphtobigwig=455=h2a80c09_1
  - unzip=6.0=h7f98852_3
  - uritemplate=4.2.0=pyhd8ed1ab_0
  - urllib3=2.5.0=pyhd8ed1ab_0
  - veracitools=0.1.3=py_0
  - wcwidth=0.2.13=pyhd8ed1ab_1
  - wget=1.21.4=hda4d442_0
  - wheel=0.45.1=pyhd8ed1ab_1
  - wrapt=1.17.2=py311h9ecbd09_0
  - xorg-fixesproto=5.0=hb9d3cd8_1003
  - xorg-inputproto=2.3.2=hb9d3cd8_1003
  - xorg-kbproto=1.0.7=hb9d3cd8_1003
  - xorg-libice=1.1.2=hb9d3cd8_0
  - xorg-libsm=1.2.6=he73a12e_0
  - xorg-libx11=1.8.9=h8ee46fc_0
  - xorg-libxau=1.0.12=hb9d3cd8_0
  - xorg-libxdmcp=1.1.5=hb9d3cd8_0
  - xorg-libxext=1.3.4=h0b41bf4_2
  - xorg-libxfixes=5.0.3=h7f98852_1004
  - xorg-libxi=1.7.10=h4bc722e_1
  - xorg-libxrender=0.9.11=hd590300_0
  - xorg-libxt=1.3.0=hd590300_1
  - xorg-libxtst=1.2.5=h4bc722e_0
  - xorg-recordproto=1.14.2=hb9d3cd8_1003
  - xorg-renderproto=0.11.1=hb9d3cd8_1003
  - xorg-xextproto=7.3.0=hb9d3cd8_1004
  - xorg-xproto=7.0.31=hb9d3cd8_1008
  - xz=5.8.1=hbcc6ac9_2
  - xz-gpl-tools=5.8.1=hbcc6ac9_2
  - xz-tools=5.8.1=hb9d3cd8_2
  - yaml=0.2.5=h280c20c_3
  - yarl=1.20.1=py311h2dc5d0c_0
  - yte=1.8.1=pyha770c72_0
  - zipp=3.23.0=pyhd8ed1ab_0
  - zlib=1.2.13=h4ab18f5_6
  - zstandard=0.23.0=py311h9ecbd09_2
  - zstd=1.5.6=ha6fb4c9_0
prefix: /work/pa-hilario/anaconda3/envs/skipper
```

### For HTseq
```bash
name: htseq-clip
channels:
  - conda-forge
  - bioconda
  - https://repo.anaconda.com/pkgs/main
  - https://repo.anaconda.com/pkgs/r
dependencies:
  - _openmp_mutex=4.5=20_gnu
  - brotli=1.2.0=hed03a55_1
  - brotli-bin=1.2.0=hb03c661_1
  - bzip2=1.0.8=hda65f42_9
  - c-ares=1.11.0=h470a237_1
  - ca-certificates=2026.1.4=hbd8a1cb_0
  - certifi=2024.8.30=pyhd8ed1ab_0
  - curl=7.71.1=he644dc0_8
  - cycler=0.11.0=pyhd8ed1ab_0
  - fonttools=4.38.0=py37h540881e_0
  - freetype=2.12.1=h267a509_2
  - htseq=2.0.3=py37hfcd875c_0
  - jpeg=9e=h0b41bf4_3
  - kiwisolver=1.4.4=py37h7cecad7_0
  - krb5=1.17.2=h926e7f8_0
  - lcms2=2.12=hddcbb42_0
  - ld_impl_linux-64=2.45.1=bootstrap_ha15bf96_1
  - libblas=3.11.0=5_h4a7cf45_openblas
  - libbrotlicommon=1.2.0=hb03c661_1
  - libbrotlidec=1.2.0=hb03c661_1
  - libbrotlienc=1.2.0=hb03c661_1
  - libcblas=3.11.0=5_h0358290_openblas
  - libcurl=7.71.1=hcdd3856_8
  - libdeflate=1.0=h14c3975_1
  - libedit=3.1.20250104=pl5321h7949ede_0
  - libev=4.33=hd590300_2
  - libffi=3.4.6=h2dba641_1
  - libgcc=15.2.0=he0feb66_18
  - libgcc-ng=15.2.0=h69a702a_18
  - libgfortran=15.2.0=h69a702a_18
  - libgfortran5=15.2.0=h68bc16d_18
  - libgomp=15.2.0=he0feb66_18
  - liblapack=3.11.0=5_h47877c9_openblas
  - libnghttp2=1.41.0=hab1572f_1
  - libnsl=2.0.1=hb9d3cd8_1
  - libopenblas=0.3.30=pthreads_h94d23a6_4
  - libpng=1.6.43=h2797004_0
  - libsqlite=3.46.0=hde9e2c9_0
  - libssh2=1.10.0=haa6b8db_3
  - libstdcxx=15.2.0=h934c35e_18
  - libstdcxx-ng=15.2.0=hdf11a46_18
  - libtiff=4.2.0=hf544144_3
  - libwebp-base=1.6.0=hd42ef1d_0
  - libzlib=1.2.13=h4ab18f5_6
  - matplotlib-base=3.5.3=py37hf395dca_2
  - munkres=1.0.7=py_1
  - ncurses=6.5=h2d0b736_3
  - numpy=1.21.6=py37h976b520_0
  - olefile=0.47=pyhd8ed1ab_0
  - openjpeg=2.4.0=hb52868f_1
  - openssl=1.1.1w=hd590300_0
  - packaging=23.2=pyhd8ed1ab_0
  - pillow=8.2.0=py37h4600e1f_1
  - pip=24.0=pyhd8ed1ab_0
  - pyparsing=3.1.4=pyhd8ed1ab_0
  - pysam=0.15.3=py37hda2845c_1
  - python=3.7.12=hb7a2778_100_cpython
  - python-dateutil=2.9.0=pyhd8ed1ab_0
  - python_abi=3.7=4_cp37m
  - readline=8.3=h853b02a_0
  - setuptools=69.0.3=pyhd8ed1ab_0
  - six=1.16.0=pyh6c4a22f_0
  - sqlite=3.46.0=h6d4b2fc_0
  - tk=8.6.13=noxft_h4845f30_101
  - typing-extensions=4.7.1=hd8ed1ab_0
  - typing_extensions=4.7.1=pyha770c72_0
  - unicodedata2=14.0.0=py37h540881e_1
  - wheel=0.42.0=pyhd8ed1ab_0
  - xz=5.2.6=h166bdaf_0
  - zlib=1.2.13=h4ab18f5_6
  - zstd=1.5.6=ha6fb4c9_0
  - pip:
      - htseq-clip==2.19.0b0
prefix: /work/pa-hilario/anaconda3/envs/htseq-clip
```

### For R
```bash
name: r_env
channels:
  - conda-forge
  - bioconda
  - https://repo.anaconda.com/pkgs/main
  - https://repo.anaconda.com/pkgs/r
dependencies:
  - _openmp_mutex=4.5=20_gnu
  - _python_abi3_support=1.0=hd8ed1ab_2
  - _r-mutex=1.0.1=anacondar_1
  - anyio=4.12.1=pyhcf101f3_0
  - argon2-cffi=25.1.0=pyhd8ed1ab_0
  - argon2-cffi-bindings=25.1.0=py314h5bd0f2a_2
  - arrow=1.4.0=pyhcf101f3_0
  - asttokens=3.0.1=pyhd8ed1ab_0
  - async-lru=2.2.0=pyhcf101f3_0
  - attrs=25.4.0=pyhcf101f3_1
  - babel=2.18.0=pyhcf101f3_0
  - backports.zstd=1.3.0=py314h680f03e_0
  - beautifulsoup4=4.14.3=pyha770c72_0
  - binutils_impl_linux-64=2.45.1=default_hfdba357_101
  - binutils_linux-64=2.45.1=default_h4852527_101
  - bleach=6.3.0=pyhcf101f3_1
  - bleach-with-css=6.3.0=hbca2aae_1
  - brotli-python=1.2.0=py314h3de4e8d_1
  - bwidget=1.10.1=ha770c72_1
  - bzip2=1.0.8=hda65f42_9
  - c-ares=1.34.6=hb03c661_0
  - ca-certificates=2026.2.25=hbd8a1cb_0
  - cached-property=1.5.2=hd8ed1ab_1
  - cached_property=1.5.2=pyha770c72_1
  - cairo=1.18.4=he90730b_1
  - certifi=2026.2.25=pyhd8ed1ab_0
  - cffi=2.0.0=py314h4a8dc5f_1
  - charset-normalizer=3.4.4=pyhd8ed1ab_0
  - comm=0.2.3=pyhe01879c_0
  - cpython=3.14.3=py314hd8ed1ab_101
  - curl=8.18.0=hcf29cc6_1
  - debugpy=1.8.20=py314h42812f9_0
  - decorator=5.2.1=pyhd8ed1ab_0
  - defusedxml=0.7.1=pyhd8ed1ab_0
  - exceptiongroup=1.3.1=pyhd8ed1ab_0
  - executing=2.2.1=pyhd8ed1ab_0
  - font-ttf-dejavu-sans-mono=2.37=hab24e00_0
  - font-ttf-inconsolata=3.000=h77eed37_0
  - font-ttf-source-code-pro=2.038=h77eed37_0
  - font-ttf-ubuntu=0.83=h77eed37_3
  - fontconfig=2.17.1=h27c8c51_0
  - fonts-conda-ecosystem=1=0
  - fonts-conda-forge=1=hc364b38_1
  - fqdn=1.5.1=pyhd8ed1ab_1
  - fribidi=1.0.16=hb03c661_0
  - gcc_impl_linux-64=15.2.0=he420e7e_18
  - gcc_linux-64=15.2.0=h862fb80_21
  - gfortran_impl_linux-64=15.2.0=h281d09f_18
  - gfortran_linux-64=15.2.0=h944c81e_21
  - graphite2=1.3.14=hecca717_2
  - gsl=2.7=he838d99_0
  - gxx_impl_linux-64=15.2.0=hda75c37_18
  - gxx_linux-64=15.2.0=h1fb793f_21
  - h11=0.16.0=pyhcf101f3_1
  - h2=4.3.0=pyhcf101f3_0
  - harfbuzz=12.3.2=h6083320_0
  - hpack=4.1.0=pyhd8ed1ab_0
  - httpcore=1.0.9=pyh29332c3_0
  - httpx=0.28.1=pyhd8ed1ab_0
  - hyperframe=6.1.0=pyhd8ed1ab_0
  - icu=78.2=h33c6efd_0
  - idna=3.11=pyhd8ed1ab_0
  - importlib-metadata=8.7.0=pyhe01879c_1
  - importlib_resources=6.5.2=pyhd8ed1ab_0
  - ipykernel=7.2.0=pyha191276_1
  - ipython=9.10.0=pyh53cf698_0
  - ipython_pygments_lexers=1.1.1=pyhd8ed1ab_0
  - isoduration=20.11.0=pyhd8ed1ab_1
  - jedi=0.19.2=pyhd8ed1ab_1
  - jinja2=3.1.6=pyhcf101f3_1
  - json5=0.13.0=pyhd8ed1ab_0
  - jsonpointer=3.0.0=pyhcf101f3_3
  - jsonschema=4.26.0=pyhcf101f3_0
  - jsonschema-specifications=2025.9.1=pyhcf101f3_0
  - jsonschema-with-format-nongpl=4.26.0=hcf101f3_0
  - jupyter-lsp=2.3.0=pyhcf101f3_0
  - jupyter_client=8.8.0=pyhcf101f3_0
  - jupyter_core=5.9.1=pyhc90fa1f_0
  - jupyter_events=0.12.0=pyhe01879c_0
  - jupyter_server=2.17.0=pyhcf101f3_0
  - jupyter_server_terminals=0.5.4=pyhcf101f3_0
  - jupyterlab=4.5.5=pyhd8ed1ab_0
  - jupyterlab_pygments=0.3.0=pyhd8ed1ab_2
  - jupyterlab_server=2.28.0=pyhcf101f3_0
  - kernel-headers_linux-64=5.14.0=he073ed8_3
  - keyutils=1.6.3=hb9d3cd8_0
  - krb5=1.22.2=ha1258a1_0
  - lark=1.3.1=pyhd8ed1ab_0
  - ld_impl_linux-64=2.45.1=default_hbd61a6d_101
  - lerc=4.0.0=h0aef613_1
  - libblas=3.11.0=5_h4a7cf45_openblas
  - libcblas=3.11.0=5_h0358290_openblas
  - libcurl=8.18.0=hcf29cc6_1
  - libdeflate=1.25=h17f619e_0
  - libedit=3.1.20250104=pl5321h7949ede_0
  - libev=4.33=hd590300_2
  - libexpat=2.7.4=hecca717_0
  - libffi=3.5.2=h3435931_0
  - libfreetype=2.14.1=ha770c72_0
  - libfreetype6=2.14.1=h73754d4_0
  - libgcc=15.2.0=he0feb66_18
  - libgcc-devel_linux-64=15.2.0=hcc6f6b0_118
  - libgcc-ng=15.2.0=h69a702a_18
  - libgfortran=15.2.0=h69a702a_18
  - libgfortran-ng=15.2.0=h69a702a_18
  - libgfortran5=15.2.0=h68bc16d_18
  - libglib=2.86.4=h6548e54_1
  - libgomp=15.2.0=he0feb66_18
  - libiconv=1.18=h3b78370_2
  - libjpeg-turbo=3.1.2=hb03c661_0
  - liblapack=3.11.0=5_h47877c9_openblas
  - liblzma=5.8.2=hb03c661_0
  - liblzma-devel=5.8.2=hb03c661_0
  - libmpdec=4.0.0=hb03c661_1
  - libnghttp2=1.67.0=had1ee68_0
  - libopenblas=0.3.30=pthreads_h94d23a6_4
  - libpng=1.6.55=h421ea60_0
  - libsanitizer=15.2.0=h90f66d4_18
  - libsodium=1.0.18=h36c2ea0_1
  - libsqlite=3.51.2=hf4e2dac_0
  - libssh2=1.11.1=hcf80075_0
  - libstdcxx=15.2.0=h934c35e_18
  - libstdcxx-devel_linux-64=15.2.0=hd446a21_118
  - libstdcxx-ng=15.2.0=hdf11a46_18
  - libtiff=4.7.1=h9d88235_1
  - libuuid=2.41.3=h5347b49_0
  - libuv=1.51.0=hb03c661_1
  - libwebp-base=1.6.0=hd42ef1d_0
  - libxcb=1.17.0=h8a09558_0
  - libxml2=2.15.1=he237659_1
  - libxml2-16=2.15.1=hca6bf5a_1
  - libzlib=1.3.1=hb9d3cd8_2
  - make=4.4.1=hb9d3cd8_2
  - markupsafe=3.0.3=pyh7db6752_0
  - matplotlib-inline=0.2.1=pyhd8ed1ab_0
  - mistune=3.2.0=pyhcf101f3_0
  - nbclient=0.10.4=pyhd8ed1ab_0
  - nbconvert-core=7.17.0=pyhcf101f3_0
  - nbformat=5.10.4=pyhd8ed1ab_1
  - ncurses=6.5=h2d0b736_3
  - nest-asyncio=1.6.0=pyhd8ed1ab_1
  - notebook=7.5.4=pyhcf101f3_0
  - notebook-shim=0.2.4=pyhd8ed1ab_1
  - openssl=3.6.1=h35e630c_1
  - overrides=7.7.0=pyhd8ed1ab_1
  - packaging=26.0=pyhcf101f3_0
  - pandoc=3.9=ha770c72_0
  - pandocfilters=1.5.0=pyhd8ed1ab_0
  - pango=1.56.4=hadf4263_0
  - parso=0.8.6=pyhcf101f3_0
  - pcre2=10.47=haa7fec5_0
  - pexpect=4.9.0=pyhd8ed1ab_1
  - pip=26.0.1=pyh145f28c_0
  - pixman=0.46.4=h54a6638_1
  - platformdirs=4.9.2=pyhcf101f3_0
  - prometheus_client=0.24.1=pyhd8ed1ab_0
  - prompt-toolkit=3.0.52=pyha770c72_0
  - psutil=7.2.2=py314h0f05182_0
  - pthread-stubs=0.4=hb9d3cd8_1002
  - ptyprocess=0.7.0=pyhd8ed1ab_1
  - pure_eval=0.2.3=pyhd8ed1ab_1
  - pycparser=2.22=pyh29332c3_1
  - pygments=2.19.2=pyhd8ed1ab_0
  - pysocks=1.7.1=pyha55dd90_7
  - python=3.14.3=h32b2ec7_101_cp314
  - python-dateutil=2.9.0.post0=pyhe01879c_2
  - python-fastjsonschema=2.21.2=pyhe01879c_0
  - python-gil=3.14.3=h4df99d1_101
  - python-json-logger=2.0.7=pyhd8ed1ab_0
  - python-tzdata=2025.3=pyhd8ed1ab_0
  - python_abi=3.14=8_cp314
  - pytz=2025.2=pyhd8ed1ab_0
  - pyyaml=6.0.3=py314h67df5f8_1
  - pyzmq=27.1.0=py312hda471dd_2
  - r-askpass=1.2.1=r45h54b55ab_1
  - r-assertthat=0.2.1=r45hc72bb7e_6
  - r-backports=1.5.0=r45h54b55ab_2
  - r-base=4.5.2=h1fbe982_4
  - r-base64enc=0.1_6=r45h54b55ab_0
  - r-bit=4.6.0=r45h54b55ab_1
  - r-bit64=4.6.0_1=r45h54b55ab_1
  - r-blob=1.3.0=r45hc72bb7e_0
  - r-boot=1.3_32=r45hc72bb7e_1
  - r-broom=1.0.12=r45hc72bb7e_0
  - r-bslib=0.10.0=r45hc72bb7e_0
  - r-cachem=1.1.0=r45h54b55ab_2
  - r-callr=3.7.6=r45hc72bb7e_2
  - r-caret=7.0_1=r45h54b55ab_0
  - r-cellranger=1.1.0=r45hc72bb7e_1008
  - r-class=7.3_23=r45h54b55ab_1
  - r-cli=3.6.5=r45h3697838_1
  - r-clipr=0.8.0=r45hc72bb7e_4
  - r-clock=0.7.4=r45h3697838_0
  - r-cluster=2.1.8.2=r45heaba542_0
  - r-codetools=0.2_20=r45hc72bb7e_2
  - r-colorspace=2.1_2=r45h54b55ab_0
  - r-commonmark=2.0.0=r45h54b55ab_1
  - r-conflicted=1.2.0=r45h785f33e_3
  - r-cpp11=0.5.3=r45h785f33e_0
  - r-crayon=1.5.3=r45hc72bb7e_2
  - r-crul=1.6.0=r45hc72bb7e_1
  - r-curl=7.0.0=r45h10955f1_1
  - r-data.table=1.17.8=r45h1c8cec4_1
  - r-dbi=1.3.0=r45hc72bb7e_0
  - r-dbplyr=2.5.2=r45hc72bb7e_0
  - r-diagram=1.6.5=r45ha770c72_4
  - r-digest=0.6.39=r45h3697838_0
  - r-dplyr=1.2.0=r45h3697838_0
  - r-dtplyr=1.3.3=r45hc72bb7e_0
  - r-e1071=1.7_17=r45h3697838_0
  - r-ellipsis=0.3.2=r45h54b55ab_4
  - r-essentials=4.5=r45hd8ed1ab_2006
  - r-evaluate=1.0.5=r45hc72bb7e_1
  - r-fansi=1.0.7=r45h54b55ab_0
  - r-farver=2.1.2=r45h3697838_2
  - r-fastmap=1.2.0=r45h3697838_2
  - r-fontawesome=0.5.3=r45hc72bb7e_1
  - r-forcats=1.0.1=r45hc72bb7e_0
  - r-foreach=1.5.2=r45hc72bb7e_4
  - r-foreign=0.8_91=r45h54b55ab_0
  - r-formatr=1.14=r45hc72bb7e_3
  - r-fs=1.6.6=r45h3697838_1
  - r-future=1.69.0=r45h785f33e_0
  - r-future.apply=1.20.2=r45hc72bb7e_0
  - r-gargle=1.6.1=r45h785f33e_0
  - r-generics=0.1.4=r45hc72bb7e_1
  - r-ggplot2=4.0.2=r45h785f33e_0
  - r-gistr=0.9.0=r45hc72bb7e_4
  - r-glmnet=4.1_10=r45ha36cffa_1
  - r-globals=0.19.0=r45hc72bb7e_0
  - r-glue=1.8.0=r45h54b55ab_1
  - r-googledrive=2.1.2=r45hc72bb7e_1
  - r-googlesheets4=1.1.2=r45h785f33e_1
  - r-gower=1.0.2=r45h54b55ab_0
  - r-gtable=0.3.6=r45hc72bb7e_1
  - r-hardhat=1.4.2=r45hc72bb7e_1
  - r-haven=2.5.5=r45h6d565e7_1
  - r-hexbin=1.28.5=r45heaba542_1
  - r-highr=0.11=r45hc72bb7e_2
  - r-hms=1.1.4=r45hc72bb7e_0
  - r-htmltools=0.5.9=r45h3697838_0
  - r-htmlwidgets=1.6.4=r45h785f33e_4
  - r-httpcode=0.3.0=r45ha770c72_5
  - r-httpuv=1.6.16=r45h6d565e7_1
  - r-httr=1.4.8=r45hc72bb7e_0
  - r-ids=1.0.1=r45hc72bb7e_5
  - r-ipred=0.9_15=r45h54b55ab_2
  - r-irdisplay=1.1=r45hd8ed1ab_4
  - r-irkernel=1.3.2=r45h785f33e_3
  - r-isoband=0.3.0=r45h3697838_0
  - r-iterators=1.0.14=r45hc72bb7e_4
  - r-jquerylib=0.1.4=r45hc72bb7e_4
  - r-jsonlite=2.0.0=r45h54b55ab_1
  - r-kernsmooth=2.23_26=r45ha0a88a1_1
  - r-knitr=1.51=r45hc72bb7e_0
  - r-labeling=0.4.3=r45hc72bb7e_2
  - r-later=1.4.7=r45h3697838_0
  - r-lattice=0.22_9=r45h54b55ab_0
  - r-lava=1.8.2=r45hc72bb7e_0
  - r-lazyeval=0.2.2=r45h54b55ab_6
  - r-lifecycle=1.0.5=r45hc72bb7e_0
  - r-listenv=0.10.0=r45hc72bb7e_0
  - r-lobstr=1.2.0=r45h3697838_0
  - r-lubridate=1.9.5=r45h54b55ab_0
  - r-magrittr=2.0.4=r45h54b55ab_0
  - r-maps=3.4.3=r45h54b55ab_1
  - r-mass=7.3_65=r45h54b55ab_0
  - r-matrix=1.7_4=r45h0e4624f_1
  - r-memoise=2.0.1=r45hc72bb7e_4
  - r-mgcv=1.9_4=r45h0e4624f_0
  - r-mime=0.13=r45h54b55ab_1
  - r-modelmetrics=1.2.2.2=r45h3697838_5
  - r-modelr=0.1.11=r45hc72bb7e_3
  - r-munsell=0.5.1=r45hc72bb7e_2
  - r-nlme=3.1_168=r45heaba542_1
  - r-nnet=7.3_20=r45h54b55ab_1
  - r-numderiv=2016.8_1.1=r45hc72bb7e_7
  - r-openssl=2.3.4=r45h50f7d53_0
  - r-otel=0.2.0=r45hc72bb7e_1
  - r-parallelly=1.46.1=r45h54b55ab_0
  - r-pbdzmq=0.3_14=r45hded8526_1
  - r-pillar=1.11.1=r45hc72bb7e_0
  - r-pkgconfig=2.0.3=r45hc72bb7e_5
  - r-plyr=1.8.9=r45h3697838_3
  - r-prettyunits=1.2.0=r45hc72bb7e_2
  - r-proc=1.19.0.1=r45h3697838_1
  - r-processx=3.8.6=r45h54b55ab_1
  - r-prodlim=2025.04.28=r45h3697838_1
  - r-progress=1.2.3=r45hc72bb7e_2
  - r-progressr=0.18.0=r45hc72bb7e_0
  - r-promises=1.5.0=r45hc72bb7e_1
  - r-proxy=0.4_29=r45h54b55ab_0
  - r-pryr=0.1.6=r45h3697838_3
  - r-ps=1.9.1=r45h54b55ab_1
  - r-purrr=1.2.1=r45h54b55ab_0
  - r-quantmod=0.4.28=r45hc72bb7e_1
  - r-r6=2.6.1=r45hc72bb7e_1
  - r-ragg=1.5.0=r45h9f1dc4d_1
  - r-randomforest=4.7_1.2=r45heaba542_1
  - r-rappdirs=0.3.4=r45h54b55ab_0
  - r-rbokeh=0.5.2=r45hc72bb7e_5
  - r-rcolorbrewer=1.1_3=r45h785f33e_4
  - r-rcpp=1.1.1=r45h3697838_0
  - r-rcppeigen=0.3.4.0.2=r45h3704496_1
  - r-readr=2.2.0=r45h3697838_0
  - r-readxl=1.4.5=r45h10e25cc_1
  - r-recipes=1.3.1=r45hc72bb7e_1
  - r-recommended=4.5=r45hd8ed1ab_1008
  - r-rematch=2.0.0=r45hc72bb7e_2
  - r-rematch2=2.1.2=r45hc72bb7e_5
  - r-repr=1.1.7=r45h785f33e_2
  - r-reprex=2.1.1=r45hc72bb7e_2
  - r-reshape2=1.4.5=r45h3697838_0
  - r-rlang=1.1.7=r45h3697838_0
  - r-rmarkdown=2.30=r45hc72bb7e_0
  - r-rpart=4.1.24=r45h54b55ab_1
  - r-rstudioapi=0.18.0=r45hc72bb7e_0
  - r-rvest=1.0.5=r45hc72bb7e_1
  - r-s7=0.2.1=r45h54b55ab_0
  - r-sass=0.4.10=r45h3697838_1
  - r-scales=1.4.0=r45hc72bb7e_1
  - r-selectr=0.5_1=r45hc72bb7e_0
  - r-shape=1.4.6.1=r45ha770c72_2
  - r-shiny=1.13.0=r45h785f33e_0
  - r-sourcetools=0.1.7_1=r45h3697838_3
  - r-sparsevctrs=0.3.6=r45h54b55ab_0
  - r-spatial=7.3_18=r45h54b55ab_1
  - r-squarem=2021.1=r45hc72bb7e_4
  - r-stringi=1.8.7=r45h3d52c89_2
  - r-stringr=1.6.0=r45h785f33e_0
  - r-survival=3.8_6=r45h54b55ab_0
  - r-sys=3.4.3=r45h54b55ab_1
  - r-systemfonts=1.3.1=r45h74f4acd_0
  - r-textshaping=1.0.4=r45h74f4acd_0
  - r-tibble=3.3.1=r45h54b55ab_0
  - r-tidyr=1.3.2=r45h3697838_0
  - r-tidyselect=1.2.1=r45hc72bb7e_2
  - r-tidyverse=2.0.0=r45h785f33e_3
  - r-timechange=0.4.0=r45h3697838_0
  - r-timedate=4052.112=r45hc72bb7e_0
  - r-tinytex=0.58=r45hc72bb7e_0
  - r-triebeard=0.4.1=r45h3697838_4
  - r-ttr=0.24.4=r45h54b55ab_2
  - r-tzdb=0.5.0=r45h3697838_2
  - r-urltools=1.7.3.1=r45h3697838_1
  - r-utf8=1.2.6=r45h54b55ab_1
  - r-uuid=1.2_2=r45h54b55ab_0
  - r-vctrs=0.7.1=r45h3697838_0
  - r-viridislite=0.4.3=r45hc72bb7e_0
  - r-vroom=1.7.0=r45h3697838_0
  - r-withr=3.0.2=r45hc72bb7e_1
  - r-xfun=0.56=r45h3697838_0
  - r-xml2=1.5.2=r45he78afff_0
  - r-xtable=1.8_8=r45hc72bb7e_0
  - r-xts=0.14.1=r45h54b55ab_1
  - r-yaml=2.3.12=r45h54b55ab_0
  - r-zoo=1.8_15=r45h54b55ab_0
  - readline=8.3=h853b02a_0
  - referencing=0.37.0=pyhcf101f3_0
  - requests=2.32.5=pyhcf101f3_1
  - rfc3339-validator=0.1.4=pyhd8ed1ab_1
  - rfc3986-validator=0.1.1=pyh9f0ad1d_0
  - rfc3987-syntax=1.1.0=pyhe01879c_1
  - rpds-py=0.30.0=py314h2e6c369_0
  - sed=4.9=h6688a6e_0
  - send2trash=2.1.0=pyha191276_1
  - setuptools=82.0.0=pyh332efcf_0
  - six=1.17.0=pyhe01879c_1
  - sniffio=1.3.1=pyhd8ed1ab_2
  - soupsieve=2.8.3=pyhd8ed1ab_0
  - stack_data=0.6.3=pyhd8ed1ab_1
  - sysroot_linux-64=2.34=h087de78_3
  - terminado=0.18.1=pyhc90fa1f_1
  - tinycss2=1.4.0=pyhd8ed1ab_0
  - tk=8.6.13=noxft_h366c992_103
  - tktable=2.10=h8d826fa_7
  - tomli=2.4.0=pyhcf101f3_0
  - tornado=6.5.4=py314h7b0bd38_0
  - traitlets=5.14.3=pyhd8ed1ab_1
  - typing-extensions=4.15.0=h396c80c_0
  - typing_extensions=4.15.0=pyhcf101f3_0
  - typing_utils=0.1.0=pyhd8ed1ab_1
  - tzdata=2025c=hc9c84f9_1
  - uri-template=1.3.0=pyhd8ed1ab_1
  - urllib3=2.6.3=pyhd8ed1ab_0
  - wcwidth=0.6.0=pyhd8ed1ab_0
  - webcolors=25.10.0=pyhd8ed1ab_0
  - webencodings=0.5.1=pyhd8ed1ab_3
  - websocket-client=1.9.0=pyhd8ed1ab_0
  - xorg-libice=1.1.2=hb9d3cd8_0
  - xorg-libsm=1.2.6=he73a12e_0
  - xorg-libx11=1.8.13=he1eb515_0
  - xorg-libxau=1.0.12=hb03c661_1
  - xorg-libxdmcp=1.1.5=hb03c661_1
  - xorg-libxext=1.3.7=hb03c661_0
  - xorg-libxrender=0.9.12=hb9d3cd8_0
  - xorg-libxt=1.3.1=hb9d3cd8_0
  - xz=5.8.2=ha02ee65_0
  - xz-gpl-tools=5.8.2=ha02ee65_0
  - xz-tools=5.8.2=hb03c661_0
  - yaml=0.2.5=h280c20c_3
  - zeromq=4.3.5=h59595ed_1
  - zipp=3.23.0=pyhcf101f3_1
  - zlib=1.3.1=hb9d3cd8_2
  - zstd=1.5.7=hb78ec9c_6
prefix: /work/pa-hilario/anaconda3/envs/r_env
```

### R packages downloaded inside R
```R
Package	Version
abind	1.4-8
askpass	1.2.1
assertthat	0.2.1
backports	1.5.0
base	4.5.2
base64enc	0.1-6
BH	1.90.0-1
Biobase	2.70.0
BiocGenerics	0.56.0
BiocManager	1.30.27
BiocParallel	1.44.0
BiocVersion	3.22.0
bit	4.6.0
bit64	4.6.0-1
blob	1.3.0
boot	1.3-32
broom	1.0.12
bslib	0.10.0
cachem	1.1.0
callr	3.7.6
caret	7.0-1
cellranger	1.1.0
class	7.3-23
cli	3.6.5
clipr	0.8.0
clock	0.7.4
cluster	2.1.8.2
codetools	0.2-20
colorspace	2.1-2
commonmark	2.0.0
compiler	4.5.2
conflicted	1.2.0
cpp11	0.5.3
crayon	1.5.3
crul	1.6.0
curl	7.0.0
data.table	1.18.2.1
datasets	4.5.2
DBI	1.3.0
dbplyr	2.5.2
DelayedArray	0.36.0
DESeq2	1.50.2
DEWSeq	1.24.0
diagram	1.6.5
digest	0.6.39
dplyr	1.2.0
dtplyr	1.3.3
e1071	1.7-17
ellipsis	0.3.2
evaluate	1.0.5
fansi	1.0.7
farver	2.1.2
fastmap	1.2.0
fdrtool	1.2.18
fontawesome	0.5.3
forcats	1.0.1
foreach	1.5.2
foreign	0.8-91
formatR	1.14
fs	1.6.6
futile.logger	1.4.9
futile.options	1.0.1
future	1.69.0
future.apply	1.20.2
gargle	1.6.1
generics	0.1.4
GenomicRanges	1.62.1
ggplot2	4.0.2
gistr	0.9.0
glmnet	4.1-10
globals	0.19.0
glue	1.8.0
googledrive	2.1.2
googlesheets4	1.1.2
gower	1.0.2
graphics	4.5.2
grDevices	4.5.2
grid	4.5.2
gtable	0.3.6
hardhat	1.4.2
haven	2.5.5
hexbin	1.28.5
highr	0.11
hms	1.1.4
htmltools	0.5.9
htmlwidgets	1.6.4
httpcode	0.3.0
httpuv	1.6.16
httr	1.4.8
ids	1.0.1
IHW	1.38.0
ipred	0.9-15
IRanges	2.44.0
IRdisplay	1.1
IRkernel	1.3.2
isoband	0.3.0
iterators	1.0.14
jquerylib	0.1.4
jsonlite	2.0.0
KernSmooth	2.23-26
knitr	1.51
labeling	0.4.3
lambda.r	1.2.4
later	1.4.7
lattice	0.22-9
lava	1.8.2
lazyeval	0.2.2
lifecycle	1.0.5
listenv	0.10.0
lobstr	1.2.0
locfit	1.5-9.12
lpsymphony	1.38.0
lubridate	1.9.5
magrittr	2.0.4
maps	3.4.3
MASS	7.3-65
Matrix	1.7-4
MatrixGenerics	1.22.0
matrixStats	1.5.0
memoise	2.0.1
methods	4.5.2
mgcv	1.9-4
mime	0.13
ModelMetrics	1.2.2.2
modelr	0.1.11
munsell	0.5.1
nlme	3.1-168
nnet	7.3-20
numDeriv	2016.8-1.1
openssl	2.3.5
otel	0.2.0
parallel	4.5.2
parallelly	1.46.1
pbdZMQ	0.3-14
pillar	1.11.1
pkgconfig	2.0.3
plyr	1.8.9
prettyunits	1.2.0
pROC	1.19.0.1
processx	3.8.6
prodlim	2025.04.28
progress	1.2.3
progressr	0.18.0
promises	1.5.0
proxy	0.4-29
pryr	0.1.6
ps	1.9.1
purrr	1.2.1
quantmod	0.4.28
R.methodsS3	1.8.2
R.oo	1.27.1
R.utils	2.13.0
R6	2.6.1
ragg	1.5.0
randomForest	4.7-1.2
rappdirs	0.3.4
rbokeh	0.5.2
RColorBrewer	1.1-3
Rcpp	1.1.1
RcppArmadillo	15.2.3-1
RcppEigen	0.3.4.0.2
readr	2.2.0
readxl	1.4.5
recipes	1.3.1
rematch	2.0.0
rematch2	2.1.2
repr	1.1.7
reprex	2.1.1
reshape2	1.4.5
rlang	1.1.7
rmarkdown	2.3
rpart	4.1.24
rstudioapi	0.18.0
rvest	1.0.5
S4Arrays	1.10.1
S4Vectors	0.48.0
S7	0.2.1
sass	0.4.10
scales	1.4.0
selectr	0.5-1
Seqinfo	1.0.0
shape	1.4.6.1
shiny	1.13.0
slam	0.1-55
snow	0.4-4
sourcetools	0.1.7-1
SparseArray	1.10.8
sparsevctrs	0.3.6
spatial	7.3-18
splines	4.5.2
SQUAREM	2021.1
stats	4.5.2
stats4	4.5.2
stringi	1.8.7
stringr	1.6.0
SummarizedExperiment	1.40.0
survival	3.8-6
sys	3.4.3
systemfonts	1.3.1
tcltk	4.5.2
textshaping	1.0.4
tibble	3.3.1
tidyr	1.3.2
tidyselect	1.2.1
tidyverse	2.0.0
timechange	0.4.0
timeDate	4052.112
tinytex	0.58
tools	4.5.2
triebeard	0.4.1
TTR	0.24.4
tzdb	0.5.0
urltools	1.7.3.1
utf8	1.2.6
utils	4.5.2
uuid	1.2-2
vctrs	0.7.1
viridisLite	0.4.3
vroom	1.7.0
withr	3.0.2
xfun	0.56
xml2	1.5.2
xtable	1.8-8
xts	0.14.1
XVector	0.50.0
yaml	2.3.12
zoo	1.8-15
```
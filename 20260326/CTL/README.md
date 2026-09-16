# CLIP-Seq analysis
- **Sample:** SNRNP70 under control conditions
- **Pipeline** Custom, with elements from Skipper (Boyle et al., Cell Genom., 2023) and seCLIP (Blue et al., Nat Protoc., 2022)
- **Library prep:** Denaturing Buffer and Denaturing Wash CLIP-Seq
- **Platform:** Illumina (PAGS)
- **Read type:** Paired-end 150bp (but here only read 1 is used)
- **Date:** 20260327

# Custom Pipeline

To run this pipeline, first, fulfill the requirements written in the section below.

Then, simply:

```bash
./00_Pipeline.sh
```

For this dataset, I reuse the genome index, mRNA annotations from prior runs, as well as the relevant public genome data.

### Therefore I skip steps 5 and 19, and genome/annotation download steps.

```bash
mkdir input output
cp ~/clipseq/20260326/OSM/input/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz ./input
cp ~/clipseq/20260326/OSM/input/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta ./input
cp ~/clipseq/20260326/OSM/input/gencode.v38.annotation.gff3.gz ./input
cp ~/clipseq/20260326/OSM/input/gencode.v38.annotation.gtf ./input
cp -r ~/clipseq/20260326/OSM/output/genome_ref ./output
cp -r ~/clipseq/20260326/OSM/output/mrna.gtf ./output
```

For each step, please see the corresponding numbered script:

01. Generate QC reports for raw reads using `FastQC`.
02. Reads were trimmed for 3' end InvRiL19 adapters using a two-pass `Cutadapt` strategy (Blue et al., Nat Protoc., 2022).
03. 5'end 10nt UMIs were extracted and placed into read headers preceded with ":"  (Boyle et al., Cell Genom., 2023) using `fastp`.
04. Generate QC reports for reads after trimming using `FastQC`.
~~05. Create a genome index based on `GRCh38_no_alt_analysis_set_GCA_000001405.15` using `STAR`.~~
06. Map reads, assigning 1 site randomly among up to 100 mappable sites (Boyle et al., Cell Genom., 2023) using `STAR`.
07. Sort and index bams using `SAMtools` (Boyle et al., Cell Genom., 2023).
08. Deduplicate reads based on UMI using `UMICollapse` (Boyle et al., Cell Genom., 2023).
09. Index deduplicated reads using `SAMtools` (Boyle et al., Cell Genom., 2023).
10. Define statistically significant peaks (based on local background) as a bed file using `CLIPper` (Blue et al., Nat Protoc., 2022).
11. Count total reads using `SAMtools` to aid in next step (Blue et al., Nat Protoc., 2022).
12. Merge overlapping peaks using custom scripts from `merge_peaks` (Blue et al., Nat Protoc., 2022).
13. ***(MAIN OUTPUT)*** Intersect peaks among replicates then blacklist mVenus peaks to create consensus, using `BEDtools`.
14. ***(MAIN OUTPUT)*** Create CPM-normalized bigwigs using `deepTools`.
15. ***(MAIN OUTPUT)*** Calculate enrichment of IP over input per consensus peak region using `BEDtools` and some `awk` commands.
16. Convert bigwigs to bedgraphs for latter analysis using `bigWigToBedGraph` from the UCSC tools suite.
17. Limit bedgraph CPM signals only to consensus peak region,s using `BEDtools`.
18. Re-convert bedgraphs to bigwigs, this time with CPM signals only in consensus peak regions, using `bedGraphToBigWig` from the UCSC tools suite.
~~19. Manually create an mRNA annotation file by filtering `gencode.v38.annotation.gtf` using some `awk` commands.~~
20. ***(MAIN OUTPUT)*** Perform metagene analysis using `deepTools`.

***Note: To save time for some commands, I parallelize the slurm job submission for multiple samples that don't need to be done sequentially. In this case, please pay attention to `#SBATCH --array=0-3`, which needs to be adjusted according to the number of samples. Technically, some other steps could also be parellelized, but I just opted for loops if they could be conveniently finished fast anyway.***

### Input requirements

Before running the pipeline, you just need to do 3 major configurations:

1. Please install the major tools required for the pipeline, as you can see from the script names:

- FastQC
- Cutadapt
- fastp
- STAR
- SAMtools
- UMICollapse
- CLIPper
- merge_peaks
- BEDtools
- deepTools
- UCSC suite: bigWigToBedGraph, bedGraphToBigWig

For currently working versions, please see dependencies listed in my different env in the section below.
Then correspond them to the env mentioned in each script.

2. Please manually edit the input/output names in each script to match your samples.

The naming pattern is typically: "RBP-Condition-rep#" (e.g. `SNRNP70-Con-rep1.`) or just "RBP-Condition" in latter steps.

Additionally, edit the file paths for your `conda`, `UMICollapse` and `merge_peaks` scripts, Slurm cluster nodes/time limit, etc.

3. In the same directory as all scripts, create an `input` directory. Then, place the following inside `input`:

- CLIP IP and input fastqs, for two replicates each.
- `GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta` from ENCODE.
- `gencode.v38.annotation.gff3.gz` from GENCODE.
- `gencode.v38.annotation.gtf` from GENCODE.

For example:
```bash
cp /project/okamura-lab/raw_sequencing_data/20260218_CLIP_NIG/aviti_run/20260129_AV244513_260129-AVITI24-FCB/SNRNP70-Con-rep1_L1_R1.fastq.gz .
cp /project/okamura-lab/raw_sequencing_data/20260218_CLIP_NIG/aviti_run/20260129_AV244513_260129-AVITI24-FCB/SNRNP70-Con-rep2_L2_R1.fastq.gz .
cp /project/okamura-lab/raw_sequencing_data/20260218_CLIP_NIG/aviti_run/20260129_AV244513_260129-AVITI24-FCB/mVenus-Con-rep5_L1_R1.fastq.gz .
cp /project/okamura-lab/raw_sequencing_data/20260218_CLIP_NIG/aviti_run/20260129_AV244513_260129-AVITI24-FCB/mVenus-Con-rep6_L2_R1.fastq.gz .
# md5:
# b5eb6dc307fe53c19d64c8fe2974f158  SNRNP70-Con-rep1_L1_R1.fastq.gz
# 36342236268856f62a1bc615281265cb  SNRNP70-Con-rep2_L2_R1.fastq.gz
# 59c5b878e0af85d3b574f7fa6bfe82f3  mVenus-Con-rep5_L1_R1.fastq.gz
# 4f9275e6333fd64838cd227f46626155  mVenus-Con-rep6_L2_R1.fastq.gz

wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
# unzipped md5: a6da8681616c05eb542f1d91606a7b2f
gunzip gunzip GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gff3.gz
# unzipped md5: bd88c2c7cb5b24815c93a1d6e7f34adc

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
# unzipped md5: 5acd8565f4dc199fea5bab1c9af854f0
gunzip gencode.v38.annotation.gtf.gz
```

# Installation

I note here instructions for the more troublesome set-ups I personally experienced.
For the rest, please install in a conda env according to typical methods (e.g. `conda install -c bioconda`, or according to their docs/github).
For convenient install of some tools (e.g. FastQC, fastp, STAR, UMICollapse, etc.), follow Skipper (Boyle et al., 2023; Xu et al., 2024).
At the end of this markdown I recorded the exact working versions I got.

### For CLIPper:

```bash
# I found from https://zhuanlan.zhihu.com/p/7046358695
git clone https://github.com/YeoLab/clipper.git
cd clipper
git log -1
# commit 8bbc3db08e6bc438edb9aa11df8ed3b8f4719367
conda config --set channel_priority flexible
conda env create -f environment3.yml
conda activate clipper3
# key step 1 (unmentioned in original clipper repo): get compiler even without root permissions
conda install -c conda-forge gxx_linux-64
conda install -c conda-forge gxx_impl_linux-64
conda install -c conda-forge gcc_linux-64
conda install -c conda-forge gcc_impl_linux-64
# key step 2 (unmentioned in original clipper repo): set compiler
cd ~/anaconda3/envs/clipper3/bin/
ln -s ~/anaconda3/envs/clipper/libexec/gcc/x86_64-conda_cos6-linux-gnu/7.3.0/gcc gcc
ln -s ~/anaconda3/envs/clipper/bin/x86_64-conda_cos6-linux-gnu-g++ g++
cd ~/clipper
pip install .
clipper --help # check
# ver 2.1.2
```

### For merge_peaks:

```bash
conda create -n merge_peaks
conda activate merge_peaks
conda install perl
conda install perl-app-cpanminus # this downgraded perl to 5.2
cpanm Statistics::Basic
cpanm Statistics::Distributions
conda install -c conda-forge r-base
cpanm Statistics::R
conda install -c bioconda idr # note: this version i got at the time only takes 9 column bed files
conda install -c bioconda -c conda-forge cwltool
conda install samtools
conda install bedtools
git clone https://github.com/YeoLab/merge_peaks.git #commit aedc0a14d4ba109ee65678a3201a52c5bb6ad473 
```

Conda environments proven to work are shown as follows:

***Note: Because I use these environments for other tasks as well, they may contain packages or tools unnecessary to this specific pipeline. Sorry.***

```bash
name: cutadapt
channels:
  - conda-forge
  - bioconda
  - defaults
  - https://repo.anaconda.com/pkgs/main
  - https://repo.anaconda.com/pkgs/r
dependencies:
  - _openmp_mutex=4.5=20_gnu
  - bzip2=1.0.8=hda65f42_9
  - ca-certificates=2026.1.4=hbd8a1cb_0
  - cffi=2.0.0=py312h460c074_1
  - cutadapt=5.2=py312h0fa9677_0
  - dnaio=1.2.2=py312hf67a6ed_0
  - icu=78.2=h33c6efd_0
  - isa-l=2.31.1=hb9d3cd8_1
  - ld_impl_linux-64=2.45.1=default_hbd61a6d_101
  - libexpat=2.7.4=hecca717_0
  - libffi=3.5.2=h3435931_0
  - libgcc=15.2.0=he0feb66_17
  - libgcc-ng=15.2.0=h69a702a_17
  - libgomp=15.2.0=he0feb66_17
  - liblzma=5.8.2=hb03c661_0
  - libnsl=2.0.1=hb9d3cd8_1
  - libsqlite=3.51.2=hf4e2dac_0
  - libstdcxx=15.2.0=h934c35e_17
  - libstdcxx-ng=15.2.0=hdf11a46_17
  - libuuid=2.41.3=h5347b49_0
  - libxcrypt=4.4.36=hd590300_1
  - libzlib=1.3.1=hb9d3cd8_2
  - ncurses=6.5=h2d0b736_3
  - openssl=3.6.1=h35e630c_1
  - packaging=26.0=pyhcf101f3_0
  - pbzip2=1.1.13=h1fcc475_2
  - pigz=2.8=h421ea60_2
  - pip=26.0.1=pyh8b19718_0
  - pycparser=2.22=pyh29332c3_1
  - python=3.12.12=hd63d673_2_cpython
  - python-isal=1.8.0=py312h4c3975b_1
  - python-zlib-ng=1.0.0=py312hc77a125_1
  - python_abi=3.12=8_cp312
  - readline=8.3=h853b02a_0
  - setuptools=82.0.0=pyh332efcf_0
  - tk=8.6.13=noxft_h366c992_103
  - tzdata=2025c=hc9c84f9_1
  - wheel=0.46.3=pyhd8ed1ab_0
  - xopen=2.0.2=pyh707e725_2
  - zlib-ng=2.3.3=hceb46e0_1
  - zstandard=0.25.0=py312h5253ce2_1
  - zstd=1.5.7=hb78ec9c_6
prefix: /work/pa-hilario/anaconda3/envs/cutadapt

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

name: clipper3
channels:
  - defaults
  - bioconda
  - conda-forge
  - https://repo.anaconda.com/pkgs/main
  - https://repo.anaconda.com/pkgs/r
dependencies:
  - _libgcc_mutex=0.1=main
  - _openmp_mutex=4.5=7_kmp_llvm
  - bedtools=2.29.2=hc088bd4_0
  - binutils_impl_linux-64=2.33.1=he1b5a44_7
  - binutils_linux-64=2.33.1=h9595d00_17
  - blas=1.0=mkl
  - bzip2=1.0.8=h7b6447c_0
  - ca-certificates=2026.2.25=hbd8a1cb_0
  - certifi=2024.8.30=pyhd8ed1ab_0
  - curl=7.71.0=hbc83047_0
  - cycler=0.10.0=py37_0
  - cython=0.29.20=py37he6710b0_0
  - dbus=1.13.16=hb2f20db_0
  - expat=2.2.9=he6710b0_2
  - fontconfig=2.13.0=h9420a91_0
  - freetype=2.10.2=h5ab3b9f_0
  - gcc_impl_linux-64=7.3.0=habb00fd_1
  - gcc_linux-64=7.3.0=h553295d_17
  - glib=2.65.0=h3eb4bd4_0
  - gst-plugins-base=1.14.0=hbbd80ab_1
  - gstreamer=1.14.0=hb31296c_0
  - gxx_impl_linux-64=7.3.0=hdf63c60_1
  - gxx_linux-64=7.3.0=h553295d_17
  - htseq=0.11.3=py37hb3f55d8_0
  - icu=58.2=he6710b0_3
  - intel-openmp=2020.1=217
  - joblib=0.15.1=py_0
  - jpeg=9b=h024ee3a_2
  - kiwisolver=1.2.0=py37hfd86e86_0
  - krb5=1.18.2=h173b8e3_0
  - ld_impl_linux-64=2.33.1=h53a641e_7
  - libcurl=7.71.0=h20c2e04_0
  - libdeflate=1.0=h14c3975_1
  - libedit=3.1.20191231=h7b6447c_0
  - libffi=3.3=he6710b0_1
  - libgcc=7.2.0=h69d50b8_2
  - libgcc-ng=9.1.0=hdf63c60_0
  - libgfortran-ng=7.3.0=hdf63c60_0
  - libpng=1.6.37=hbc83047_0
  - libssh2=1.9.0=h1ba5d50_1
  - libstdcxx-ng=9.1.0=hdf63c60_0
  - libuuid=1.0.3=h1bed415_2
  - libxcb=1.13=h1bed415_1
  - libxml2=2.9.10=he19cac6_1
  - llvm-openmp=20.1.8=h4922eb0_0
  - matplotlib=3.2.2=0
  - matplotlib-base=3.2.2=py37hef1b27d_0
  - mkl=2020.1=217
  - mkl-service=2.3.0=py37he904b0f_0
  - mkl_fft=1.1.0=py37h23d657b_0
  - mkl_random=1.1.1=py37h0573a6f_0
  - ncurses=6.2=he6710b0_1
  - numpy=1.18.5=py37ha1c710e_0
  - numpy-base=1.18.5=py37hde5b4d6_0
  - openssl=1.1.1h=h516909a_0
  - pandas=1.0.5=py37h0573a6f_0
  - pcre=8.44=he6710b0_0
  - pip=20.1.1=py37_1
  - pybedtools=0.8.1=py37h4ef193e_1
  - pyparsing=2.4.7=py_0
  - pyqt=5.9.2=py37h05f1152_2
  - pysam=0.15.3=py37hda2845c_1
  - python=3.7.7=hcff3b4d_5
  - python-dateutil=2.8.1=py_0
  - python_abi=3.7=1_cp37m
  - pytz=2020.1=py_0
  - qt=5.9.7=h5867ecd_1
  - readline=8.0=h7b6447c_0
  - samtools=1.7=1
  - scikit-learn=0.23.1=py37h423224d_0
  - scipy=1.5.0=py37h0b6359f_0
  - setuptools=47.3.1=py37_0
  - sip=4.19.8=py37hf484d3e_0
  - six=1.15.0=py_0
  - sqlite=3.32.3=h62c20be_0
  - threadpoolctl=2.1.0=pyh5ca1d4c_0
  - tk=8.6.10=hbc83047_0
  - tornado=6.0.4=py37h7b6447c_1
  - wheel=0.34.2=py37_0
  - xz=5.2.5=h7b6447c_0
  - zlib=1.2.11=h7b6447c_3
  - pip:
      - clipper==2.1.2
prefix: /work/pa-hilario/anaconda3/envs/clipper3

name: deeptools
channels:
  - bioconda
  - conda-forge
  - https://repo.anaconda.com/pkgs/main
  - https://repo.anaconda.com/pkgs/r
dependencies:
  - _openmp_mutex=4.5=20_gnu
  - alabaster=1.0.0=pyhd8ed1ab_1
  - babel=2.18.0=pyhcf101f3_0
  - backports.zstd=1.3.0=py312h90b7ffd_0
  - brotli=1.2.0=hed03a55_1
  - brotli-bin=1.2.0=hb03c661_1
  - brotli-python=1.2.0=py312hdb49522_1
  - bzip2=1.0.8=hda65f42_8
  - c-ares=1.34.6=hb03c661_0
  - ca-certificates=2026.1.4=hbd8a1cb_0
  - certifi=2026.1.4=pyhd8ed1ab_0
  - charset-normalizer=3.4.4=pyhd8ed1ab_0
  - colorama=0.4.6=pyhd8ed1ab_1
  - contourpy=1.3.3=py312h0a2e395_4
  - cycler=0.12.1=pyhcf101f3_2
  - deeptools=3.5.6=pyhdfd78af_0
  - deeptoolsintervals=0.1.9=py312ha9c1134_11
  - docutils=0.22.4=pyhd8ed1ab_0
  - fonttools=4.61.1=py312h8a5da7c_0
  - freetype=2.14.1=ha770c72_0
  - h2=4.3.0=pyhcf101f3_0
  - hpack=4.1.0=pyhd8ed1ab_0
  - hyperframe=6.1.0=pyhd8ed1ab_0
  - icu=78.2=h33c6efd_0
  - idna=3.11=pyhd8ed1ab_0
  - imagesize=1.4.1=pyhd8ed1ab_0
  - importlib-metadata=8.7.0=pyhe01879c_1
  - jinja2=3.1.6=pyhcf101f3_1
  - keyutils=1.6.3=hb9d3cd8_0
  - kiwisolver=1.4.9=py312h0a2e395_2
  - krb5=1.22.2=ha1258a1_0
  - lcms2=2.18=h0c24ade_0
  - ld_impl_linux-64=2.45.1=default_hbd61a6d_101
  - lerc=4.0.0=h0aef613_1
  - libblas=3.11.0=5_h4a7cf45_openblas
  - libbrotlicommon=1.2.0=hb03c661_1
  - libbrotlidec=1.2.0=hb03c661_1
  - libbrotlienc=1.2.0=hb03c661_1
  - libcblas=3.11.0=5_h0358290_openblas
  - libcurl=8.18.0=hcf29cc6_1
  - libdeflate=1.25=h17f619e_0
  - libedit=3.1.20250104=pl5321h7949ede_0
  - libev=4.33=hd590300_2
  - libexpat=2.7.3=hecca717_0
  - libffi=3.5.2=h3435931_0
  - libfreetype=2.14.1=ha770c72_0
  - libfreetype6=2.14.1=h73754d4_0
  - libgcc=15.2.0=he0feb66_17
  - libgcc-ng=15.2.0=h69a702a_17
  - libgfortran=15.2.0=h69a702a_17
  - libgfortran5=15.2.0=h68bc16d_17
  - libgomp=15.2.0=he0feb66_17
  - libjpeg-turbo=3.1.2=hb03c661_0
  - liblapack=3.11.0=5_h47877c9_openblas
  - liblzma=5.8.2=hb03c661_0
  - libnghttp2=1.67.0=had1ee68_0
  - libnsl=2.0.1=hb9d3cd8_1
  - libopenblas=0.3.30=pthreads_h94d23a6_4
  - libpng=1.6.55=h421ea60_0
  - libsqlite=3.51.2=hf4e2dac_0
  - libssh2=1.11.1=hcf80075_0
  - libstdcxx=15.2.0=h934c35e_17
  - libstdcxx-ng=15.2.0=hdf11a46_17
  - libtiff=4.7.1=h9d88235_1
  - libuuid=2.41.3=h5347b49_0
  - libwebp-base=1.6.0=hd42ef1d_0
  - libxcb=1.17.0=h8a09558_0
  - libxcrypt=4.4.36=hd590300_1
  - libzlib=1.3.1=hb9d3cd8_2
  - markupsafe=3.0.3=py312h8a5da7c_0
  - matplotlib-base=3.10.8=py312he3d6523_0
  - munkres=1.1.4=pyhd8ed1ab_1
  - narwhals=2.16.0=pyhcf101f3_0
  - ncurses=6.5=h2d0b736_3
  - numpy=2.4.2=py312h33ff503_1
  - numpydoc=1.10.0=pyhcf101f3_0
  - openjpeg=2.5.4=h55fea9a_0
  - openssl=3.6.1=h35e630c_1
  - packaging=26.0=pyhcf101f3_0
  - pillow=12.1.1=py312h50c33e8_0
  - pip=26.0.1=pyh8b19718_0
  - plotly=6.5.2=pyhd8ed1ab_0
  - pthread-stubs=0.4=hb9d3cd8_1002
  - py2bit=0.3.3=py312h0fa9677_1
  - pybigwig=0.3.25=py312h83fe50f_0
  - pygments=2.19.2=pyhd8ed1ab_0
  - pyparsing=3.3.2=pyhcf101f3_0
  - pysam=0.23.3=py312h8f9e533_2
  - pysocks=1.7.1=pyha55dd90_7
  - python=3.12.12=hd63d673_2_cpython
  - python-dateutil=2.9.0.post0=pyhe01879c_2
  - python_abi=3.12=8_cp312
  - pytz=2025.2=pyhd8ed1ab_0
  - qhull=2020.2=h434a139_5
  - readline=8.3=h853b02a_0
  - requests=2.32.5=pyhcf101f3_1
  - roman-numerals=4.1.0=pyhd8ed1ab_0
  - scipy=1.17.0=py312h54fa4ab_1
  - setuptools=82.0.0=pyh332efcf_0
  - six=1.17.0=pyhe01879c_1
  - snowballstemmer=3.0.1=pyhd8ed1ab_0
  - sphinx=9.1.0=pyhd8ed1ab_0
  - sphinxcontrib-applehelp=2.0.0=pyhd8ed1ab_1
  - sphinxcontrib-devhelp=2.0.0=pyhd8ed1ab_1
  - sphinxcontrib-htmlhelp=2.1.0=pyhd8ed1ab_1
  - sphinxcontrib-jsmath=1.0.1=pyhd8ed1ab_1
  - sphinxcontrib-qthelp=2.0.0=pyhd8ed1ab_1
  - sphinxcontrib-serializinghtml=1.1.10=pyhd8ed1ab_1
  - tk=8.6.13=noxft_h366c992_103
  - tomli=2.4.0=pyhcf101f3_0
  - tzdata=2025c=hc9c84f9_1
  - unicodedata2=17.0.1=py312h4c3975b_0
  - urllib3=2.6.3=pyhd8ed1ab_0
  - wheel=0.46.3=pyhd8ed1ab_0
  - xorg-libxau=1.0.12=hb03c661_1
  - xorg-libxdmcp=1.1.5=hb03c661_1
  - zipp=3.23.0=pyhcf101f3_1
  - zlib=1.3.1=hb9d3cd8_2
  - zlib-ng=2.3.3=hceb46e0_1
  - zstd=1.5.7=hb78ec9c_6
prefix: /work/pa-hilario/anaconda3/envs/deeptools
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

name: merge_peaks
channels:
  - conda-forge
  - bioconda
  - defaults
  - https://repo.anaconda.com/pkgs/main
  - https://repo.anaconda.com/pkgs/r
dependencies:
  - _openmp_mutex=4.5=20_gnu
  - _r-mutex=1.0.1=anacondar_1
  - bedtools=2.31.1=h13024bc_3
  - binutils_impl_linux-64=2.45.1=default_hfdba357_101
  - blas=2.305=openblas
  - blas-devel=3.11.0=5_h1ea3ea9_openblas
  - brotli=1.2.0=hed03a55_1
  - brotli-bin=1.2.0=hb03c661_1
  - bwidget=1.10.1=ha770c72_1
  - bzip2=1.0.8=hda65f42_9
  - c-ares=1.34.6=hb03c661_0
  - ca-certificates=2026.2.25=hbd8a1cb_0
  - cairo=1.18.4=h3394656_0
  - certifi=2025.8.3=pyhd8ed1ab_0
  - curl=8.19.0=hcf29cc6_0
  - cycler=0.12.1=pyhd8ed1ab_1
  - font-ttf-dejavu-sans-mono=2.37=hab24e00_0
  - font-ttf-inconsolata=3.000=h77eed37_0
  - font-ttf-source-code-pro=2.038=h77eed37_0
  - font-ttf-ubuntu=0.83=h77eed37_3
  - fontconfig=2.17.1=h27c8c51_0
  - fonts-conda-ecosystem=1=0
  - fonts-conda-forge=1=hc364b38_1
  - fonttools=4.59.1=py39heb7d2ae_0
  - freetype=2.14.2=ha770c72_0
  - fribidi=1.0.16=hb03c661_0
  - gcc_impl_linux-64=15.2.0=he420e7e_18
  - gfortran_impl_linux-64=15.2.0=h281d09f_18
  - graphite2=1.3.14=hecca717_2
  - gsl=2.7=he838d99_0
  - gxx_impl_linux-64=15.2.0=hda75c37_18
  - harfbuzz=11.2.1=h3beb420_0
  - htslib=1.23=h566b1c6_0
  - icu=75.1=he02047a_0
  - idr=2.0.4.2=py39h031d066_12
  - kernel-headers_linux-64=5.14.0=he073ed8_3
  - keyutils=1.6.3=hb9d3cd8_0
  - kiwisolver=1.4.7=py39h74842e3_0
  - krb5=1.22.2=ha1258a1_0
  - lcms2=2.17=h717163a_0
  - ld_impl_linux-64=2.45.1=default_hbd61a6d_101
  - lerc=4.1.0=hdb68285_0
  - libblas=3.11.0=5_h4a7cf45_openblas
  - libbrotlicommon=1.2.0=hb03c661_1
  - libbrotlidec=1.2.0=hb03c661_1
  - libbrotlienc=1.2.0=hb03c661_1
  - libcblas=3.11.0=5_h0358290_openblas
  - libcurl=8.19.0=hcf29cc6_0
  - libdeflate=1.22=hb9d3cd8_0
  - libedit=3.1.20250104=pl5321h7949ede_0
  - libev=4.33=hd590300_2
  - libexpat=2.7.4=hecca717_0
  - libffi=3.4.6=h2dba641_1
  - libfreetype=2.14.2=ha770c72_0
  - libfreetype6=2.14.2=h73754d4_0
  - libgcc=15.2.0=he0feb66_18
  - libgcc-devel_linux-64=15.2.0=hcc6f6b0_118
  - libgcc-ng=15.2.0=h69a702a_18
  - libgfortran=15.2.0=h69a702a_18
  - libgfortran-ng=15.2.0=h69a702a_18
  - libgfortran5=15.2.0=h68bc16d_18
  - libglib=2.84.1=h2ff4ddf_0
  - libgomp=15.2.0=he0feb66_18
  - libiconv=1.18=h3b78370_2
  - libjpeg-turbo=3.1.2=hb03c661_0
  - liblapack=3.11.0=5_h47877c9_openblas
  - liblapacke=3.11.0=5_h6ae95b6_openblas
  - liblzma=5.8.2=hb03c661_0
  - libnghttp2=1.67.0=had1ee68_0
  - libnsl=2.0.1=hb9d3cd8_1
  - libopenblas=0.3.30=pthreads_h94d23a6_4
  - libpng=1.6.55=h421ea60_0
  - libsanitizer=15.2.0=h90f66d4_18
  - libsqlite=3.52.0=h0c1763c_0
  - libssh2=1.11.1=hcf80075_0
  - libstdcxx=15.2.0=h934c35e_18
  - libstdcxx-devel_linux-64=15.2.0=hd446a21_118
  - libstdcxx-ng=15.2.0=hdf11a46_18
  - libtiff=4.7.0=hc4654cb_2
  - libuuid=2.41.3=h5347b49_0
  - libwebp-base=1.6.0=hd42ef1d_0
  - libxcb=1.17.0=h8a09558_0
  - libxcrypt=4.4.36=hd590300_1
  - libzlib=1.3.1=hb9d3cd8_2
  - make=4.4.1=hb9d3cd8_2
  - matplotlib-base=3.5.3=py39h19d6b11_2
  - munkres=1.0.7=py_1
  - ncurses=6.5=h2d0b736_3
  - numpy=1.19.5=py39hd249d9e_3
  - openblas=0.3.30=pthreads_h6ec200e_4
  - openjpeg=2.5.3=h55fea9a_1
  - openssl=3.6.1=h35e630c_1
  - packaging=26.0=pyhcf101f3_0
  - pango=1.56.3=h9ac818e_1
  - pcre2=10.44=hc749103_2
  - perl=5.26.2=h36c2ea0_1008
  - perl-app-cpanminus=1.7044=pl526_1
  - pillow=11.3.0=py39h15c0740_0
  - pip=25.2=pyh8b19718_0
  - pixman=0.46.4=h54a6638_1
  - pthread-stubs=0.4=hb9d3cd8_1002
  - pyparsing=3.2.3=pyhe01879c_2
  - python=3.9.23=hc30ae73_0_cpython
  - python-dateutil=2.9.0.post0=pyhe01879c_2
  - python_abi=3.9=8_cp39
  - r-base=4.4.2=hbab086e_1
  - readline=8.3=h853b02a_0
  - samtools=1.23=h96c455f_0
  - scipy=1.9.3=py39h32ae08f_2
  - sed=4.9=h6688a6e_0
  - setuptools=80.9.0=pyhff2d567_0
  - six=1.17.0=pyhe01879c_1
  - sysroot_linux-64=2.34=h087de78_3
  - tk=8.6.13=noxft_h366c992_103
  - tktable=2.10=h8d826fa_7
  - tzdata=2025c=hc9c84f9_1
  - unicodedata2=16.0.0=py39h8cd3c5a_0
  - wheel=0.45.1=pyhd8ed1ab_1
  - xorg-libice=1.1.2=hb9d3cd8_0
  - xorg-libsm=1.2.6=he73a12e_0
  - xorg-libx11=1.8.13=he1eb515_0
  - xorg-libxau=1.0.12=hb03c661_1
  - xorg-libxdmcp=1.1.5=hb03c661_1
  - xorg-libxext=1.3.7=hb03c661_0
  - xorg-libxrender=0.9.12=hb9d3cd8_0
  - xorg-libxt=1.3.1=hb9d3cd8_0
  - zstd=1.5.7=hb78ec9c_6
  - pip:
      - argcomplete==3.6.3
      - cachecontrol==0.14.3
      - charset-normalizer==3.4.6
      - coloredlogs==15.0.1
      - cwl-upgrader==1.2.12
      - cwl-utils==0.40
      - cwlref-runner==1.0
      - cwltool==3.1.20251031082601
      - filelock==3.19.1
      - humanfriendly==10.0
      - idna==3.11
      - isodate==0.7.2
      - lxml==6.0.2
      - markdown-it-py==3.0.0
      - mdurl==0.1.2
      - mistune==3.1.4
      - msgpack==1.1.2
      - mypy-extensions==1.1.0
      - networkx==3.2.1
      - prov==1.5.1
      - psutil==7.2.2
      - pydot==4.0.1
      - pygments==2.19.2
      - rdflib==7.4.0
      - requests==2.32.5
      - rich==14.3.3
      - rich-argparse==1.7.2
      - ruamel-yaml==0.18.17
      - ruamel-yaml-clib==0.2.15
      - schema-salad==8.9.20251102115403
      - spython==0.3.14
      - typing-extensions==4.15.0
      - urllib3==2.6.3
prefix: /work/pa-hilario/anaconda3/envs/merge_peaks

name: ucsc_tools
channels:
  - conda-forge
  - bioconda
  - defaults
  - https://repo.anaconda.com/pkgs/main
  - https://repo.anaconda.com/pkgs/r
dependencies:
  - _openmp_mutex=4.5=20_gnu
  - bedtools=2.31.1=h13024bc_3
  - bzip2=1.0.8=hda65f42_9
  - ca-certificates=2026.2.25=hbd8a1cb_0
  - libgcc=15.2.0=he0feb66_18
  - libgcc-ng=15.2.0=h69a702a_18
  - libgomp=15.2.0=he0feb66_18
  - libiconv=1.18=h3b78370_2
  - liblzma=5.8.2=hb03c661_0
  - libopenssl-static=3.6.1=hb03c661_1
  - libpng=1.6.55=h421ea60_0
  - libstdcxx=15.2.0=h934c35e_18
  - libstdcxx-ng=15.2.0=hdf11a46_18
  - libuuid=2.41.3=h5347b49_0
  - libzlib=1.3.1=hb9d3cd8_2
  - mysql-connector-c=6.1.11=h659d440_1008
  - openssl=3.6.1=h35e630c_1
  - ucsc-bedgraphtobigwig=482=hdc0a859_0
  - ucsc-bigwiginfo=482=h0b57e2e_0
  - ucsc-bigwigsummary=482=h0b57e2e_0
  - ucsc-bigwigtobedgraph=482=h0b57e2e_0
prefix: /work/pa-hilario/anaconda3/envs/ucsc_tools

Perl installations:
  Mon Mar 16 20:19:17 2026: "Module" Number::Format
    *   "installed into:
        /work/pa-hilario/anaconda3/envs/merge_peaks/lib/site_perl/5.26.2"

    *   "LINKTYPE: dynamic"

    *   "VERSION: 1.76"

    *   "EXE_FILES: "

  Mon Mar 16 20:19:21 2026: "Module" Statistics::Basic
    *   "installed into:
        /work/pa-hilario/anaconda3/envs/merge_peaks/lib/site_perl/5.26.2"

    *   "LINKTYPE: dynamic"

    *   "VERSION: 1.6611"

    *   "EXE_FILES: "

  Mon Mar 16 20:19:32 2026: "Module" Statistics::Distributions
    *   "installed into:
        /work/pa-hilario/anaconda3/envs/merge_peaks/lib/site_perl/5.26.2"

    *   "LINKTYPE: dynamic"

    *   "VERSION: 1.02"

    *   "EXE_FILES: "

  Mon Mar 16 20:19:52 2026: "Module" Regexp::Common
    *   "installed into:
        /work/pa-hilario/anaconda3/envs/merge_peaks/lib/site_perl/5.26.2"

    *   "LINKTYPE: dynamic"

    *   "VERSION: 2024080801"

    *   "EXE_FILES: "

  Mon Mar 16 20:20:11 2026: "Module" IPC::Run
    *   "installed into:
        /work/pa-hilario/anaconda3/envs/merge_peaks/lib/site_perl/5.26.2"

    *   "LINKTYPE: dynamic"

    *   "VERSION: 20250809.0"

    *   "EXE_FILES: "

  Mon Mar 16 20:21:53 2026: "Module" Statistics::R
    *   "installed into:
        /work/pa-hilario/anaconda3/envs/merge_peaks/lib/site_perl/5.26.2"

    *   "LINKTYPE: dynamic"

    *   "VERSION: 0.34"

    *   "EXE_FILES: "
```

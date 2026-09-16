# /work/pa-hilario/clipseq
### cc25dev0.naist.jp (Date: 20260508)

20251024:
dFP-CLIP libraries (SNRNP70 OSM,CTL) (n=1) with pseudoreplicates to complete Skipper pipeline.

20260218:
ddCLIP libraries (SNRNP70 OSM,CTL,OXI) (n=2) with Skipper read processing but incomplete pipeline due to lack of any window passing filters.
ddCLIP libraries (SNRNP70 OSM,CTL,OXI) (n=2) processed with Skipper (20260218) for DEWSeq analysis (done on 20260227).

20260326:
ddCLIP libraries (SNRNP70 OSM,CTL,OXI) (n=2) processed with custom pipeline: mRNA annotation **(!) without unique transcripts per gene (!)** (my custom), Cutadapt (eCLIP/seCLIP), 1 site per multimapping read (Skipper), UMI deduplication by directional adjacency (Skipper), Peak-calling (CLIPper), Consensus peak regions (my custom), RPM normalization per consensus peak (my custom), read QC (Skipper), metagene analysis (my custom).
ddCLIP libraries (SNRNP70 OSM,CTL) (n=2) processed allowing multimapping (20260326) for repeat analysis using RepeatMasker (my custom) (done on 20260423).

20260421:
ddCLIP libraries (SNRNP70 OSM,CTL) (n=2) processed with custom pipeline: mRNA annotation (RBPBench, extract-transcript-regions, my custom), Cutadapt (eCLIP/seCLIP), unique mapping (eCLIP/seCLIP), UMI deduplication by directional adjacency (Skipper), RPM normalization (my custom), read QC (Skipper).
ddCLIP libraries (SNRNP70 OSM,CTL) (n=2) processed for unique mapping (20260421) for metagene analysis using deepTools (my custom) (done on 20260423).

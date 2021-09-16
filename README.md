# ERV cancer enhancers

### Scripts and files used in:

Ivancevic A, Simpson D, Chuong E (2021) "Endogenous retroviruses mediate transcriptional rewiring in response to oncogenic signaling in colorectal cancer" JOURNAL, ISSUE: PAGES. **[add biorxiv link]**

### GEO links:

### UCSC sessions:
Add a few of my UCSC sessions

### Programs used:
- GIGGLE v0.6.3 (https://github.com/ryanlayer/giggle)
- bedtools v2.28.0 (http://bedtools.readthedocs.io/en/latest/)
- deepTools v3.0.1 (https://deeptools.readthedocs.io/en/develop/index.html)
- samtools v1.10 (http://www.htslib.org/)
- BBDuk/BBMap v38.05 (https://jgi.doe.gov/data-and-tools/bbtools/)
- FastQC v0.11.8 (https://github.com/s-andrews/FastQC)
- MultiQC v1.7 (https://github.com/ewels/MultiQC)
- hisat2 v2.1.0 (https://github.com/DaehwanKimLab/hisat2)
- subread/featureCounts v1.6.2 (http://subread.sourceforge.net/)
- TEtranscripts v2.1.4 (https://github.com/mhammell-laboratory/TEtranscripts)
- MACS2 v2.1.1 (https://pypi.org/project/MACS2/)
- BWA v0.7.15 (https://github.com/lh3/bwa)
- DESeq2 v1.32.0 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- fastq-dump v2.10.5 from SRA Toolkit (https://hpc.nih.gov/apps/sratoolkit.html)
- MEME Suite v5.4.1 (https://meme-suite.org/meme/)
- vsearch v2.14.1 (https://github.com/torognes/vsearch)
- MUSCLE v3.8.1551 (https://www.drive5.com/muscle/)
- Jalview v2.11.1.4 (https://www.jalview.org/)
- singularity v3.1.1 (https://github.com/hpcng/singularity)

### Public databases:
- The Cancer Genome Atlas (TCGA), via the Genomic Data Commons (https://gdc.cancer.gov/)
- CEMT epigenomes (http://www.epigenomes.ca/data-release/)
- Cistrome (http://cistrome.org/)
- Roadmap Epigenomics Project (http://www.roadmapepigenomics.org/)
- ENCODE (https://www.encodeproject.org/)
- Dfam (https://dfam.org/home)
- Human Endogenous Retrovirus Database (https://herv.img.cas.cz/)
- UCSC (https://genome.ucsc.edu/)
- GENCODE Release 34 (https://www.gencodegenes.org/human/)
- JASPAR (http://jaspar.genereg.net/)
- gnomAD (https://gnomad.broadinstitute.org/)

## Pan-cancer epigenomic analysis of TE activity

### 1. Identify cancer-specific regulatory regions

Predicted regulatory regions in healthy adult tissues were defined using Roadmap categories 1_TssA, 6_EnhG & 7_Enh (see [Roadmap_healthy_adult_tissues_list.txt](1_pancancer_epigenomics/Roadmap_healthy_adult_tissues_list.txt)). Predicted regulatory regions for 21 different cancer types were obtained from TCGA ATACseq peaks (see [TCGA_cancers_list.txt](1_pancancer_epigenomics/TCGA_cancers_list.txt)). Cancer-specific regulatory regions were identified by subtracting "healthy" regulatory regions from each cancer peak set.

### 2. Test for family-level TE enrichment
GIGGLE was used to create a database of all TE families in the hg38 human genome, based on Dfam annotation. Cancer-specific regulatory regions were then searched against the repeat database (see [find_enriched_TEs.sh](1_pancancer_epigenomics/find_enriched_TEs.sh)). Results were ranked from highest to lowest enrichment score, and filtered based on odds ratio, score, and the number of overlaps ([SuppTable1_TCGA_giggle_results_top23TEs.tab](1_pancancer_epigenomics/SuppTable1_TCGA_giggle_results_top23TEs.tab)). Cancer-TE associations were visualized as bubble and volcano plots (e.g. [Fig1_bubbles.py](1_pancancer_epigenomics/Fig1_bubbles.py), [Fig1_volcano.R](1_pancancer_epigenomics/Fig1_volcano.R)). 

### 3. Element-level TE analysis

hg38 genome coordinates of LTR10 elements were obtained from Dfam. LTR10A and LTR10F elements were merged (2kb window) and used to generate signal heatmaps over TCGA tumor ATACseq (e.g. [Fig1_deeptools_atacseq.sbatch](1_pancancer_epigenomics/Fig1_deeptools_atacseq.sbatch)). Similar heatmaps were plotted using other LTR10 subfamilies (e.g. [deeptools_atacseq_allLTR10.sbatch](1_pancancer_epigenomics/deeptools_atacseq_allLTR10.sbatch)) and other TEs. 

We further assessed TE signal in public ChIPseq datasets from Cistrome, particularly looking at histone marks indicative of enhancer activity in cancer cell lines (e.g. HCT116 H3K27ac ChIPseq & H3K4me1 ChIPseq; see [Fig1_deeptools_hct116.sbatch](1_pancancer_epigenomics/Fig1_deeptools_hct116.sbatch)). 

## Regulatory activity of LTR10 in colorectal cancer

Fig2a: overall heatmap and deeptools signal heatmap of CEMT patient chipseq

2b: Roadmap category giggle scores

Supp table 2: roadmap giggle results

2c: barchart top cistrome tfs (activators)

2d: tf activators deeptools heatmap

2e: barchart cistrome tf repressors 

2f: tf repressors deeptools heatmap

2g: python script coolbox screenshot

2h: meme results and input files e.g. primary set H3K27ac-marked LTR10s vs control set non-H3K27ac-marked LTR10s (with shell script process of how generated those files). Also deeptools motif plot. And fimo script to get FOSL1 motif coordinates genome-wide. 

## Evolutionary history of LTR10 elements

2i: muscle/gblocks/fasttree scripts for LTR10 consensus tree

Fig 2 Supp: pca plot, alignment of all LTR10 sequences, etc. 

## Control of LTR10 activity by AP1/MAPK signaling

Started with RNAseq fastq files (Cobi_24hr x2, TNF_24hr x2, Untreated_24hr x2). 

**RNAseq files [add GEO link here]:** 
- Cobi_24_1_S54_R1.fastq.gz
- Cobi_24_1_S54_R2.fastq.gz
- Cobi_24_2_S55_R1.fastq.gz
- Cobi_24_2_S55_R2.fastq.gz
- TNF_24_1_S58_R1.fastq.gz
- TNF_24_1_S58_R2.fastq.gz
- TNF_24_2_S59_R1.fastq.gz
- TNF_24_2_S59_R2.fastq.gz
- UT_24_1_S52_R1.fastq.gz 
- UT_24_1_S52_R2.fastq.gz
- UT_24_2_S53_R1.fastq.gz
- UT_24_2_S53_R2.fastq.gz

**RNAseq Workflow:**
1) [bbduk_PE.sbatch](rnaseq/bbduk_PE.sbatch)
2) [fastqc.sbatch](rnaseq/fastqc.sbatch)
3) [multiqc.sbatch](rnaseq/multiqc.sbatch) 
4) [hisat2_PE.sbatch](rnaseq/hisat2_PE.sbatch)
5) [bam_to_bw.sbatch](rnaseq/bam_to_bw.sbatch)
6) [feature_counts.sbatch](rnaseq/feature_counts.sbatch)
7) [add_gene_names.sh](rnaseq/add_gene_names.sh), using [gencode34_geneid_genename.txt](rnaseq/gencode34_geneid_genename.txt)
8) [deseq2_genes.R](rnaseq/deseq2_genes.R)

**For TE-transcripts, re-aligned bams to allow multiple alignments per read:**
1) [hisat2_PE_k100.sbatch](rnaseq/hisat2_PE_k100.sbatch)
2) [te_transcripts.sbatch](rnaseq/te_transcripts.sbatch)
3) [extract_TEs.sh](rnaseq/extract_TEs.sh)
4) [deseq2_TEs.R](rnaseq/deseq2_TEs.R)

**CUT&RUN files [add GEO link here]:** 
- List them here
- H3K27ac CUT&RUN Untreated (2 reps, second rep was resequenced but I guess we just call it Rep 2?), Cobi-treated (x2) and TNF-treated (x2).
- Also have the FOSL1, JUN, H3K9me3, POL2 CUT&RUN, but not sure if we're including those in this paper?
- Only list the ones that we're going to add to GEO, e.g. maybe only the H3K27ac for now

**CUT&RUN Workflow:**

1) [setup_workspace.sbatch](cutnrun/setup_workspace.sbatch)
2) [bbduk_PE.sbatch](cutnrun/bbduk_PE.sbatch)
3) [fastqc.sbatch](cutnrun/fastqc.sbatch)
4) [multiqc.sbatch](cutnrun/multiqc.sbatch) 
5) [bwa_PE.sbatch](cutnrun/bwa_PE.sbatch)
6) [subset_by_fragment_size.sbatch](cutnrun/subset_by_fragment_size.sbatch) (for transcription factors only) 
7) [get_fragment_size.sbatch](cutnrun/get_fragment_size.sbatch)



[bam_to_bw.sbatch](cutnrun/bam_to_bw.sbatch)
peak calling without igg
[convert_bam_to_fragment_bdg.sbatch](cutnrun/convert_bam_to_fragment_bdg.sbatch)
[bgzip_and_giggle.sbatch](cutnrun/bgzip_and_giggle.sbatch)
[calculate_frip_score.sbatch](cutnrun/calculate_frip_score.sbatch)
[convert_macs2_bdg_to_bigwig.sbatch](cutnrun/convert_macs2_bdg_to_bigwig.sbatch)
[deeptools_heatmap_from_gencode_bed.sbatch](cutnrun/deeptools_heatmap_from_gencode_bed.sbatchh)
[deeptools_heatmap_from_narrowPeak.sbatch](cutnrun/deeptools_heatmap_from_narrowPeak.sbatch)
[get_fragment_size.sbatch](cutnrun/get_fragment_size.sbatch)
[merge_peak_files.sbatch](cutnrun/merge_peak_files.sbatch)



### 1. MAPK-treated HCT116 RNAseq, TE transcripts analysis

- Raw counts table of TEs for Cobimetinib vs Ctrl: [raw_counts_tetranscripts_cobi_24hr.tab](4_mapk_signaling/raw_counts_tetranscripts_cobi_24hr.tab)
- Raw counts table of TEs for TNF-alpha vs Ctrl: [raw_counts_tetranscripts_tnf_24hr.tab](4_mapk_signaling/raw_counts_tetranscripts_tnf_24hr.tab)

- Normalized counts table of TEs for Cobimetinib vs Ctrl: [normalized_counts_tetranscripts_cobi_24hr.tab](4_mapk_signaling/normalized_counts_tetranscripts_cobi_24hr.tab)
- Normalized counts table of TEs for TNF-alpha vs Ctrl: [normalized_counts_tetranscripts_tnf_24hr.tab](4_mapk_signaling/normalized_counts_tetranscripts_tnf_24hr.tab)

- DEseq2 results table of TEs for Cobimetinib vs Ctrl: [deseq2_results_tetranscripts_cobi_vs_ctrl.tab](4_mapk_signaling/deseq2_results_tetranscripts_cobi_vs_ctrl.tab)
- DEseq2 results table of TEs for TNF-alpha vs Ctrl: [deseq2_results_tetranscripts_tnf_vs_ctrl.tab](4_mapk_signaling/deseq2_results_tetranscripts_tnf_vs_ctrl.tab)

- MA plot of TEs for Cobimetinib vs Ctrl: [Fig3_MAplot_TEs_cobi_vs_ctrl.R](4_mapk_signaling/Fig3_MAplot_TEs_cobi_vs_ctrl.R), [Fig3_MAplot_TEs_cobi_vs_ctrl.pdf](4_mapk_signaling/Fig3_MAplot_TEs_cobi_vs_ctrl.pdf)
- MA plot of TEs for TNF-alpha vs Ctrl: [Fig3_MAplot_TEs_tnf_vs_ctrl.R](4_mapk_signaling/Fig3_MAplot_TEs_tnf_vs_ctrl.R), [Fig3_MAplot_TEs_tnf_vs_ctrl.pdf](4_mapk_signaling/Fig3_MAplot_TEs_tnf_vs_ctrl.pdf)

### 2. MAPK-treated HCT116 RNAseq, gene analysis

- Raw counts table of genes for all samples: [raw_counts_genes_24hr.tab](4_mapk_signaling/raw_counts_genes_24hr.tab)

- Normalized counts table of genes for all samples: [normalized_counts_genes_24hr.tab](4_mapk_signaling/normalized_counts_genes_24hr.tab)

- DEseq2 results table of genes for Cobimetinib vs Ctrl: [deseq2_results_genes_cobi_vs_ctrl.tab](4_mapk_signaling/deseq2_results_genes_cobi_vs_ctrl.tab)
- DEseq2 results table of genes for TNF-alpha vs Ctrl: [deseq2_results_genes_tnf_vs_ctrl.tab](4_mapk_signaling/deseq2_results_genes_tnf_vs_ctrl.tab)

- MA plot of genes for Cobimetinib vs Ctrl: [MAplot_genes_cobi_vs_ctrl.R](4_mapk_signaling/MAplot_genes_cobi_vs_ctrl.R), [MAplot_genes_cobi_vs_ctrl.pdf](4_mapk_signaling/MAplot_genes_cobi_vs_ctrl.pdf)
- MA plot of genes for TNF-alpha vs Ctrl: [MAplot_genes_tnf_vs_ctrl.R](4_mapk_signaling/MAplot_genes_tnf_vs_ctrl.R), [MAplot_genes_tnf_vs_ctrl.pdf](4_mapk_signaling/MAplot_genes_tnf_vs_ctrl.pdf)

- Merged table of Cobimetinib and TNF-alpha results: [merge_genes_cobi_and_tnf.R](4_mapk_signaling/merge_genes_cobi_and_tnf.R), [merge_genes_cobi_and_tnf.tab](4_mapk_signaling/merge_genes_cobi_and_tnf.tab)

### 3. MAPK-treated HCT116 H3K27ac CUT&RUN

E.g. MA plots and giggle enrichment volcano showing LTR10A/F

Analysis scripts/workflow for Cut&Run.
Then the Cut&Run scripts. 
3d: giggle volcano plot of treated C&R 

### 4. Enhancer-gene predictions (combining MAPK-treated HCT116 RNAseq, HCT116 H3K27ac CUT&RUN and Activity-by-Contact model)

**Activity by Contact model predictions** (https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/master/README.md).

Followed the GitHub instructions to make enhancer-gene predictions in HCT116 cells, using as input:
1) Publicly available HCT116 ATACseq, GEO accession GSM3593802 (https://www.ncbi.nlm.nih.gov/sra/?term=SRR8544480)
2) In-house HCT116 H3K27ac CUT&RUN, GEO accession DUMMY (LINK)
3) Average HiC file provided by the ABC model

**Workflow:**
1) Align ATACseq to hg19 genome: [bwa_SE_atacseq_hg19.sbatch](4_mapk_signaling/abc_model/bwa_SE_atacseq_hg19.sbatch)
2) Index bam file: [bam_index.sbatch](4_mapk_signaling/abc_model/bam_index.sbatch) 
3) Call peaks using the ATACseq file: [call_peaks_atacseq.sbatch:](4_mapk_signaling/abc_model/call_peaks_atacseq.sbatch) 
4) Sort the narrowPeak file: [sort_narrowpeak.sh](4_mapk_signaling/abc_model/sort_narrowpeak.sh) 
5) Call candidate regions: [call_candidate_regions.sbatch](4_mapk_signaling/abc_model/call_candidate_regions.sbatch) 
6) Align H3K27ac CUT&RUN to hg19 genome: [bwa_PE_cutnrun_hg19.sbatch](4_mapk_signaling/abc_model/bwa_PE_cutnrun_hg19.sbatch)
7) Find ABC neighbourhoods by combining ATACseq and H3K27ac CUT&RUN: [find_neighborhoods.sbatch](4_mapk_signaling/abc_model/find_neighborhoods.sbatch)
8) Predict enhancer regions by incorporating HiC: [predict_enhancers.sbatch](4_mapk_signaling/abc_model/predict_enhancers.sbatch)
9) Liftover to hg38 genome coordinates: [liftover.sh](4_mapk_signaling/abc_model/liftover.sh)
10) Intersect predicted enhancers with H3K27ac-marked LTR10A/F elements: [intersect_and_filter.sh](4_mapk_signaling/abc_model/intersect_and_filter.sh)

**Output:** [EnhancerPredictionsAllPutativehg38_intersect_LTR10AF_scoreOver0.001.bed](4_mapk_signaling/abc_model/EnhancerPredictionsAllPutativehg38_intersect_LTR10AF_scoreOver0.001.bed)

Intersect MAPK genes with ABC model: shell script showing process to get to final list of gene/LTR10 enhancer. Table of final list.

3e: Redo Log2FC plot to include final gene candidates as bigger bubbles: [Fig3_Log2FCplot_genes_treatments.R](4_mapk_signaling/Fig3_Log2FCplot_genes_treatments.R)

3g: R script pheatmap of genes with treatments. Also R script pheatmap of H3K27ac Cut&Run LTR10 elements

## CRISPR silencing & deletion of LTR10 elements

Analysis scripts for RNAseq, like above. 

Python coolbox scripts and R scripts for deseq2/genome distance plots. Need to redo with lfcShrink. Also python scripts for bargraphs. And MA plots for e.g. ATG12 LTR10F. 

## Human variation at AP1 sites within LTR10 elements

Python scripts for AP1 count histogram/correlation plots. 

Also python/jupyter scripts for gnomad deletions in human popn. 

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
- FIMO v5.1.0 (https://meme-suite.org/meme/doc/fimo.html?man_type=web)
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
- Dfam 2.0 (https://dfam.org/home)
- Human Endogenous Retrovirus Database (https://herv.img.cas.cz/)
- UCSC (https://genome.ucsc.edu/)
- GENCODE Release 34 (https://www.gencodegenes.org/human/)
- JASPAR (http://jaspar.genereg.net/)
- gnomAD (https://gnomad.broadinstitute.org/)

## Pan-cancer epigenomic analysis of TE activity

### 1. Identify cancer-specific regulatory regions

Predicted regulatory regions in healthy adult tissues were obtained from Roadmap (https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all_hg38lift.mnemonics.bedFiles.tgz), using categories 1_TssA, 6_EnhG & 7_Enh (see [Roadmap_healthy_adult_tissues_list.txt](1_pancancer_epigenomics/Roadmap_healthy_adult_tissues_list.txt)). Predicted regulatory regions for 21 different cancer types were obtained from TCGA ATACseq peaks (see [TCGA_cancers_list.txt](1_pancancer_epigenomics/TCGA_cancers_list.txt)). Cancer-specific regulatory regions were identified by subtracting "healthy" regulatory regions from each cancer peak set (see [find_enriched_TEs.sh](1_pancancer_epigenomics/find_enriched_TEs.sh)).

### 2. Test for family-level TE enrichment
GIGGLE was used to create a database of all TE families in the hg38 human genome, based on Dfam annotation. Cancer-specific regulatory regions were then searched against the repeat database ([find_enriched_TEs.sh](1_pancancer_epigenomics/find_enriched_TEs.sh)). Results were ranked from highest to lowest enrichment score, and filtered based on odds ratio, score, and the number of overlaps ([SuppTable1_TCGA_giggle_results_top23TEs.tab](1_pancancer_epigenomics/SuppTable1_TCGA_giggle_results_top23TEs.tab)). Cancer-TE associations were visualized as bubble and volcano plots (e.g. [Fig1_bubbles.py](1_pancancer_epigenomics/Fig1_bubbles.py), [Fig1_volcano.R](1_pancancer_epigenomics/Fig1_volcano.R)). 

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
- Cobi_24_1_S54_R1.fastq.gz, Cobi_24_1_S54_R2.fastq.gz
- Cobi_24_2_S55_R1.fastq.gz, Cobi_24_2_S55_R2.fastq.gz
- TNF_24_1_S58_R1.fastq.gz, TNF_24_1_S58_R2.fastq.gz
- TNF_24_2_S59_R1.fastq.gz, TNF_24_2_S59_R2.fastq.gz
- UT_24_1_S52_R1.fastq.gz, UT_24_1_S52_R2.fastq.gz
- UT_24_2_S53_R1.fastq.gz, UT_24_2_S53_R2.fastq.gz

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
Temporarily put them in /scratch/Users/ativ2716/data/2_files_for_GEO_upload/MAPK_treated_H3K27ac_CutnRun

Rep1: 210205_A00405_0345_BHT2K3DSXY

Rep2: 210409_A00405_0379_AH2YHCDSX2

Rep2 resequenced (untreated H3K27ac only): 210723_A00405_0433_AHHYFCDSX2

- Untreated rep 1: Untreated_H3K27ac_rep1_R1.fastq.gz, Untreated_H3K27ac_rep1_R2.fastq.gz
- Untreated rep 2: Untreated_H3K27ac_rep2_R1.fastq.gz, Untreated_H3K27ac_rep2_R2.fastq.gz
- Untreated rep 2, resequenced (merged with above at the bam stage): Untreated_H3K27ac_rep2_reseq_R1.fastq.gz, Untreated_H3K27ac_rep2_reseq_R2.fastq.gz
- Cobi-treated rep1: Cobimetinib_H3K27ac_rep1_R1.fastq.gz, Cobimetinib_H3K27ac_rep1_R2.fastq.gz
- Cobi-treated rep2: Cobimetinib_H3K27ac_rep2_R1.fastq.gz, Cobimetinib_H3K27ac_rep2_R2.fastq.gz
- Tnf-treated rep1: TNFalpha_H3K27ac_rep1_R1.fastq.gz, TNFalpha_H3K27ac_rep1_R2.fastq.gz
- Tnf-treated rep2: TNFalpha_H3K27ac_rep2_R1.fastq.gz, TNFalpha_H3K27ac_rep2_R2.fastq.gz
- Untreated IgG control rep1: Untreated_IgG_rep1_R1.fastq.gz, Untreated_IgG_rep1_R2.fastq.gz
- Untreated IgG control rep2: Untreated_IgG_rep2_R1.fastq.gz, Untreated_IgG_rep2_R2.fastq.gz  
- Cobi-treated IgG control rep1: Cobimetinib_IgG_rep1_R1.fastq.gz, Cobimetinib_IgG_rep1_R2.fastq.gz
- Cobi-treated IgG control rep2: Cobimetinib_IgG_rep2_R1.fastq.gz, Cobimetinib_IgG_rep2_R2.fastq.gz 
- Tnf-treated IgG control rep1: TNFalpha_IgG_rep1_R1.fastq.gz, TNFalpha_IgG_rep1_R2.fastq.gz
- Tnf-treated IgG control rep2: TNFalpha_IgG_rep2_R1.fastq.gz, TNFalpha_IgG_rep2_R2.fastq.gz

**CUT&RUN Workflow:**

1) [setup_workspace.sbatch](cutnrun/setup_workspace.sbatch)
2) [bbduk_PE.sbatch](cutnrun/bbduk_PE.sbatch)
3) [fastqc.sbatch](cutnrun/fastqc.sbatch)
4) [multiqc.sbatch](cutnrun/multiqc.sbatch) 
5) [bwa_PE.sbatch](cutnrun/bwa_PE.sbatch)
6) [subset_by_fragment_size.sbatch](cutnrun/subset_by_fragment_size.sbatch) (for transcription factors only) 
7) [get_fragment_size.sbatch](cutnrun/get_fragment_size.sbatch)
8) [call_peaks_with_macs2_PEmode_noIgGControl.sbatch](cutnrun/call_peaks_with_macs2_PEmode_noIgGControl.sbatch)
10) [call_peaks_with_macs2_SEmode_noIgGControl.sbatch](cutnrun/call_peaks_with_macs2_SEmode_noIgGControl.sbatch)
11) [merge_peak_files.sbatch](cutnrun/merge_peak_files.sbatch)
12) [convert_macs2_bdg_to_bigwig.sbatch](cutnrun/convert_macs2_bdg_to_bigwig.sbatch)
13) [calculate_frip_score.sbatch](cutnrun/calculate_frip_score.sbatch)
14) [deeptools_heatmap_from_gencode_bed.sbatch](cutnrun/deeptools_heatmap_from_gencode_bed.sbatch)
15) [bgzip_and_giggle.sbatch](cutnrun/bgzip_and_giggle.sbatch)
16) [extract_top_peaks.sh](cutnrun/extract_top_peaks.sh)
17) [compute_bam_count_table.sbatch](cutnrun/compute_bam_count_table.sbatch)
18) [deseq2_CnR.R](cutnrun/deseq2_CnR.R)
19) [giggle_de_regions.sh](cutnrun/giggle_de_regions.sh)

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

Peak files:
- [UT_K27ac_CnR_rep1_peaks.bed](4_mapk_signaling/UT_K27ac_CnR_rep1_peaks.bed)
- [UT_K27ac_CnR_rep2_peaks.bed](4_mapk_signaling/UT_K27ac_CnR_rep2_peaks.bed)
- [Cobi_K27ac_CnR_rep1_peaks.bed](4_mapk_signaling/Cobi_K27ac_CnR_rep1_peaks.bed)
- [Cobi_K27ac_CnR_rep2_peaks.bed](4_mapk_signaling/Cobi_K27ac_CnR_rep2_peaks.bed)
- [TNF_K27ac_CnR_rep1_peaks.bed](4_mapk_signaling/TNF_K27ac_CnR_rep1_peaks.bed)
- [TNF_K27ac_CnR_rep2_peaks.bed](4_mapk_signaling/TNF_K27ac_CnR_rep2_peaks.bed)

The top 20k peaks were extracted from each sample (based on macs2 peak score) and merged: [ALL_K27ac_CnR_top20k_peaks.bed](4_mapk_signaling/ALL_K27ac_CnR_top20k_peaks.bed)

Deeptools: 
[deeptools_H3K27ac_CnR_over_LTR10AF.sbatch](4_mapk_signaling/deeptools_H3K27ac_CnR_over_LTR10AF.sbatch), [deeptools_H3K27ac_CnR_over_LTR10AF_CnR.pdf](4_mapk_signaling/deeptools_H3K27ac_CnR_over_LTR10AF_CnR.pdf)

- Raw counts table of H3K27ac-marked regions for all samples: [raw_counts_H3K27ac_CnR.tab](4_mapk_signaling/raw_counts_H3K27ac_CnR.tab)
- Normalized counts table of H3K27ac-marked regions for all samples: [normalized_counts_H3K27ac_CnR.tab](4_mapk_signaling/normalized_counts_H3K27ac_CnR.tab)
- DEseq2 results table of H3K27ac-marked regions for Cobimetinib vs Ctrl: [deseq2_results_H3K27ac_CnR_cobi_vs_ctrl.tab](4_mapk_signaling/deseq2_results_H3K27ac_CnR_cobi_vs_ctrl.tab)
- DEseq2 results table of H3K27ac-marked regions for TNF-alpha vs Ctrl: [deseq2_results_H3K27ac_CnR_tnf_vs_ctrl.tab](4_mapk_signaling/deseq2_results_H3K27ac_CnR_tnf_vs_ctrl.tab)
- MA plot of H3K27ac-marked regions for Cobimetinib vs Ctrl: [MAplot_H3K27ac_CnR_cobi_vs_ctrl.R](4_mapk_signaling/MAplot_H3K27ac_CnR_cobi_vs_ctrl.R), [MAplot_H3K27ac_CnR_cobi_vs_ctrl.pdf](4_mapk_signaling/MAplot_H3K27ac_CnR_cobi_vs_ctrl.pdf)
- MA plot of H3K27ac-marked regions for TNF-alpha vs Ctrl: [MAplot_H3K27ac_CnR_tnf_vs_ctrl.R](4_mapk_signaling/MAplot_H3K27ac_CnR_tnf_vs_ctrl.R), [MAplot_H3K27ac_CnR_tnf_vs_ctrl.pdf](4_mapk_signaling/MAplot_H3K27ac_CnR_tnf_vs_ctrl.pdf)

Giggle tables: 
- [DE_regions_H3K27ac_CnR_cobimetinib_giggle_results.tab](4_mapk_signaling/DE_regions_H3K27ac_CnR_cobimetinib_giggle_results.tab)
- [DE_regions_H3K27ac_CnR_tnfalpha_giggle_results.tab](4_mapk_signaling/DE_regions_H3K27ac_CnR_tnfalpha_giggle_results.tab)

Giggle plots:
- Volcano plot of giggle enrichment for H3K27ac-marked DE regions for Cobimetinib vs Ctrl: [giggle_volcano_H3K27ac_CnR_DE_cobi_vs_ctrl.R](4_mapk_signaling/giggle_volcano_H3K27ac_CnR_DE_cobi_vs_ctrl.R), [giggle_volcano_H3K27ac_CnR_DE_cobi_vs_ctrl.pdf](4_mapk_signaling/giggle_volcano_H3K27ac_CnR_DE_cobi_vs_ctrl.pdf)
- Volcano plot of giggle enrichment for H3K27ac-marked DE regions for TNF-alpha vs Ctrl: [giggle_volcano_H3K27ac_CnR_DE_tnf_vs_ctrl.R](4_mapk_signaling/giggle_volcano_H3K27ac_CnR_DE_tnf_vs_ctrl.R), [giggle_volcano_H3K27ac_CnR_DE_tnf_vs_ctrl.pdf](4_mapk_signaling/giggle_volcano_H3K27ac_CnR_DE_tnf_vs_ctrl.pdf)

### 4. Enhancer-gene predictions 

Combining MAPK-treated HCT116 RNAseq, HCT116 H3K27ac CUT&RUN and Activity-by-Contact model.

**4a) Activity by Contact model predictions** (https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/master/README.md).

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

**4b) Combine ABC model predictions with MAPK treatment results**

Finally, prioritize LTR10 predicted enhancers based on proximity to MAPK-regulated genes: [prioritize_candidates.sh](4_mapk_signaling/prioritize_candidates.sh), using [gencode34_genename_genecoords.txt](rnaseq/gencode34_genename_genecoords.txt)

**Output:** [LTR10_predicted_enhancer_gene_pairs.bed](4_mapk_signaling/LTR10_predicted_enhancer_gene_pairs.bed)

Plot the candidates! 

Gene heatmap across treatments:
- [74_genes_normalized_counts_rnaseq.tab](4_mapk_signaling/74_genes_normalized_counts_rnaseq.tab)
- [74_genes_rnaseq_heatmap.R](4_mapk_signaling/74_genes_rnaseq_heatmap.R)
- [74_genes_rnaseq_heatmap.pdf](4_mapk_signaling/74_genes_rnaseq_heatmap.pdf)

LTR10 heatmap across treatments:
- [57_LTR10_enhancers_normalized_counts_h3K27ac_cutnrun.tab](4_mapk_signaling/57_LTR10_enhancers_normalized_counts_h3K27ac_cutnrun.tab)
- [57_LTR10_enhancers_k27acCnR_heatmap.R](4_mapk_signaling/57_LTR10_enhancers_k27acCnR_heatmap.R)
- [57_LTR10_enhancers_k27acCnR_heatmap.pdf](4_mapk_signaling/57_LTR10_enhancers_k27acCnR_heatmap.pdf)

Log2FC plot with the final 74 gene candidates as bigger bubbles: 
- [merge_genes_cobi_and_tnf_annotated.tab](4_mapk_signaling/merge_genes_cobi_and_tnf_annotated.tab)
- [Log2FCplot_genes_treatments.R](4_mapk_signaling/Log2FCplot_genes_treatments.R)
- [Log2FCplot_genes_treatments.pdf](4_mapk_signaling/Log2FCplot_genes_treatments.pdf)


## CRISPR silencing & deletion of LTR10 elements

Just give one example e.g. ATG12 

1) Coolbox python script

2) DEseq2 analysis of CRISPRi candidates vs GFP control

Processed each sample using the same RNAseq workflow as above.

- Raw feature counts table of all CRISPRi candidates: [feature_counts_crispri.tab](5_crispr_results/feature_counts_crispri.tab)
- Normalized counts table of all CRISPRi candidates: [normalized_counts_crispri.tab](5_crispr_results/normalized_counts_crispri.tab)
- PCA of all CRISPRi and control samples: [deseq2_and_pca_crispri.R](5_crispr_results/deseq2_and_pca_crispri.R), [pca_crispri.pdf](5_crispr_results/pca_crispri.pdf)
- Count plot bargraphs (python scripts) showing norm counts e.g. for genes ATG12 and AP3S1. 
- DEseq2 results table for all comparisons:
1. ATG12_LTR10i vs Ctrl: [deseq2_results_genes_atg12ltr10i_vs_ctrl.tab](5_crispr_results/deseq2_results_genes_atg12ltr10i_vs_ctrl.tab)
2. XRCC4_LTR10i vs Ctrl: [deseq2_results_genes_xrcc4ltr10i_vs_ctrl.tab](5_crispr_results/deseq2_results_genes_xrcc4ltr10i_vs_ctrl.tab)
3. MEF2D_LTR10i vs Ctrl: [deseq2_results_genes_mef2dltr10i_vs_ctrl.tab](5_crispr_results/deseq2_results_genes_mef2dltr10i_vs_ctrl.tab)
4. FGF2_LTR10i vs Ctrl: [deseq2_results_genes_fgf2ltr10i_vs_ctrl.tab](5_crispr_results/deseq2_results_genes_fgf2ltr10i_vs_ctrl.tab)
5. MCPH1_LTR10i vs Ctrl: [deseq2_results_genes_mcph1ltr10i_vs_ctrl.tab](5_crispr_results/deseq2_results_genes_mcph1ltr10i_vs_ctrl.tab)
6. ATG12_TSSi vs Ctrl: [deseq2_results_genes_atg12i_vs_ctrl.tab](5_crispr_results/deseq2_results_genes_atg12i_vs_ctrl.tab)
7. FOSL1_TSSi vs Ctrl: [deseq2_results_genes_fosli_vs_ctrl.tab](5_crispr_results/deseq2_results_genes_fosli_vs_ctrl.tab)
8. JUN_TSSi vs Ctrl: [deseq2_results_genes_cjuni_vs_ctrl.tab](5_crispr_results/deseq2_results_genes_atg12ltr10i_vs_ctrl.tab)

- MA plots for all comparisons (template R script, then separate pdfs):
1. ATG12_LTR10i vs Ctrl: [MAplot_genes_atg12ltr10i_vs_ctrl.pdf](5_crispr_results/MAplot_genes_atg12ltr10i_vs_ctrl.pdf)
2. XRCC4_LTR10i vs Ctrl: [MAplot_genes_xrcc4ltr10i_vs_ctrl.pdf](5_crispr_results/MAplot_genes_xrcc4ltr10i_vs_ctrl.pdf)
3. MEF2D_LTR10i vs Ctrl: [MAplot_genes_mef2dltr10i_vs_ctrl.pdf](5_crispr_results/MAplot_genes_mef2dltr10i_vs_ctrl.pdf)
4. FGF2_LTR10i vs Ctrl: [MAplot_genes_fgf2ltr10i_vs_ctrl.pdf](5_crispr_results/MAplot_genes_fgf2ltr10i_vs_ctrl.pdf)
5. MCPH1_LTR10i vs Ctrl: [MAplot_genes_mcph1ltr10i_vs_ctrl.pdf](5_crispr_results/MAplot_genes_mcph1ltr10i_vs_ctrl.pdf)
6. ATG12_TSSi vs Ctrl: [MAplot_genes_atg12i_vs_ctrl.pdf](5_crispr_results/MAplot_genes_atg12i_vs_ctrl.pdf)
7. FOSL1_TSSi vs Ctrl: [MAplot_genes_fosli_vs_ctrl.pdf](5_crispr_results/MAplot_genes_fosli_vs_ctrl.pdf)
8. JUN_TSSi vs Ctrl: [MAplot_genes_cjuni_vs_ctrl.pdf](5_crispr_results/MAplot_genes_atg12ltr10i_vs_ctrl.pdf)

- Manhattan-like genomic distance plots for all candidates (template R script, then separate pdfs):
1. ATG12_LTR10i vs Ctrl: [chrplot_genes_atg12ltr10i_vs_ctrl.pdf](5_crispr_results/chrplot_genes_atg12ltr10i_vs_ctrl.pdf)
2. XRCC4_LTR10i vs Ctrl: [chrplot_genes_xrcc4ltr10i_vs_ctrl.pdf](5_crispr_results/chrplot_genes_xrcc4ltr10i_vs_ctrl.pdf)
3. MEF2D_LTR10i vs Ctrl: [chrplot_genes_mef2dltr10i_vs_ctrl.pdf](5_crispr_results/chrplot_genes_mef2dltr10i_vs_ctrl.pdf)
4. FGF2_LTR10i vs Ctrl: [chrplot_genes_fgf2ltr10i_vs_ctrl.pdf](5_crispr_results/chrplot_genes_fgf2ltr10i_vs_ctrl.pdf)
5. MCPH1_LTR10i vs Ctrl: [chrplot_genes_mcph1ltr10i_vs_ctrl.pdf](5_crispr_results/chrplot_genes_mcph1ltr10i_vs_ctrl.pdf)
6. ATG12_TSSi vs Ctrl: [chrplot_genes_atg12i_vs_ctrl.pdf](5_crispr_results/chrplot_genes_atg12i_vs_ctrl.pdf)
7. FOSL1_TSSi vs Ctrl: [chrplot_genes_fosli_vs_ctrl.pdf](5_crispr_results/chrplot_genes_fosli_vs_ctrl.pdf)
8. JUN_TSSi vs Ctrl: [chrplot_genes_cjuni_vs_ctrl.pdf](5_crispr_results/chrplot_genes_atg12ltr10i_vs_ctrl.pdf)

3) DEseq2 analysis of CRISPR-KO candidates vs untreated HCT116 cells
- Mainly the feature counts, deseq2 output tables, and plots. 
KDM6A KO: don't need to provide all the scripts but specify which samples we used as controls (will need to upload those to GEO). Put them all in a dir for easy access. Maybe include the date sequenced? Same for all the other crispri candidates, maybe make a table of all the samples that will need to be uploaded to GEO. 

## Human variation at AP1 sites within LTR10 elements

Python scripts for AP1 count histogram/correlation plots. 

Also python/jupyter scripts for gnomad deletions in human popn. 

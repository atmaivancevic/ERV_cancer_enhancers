# ERV cancer enhancers

### Scripts and files used in:

Ivancevic A, Simpson D, Chuong E (2021) "Endogenous retroviruses mediate transcriptional rewiring in response to oncogenic signaling in colorectal cancer" JOURNAL, ISSUE: PAGES. **[add biorxiv link]**

### GEO links:

### UCSC sessions:
Add a few of my UCSC sessions e.g. for each section

### Programs used:
- GIGGLE v0.6.3 (https://github.com/ryanlayer/giggle)
- bedtools v2.28.0 (http://bedtools.readthedocs.io/en/latest/)
- deepTools v3.0.1 (https://deeptools.readthedocs.io/en/develop/index.html)
- samtools v1.10 (http://www.htslib.org/)
- BBDuk/BBMap v38.05 (https://jgi.doe.gov/data-and-tools/bbtools/)
- FastQC v0.11.8 (https://github.com/s-andrews/FastQC)
- MultiQC (https://github.com/ewels/MultiQC)
- hisat2 v2.1.0 (https://github.com/DaehwanKimLab/hisat2)
- subread/featureCounts v1.6.2 (http://subread.sourceforge.net/)
- TEtranscripts v2.1.4 (https://github.com/mhammell-laboratory/TEtranscripts)
- MACS2 (https://pypi.org/project/MACS2/)
- BWA (https://github.com/lh3/bwa)
- DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- SRA Toolkit (https://hpc.nih.gov/apps/sratoolkit.html)
- MEME Suite (https://meme-suite.org/meme/)
- VSEARCH (https://github.com/torognes/vsearch)
- MUSCLE (https://www.drive5.com/muscle/)
- Jalview (https://www.jalview.org/)
- Gblocks (https://ls23l.lscore.ucla.edu/MakeTree/documentation/gblocks.html)
- FastTree (http://www.microbesonline.org/fasttree/)
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

#### 1. Identify cancer-specific regulatory regions

Predicted regulatory regions in healthy adult tissues were defined using Roadmap categories 1_TssA, 6_EnhG & 7_Enh (see [Roadmap_healthy_adult_tissues_list.txt](pancancer_epigenomic_analysis/Roadmap_healthy_adult_tissues_list.txt)). Predicted regulatory regions for 21 different cancer types were obtained from TCGA ATACseq peaks (see [TCGA_cancers_list.txt](pancancer_epigenomic_analysis/TCGA_cancers_list.txt)). Cancer-specific regulatory regions were identified by subtracting "healthy" regulatory regions from each cancer peak set.

#### 2. Test for family-level TE enrichment
GIGGLE was used to create a database of all TE families in the hg38 human genome, based on Dfam annotation. Cancer-specific regulatory regions were then searched against the repeat database (see [find_enriched_TEs.sh](pancancer_epigenomic_analysis/find_enriched_TEs.sh)). Results were ranked from highest to lowest enrichment score, and filtered based on odds ratio, score, and the number of overlaps ([SuppTable1_TCGA_giggle_results_top23TEs.tab](pancancer_epigenomic_analysis/SuppTable1_TCGA_giggle_results_top23TEs.tab)). Cancer-TE associations were visualized as bubble and volcano plots (e.g. [Fig1_bubbles.py](pancancer_epigenomic_analysis/Fig1_bubbles.py), [Fig1_volcano.R](pancancer_epigenomic_analysis/Fig1_volcano.R)). 

#### 3. Element-level TE analysis

hg38 genome coordinates of LTR10 elements were obtained from Dfam. LTR10A and LTR10F elements were merged (2kb window) and used to generate signal heatmaps over TCGA tumor ATACseq (e.g. [Fig1_deeptools_atacseq.sbatch](pancancer_epigenomic_analysis/Fig1_deeptools_atacseq.sbatch)). Similar heatmaps were plotted using other LTR10 subfamilies (e.g. [deeptools_atacseq_allLTR10.sbatch](pancancer_epigenomic_analysis/deeptools_atacseq_allLTR10.sbatch)) and other TEs. 

We further assessed TE signal in public ChIPseq datasets from Cistrome, particularly looking at histone marks indicative of enhancer activity in cancer cell lines (e.g. HCT116 H3K27ac ChIPseq & H3K4me1 ChIPseq; see [Fig1_deeptools_hct116.sbatch](pancancer_epigenomic_analysis/Fig1_deeptools_hct116.sbatch)). 

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

2i: muscle/gblocks/fasttree scripts for LTR10 consensus tree

## Control of LTR10 activity by AP1/MAPK signaling

Started with RNAseq fastq files (Cobi_24hr x2, TNF_24hr x2, Untreated_24hr x2). 

**RNAseq files [will be uploaded to GEO, add link here]:** 
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

**Workflow:**
1) [bbduk_PE.sbatch](rnaseq/bbduk_PE.sbatch)
2) [fastqc.sbatch](rnaseq/fastqc.sbatch)
3) [multiqc.sbatch](rnaseq/multiqc.sbatch) 
4) [hisat2_PE.sbatch](rnaseq/hisat2_PE.sbatch)
5) [bam_to_bw.sbatch](rnaseq/bam_to_bw.sbatch)
6) [feature_counts.sbatch](rnaseq/feature_counts.sbatch)
7) [add_gene_names.sh](rnaseq/add_gene_names.sh), using [gencode34_geneid_genename.txt](rnaseq/gencode34_geneid_genename.txt)
8) [deseq2_genes.R](rnaseq/deseq2_genes.R)

**For TE-transcripts:**
1) [hisat2_PE_k100.sbatch](rnaseq/hisat2_PE_k100.sbatch)
2) [te_transcripts.sbatch](rnaseq/te_transcripts.sbatch)
3) [extract_TEs_only.sh]
4) [deseq2_TEs.R](rnaseq/deseq2_TEs.R)

Analysis scripts for Cut&Run and RNAseq (slurm scripts would be everything from fastq, peak calling with nacs2, etc; R scripts would be deseq2 and plottting)

3a: MA plot Cobi. Redo with lfcShrink. 

3b: MA plot TNF-alpha, also need to redo. Only need to provide script for one of them (since it's the same with diff input table). Upload input feature counts tables (raw counts), normalized count tables and Deseq2 results tables for each case. 

3d: giggle volcano plot of treated C&R 

3e: R script for Log2FC plot (genes and TEs)

3f: shell script showing process to get to final list of gene/LTR10 enhancer. Table of final list.

3g: R script pheatmap of genes with treatments. Also R script pheatmap of H3K27ac Cut&Run LTR10 elements

## CRISPR silencing & deletion of LTR10 elements

Analysis scripts for RNAseq, like above. 

Python coolbox scripts and R scripts for deseq2/genome distance plots. Need to redo with lfcShrink. Also python scripts for bargraphs. And MA plots for e.g. ATG12 LTR10F. 

## Human variation at AP1 sites within LTR10 elements

Python scripts for AP1 count histogram/correlation plots. 

Also python/jupyter scripts for gnomad deletions in human popn. 

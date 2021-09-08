# ERV cancer enhancers

### Scripts and files used in:

Ivancevic A, Simpson D, Chuong E (2021) "Endogenous retroviruses mediate transcriptional rewiring in response to oncogenic signaling in colorectal cancer" JOURNAL, ISSUE: PAGES. **[add biorxiv link]**

### Programs used (add version numbers!!):
- GIGGLE v0.6.3 (https://github.com/ryanlayer/giggle)
- bedtools v2.28.0 (http://bedtools.readthedocs.io/en/latest/)
- deepTools (https://deeptools.readthedocs.io/en/develop/index.html)
- MACS2 (https://pypi.org/project/MACS2/)
- BWA (https://github.com/lh3/bwa)
- Samtools (http://www.htslib.org/)
- hisat2 (https://github.com/DaehwanKimLab/hisat2)
- BBMap/BBduk (https://jgi.doe.gov/data-and-tools/bbtools/)
- FastQC (https://github.com/s-andrews/FastQC)
- MultiQC (https://github.com/ewels/MultiQC)
- DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- SRA Toolkit (https://hpc.nih.gov/apps/sratoolkit.html)
- MEME Suite (https://meme-suite.org/meme/)
- VSEARCH (https://github.com/torognes/vsearch)
- MUSCLE (https://www.drive5.com/muscle/)
- Jalview (https://www.jalview.org/)
- Gblocks (https://ls23l.lscore.ucla.edu/MakeTree/documentation/gblocks.html)
- FastTree (http://www.microbesonline.org/fasttree/)

### Public databases:
- The Cancer Genome Atlas (TCGA), via the Genome Data Commons (https://gdc.cancer.gov/)
- CEMT epigenomes (http://www.epigenomes.ca/data-release/)
- Cistrome (http://cistrome.org/)
- Roadmap Epigenomics Project (http://www.roadmapepigenomics.org/)
- ENCODE (https://www.encodeproject.org/)
- Dfam (https://dfam.org/home)
- Human Endogenous Retrovirus Database (https://herv.img.cas.cz/)
- UCSC (https://genome.ucsc.edu/)
- JASPAR (http://jaspar.genereg.net/)
- gnomAD (https://gnomad.broadinstitute.org/)

## Pan-cancer epigenomic analysis of TE activity

#### 1. Identify cancer-specific regulatory regions

Predicted regulatory regions in healthy adult tissues were defined using Roadmap categories 1_TssA, 6_EnhG & 7_Enh (see [Roadmap_healthy_adult_tissues_list.txt](pancancer_epigenomic_analysis/Roadmap_healthy_adult_tissues_list.txt)). Predicted regulatory regions for 21 different cancer types were obtained from TCGA ATACseq peaks (see [TCGA_cancers_list.txt](pancancer_epigenomic_analysis/TCGA_cancers_list.txt)). Cancer-specific regulatory regions were identified by subtracting "healthy" regulatory regions from each cancer peak set.

#### 2. Test for family-level TE enrichment
GIGGLE was used to create a database of all TE families in the hg38 human genome, based on Dfam annotation. Cancer-specific regulatory regions were then searched against the repeat database (see [find_enriched_TEs.sh](pancancer_epigenomic_analysis/find_enriched_TEs.sh)). Results were ranked from highest to lowest enrichment score, and filtered based on odds ratio, score, and the number of overlaps ([SuppTable1_TCGA_giggle_results_top23TEs.tab](pancancer_epigenomic_analysis/SuppTable1_TCGA_giggle_results_top23TEs.tab)). Cancer-TE associations were visualized as bubble and volcano plots (e.g. [Fig1_bubbles.py](pancancer_epigenomic_analysis/Fig1_bubbles.py), [Fig1_volcano.R](pancancer_epigenomic_analysis/Fig1_volcano.R)). 

#### 3. Element-level TE analysis

hg38 genome coordinates of LTR10 elements were obtained from Dfam. LTR10A and LTR10F elements were merged (within 2kb) and used to generate signal heatmaps over TCGA tumor ATACseq (e.g. [Fig1_deeptools_atacseq.sbatch](pancancer_epigenomic_analysis/Fig1_deeptools_atacseq.sbatch)). Similar heatmaps were plotted using other LTR10 subfamilies (e.g. [deeptools_atacseq_allLTR10.sbatch](pancancer_epigenomic_analysis/deeptools_atacseq_allLTR10.sbatch)) and other TEs. 

We further assessed LTR10A/F signal in public ChIPseq datasets from Cistrome, particularly looking at histone marks from colorectal cancer cell lines (e.g. HCT116 H3K27ac ChIPseq & H3K4me1 ChIPseq; see [Fig1_deeptools_hct116.sbatch](pancancer_epigenomic_analysis/Fig1_deeptools_hct116.sbatch)). 

## Regulatory activity of LTR10 in colorectal cancer
 

## Control of LTR10 activity by MAPK signaling

## CRISPR silencing & deletion of LTR10 elements

## Human variation at AP1 sites within LTR10 elements

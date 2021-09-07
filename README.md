# ERV cancer enhancers

### Scripts and files used in:

Ivancevic A, Simpson D, Chuong E (2021) "Endogenous retroviruses mediate transcriptional rewiring in response to oncogenic signaling in colorectal cancer" JOURNAL, ISSUE: PAGES. **[add link to paper]**

### Programs used:
- GIGGLE (https://github.com/ryanlayer/giggle)
- bedtools (http://bedtools.readthedocs.io/en/latest/)
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

#### 2. Test for enrichment of TEs in cancer-specific regulatory regions
GIGGLE was used to create a searchable database of all TE families in the human genome, based on Dfam annotation. Cancer-specific regulatory regions were then searched against the repeat database (see [find_enriched_TEs.sh](pancancer_epigenomic_analysis/find_enriched_TEs.sh)). Results were ranked from highest to lowest enrichment score ([tcgaMinusRoadmap_vs_TEs_rankedByScore.tab](pancancer_epigenomic_analysis/tcgaMinusRoadmap_vs_TEs_rankedByScore.tab)). The top cancer-TE associations were shown in a bubble plot (Jupyter script FILE), and for each cancer type, the top enriched TEs were shown as volcano plots (example Rscript FILE). 

## Regulatory activity of LTR10 in colorectal cancer

## Control of LTR10 activity by MAPK signaling

## CRISPR silencing & deletion of LTR10 elements

## Human variation at AP1 sites within LTR10 elements

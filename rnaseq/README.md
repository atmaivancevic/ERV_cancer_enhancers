A typical RNAseq workflow might look like this:

1) bbduk_PE.sbatch
2) fastqc.sbatch
3) multiqc.sbatch
4) hisat2_PE.sbatch
5) bam_to_bw.sbatch
6) feature_counts.sbatch
7) add_gene_names.sh, using gencode34_geneid_genename.txt
8) deseq2_genes.R

For TE-transcripts, re-align bams to allow multiple alignments per read:

9) hisat2_PE_k100.sbatch
10) te_transcripts.sbatch
11) extract_TEs.sh
12) deseq2_TEs.R

#!/bin/bash

module load bedtools

bedtools sort -i HCT116_ATACseq_SRR8544480_peaks.narrowPeak \
> HCT116_ATACseq_SRR8544480_peaks.narrowPeak.sorted

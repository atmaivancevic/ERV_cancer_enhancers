# ERV cancer enhancers

Scripts and files used in:

Ivancevic A, Simpson D, Chuong E (2021) Endogenous retroviruses mediate transcriptional rewiring in response to oncogenic signaling in colorectal cancer. JOURNAL, ISSUE: PAGES.

#### Recommended programs
- BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- SAMtools (http://www.htslib.org/)
- BEDtools (http://bedtools.readthedocs.io/en/latest/)
- CENSOR, which requires wu-blast and bioperl (https://girinst.org/downloads/software/censor/)
- USEARCH (https://www.drive5.com/usearch/)
- VSEARCH (https://github.com/torognes/vsearch)
- MUSCLE (https://www.drive5.com/muscle/)

#### Optional
- LASTZ (https://www.bx.psu.edu/~rsharris/lastz/)
- SiLiX (http://lbbe.univ-lyon1.fr/-SiLiX-?lang=en)
- RepeatMasker (http://www.repeatmasker.org/)
- Gblocks (https://ls23l.lscore.ucla.edu/MakeTree/documentation/gblocks.html)
- HMMer (http://hmmer.org/)
- FastTree (http://www.microbesonline.org/fasttree/)

#### Prerequisites
- Some level of familiarity with computers and queuing systems (SLURM)

#### Scripts and test genome
A test genome (fungus *Yarrowia lipolytica*) has been placed in [genomes/fungi/Yarrowia.lipolytica](genomes/fungi/Yarrowia.lipolytica). We recommend trying out the workflow below on this genome first. Intermediate files for each step can be found in [results/L1](results/L1), to help with troubleshooting. Analysis scripts are provided in [scripts](scripts). LaTeX files for the Additional Files (Supp Tables and Figures) are provided in [latex](latex).

## Workflow

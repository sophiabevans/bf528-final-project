A partial reproduction of the study Transcriptional Reversion of Cardiac Myocyte Fate During Mammalian Cardiac Regeneration [https://pubmed.ncbi.nlm.nih.gov/25477501/], performing Programmer and Analyst tasks from https://bf528.readthedocs.io/en/latest/content/projects/project_2_rnaseq_1/project_2_rnaseq_1.html.

Directory scripts/ contains all scripts used in this analysis, with output from scripts in scripts/output. Large output files are not included in this directory, and are on the scc in a separate folder. Directory plots/ contains plots and tables used in the report. The report can be found in the parent directory (Bevans_528_final.pdf).

To replicate this portion of the analysis, run scripts in this order:
run_tophat.qsub -> geneBody.qsub -> inner_dist.qsub -> bam_stat.qsub -> run_cufflinks.qsub -> run_cuffdiff.qsub -> final.R

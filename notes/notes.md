# Cost and Time Estimates for Exhaustive Screening

- # of CTD2 compounds = 546
- # of COSMIC compounds = 1390
- Assume best 3 of 6 input dataset types like gene expression, pathway expression, rna-seq, exome seq, etc.
- Assume 3x COSMIC training times with 400 cell lines
- Assume COSMIC training time of 1500 seconds per core
- 1600 compounds x 3 datasets x 3 for 1000 cell lines x 1500 seconds = 21,600,000 seconds for one core = 250 days
- Numbers around what a high number of compounds is: http://news.mit.edu/2016/faster-gene-expression-profiling-drug-discovery-0128
# TAAIgM
The scripts related to the manuscript Patterns of IgM Binding to Tumor Associated Antigen Peptides Correlate with the Type of Brain Tumors

Description:
rawDataProc.R - Gpr files reading, summarizing, cleanin and normalization of data. 
vn1 - Data - the matrix of cleaned and normalized data produced by rawDataProc.R
TAAnlz.R - Main analysis script using vn1 as a input.
SVMtest.R - A script producing the SVM models for the leave one out validation.
mapll.csv - MAP of the gpr files.
pep_info0.csv - Biological information about the peptides.

All the remaining functions are called from the described files and have descriptons in them.

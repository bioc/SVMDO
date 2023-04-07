# First of all, assign your dataset to a variable (dataset has to be in txt form)

# Assume that, variable name is "x", tab-separated and first line is header

# x<-read.table("dataset_name.txt",sep="\t",header=TRUE)

# To create SummarizedExperiment object, SummarizedExperiment() function from the
# SummarizedExperiment Bioconductor package is used

# Assume that a variable named as "sum_data" is used for SummarizedExperiment object creation

# library(SummarizedExperiment)
# sum_data<-SummarizedExperiment(x)

# To save your object into a rda file, use saveRDS() base function

# Assume that rda file is named as prep_sum_data

# saveRDS(sum_data,file="prep_sum_data")



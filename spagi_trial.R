library(data.table)
#> Warning: package 'data.table' was built under R version 3.3.3
library(spagi)
example(spagi)

# Pre-process the query data

# The query data is already in CPM and log2 normalized form. 
# Here, we will use the expression cutoff as 1.8 of the query data.

ROR1.processed.data <- preprocess_querydata(cell.tissue.data = ROR1.data, exp.cutoff.th = 1.8)

# Identify active pathway paths of the processed query data

# Here we will use the background pathway.path data
# to get the active pathway paths of the processed query data.

ROR1.active.pathway <- identify_active_pathway_path(pathway.path = pathway.path, processed.query.data = ROR1.processed.data)


ROR1.active.pathway2 <- identify_active_pathway_path(pathway.path = pathway.path, processed.query.data = ROR1.data)

# To save the active pathway paths in csv format

# After getting the active pathway paths, you can save the data in csv file format
# where each row will denote a path and all the paths starting from the same source comprises of a pathway of the source protein.

lapply(unlist(ROR1.active.pathway$ROR1_LEC, recursive = F, use.names = F), write, "ROR1.active.pathway.csv", append=T, ncolumns=10)

# Get active pathway ranking metric (i.e., activity score and number of downstream transcription factors)

# Here we will use the ROR1.active.pathway and ROR1.processed.data data sets to get the acitve pathway ranking metric. Also we will use a high expression threshold (here we will use 7) for the processed query data.

ROR1.active.pathway.ranking.metric <- get_pathway_ranking_metric(active.pathway.path = ROR1.active.pathway, processed.query.data = ROR1.processed.data, high.exp.th = 7)

# Plot the ranking metric result in a 2D plane

# After getting the active pathway ranking metric result you can display them in your preferred format.
# Here we will plot them in a 2D plane where x-axis denotes the number of downstream transcription factors and 
# y-axis denotes the activity score for each pathway.

display_pathway_ranking_metric(pathway.ranking.metric = ROR1.active.pathway.ranking.metric)
#> [1] "ROR1_LEC -- result plotting done!!"
#To separate the top ranked pathways we can do this
abline(v=45, h=0.2, lty=2, col="black")





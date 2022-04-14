#options(warn=-1)


library(rnetcarto)
library(igraph)


get_cartography_and_save_result <- function(g, folder){

  adj <- get.adjacency(g,sparse=FALSE)
  result <- netcarto(adj)[[1]]
  rownames(result) <- result$name

  result <- result[,2:5]
  write.csv(result, file = paste0(folder, "/base_natwork_analysis.csv"))

  #return()
}




## main

folder <- commandArgs(trailingOnly = T)[1]


d <- read.csv(paste0(folder, "/linkList.csv"))
g <- graph.data.frame(d[1:2], directed = T)
E(g)$weight <- d[[3]]
get_cartography_and_save_result(g, folder)

message("finished")

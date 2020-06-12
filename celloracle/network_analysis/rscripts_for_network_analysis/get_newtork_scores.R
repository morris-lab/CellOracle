library(igraph)
library(linkcomm)

library(rnetcarto)

#Functional Cartography of Complex Networks
category_analysis <- function(g, data){
  # data: result of community detection

  mem <- data$membership
  num_mod <- max(mem)
  num_nodes <- vcount(g)
  deg <- degree(g)

  # Calculation of the within-module degree
  z_score <- numeric(num_nodes)
  for(s in 1:num_mod){
    v_seq <- subset(V(g),mem==s)
    g_sub <- delete.vertices(g,subset(V(g),mem!=s))
    deg_sub <- degree(g_sub)
    z_score[v_seq] <- (deg_sub - mean(deg_sub)) / sd(deg_sub)
  }

  # Calculation of the participation coefficient
  participation_coeff <- numeric(num_nodes) + 1.0
  for(i in 1:num_nodes){
    mem_nei <- mem[neighbors(g,i)]
    for(s in 1:num_mod){
      deg_to_s <- length(subset(mem_nei,mem_nei == s))
      participation_coeff[[i]] <- participation_coeff[[i]] - (deg_to_s / deg[[i]]) ** 2.0
    }
  }

  # Classification
  role <- numeric(num_nodes)

  result <- as.data.frame(cbind(z_score,participation_coeff))
  names(result) <- c("within_module_degree","participation_coefficient")

  # R1: Ultra-peripheral nodes
  v_seq <- which(result$within_module_degree<2.5 & result$participation_coeff<0.05)
  role[v_seq] <- "R1: Ultra-peripheral nodes"

  # R2: Peripheral nodes
  v_seq <- which(result$within_module_degree<2.5 & result$participation_coeff>=0.05 & result$participation_coeff<0.625)
  role[v_seq] <- "R2: Peripheral nodes"

  # R3: Non-hub connectors
  v_seq <- which(result$within_module_degree<2.5 & result$participation_coeff>=0.625 & result$participation_coeff<0.8)
  role[v_seq] <- "R3: Non-hub connectors"

  # R4: Non-hub kinless nodes
  v_seq <- which(result$within_module_degree<2.5 & result$participation_coeff>=0.8)
  role[v_seq] <- "R4: Non-hub kinless nodes"

  # R5: Provincial hubs
  v_seq <- which(result$within_module_degree>=2.5 & result$participation_coeff<0.3)
  role[v_seq] <- "R5: Provincial hubs"

  # R6: Connector hubs
  v_seq <- which(result$within_module_degree>=2.5 & result$participation_coeff>=0.3 & result$participation_coeff<0.75)
  role[v_seq] <- "R6: Connector hubs"

  # R7: Kinless hubs
  v_seq <- which(result$within_module_degree>=2.5 & result$participation_coeff>=0.75)
  role[v_seq] <- "R7: Kinless hubs"

  # 結果をまとめる
  result <- cbind(result,role)
  rownames(result) <- names(degree(g))
  return(result)

}

get_scores_cartography <- function(g){
  adj <- get.adjacency(g,sparse=FALSE)
  score <- netcarto(adj)[[1]]
  rownames(score) <- score$name
  return(score)
}

calculateNetworkScores <- function(g, folder){

  ### degree

  res <- degree(g, mode = "all")
  res <- data.frame(res)
  colnames(res) <- c("degree_all")

  res["degree_in"] <- degree(g, mode = "in")
  res["degree_out"] <- degree(g, mode = "out")
  ### clustering coefficient
  # clustering coefficient
  res["clustering_coefficient"] <- transitivity(g,type="local",isolates="zero")
  # clustering coefficient for weighted network
  res["clustering_coefficient_weighted"] <- transitivity(g,vids=V(g),type="weighted",isolates="zero")

  ### centrality
  # degree centrality
  res["degree_centrality_all"] <- res["degree_all"] / (vcount(g) - 1)
  res["degree_centrality_in"] <- res["degree_in"] / (vcount(g) - 1)
  res["degree_centrality_out"] <- res["degree_out"] / (vcount(g) - 1)
  # Betweenness centrality
  res["betweenness_centrality"] <- betweenness(g)
  # Closeness centrality
  res["closeness_centrality"] <- closeness(g)
  # eigenvector centrality
  res["eigenvector_centrality"] <- evcent(g)$vector
  # page rank
  res["page_rank"] <-  page.rank(g)$vector

  res["assortative_coefficient"] <- assortativity.degree(g)

  res["average_path_length"] <- average.path.length(g)

  ## community detection (non-overlapping)


  # community based on Edge betweenness (Girvan-Newman algorithm) THIS IS DEPRECATED NOW
  #data <- edge.betweenness.community(g)
  #res["community_edge_betweenness"] <- data$membership

  # random walk
  data <- walktrap.community(g, modularity=TRUE)

  res["community_random_walk"] <- data$membership
  # eigenvector THIS IS DEPRECATED NOW
  #suppressWarnings(data <- leading.eigenvector.community(g,options=list(maxiter=1000000, ncv=5)))
  #res["community_eigenvector"] <- data$membership

  #result <- category_analysis(g, data_)

  result <- get_scores_cartography(g)
  result <- result[rownames(res),2:5]

  result <- cbind(res, result)

  write.csv(result, file = paste0(folder, "/base_natwork_analysis.csv"))

  #return()
}




## main

folder <- commandArgs(trailingOnly = T)[1]


d <- read.csv(paste0(folder, "/linkList.csv"))
g <- graph.data.frame(d[1:2], directed = T)
E(g)$weight <- d[[3]]

calculateNetworkScores(g, folder)

message("finished")

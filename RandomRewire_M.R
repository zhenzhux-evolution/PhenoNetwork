rm(list=ls())
library(foreach)
library(dplyr)
library(psych)
library(igraph)
require(doParallel)
num_cores <- detectCores()
doParallel::registerDoParallel(cores = num_cores)

setwd("/home/zhenzhux/pheno")
output_directory = "/home/zhenzhux/pheno/Res_202506"

phenotype <- read.table("./Dorey_2024.txt", header = TRUE)

Traits_2024 <- phenotype |>
  select(1:7, 12:21, 27:40, 44, 47, 49, 50, 51:64)
colnames(Traits_2024[8:49])

colnames(Traits_2024)
Traits_2024 |>
  group_by(Treatment, Generation, Replicate) |>
  summarise(n= n())

Traits_2024 <- Traits_2024 |> mutate(Treatment = recode(Treatment,
                                                        "THB"="ABC",  
                                                        "THH"="aBC",  
                                                        "TNHB"="AbC", 
                                                        "TNHH"="abC",
                                                        "LHB"="ABc",  
                                                        "LHH"="aBc",   
                                                        "LNHB"="Abc", 
                                                        "LNHH"="abc"  
))


get_corr_matrix <- function(data) {
  res <- corr.test(data, use = "pairwise", method = "spearman", adjust = "none", ci = FALSE)
  r <- res$r
  p <- res$p
  p_adj <- p.adjust(p[upper.tri(p)], method = "BH")
  p_adj_matrix <- matrix(0, nrow = nrow(p), ncol = ncol(p))
  p_adj_matrix[upper.tri(p_adj_matrix)] <- p_adj
  p_adj_matrix[lower.tri(p_adj_matrix)] <- t(p_adj_matrix)[lower.tri(p_adj_matrix)]
  r[p_adj_matrix > 0.05] <- 0
  diag(r) <- 0
  r
}

# select the best M and membership through 100 runs of Louvain cluster
run_best_louvain <- function(graph, n_runs = 100) {
  best_modularity <- -Inf
  best_result <- NULL
  
  for (i in 1:n_runs) {
    result <- cluster_louvain(graph)
    mod <- modularity(result)
    if (mod > best_modularity) {
      best_modularity <- mod
      best_result <- result
    }
  }
  
  return(list(
    clustering = best_result,
    modularity = round(best_modularity, 4),
    membership = membership(best_result)
  ))
}

# single network stats
Network_metrics <- data.frame()

for (Trt in unique(Traits_2024$Treatment)) { 
  SingNet <- Traits_2024 |> filter(Treatment == Trt) %>% select(8:49)
  SingNet_CorrMatrix <- get_corr_matrix(SingNet)
  SingNet_r <- ifelse(SingNet_CorrMatrix != 0, 1, 0) 
  
  diag(SingNet_r) <- 0  
  isSymmetric(SingNet_r) 
  sum(SingNet_r[upper.tri(SingNet_r)] != 0)
  
  SingNet_graph <- graph_from_adjacency_matrix(SingNet_r, mode = "undirected", weighted = NULL, diag = FALSE)
  
  # Apply Louvain algorithm to two generations in each trt
  SingNet_Louvain <- run_best_louvain(SingNet_graph, n_runs = 100)
  SingleNet_Modu_Louvain <- round(SingNet_Louvain$modularity, 4); SingleNet_Modu_Louvain
  memb_SingleNet <- SingNet_Louvain$membership;unique(memb_SingleNet)
  
  Node_num_SingleNet <- vcount(SingNet_graph)
  edge_num_SingleNet <- ecount(SingNet_graph);
  node_Mdegree_SingleNet <- mean(degree(SingNet_graph))
  betweenModu <- sum(crossing(SingNet_Louvain, SingNet_graph))
  withinModu <- sum(ecount(SingNet_graph) - betweenModu)
  
  
  sum(SingNet_r[upper.tri(SingNet_r)] != 0)
  sum(SingNet_r[upper.tri(SingNet_r)] == 0)
  num_pair_SingleNet <- sum(SingNet_r[upper.tri(SingNet_r)] != 0)
  num_modules_SingleNet <- length(unique(memb_SingleNet))
  module_sizes_SingleNet <- mean(sizes(SingNet_Louvain))
  density_ratio_SingleNet <- sum(!crossing(SingNet_Louvain, SingNet_graph))/sum(crossing(SingNet_Louvain, SingNet_graph))
  
  
  
  # Simulate random networks
  num_simu = 1000
  random_modu_SingNet <- numeric(num_simu)
  for (i in 1:num_simu) {
    #random_graph_SingleNet <- sample_gnm(Node_num_SingleNet, edge_num_SingleNet, directed = FALSE)
    random_graph_SingleNet <- rewire(SingNet_graph, with = keeping_degseq(niter = ecount(SingNet_graph) * 10))
    random_clusters_SingleNet <- run_best_louvain(random_graph_SingleNet, n_runs = 100)
    random_modu_SingNet[i] <- round(random_clusters_SingleNet$modularity, 4) 
  }
  
  p_value_SingleNet <- sum(random_modu_SingNet >= SingleNet_Modu_Louvain)/num_simu
  Modu_mean_SingleNet <- mean(random_modu_SingNet); Modu_sd_SingleNet <- sd(random_modu_SingNet)
  z_score_SingleNet <- (SingleNet_Modu_Louvain - Modu_mean_SingleNet) / Modu_sd_SingleNet
  
  
  Network_metrics <- rbind(Network_metrics, data.frame(
    Treatment = Trt, 
    numNodes = Node_num_SingleNet,
    numPairs = num_pair_SingleNet,
    betweenModu = betweenModu, 
    withinModu = withinModu,
    MeanDegree = round(node_Mdegree_SingleNet,4),
    Num_module = num_modules_SingleNet,
    MeanNode_modu = round(module_sizes_SingleNet,4),
    densityRatio = round(density_ratio_SingleNet, 4),
    Modularity = SingleNet_Modu_Louvain, 
    mean_modu = round(Modu_mean_SingleNet, 4),
    sdModu = round(Modu_sd_SingleNet, 4),
    Modu_pvalue = round(p_value_SingleNet,4),padj=NA,
    Zscore_modu = round(z_score_SingleNet,4)
  ))
  
}
Network_metrics$padj <- p.adjust(Network_metrics$Modu_pvalue, "BH")
write.table(Network_metrics, file.path(output_directory,"Network_q01_metrics_G1G10.txt"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")



# edge number 
memb_SingleNet <- SingNet_Louvain$membership;unique(memb_SingleNet)
between_edges <- E(SingNet_graph)[memb_SingleNet[ends(SingNet_graph, E(SingNet_graph))[,1]] != memb_SingleNet[ends(SingNet_graph, E(SingNet_graph))[,2]]]
length(between_edges)

ends(SingNet_graph, E(SingNet_graph))


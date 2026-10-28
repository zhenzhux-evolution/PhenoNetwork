#!/usr/bin/env Rscript
rm(list=ls())
library(foreach)
library(dplyr)
library(psych)
library(igraph)
library(aricode)
library(ggplot2)
library(patchwork)
require(doParallel)
num_cores <- detectCores()
doParallel::registerDoParallel(cores = 20)

setwd("/home/zhenzhux/pheno")
output_directory = "/home/zhenzhux/pheno"
phenotype <- read.table("./Dorey_2024.txt", header = TRUE)

Traits_2024 <- phenotype |>
  select(1:7, 12:21, 27:40, 51:58, 44, 47, 49, 50, 59:64)
colnames(Traits_2024[8:49])

head(Traits_2024)
unique(Traits_2024$Treatment)

# rename treatments
Traits_2024 <- Traits_2024 |> mutate(Treatment = recode(Treatment, "G1"="G1", 
                                                        "THB"="ABC",  
                                                        "THH"="aBC",  
                                                        "TNHB"="AbC", 
                                                        "TNHH"="abC",
                                                        "LHB"="ABc",  
                                                        "LHH"="aBc",   
                                                        "LNHB"="Abc", 
                                                        "LNHH"="abc"  
))

# prepare binary matrix based on Spearman's correlation
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



Treatments <- setdiff(Traits_2024$Treatment, "G0")
Trt_pairs <- combn(Treatments, 2, simplify = FALSE)


perm_results_all <- data.frame()


G0 <- Traits_2024 |> filter(Treatment == "G0") %>% select(8:49)
G0_r <- get_corr_matrix(G0)
G0_matrix <- ifelse(G0_r != 0, 1, 0) 
G0_upper <- G0_matrix[upper.tri(G0_matrix)];sum(G0_upper)
G0_upper_idx <- which(upper.tri(G0_matrix), arr.ind = TRUE)
G0_idx_only <- which((G0_upper == 1)); 


perm_results_all <- foreach(Trt = Trt_pairs, .combine = 'rbind') %do% {
  Trt1 <- Trt[1]; Trt2 <- Trt[2]
  group1 <- Traits_2024 |> filter(Treatment == Trt1) %>% select(8:49)
  group2 <- Traits_2024 |> filter(Treatment == Trt2) %>% select(8:49)
  
  # Compute matrices
  group1_r <- get_corr_matrix(group1)
  group2_r <- get_corr_matrix(group2)
  
  # Binary 0/1 matrices
  group1_r_bin <- ifelse(group1_r != 0, 1, 0)
  group2_r_bin <- ifelse(group2_r != 0, 1, 0)
  
  # Symmetry check
  isSymmetric(group1_r_bin)
  isSymmetric(group2_r_bin)
  
  #network construction 
  group1_graph <- graph_from_adjacency_matrix(group1_r_bin, mode = "undirected", weighted = NULL, diag = FALSE)
  group2_graph <- graph_from_adjacency_matrix(group2_r_bin, mode = "undirected", weighted = NULL, diag = FALSE)
  
  group1_result <- run_best_louvain(group1_graph, n_runs = 100)
  group1_opt <- group1_result$clustering
  group1_Modu_opt <- group1_result$modularity
  memb_group1_opt <- group1_result$membership
  
  group2_result <- run_best_louvain(group2_graph, n_runs = 100)
  group2_opt <- group2_result$clustering
  group2_Modu_opt <- group2_result$modularity
  memb_group2_opt <- group2_result$membership
  
  Node_num_group1 <- vcount(group1_graph); edge_num_group1 <- ecount(group1_graph);
  Node_num_group2 <- vcount(group2_graph); edge_num_group2 <- ecount(group2_graph);
  
  
  # Compare mudular structure
  nmi_obs <- compare(group1_opt, group2_opt, method = "nmi")
  moduDiff_obs <- abs(group2_Modu_opt - group1_Modu_opt)
  
  group1_upper <- group1_r_bin[upper.tri(group1_r_bin)];sum(group1_upper)
  group2_upper <- group2_r_bin[upper.tri(group2_r_bin)];sum(group2_upper)
  upper_idx <- which(upper.tri(group1_r_bin), arr.ind = TRUE)
  #group1_idx <- which((group1_upper == 1)); group2_idx <- which((group2_upper == 1)); 
  
  # Function to add symmetric edge
  add_edge <- function(mat, idx_list) {
    for (i in idx_list) {
      row <- upper_idx[i, 1]
      col <- upper_idx[i, 2]
      mat[row, col] <- 1
      mat[col, row] <- 1
    }
    return(mat)
  }  
  
  
  group1_upper <- group1_r_bin[upper.tri(group1_r_bin)];sum(group1_upper)
  group1_sharedG0_idx_only <- which((G0_upper == 1) & (group1_upper == 1))
  group1_n_sharedG0 <- length(group1_sharedG0_idx_only)
  group1_n_gain <- sum((G0_upper == 0) & (group1_upper == 1));
  group1_n_loss <- sum((G0_upper == 1) & (group1_upper == 0))
  
  group2_upper <- group2_r_bin[upper.tri(group2_r_bin)];sum(group2_upper)
  group2_sharedG0_idx_only <- which((G0_upper == 1) & (group2_upper == 1))
  group2_n_sharedG0 <- length(group2_sharedG0_idx_only)
  group2_n_gain <- sum((G0_upper == 0) & (group2_upper == 1)) 
  group2_n_loss <- sum((G0_upper == 1) & (group2_upper == 0))
  
  
  # unique edge number simulations 
  permu_num <- 1000
  
  perm_results <- foreach(perm = 1:permu_num, .combine = 'rbind') %dopar% {
    # randomly sample unique edge number from group1_only
    group1_loss_idx_sim <- sample(G0_idx_only, group1_n_loss, replace = FALSE);
    length(group1_loss_idx_sim)
    group1_after_loss_sim <- setdiff(G0_idx_only, group1_loss_idx_sim);
    length(group1_after_loss_sim) == group1_n_sharedG0
    sum(length(group1_after_loss_sim) + sum(group1_n_gain)) == sum(group1_upper)
    group1_all_possible_idx <- which(G0_upper == 0) 
    group1_gain_idx_sim <- sample(group1_all_possible_idx, group1_n_gain, replace = FALSE)
    group1_matrix_simu <- matrix(0, nrow = nrow(group1_r_bin), ncol = ncol(group1_r_bin))
    group1_matrix_simu <- add_edge(group1_matrix_simu, group1_after_loss_sim)
    group1_matrix_simu <- add_edge(group1_matrix_simu, group1_gain_idx_sim)
    
    group2_loss_idx_sim <- sample(G0_idx_only, group2_n_loss, replace = FALSE);length(group2_loss_idx_sim)
    group2_after_loss_sim <- setdiff(G0_idx_only, group2_loss_idx_sim);
    length(group2_after_loss_sim) == group2_n_sharedG0
    sum(length(group2_after_loss_sim) + sum(group2_n_gain)) == sum(group2_upper)
    group2_all_possible_idx <- which(G0_upper == 0) 
    group2_gain_idx_sim <- sample(group2_all_possible_idx, group2_n_gain, replace = FALSE)
    group2_matrix_simu <- matrix(0, nrow = nrow(group2_r_bin), ncol = ncol(group2_r_bin))
    group2_matrix_simu <- add_edge(group2_matrix_simu, group2_after_loss_sim)
    group2_matrix_simu <- add_edge(group2_matrix_simu, group2_gain_idx_sim)
    
    
    group1_perm_graph <- graph_from_adjacency_matrix(group1_matrix_simu, mode = "undirected", weighted = NULL, diag = FALSE)
    group2_perm_graph <- graph_from_adjacency_matrix(group2_matrix_simu, mode = "undirected", weighted = NULL, diag = FALSE)
    
    # Apply Best Louvain algorithm
    group1_perm_result <- run_best_louvain(group1_perm_graph, n_runs = 100)
    group1_perm_opt <- group1_perm_result$clustering
    group1_perm_Modu_opt <- group1_perm_result$modularity
    memb_group1_perm_opt <- group1_perm_result$membership
    
    group2_perm_result <- run_best_louvain(group2_perm_graph, n_runs = 100)
    group2_perm_opt <- group2_perm_result$clustering
    group2_perm_Modu_opt <- group2_perm_result$modularity
    memb_group2_perm_opt <- group2_perm_result$membership
    
    # Compare permuted networks
    nmi_perm <- compare(group1_perm_opt, group2_perm_opt, method = "nmi")
    moduDiff_perm <- abs(group2_perm_Modu_opt - group1_perm_Modu_opt)
    
    return(data.frame(
      Permutation = perm,
      Trt_perm_1 = Trt1, 
      Trt_perm_2 = Trt2, 
      Perm_n1_sample = nrow(group1),
      Perm_n2_sample = nrow(group2),
      Modu_g1_perm = group1_perm_Modu_opt,
      Modu_g2_perm = group2_perm_Modu_opt,
      moduDiff = round(moduDiff_perm, 3),
      nmi = round(nmi_perm, 3)
    ))
  }
  
  
  
  group1_mean_modu_perm <- mean(perm_results$Modu_g1_perm)  
  group1_sdModu_perm <- sd(perm_results$Modu_g1_perm)
  group1_p_modu <- sum(abs(perm_results$Modu_g1_perm-group1_mean_modu_perm) >= abs(group1_Modu_opt-group1_mean_modu_perm)) / permu_num
  
  
  group2_mean_modu_perm <- mean(perm_results$Modu_g2_perm)  
  group2_sdModu_perm <- sd(perm_results$Modu_g2_perm)
  group2_p_modu <- sum(abs(perm_results$Modu_g2_perm-group2_mean_modu_perm) >= abs(group2_Modu_opt-group2_mean_modu_perm)) / permu_num
  
  
  mean_moduDiff <- mean(perm_results$moduDiff);sd_moduDiff <- sd(perm_results$moduDiff)
  p_moduDiff <- sum(abs(perm_results$moduDiff-mean_moduDiff) >= abs(moduDiff_obs-mean_moduDiff)) / permu_num
  
  
  meanNMI <- mean(perm_results$nmi);sdNMI <- sd(perm_results$nmi)
  p_nmi <- sum(abs(perm_results$nmi-meanNMI) >= abs(nmi_obs-meanNMI)) / permu_num
  
  
  return(data.frame(
    Treatment_1 = Trt1,
    Treatment_2 = Trt2,
    Trt1_Modularity = round(group1_Modu_opt,3),
    Trt2_Modularity = round(group2_Modu_opt,3),
    deltaM_obs = round(moduDiff_obs, 3),
    deltaM_p = round(p_moduDiff, 3), moduDiff_p_adj = NA,
    deltaM_Expectchance_mean = round(mean_moduDiff, 3),
    deltaM_Expectchance_sd = round(sd_moduDiff, 3),
    nmi = round(nmi_obs, 3),
    nmi_p = round(p_nmi, 3),nmi_p_adj = NA, 
    NMI_Expectchance_mean = round(meanNMI, 3),
    NMI_Expectchance_sd = round(sdNMI, 3)
  ))
}


perm_results_all$deltaM_p_adj <- round(p.adjust(perm_results_all$deltaM_p, "BH"), 3)
perm_results_all$nmi_p_adj <- round(p.adjust(perm_results_all$nmi_p, "BH"), 3)

write.table(perm_results_all,
            file.path(output_directory, "./edgeSimu_Trtselection.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)   





# this is to test the observed NMI between G0 and G10 networks was upward biased significantly across 8 pairs
binom.test(8, 8, p = 0.5, alternative = "two.sided")

# this is to test the observed NMI between treatment networks was upward biased significantly across 28 pairs
binom.test(28, 28, p = 0.5, alternative = "two.sided")



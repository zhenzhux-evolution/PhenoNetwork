#!/usr/bin/env Rscript

rm(list=ls())
library(foreach)
library(tidyr)
library(dplyr)
library(tibble)
library(psych)
library(igraph)
library(ggplot2)
library(patchwork)
require(doParallel)
registerDoParallel(cores = 64)

setwd("~/Downloads/Bioinfo/data_analysis/Network/dataSets/Dorey_2024_NC/2024_NC/NetworkTest/SingleNet")
output_directory ="~/Documents/demon/UZH/DataSci_UM/network_genomics/Script4Publication"

phenotype <- read.table("~/Downloads/Bioinfo/data_analysis/Network/dataSets/Dorey_2024_NC/2024_NC/Dorey_2024.txt", header = TRUE)
#colnames(phenotype)
Traits_2024 <- phenotype |>
  select(1:7, 12:21, 27:40, 51:58, 44, 47, 49, 50, 59:64)

Traits_2024 <- Traits_2024 |> mutate(Treatment = recode(Treatment, 
                                                        "G1"="G0",    
                                                        "THB"="ABC",  
                                                        "THH"="aBC",  
                                                        "TNHB"="AbC", 
                                                        "TNHH"="abC",
                                                        "LHB"="ABc",  
                                                        "LHH"="aBc",   
                                                        "LNHB"="Abc", 
                                                        "LNHH"="abc",  
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

build_graph <- function(data) {
  bin_matrix <- ifelse(data != 0, 1, 0)
  graph_from_adjacency_matrix(bin_matrix, mode = "undirected", weighted = NULL, diag = FALSE)
}


Treatments <- c("G0", "ABC", "aBC", "AbC", "abC", "ABc", "aBc", "Abc", "abc")
#unique(Traits_2024$Treatment);
n_iter <- 100
Trt_results_all <- foreach(i = seq_along(Treatments), .combine = 'rbind') %do% {
  Trt <- Treatments[i]
  
  group1 <- Traits_2024 %>% filter(Treatment == Trt) %>% select(8:49)
  corr_matrix <- get_corr_matrix(group1) 
  graph <- build_graph(corr_matrix)
  obs_cluster <- cluster_optimal(graph, weights = NULL)
  obs_modu <- obs_cluster$modularity
  obs_result <- data.frame(
    Treatment = Trt,
    Modu = round(obs_modu, 4),
    type = "obs"
  )
  
  perm_results <- foreach(i = 1:n_iter, .combine = rbind) %dopar% {
    sim_cluster <- cluster_optimal(graph, weights = NULL)
    sim_modu <- sim_cluster$modularity
    nmi <- compare(obs_cluster, sim_cluster)
    
    data.frame(
      Treatment = Trt, 
      Modu = round(sim_modu, 4), 
      type = "perm"
    )
    
  }
  rbind(obs_result, perm_results)
}

obs_mod_df <- Trt_results_all %>% filter(type == "obs")
obs_mod_values <- deframe(obs_mod_df %>% select(Treatment, Modu))

plot_all_modularity <- function(results_df, obs_mod_values) {
  Treatments <- unique(results_df$Treatment)
  plots <- lapply(Treatments, function(Trt) {
    df_trt <- filter(results_df, Treatment == Trt & type == "perm")
    obs_val <- obs_mod_values[Trt]
    max_val <- max(df_trt$Modu)
    maxNumber <- sum(abs(df_trt$Modu - max_val) < 1e-6)
    
    ggplot(df_trt, aes(x = Modu)) +
      geom_histogram(bins = 30, fill = "grey80", color = "black") +
      geom_vline(xintercept = obs_val, color = "blue", linetype = "dashed", linewidth = 1) +
      geom_vline(xintercept = max_val, color = "red", linetype = "dotted", linewidth = 1) +
      labs(title = paste0(Trt, "\nObs: ", round(obs_val, 3), ", Max: ", round(max_val, 3)),
           x = "Modularity", y = "Count") +
      theme_minimal()
  })
  
  wrap_plots(plots, ncol = 3) + plot_annotation(title = "Optimal modularity distribution")
}

#pdf(file = file.path(output_directory, "opt_Modu_distribution.pdf"), width = 12, height = 7)
png(file = file.path(output_directory, "opt_Modu_distribution_100.png"), width = 12, height = 7, 
    units = "in", res = 300)
plot <- plot_all_modularity(Trt_results_all, obs_mod_values)
print(plot)
dev.off()    




results <- list()
obs_mod_values <- list()
for (Trt in Treatments) {
  
  group1 <- Traits_2024 %>% filter(Treatment == Trt) %>% select(8:49)
  
  corr_matrix <- get_corr_matrix(group1)
  graph <- build_graph(corr_matrix)
  
  # Observed clustering results
  obs_cluster <- cluster_louvain(graph)
  obs_modu <- modularity(obs_cluster)
  obs_mod_values[[Trt]] <- obs_modu
  obs_members <- factor(membership(obs_cluster))
  
  # Parallel bootstrap
  group_results <- foreach(i = 1:n_iter, .combine = rbind) %dopar% {
    sim_cluster <- cluster_louvain(graph)
    sim_modu <- modularity(sim_cluster)
    sim_members <- factor(membership(sim_cluster))
    
    data.frame(
      Treatment = Trt,
      Louvain_mod = round(sim_modu, 6)
    )
  }
  results[[Trt]] <- group_results 
  cat("Finished group:", Trt, "\n")
}


results_combined <- bind_rows(results)
head(results_combined)

plot_all_modularity <- function(results_df, obs_mod_values) {
  Treatments <- unique(results_df$Treatment)
  plots <- lapply(Treatments, function(trt) {
    df_trt <- filter(results_df, Treatment == trt)
    obs_val <- round(obs_mod_values[[trt]], 6)
    max_val <- round(max(df_trt$Louvain_mod), 6)
    maxNumber <- sum(abs(df_trt$Louvain_mod - max_val) < 1e-6)
    
    ggplot(df_trt, aes(x = Louvain_mod)) +
      geom_histogram(bins = 30, fill = "grey80", color = "black") +
      geom_vline(xintercept = obs_val, color = "blue", linetype = "dashed", linewidth = 1) +
      geom_vline(xintercept = max_val, color = "red", linetype = "dotted", linewidth = 1) +
      labs(title = paste0(trt, "\nObs: ", round(obs_val, 3), ", Max: ", round(max_val, 3)),
           x = "Modularity", y = "Count") +
      theme_minimal()
  })
  
  wrap_plots(plots, ncol = 3) + plot_annotation(title = "Louvain modularity distribution")
}

#pdf(file = file.path(output_directory, "Louvain_Modu_distribution_100test.pdf"), width = 12, height = 7)
png(file = file.path(output_directory, "Louvain_Modu_distribution_100.png"), width = 12, height = 7, 
    units = "in", res = 300)
plot_all_modularity(results_combined, obs_mod_values)
dev.off()





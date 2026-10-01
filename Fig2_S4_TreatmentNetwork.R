#!/usr/bin/env Rscript

# 2026-03-08 21:14:52 #### modularity in M
rm(list=ls())
library(foreach)
library(dplyr)
library(psych)
library(igraph)

setwd("~/Downloads/Bioinfo/data_analysis/Network/dataSets/Dorey_2024_NC/2024_NC/NetworkTest/SingleNet")
output_directory ="/Users/xiaodemon/Documents/demon/UZH/DataSci_UM/network_genomics/Script4Publication/"

phenotype <- read.table("~/Downloads/Bioinfo/data_analysis/Network/dataSets/Dorey_2024_NC/2024_NC/Dorey_2024.txt", header = TRUE)
#colnames(phenotype)
Traits_2024 <- phenotype |>
  select(1:7, 12:21, 27:40, 51:58, 44, 47, 49, 50, 59:64)
colnames(Traits_2024[8:49])

colnames(Traits_2024)
Traits_2024 |>
  group_by(Treatment, Generation, Replicate) |>
  summarise(n= n())

pop_size <- Traits_2024 |>
  group_by(Treatment, Generation) |>
  summarise(n= n())
#write.table(pop_size, file.path(output_directory,"pop_size.txt"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

Traits_2024 <- Traits_2024 |> mutate(Treatment = recode(Treatment, 
                                                        "G1"="G1",    
                                                        "THB"="ABC",  
                                                        "THH"="aBC",  
                                                        "TNHB"="AbC", 
                                                        "TNHH"="abC",
                                                        "LHB"="ABc",  
                                                        "LHH"="aBc",   
                                                        "LNHB"="Abc", 
                                                        "LNHH"="abc",  
))

unique(Traits_2024$Treatment)
# create correlation matrix
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

# Get unique treatments and generations
all_treatments <- c("G1", "ABC", "aBC", "AbC", "abC", "ABc", "aBc", "Abc", "abc")
subplot_labels <- c("A", "B", "C","D", "E", "F", "G", "H", "I")

layout_matrix <- matrix(
  c(1, 2, 3,
    4, 5, 6,
    7, 8, 9), 
  nrow = 3, byrow = TRUE
)


#pdf(file = file.path(output_directory, "all_Trt_networks_BLv.pdf"), 
#    width = 11, height = 10, bg = "white")
png(filename = file.path(output_directory, "all_Trt_networks_BLv.png"),
    width = 11, height = 10, units = "in", res = 300, bg = "white")

# Set up layout
layout(layout_matrix)
par(mar = c(1, 1, 2.5, 1) + 0.1, oma = c(0.5, 0.2, 1, 0.5))

# Define colors
category_colors <- c("A" = '#FF8B8B', "B" = '#00D7FF', "C" = "#F7FD04", "D" = '#B5FE83')

# Loop through treatments (first one in the center)
for (i in seq_along(all_treatments)) {
  
  TrtSingle <- Traits_2024 |> 
    filter(Treatment == all_treatments[i]) |> 
    select(8:49)
  
  if (nrow(TrtSingle) > 1) {
    # Calculate correlation matrix
    TrtSingle_r <- get_corr_matrix(TrtSingle)
    TrtSingle_01matrix <- ifelse(TrtSingle_r != 0, 1, 0) 
    
    # Build network graph
    TrtSingle_graph <- graph_from_adjacency_matrix(TrtSingle_01matrix, mode = "undirected", weighted = NULL, diag = FALSE)
    
    # Best Louvain clustering
    TrtSingle_result <- run_best_louvain(TrtSingle_graph, n_runs = 100)
    TrtSingle_Louvain <- TrtSingle_result$clustering
    TrtSingle_Modu_Louvain <- round(TrtSingle_result$modularity, 3)
    
    # Set node colors by category
    V(TrtSingle_graph)$category <- c(rep("A", 10), rep("B", 14), rep("C", 8), rep("D", 10))
    #category_colors <- c("A" = '#FF8B8B', "B" = '#00D7FF', "C" = "#F7FD04")
    V(TrtSingle_graph)$color <- category_colors[V(TrtSingle_graph)$category]
    
    
    set.seed(37)
    plot(TrtSingle_Louvain, TrtSingle_graph, 
         mark.col = NA, mark.border = "grey", mark.lwd = 2,  
         vertex.label = NA, 
         vertex.size = 9, 
         vertex.label.cex = 1.05, 
         vertex.frame.color = V(TrtSingle_graph)$color, 
         vertex.frame.width = 5, 
         shapes(shape = "circle"))
    
    if (i == 1) {
      title(main = bquote(.(as.character(all_treatments[i])) *": " * 
                            italic(M) == .(TrtSingle_Modu_Louvain)), cex.main = 2, col.main = "red")
    } else {
      title(main = bquote("G10 " * "(" * .(as.character(all_treatments[i])) * ")" *
                            ": " * italic(M) == .(TrtSingle_Modu_Louvain)), cex.main = 2, col.main = "black")
    }
    
    mtext(subplot_labels[i], side = 3, line = 0, adj = 0, cex = 2, font = 2)
  }
}
#subplot_labels
# Reset layout and close device
layout(1)
dev.off()




# with numbers 

#pdf(file = file.path(output_directory, "all_Trt_networks_BLvSeq.pdf"), 
#    width = 11, height = 10, bg = "white")
png(filename = file.path(output_directory, "all_Trt_networks_BLvSeq.png"),
    width = 11, height = 10, units = "in", res = 300, bg = "white")

# Set up layout
layout(layout_matrix)
par(mar = c(1, 1, 2.5, 1) + 0.1, oma = c(0.5, 0.2, 1, 0.5))

# Define colors
category_colors <- c("A" = '#FF8B8B', "B" = '#00D7FF', "C" = "#F7FD04", "D" = '#B5FE83')

# Loop through treatments (first one in the center)
for (i in seq_along(all_treatments)) {
  
  TrtSingle <- Traits_2024 |> 
    filter(Treatment == all_treatments[i]) |> 
    select(8:49)
  
  if (nrow(TrtSingle) > 1) {
    # Calculate correlation matrix
    TrtSingle_r <- get_corr_matrix(TrtSingle)
    TrtSingle_01matrix <- ifelse(TrtSingle_r != 0, 1, 0) 
    
    # Build network graph
    TrtSingle_graph <- graph_from_adjacency_matrix(TrtSingle_01matrix, mode = "undirected", weighted = NULL, diag = FALSE)
    
    # Best Louvain clustering
    TrtSingle_result <- run_best_louvain(TrtSingle_graph, n_runs = 100)
    TrtSingle_Louvain <- TrtSingle_result$clustering
    TrtSingle_Modu_Louvain <- round(TrtSingle_result$modularity, 3)
    

    # Set node colors by category
    V(TrtSingle_graph)$category <- c(rep("A", 10), rep("B", 14), rep("C", 8), rep("D", 10))
    V(TrtSingle_graph)$color <- category_colors[V(TrtSingle_graph)$category]
    set.seed(37)
    #l <- layout_with_fr(TrtSingle_graph)
    plot(TrtSingle_Louvain, TrtSingle_graph, 
         mark.col = NA, mark.border = "grey", mark.lwd = 2,  
         vertex.label = c(1:42), 
         vertex.size = 15, 
         vertex.label.cex = 0.75, 
         vertex.frame.color = V(TrtSingle_graph)$color, 
         vertex.frame.width = 2.5, 
         shapes(shape = "circle"))
    
    if (i == 1) {
      title(main = bquote(.(as.character(all_treatments[i])) *": " * 
                            italic(M) == .(TrtSingle_Modu_Louvain)), cex.main = 2, col.main = "red")
    } else {
      title(main = bquote("G10 " * "(" * .(as.character(all_treatments[i])) * ")" *
                            ": " * italic(M) == .(TrtSingle_Modu_Louvain)), cex.main = 2, col.main = "black")
    }
    
    mtext(subplot_labels[i], side = 3, line = 0, adj = 0, cex = 2, font = 2)
  }
}
#subplot_labels
# Reset layout and close device
layout(1)
dev.off()




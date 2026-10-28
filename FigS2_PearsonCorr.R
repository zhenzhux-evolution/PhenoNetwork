#!/usr/bin/env Rscript

# 2026-03-05 16:01:55
rm(list=ls())
library(foreach)
library(dplyr)
library(psych)
library(igraph)
library(bestNormalize)
library(ggplot2)
library(tidyr)


setwd("/home/zhenzhux/pheno")
output_directory = "/home/zhenzhux/pheno"
phenotype <- read.table("./Dorey_2024.txt", header = TRUE)


#colnames(phenotype)
Traits_2024 <- phenotype |>
  select(1:7, 12:21, 27:40, 51:58, 44, 47, 49, 50, 59:64)
colnames(Traits_2024[8:49])

colnames(Traits_2024)
Traits_2024 |>
  group_by(Treatment, Generation, Replicate) |>
  summarise(n= n())

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

Traits_2024_Norm <- as.data.frame(apply(Traits_2024[8:49], 2, function(x) bestNormalize(x, loo = TRUE, na.rm = TRUE)$x.t))
Traits_2024_Norm_new <- bind_cols(Traits_2024[1:7], Traits_2024_Norm)



colnames(Traits_2024)[8:49] <- paste0("Trait", sprintf("%02d", 1:42))
hist(Traits_2024_Norm_new$Glucoraphanin)
Traits_long <- Traits_2024 %>% pivot_longer(
  cols = 8:49, names_to = "traits", values_to = "value"
)
plot <- ggplot(Traits_long, aes(value)) + 
  geom_histogram(bins = 30, fill= "grey90", color = "black") + 
  facet_wrap(~traits, scales = "free", ncol = 7) + 
  labs(x = "", y = "Count") +
  theme_bw()
ggsave("Trait_distribution.png", path = output_directory, plot = plot, width = 14, height = 10, dpi = 300)



get_corr_matrix <- function(data) {
  res <- corr.test(data, use = "pairwise", method = "pearson", adjust = "none", ci = FALSE)
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
all_treatments <- c("G0", "ABC", "aBC", "AbC", "abC", "ABc", "aBc", "Abc", "abc")
subplot_labels <- c("A", "B", "C","D", "E", "F", "G", "H", "I")

layout_matrix <- matrix(
  c(1, 2, 3,
    4, 5, 6,
    7, 8, 9), 
  nrow = 3, byrow = TRUE
)
# with numbers 
pdf(file = file.path(output_directory, "CorrelationTest_M_Seq.pdf"), 
    width = 11, height = 10, bg = "white")

# Set up layout
layout(layout_matrix)
par(mar = c(1, 1, 2.5, 1) + 0.1, oma = c(0.5, 0.2, 1, 0.5))

# Define colors
category_colors <- c("A" = '#FF8B8B', "B" = '#00D7FF', "C" = "#F7FD04", "D" = '#B5FE83')

# Loop through treatments (first one in the center)
for (i in seq_along(all_treatments)) {
  
  TrtSingle <- Traits_2024_Norm_new |> 
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
         vertex.label.cex = 0.65, 
         vertex.frame.color = V(TrtSingle_graph)$color, 
         vertex.frame.width = 3, 
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

layout(1)
dev.off()




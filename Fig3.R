
###########
######### Fig. 3. simulation and a simulated plot ########
###########
rm(list=ls())
library(dplyr)
library(psych)
library(igraph)
require(doParallel)
num_cores <- detectCores()
doParallel::registerDoParallel(cores = num_cores)

setwd("/home/zhenzhux/pheno")
output_directory = "/home/zhenzhux/pheno"
phenotype <- read.table("./Dorey_2024.txt", header = TRUE)


#colnames(phenotype)
Traits_2024 <- phenotype |>
  select(1:7, 12:21, 27:40, 51:58, 44, 47, 49, 50, 59:64)
colnames(Traits_2024[8:49])

AbC_simuModu <- read.table("./edgeSimu_Louvain_gainloss_AbC.txt", header=TRUE)

unique(Traits_2024$Treatment);unique(Traits_2024$Generation)

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


G1 <- Traits_2024 |> filter(Treatment == "G1" & Generation == "G1") %>% select(8:49)
G10 <- Traits_2024 |> filter(Treatment == "TNHB" & Generation == "G10") %>% select(8:49)
# Compute matrices
group1_r <- get_corr_matrix(G1)
group10_r <- get_corr_matrix(G10)

# Binary 0/1 matrices
group1_r_bin <- ifelse(group1_r != 0, 1, 0)
group10_r_bin <- ifelse(group10_r != 0, 1, 0)

isSymmetric.matrix(group1_r_bin)
isSymmetric.matrix(group10_r_bin)

#network construction 
G1_graph <- graph_from_adjacency_matrix(group1_r_bin, mode = "undirected", weighted = NULL, diag = FALSE)
G10_graph <- graph_from_adjacency_matrix(group10_r_bin, mode = "undirected", weighted = NULL, diag = FALSE)

G1_result <- run_best_louvain(G1_graph, n_runs = 100)
G1_opt <- G1_result$clustering
G1_Modu_opt <- round(G1_result$modularity,3)

G10_result <- run_best_louvain(G10_graph, n_runs = 100)
G10_opt <- G10_result$clustering
G10_Modu_opt <- round(G10_result$modularity,3)


# M simulations
G1_upper <- group1_r_bin[upper.tri(group1_r_bin)];sum(G1_upper)
G10_upper <- group10_r_bin[upper.tri(group10_r_bin)];sum(G10_upper)
upper_idx <- which(upper.tri(group1_r_bin), arr.ind = TRUE)
G1_idx <- which((G1_upper == 1)); G10_idx <- which((G10_upper == 1)); 
shared_idx <- which((G1_upper == 1) & (G10_upper == 1))      
G1_only_idx <- which((G1_upper == 1) & (G10_upper == 0))
G10_only_idx <- which((G1_upper == 0) & (G10_upper == 1))

n_shared <- length(shared_idx)
n_gain <- length(G10_only_idx);#sum((G1_upper == 0) & (G10_upper == 1))
n_loss <- sum((G1_upper == 1) & (G10_upper == 0))

# Randomly sample loss edges from G1-only
loss_idx_sim <- sample(G1_idx, n_loss, replace = FALSE);length(loss_idx_sim)
after_loss_sim <- setdiff(G1_idx, loss_idx_sim);
length(after_loss_sim) == n_shared
#loss_idx_sim %in% G1_idx

# Randomly sample gain edges from G1-zero 
all_possible_idx <- which(G1_upper == 0) #& (G10_upper == 0)
gain_idx_sim <- sample(all_possible_idx, n_gain, replace = FALSE)

Matrix_simu <- matrix(0, nrow = nrow(group10_r_bin), ncol = ncol(group10_r_bin)) 
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


Matrix_simu <- add_edge(Matrix_simu, after_loss_sim)
Matrix_simu <- add_edge(Matrix_simu, gain_idx_sim)

rownames(Matrix_simu) <- rownames(group10_r_bin)
colnames(Matrix_simu) <- colnames(group10_r_bin)
isSymmetric(Matrix_simu)
sum(Matrix_simu[upper.tri(Matrix_simu)]) == sum(G10_upper)

simu_graph <- graph_from_adjacency_matrix(Matrix_simu, mode = "undirected", weighted = NULL, diag = FALSE)

simu_result <- run_best_louvain(simu_graph, n_runs = 100)
simu_opt <- simu_result$clustering
simu_Modu_opt <- round(simu_result$modularity,3)
memb_simu <- simu_result$membership
simu_Modu_opt 

# color of G1 and G10
category_colors <- c("A" = '#FF8B8B', "B" = '#00D7FF', "C" = "#F7FD04", "D" = "#90EE90")

V(G1_graph)$category <- c(rep("A", 10), rep("B", 14), rep("C", 8), rep("D", 10))
V(G1_graph)$color <- category_colors[V(G1_graph)$category]
V(G10_graph)$category <- c(rep("A", 10), rep("B", 14), rep("C", 8), rep("D", 10))
V(G10_graph)$color <- category_colors[V(G10_graph)$category]


# color of simulated matrix 
V(simu_graph)$category <- c(rep("A", 10), rep("B", 14), rep("C", 8), rep("D", 10))
V(simu_graph)$color <- category_colors[V(simu_graph)$category]




pdf(file = file.path(output_directory, paste0("G1_AbC", "_Simu.pdf")), 
    width = 11, height = 11, bg = "white")
par(mfrow = c(2, 2),  # 2x2 grid
    mar = c(4, 4, 4, 0) + 0.1,  # Margins: bottom, left, top, right
    oma = c(0, 0, 0.5, 0))  # Outer margins for spacing
# Plot G1
set.seed(37) #33, 40
plot(G1_opt, G1_graph, 
     mark.col = NA, mark.border = "grey", mark.lwd = 2, 
     vertex.label = c(1:42), 
     vertex.size = 15, 
     vertex.label.cex = 1.05, 
     vertex.frame.color = V(G1_graph)$color, 
     vertex.frame.width = 5, 
     cex.main = 5, 
     shapes(shape = "circle"))
title(main = bquote("G1: " * italic(M)== .(G1_Modu_opt)), cex.main = 2, line = 0.5)  
mtext("A", side = 3, line = 0, adj = 0, cex = 2, font = 2) 

# Plot G10
plot(G10_opt, G10_graph, 
     mark.col = NA, mark.border = "grey", mark.lwd = 2,  
     vertex.label = c(1:42), 
     vertex.size = 15, 
     vertex.label.cex = 1.05, 
     vertex.frame.color = V(G10_graph)$color, 
     vertex.frame.width = 5,
     shapes(shape = "circle"))
title(main = bquote("G10 (AbC): " * italic(M)== .(G10_Modu_opt)), cex.main = 2, line = 0.5)  
mtext("B", side = 3, line = 0, adj = 0, cex = 2, font = 2)

# Plot simulation
plot(simu_opt, simu_graph, 
     mark.col = NA, mark.border = "grey", mark.lwd = 2,  
     vertex.label = c(1:42), 
     vertex.size = 15, 
     vertex.label.cex = 1.05, 
     vertex.frame.color = V(simu_graph)$color, 
     vertex.frame.width = 5, 
     shapes(shape = "circle"))
title(main = bquote(atop("A simulated PCN at G10: ", italic(M) == .(simu_Modu_opt))), cex.main = 2)  
mtext("C", side = 3, line = 0, adj = 0, cex = 2, font = 2) 

# Histogram for simulated modularity distribution
hist(AbC_simuModu$Modu_g2_perm,
     col = gray.colors(5, alpha = 0.1), # light fill to match node colors
     border = "black",                  # clean edges
     main = "Simulated             Oberseved",
     xlab = bquote(italic(M)),
     ylab = "Frequency",
     xlim =c(0.25,0.75),
     cex.main = 2, font.main =1,                    # match title size from network plot
     cex.lab = 1.3,                      # axis label size
     cex.axis = 1.1                      # tick label size
)
abline(v = G10_Modu_opt, col = "red", lwd = 2, lty = 2)
mtext("D", side = 3, line = 0, adj = 0, cex = 2, font = 2)

par(xpd = NA)
legend(x = 0.04, y = 275,
       legend = c("Leaf glucosinolate (1-10)", "Floral volatile (11-24)", 
                  "Floral organ morphology (25-32)", "Non-floral morphology (33-42)"), 
       col = unique(V(G1_graph)$color), 
       pch = 21, pt.bg = unique(V(G1_graph)$color), 
       pt.cex = 2, cex = 1.2, box.lty = 1,
       title = "Trait category")
dev.off()


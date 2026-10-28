rm(list=ls())
library(data.table)
library(dplyr)
library(purrr)
library(psych)

setwd("~/Documents/demon/UZH/DataSci_UM/network_genomics/all_pheno/gwas/PolygScore/polyScore42/trt_trait_PGS/")
out_dir="~/Documents/demon/UZH/DataSci_UM/network_genomics/all_pheno/gwas/PolygScore/polyScore42/trt_trait_PGS/"
out_dir="~/Documents/demon/UZH/DataSci_UM/network_genomics/Script4Publication"


# Treatments to compare against G1_9
treatments <- c("ABC_1", "aBC_2", "AbC_3", "abC_4",
                "ABc_5", "aBc_6", "Abc_7", "abc_8")

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

base_path <- "~/Documents/demon/UZH/DataSci_UM/network_genomics/all_pheno/gwas/PolygScore/polyScore42/trt_trait_PGS"
phenotype <- read.table(
  "~/Documents/demon/UZH/DataSci_UM/network_genomics/all_pheno/gwas/PolygScore/allTrait.pheno",
  header = TRUE
)

################################
## Load G1_9 once
## -----------------------------
path_G1_9 <- file.path(base_path, "G1_9_gwas_LD0.2")
files_G1_9 <- list.files(
  path = path_G1_9, 
  pattern = "_G1_9_LD0\\.2_PGS\\.summary$",
  full.names = TRUE
)

G1_9_LD0.2 <- do.call(rbind, lapply(files_G1_9, function(f) {
  df <- read.table(f, header = TRUE)
  df$Trait <- sub("_G1_9_LD0\\.2_PGS\\.summary$", "", basename(f))
  df
}))

sig_G1_9_LD0.2 <- G1_9_LD0.2 %>% filter(P < 0.05)
sig_traits_G1_9 <- sig_G1_9_LD0.2$Trait


G1_9_sig_traits <- map(sig_traits_G1_9, function(f) {
  read.table(file.path(path_G1_9,
                       paste0(f, "_G1_9_LD0.2_PGS.best")),
             header = TRUE) %>%
    select(FID, IID, PRS) %>%
    rename(!!f := PRS)
}) %>%
  reduce(full_join, by = c("FID", "IID"))

phenotype_raw_G1_9 <- phenotype %>%
  semi_join(G1_9_sig_traits, by = c("FID", "IID")) %>%
  select(FID, IID, all_of(sig_traits_G1_9))


## -----------------------------
## Loop over all treatments
## -----------------------------
results <- list()
perm_list <- list()
for(trt in treatments){
  
  message("Processing ", trt)
  
  ## Summary files
  trt_path <- file.path(base_path, paste0(trt, "_gwas_LD0.2"))
  
  files <- list.files(
    path = trt_path,
    pattern = paste0("_", trt, "_LD0\\.2_PGS\\.summary$"),
    full.names = TRUE
  )
  
  trt_summary <- map_dfr(files, function(f){
    df <- read.table(f, header = TRUE)
    df$Trait <- sub(
      paste0("_", trt, "_LD0\\.2_PGS\\.summary$"),
      "",
      basename(f)
    )
    df
  })
  
  sig_traits <- trt_summary %>%
    filter(P < 0.05) %>%
    pull(Trait)
  
  ## Predicted PRS
  trt_sig_traits <- map(sig_traits, function(f){
    read.table(file.path(trt_path,
                         paste0(f, "_", trt, "_LD0.2_PGS.best")),
               header = TRUE) %>%
      select(FID, IID, PRS) %>%
      rename(!!f := PRS)
  }) %>%
    reduce(full_join, by = c("FID", "IID"))
  
  ## Raw phenotype
  phenotype_raw_trt <- phenotype %>%
    semi_join(trt_sig_traits, by = c("FID","IID")) %>%
    select(FID, IID, all_of(sig_traits))
  
  ## Shared traits
  shared_traits <- intersect(sig_traits_G1_9, sig_traits)
  
  ## Skip if too few traits
  if(length(shared_traits) < 9){
    warning(trt, ": fewer than 9 shared traits.")
    next
  }
  
  ## Reduced matrices
  G1_pred <- G1_9_sig_traits %>%
    select(FID, IID, all_of(shared_traits)) %>%
    arrange(FID, IID)
  
  TRT_pred <- trt_sig_traits %>%
    select(FID, IID, all_of(shared_traits)) %>%
    arrange(FID, IID)
  
  G1_raw <- phenotype_raw_G1_9 %>%
    select(FID, IID, all_of(shared_traits)) %>%
    arrange(FID, IID)
  
  TRT_raw <- phenotype_raw_trt %>%
    select(FID, IID, all_of(shared_traits)) %>%
    arrange(FID, IID)
  
  ## Correlation matrices
  r_G1_pred <- get_corr_matrix(G1_pred[, shared_traits])
  r_TRT_pred <- get_corr_matrix(TRT_pred[, shared_traits])
  r_pred_change <- r_TRT_pred - r_G1_pred
  
  
  r_G1_raw <- get_corr_matrix(G1_raw[, shared_traits])
  r_TRT_raw <- get_corr_matrix(TRT_raw[, shared_traits])
  r_raw_change <- r_TRT_raw - r_G1_raw
  
  ## Binary
  b_G1_pred <- (r_G1_pred != 0)
  b_TRT_pred <- (r_TRT_pred != 0)
  
  b_G1_raw <- (r_G1_raw != 0)
  b_TRT_raw <- (r_TRT_raw != 0)
  
  ## Gain/Loss
  gain_pred <- !b_G1_pred & b_TRT_pred
  gain_raw  <- !b_G1_raw & b_TRT_raw
  gain_pred_idx <- which(gain_pred & upper.tri(gain_pred), arr.ind = TRUE)
  gain_raw_idx <- which(gain_raw & upper.tri(gain_raw), arr.ind = TRUE)
  
  loss_pred <- b_G1_pred & !b_TRT_pred
  loss_raw  <- b_G1_raw & !b_TRT_raw
  loss_pred_idx <- which(loss_pred & upper.tri(loss_pred), arr.ind = TRUE)
  loss_raw_idx <- which(loss_raw & upper.tri(loss_raw), arr.ind = TRUE)
  
  ## Predicted edge gain and loss in one matrix
  edge_pred <- matrix(
    0,
    nrow = nrow(b_G1_pred),
    ncol = ncol(b_G1_pred),
    dimnames = dimnames(b_G1_pred)
  )
  edge_pred[gain_pred] <- 1
  edge_pred[loss_pred] <- -1
  
  ## Raw edge status
  edge_raw <- matrix(
    0,
    nrow = nrow(b_G1_raw),
    ncol = ncol(b_G1_raw),
    dimnames = dimnames(b_G1_raw)
  )
  
  edge_raw[gain_raw] <- 1
  edge_raw[loss_raw] <- -1
  
  write.table(edge_pred, file = file.path(out_dir, paste0(trt, "_pred_change_signs.txt")),
              sep = "\t",quote = FALSE, col.names = NA)
  write.table(edge_raw, file = file.path(out_dir, paste0(trt, "_raw_change_signs.txt")),
              sep = "\t", quote = FALSE, col.names = NA)
  
  
  
  upper <- upper.tri(edge_pred)
  pred_vec <- edge_pred[upper]; raw_vec <- edge_raw[upper]
  
  ## Observed number of shared gain/loss
  observed <- sum(pred_vec == raw_vec & pred_vec != 0)
  
  
  perm_list[[trt]] <- list(
    pred_vec = pred_vec, 
    raw_vec = raw_vec, 
    observed = observed
  )
  
}

# ## Permutation test
# perm <- replicate(1000, {
#   ## Randomly permute the raw edge labels
#   shuffled <- sample(raw_vec)

#   ## Count total shared gain/loss
#   sum(pred_vec == shuffled & pred_vec != 0)
# })
# p_value <- mean(perm >= observed)
# mean_perm <- mean(perm)
# sd_perm <- sd(perm)

#   results[[trt]] <- list(
#     shared_traits = shared_traits,
#     observed = observed,
#     p_value = p_value,
#     mean_perm_shared = mean_perm,
#     sd_perm_shared = sd_perm
#   )


observed_sum <- sum(sapply(perm_list, function(x) x$observed))
perm_sum <- replicate(100000, {
  
  total <- 0
  
  for (trt in names(perm_list)) {
    
    pred_vec <- perm_list[[trt]]$pred_vec
    raw_vec  <- perm_list[[trt]]$raw_vec
    
    shuffled <- sample(raw_vec)
    
    total <- total + sum(pred_vec == shuffled & pred_vec != 0)
  }
  
  total
})
p_value_sum <- mean(perm_sum >= observed_sum);p_value_sum
mean_perm_sum <- mean(perm_sum);mean_perm_sum
sd_perm_sum <- sd(perm_sum);sd_perm_sum

png(filename = file.path(out_dir, "Perm_concordancy_100k.png"),
    width = 3, height = 4, units = "in", res = 600)
hist(perm_sum, main = "", xlab = "# of concordant changes")
abline(v = observed_sum, col = "blue", lty = 2, lwd = 2)
dev.off()




summary_df <- bind_rows(lapply(names(results), function(trt){
  
  r <- results[[trt]]
  
  data.frame(
    Treatment = trt,
    n_shared_traits = length(r$shared_traits),
    observed_match = r$observed,
    expected = r$mean_perm_shared,
    sd = round(r$sd_perm_shared,3),
    p_value = r$p_value
  )
}))
write.table(
  summary_df,
  "permutation_test_summary.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)




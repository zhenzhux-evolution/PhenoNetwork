
# gwas output 
#gemma results visualization
rm(list=ls())
library(ggplot2)
library(qqman)
library(ggpubr)
library(data.table)
library(grid)
library(ggplotify)

setwd("/Users/UZH/DataSci_UM/network_genomics/all_pheno")
out_dir="/Users/UZH/DataSci_UM/network_genomics"

n_traits <- 28

trait_name <- c("Trait33", "Trait34", "Trait37", "Trait38", "Trait39", "Trait40", 
                "Trait25", "Trait11", "Trait13", "Trait15", "Trait16", "Trait17", 
                "Trait18", "Trait19", "Trait20", "Trait22", "Trait23", "Trait24")

foreach(i = 1:n_traits) %dopar% {
  trait_name <- trait_name[i]
  file_path <- file.path(".", paste0(trait_name, ".assoc.txt"))
  
  data <- fread(file_path)
  data[, chr := as.integer(sub("^A0*([0-9]+):.*", "\\1", rs))]
  data$fdr <- p.adjust(data$p_wald, method = 'BH')
  
  # significant snps
  sig_data <- data[fdr < 0.05]
  if (nrow(sig_data) > 0) {
    write.table(sig_data, file.path(out_dir, paste0("sig_", trait_name, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  data_plot <- data[fdr < 0.99]
  highlight_snps <- sig_data$rs 
  
  # manhatten plot
  plot_file <- file.path(out_dir, paste0(trait_name, ".tiff"))
  tiff(plot_file, res = 300, width = 15, height = 10, units = "cm")
  manhattan(
    data_plot, chr = 'chr', bp="ps", snp="rs", p='p_wald', logp=TRUE, ylab="-log10(p)", 
    col = c("grey30", "grey80"), 
    suggestiveline=FALSE, 
    genomewideline=FALSE,
    highlight = highlight_snps,
    cex =0.1, 
    main=paste0(trait_name) 
  )
  dev.off()
  
  # QQ plot 
  qq_file <- file.path(out_dir, paste0(trait_name, "_qq",".tiff"))
  tiff(qq_file, res = 300, width = 10, height = 10, units = "cm")
  qq(data$p_wald, main=paste0("QQ plot of", trait_name), cex =0.2)
  dev.off()
}



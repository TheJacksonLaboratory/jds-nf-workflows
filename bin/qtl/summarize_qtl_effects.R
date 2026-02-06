#!/usr/bin/env Rscript

################################################################################
# Summarize QTL effects
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20251113
################################################################################

library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# Get project ID
id <- args[1]

# Wrangle genome coordinates
primary_chrom_bed <- args[2]
primary_chrom_bed <- read.table(primary_chrom_bed)
colnames(primary_chrom_bed) <- c("chrom","start","stop")
chrom_levels <- c(as.character(1:19), "X")
primary_chrom_bed <- primary_chrom_bed %>%
  dplyr::filter(chrom %in% chrom_levels) %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_levels))

# Read in peak files
peak_files <- list.files(pattern = "peaks.csv")

# Bind all the peaks together and reduce
all_peaks <- Reduce(dplyr::bind_rows,lapply(peak_files, function(x) read.csv(x) %>%
                                              dplyr::mutate(chrom = factor(as.character(chrom), levels = chrom_levels))))
peaks_viewer <- all_peaks %>%
  dplyr::select(phenotype, chrom, pos_peak, LOD, LETTERS[1:8]) %>%
  dplyr::mutate(marker.id = paste(chrom,as.integer(pos_peak*1e6),sep = "_")) %>%
  dplyr::select(phenotype, marker.id, LOD, LETTERS[1:8])
write.csv(peaks_viewer, paste(id,"peaks_viewer_file.csv", sep = "_"), row.names = FALSE, quote = FALSE)
write.csv(all_peaks, paste(id,"all_peaks_file.csv", sep = "_"), row.names = FALSE, quote = FALSE)

# QTL heatmap
plot_peaks <- all_peaks %>%
  dplyr::distinct(chrom, phenotype) %>%
  dplyr::left_join(., primary_chrom_bed) %>%
  dplyr::left_join(., all_peaks)

hm <- ggplot(plot_peaks) +
  theme_bw() +
  geom_point(aes(x = pos_peak, y = phenotype, colour = LOD), na.rm = TRUE) +
  geom_segment(aes(x = pos_start, y = phenotype, xend = pos_end, yend = phenotype, colour = LOD)) + 
  geom_segment(aes(x = 0, y = phenotype, xend = stop/1e6, yend = phenotype), size = 2.5, alpha = 0) + 
  facet_grid(. ~ chrom, scales = "free_x") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,200,50)) + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank()) + 
  labs(x = "Genomic position (Mb)", y = "Trait")
ggsave(plot = hm, filename = paste(id,"qtl_effects_heatmap.png", sep = "_"), width = 10, height = 7)

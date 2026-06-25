library(dplyr)

setwd("/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/compare_genomes/TEs/")

# h2 <- read_tsv("austro_to_JGV016_h2_extract_austroID_pIdent_austroAlignedLength.txt", col_names = F)
# colnames(h2) <- c('id', 'pIdent', 'alnLength')
# 
# maxAln_h2 <- h2 %>%
#   group_by(id) %>%
#   filter(alnLength == max(alnLength))
# 
# mean(maxAln_h1$pIdent)
# # 94.35
# mean(maxAln_h1$alnLength)
# # 8394.24
# length(unique(maxAln_h1$id))
# # 9606
# 
# mean(maxAln_h2$pIdent)
# # 95.51
# mean(maxAln_h2$alnLength)
# # 8929.589
# length(unique(maxAln_h2$id))
# # 9605


#hist(maxAln_h1$pIdent, breaks = 50)
#hist(maxAln_h2$pIdent, breaks = 50)





# library
library(ggplot2)
library(viridis)

meta <- read.table("TE_content_metadata.tsv", header = T)
diff <- data.frame(accession=meta$accession,
                   genome=meta$genomeSize1-meta$genomeSize2,
                   TE=meta$TEsize1-meta$TEsize2,
                   gene=meta$geneSize1-meta$geneSize2,
                   Tyba=meta$TybaSize1-meta$TybaSize2)

wide_to_long <- function(df) {
  data.frame(
    accession = rep(df[[1]], 2),
    group = c(rep("TE", nrow(df)),
              rep("other", nrow(df))),
    value = c(df[[3]], df[[4]])
  )
}


# plot
ggplot(wide_to_long(diff), aes(fill=group, y=value, x=accession)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_viridis(discrete=TRUE, name="") +
  theme_bw() +
  ylab("Length (bp)") + 
  xlab("Accession")

ggsave("TE_content_difference_barplot.pdf", width = 8, height = 4)

# 
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# Read the data
data <- read.table("TE_content_metadata.tsv", header = TRUE)

# Calculate differences and contributions
data_diff <- data %>%
  mutate(
    genome_diff = genomeSize1 - genomeSize2,
    TE_diff = TEsize1 - TEsize2,
    gene_diff = geneSize1 - geneSize2,
    Tyba_diff = TybaSize1 - TybaSize2,
    
    # Calculate contribution of each element to genome size difference
    TE_contribution = TE_diff / genome_diff,
    gene_contribution = gene_diff / genome_diff,
    Tyba_contribution = Tyba_diff / genome_diff,
    
    # Calculate "other" category (elements not tracked)
    other_diff = genome_diff - (TE_diff + gene_diff + Tyba_diff),
    other_contribution = other_diff / genome_diff
  )

pdf("visual_test.pdf", 8, 5)

############# plot 1 ##############
# visualization: Show absolute differences
abs_plot_data <- data_diff %>%
  select(accession, TE_diff, gene_diff, Tyba_diff, other_diff) %>%
  pivot_longer(cols = -accession, names_to = "element", values_to = "difference") %>%
  mutate(
    element = case_when(
      element == "TE_diff" ~ "Transposable Elements",
      element == "gene_diff" ~ "Genes", 
      element == "Tyba_diff" ~ "Tyba Repeats",
      element == "other_diff" ~ "Other Elements"
    ),
    element = factor(element, levels = c("Transposable Elements", "Genes", "Tyba Repeats", "Other Elements"))
  )

# Absolute differences plot
#ggplot(abs_plot_data, aes(x = difference/1000000, y = fct_reorder(accession, difference, .fun = sum), fill = element)) +
ggplot(abs_plot_data, aes(x = difference/1000000, y = accession, fill = element)) +
  geom_col(position = "stack") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Transposable Elements" = "orange",
                               "Genes" = "steelblue",
                               "Tyba Repeats" = "darkgreen", 
                               "Other Elements" = "#999999")) +
  labs(
    title = "Absolute Size Differences Between Haplotypes by Genomic Element",
    subtitle = "Positive values = larger in haplotype 1, Negative values = larger in haplotype 2",
    x = "Size Difference (Mbp)",
    y = "Accession",
    fill = "Genomic Element"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15),
    plot.subtitle = element_text(size = 10, color = "#3E3E3E"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom"
  )



############# plot 2 ##############
# Prepare data for plotting
plot_data <- data_diff %>%
  select(accession, TE_contribution, gene_contribution, Tyba_contribution, other_contribution) %>%
  pivot_longer(cols = -accession, names_to = "element", values_to = "contribution") %>%
  mutate(
    element = case_when(
      element == "TE_contribution" ~ "Transposable Elements",
      element == "gene_contribution" ~ "Genes",
      element == "Tyba_contribution" ~ "Tyba Repeats",
      element == "other_contribution" ~ "Other Elements"
    ),
    element = factor(element, levels = c("Transposable Elements", "Genes", "Tyba Repeats", "Other Elements"))
  )

# Create the visualization
#ggplot(plot_data, aes(x = contribution, y = fct_reorder(accession, contribution, .fun = sum), fill = element)) +
ggplot(plot_data, aes(x = contribution, y = accession, fill = element)) +
  geom_col(position = "stack") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  scale_x_continuous(labels = scales::percent, 
                     breaks = seq(0, 1, 0.2),
                     limits = c(-0.1, 1.1)) +
  scale_fill_manual(values = c("Transposable Elements" = "orange",
                              "Genes" = "steelblue", 
                              "Tyba Repeats" = "darkgreen",
                              "Other Elements" = "#999999")) +
  labs(
    title = "Contribution of Genomic Elements to Haplotype Size Differences",
    subtitle = "Positive values = same as genome size difference, Negative values = opposite to genome size difference",
    x = "Contribution to Genome Size Difference (%)",
    y = "Accession",
    fill = "Genomic Element"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 15),
    plot.subtitle = element_text(size = 10, color = "#3E3E3E")
  )
dev.off()

# Print summary statistics
cat("Summary of contributions to genome size differences:\n")
summary_stats <- plot_data %>%
  group_by(element) %>%
  summarise(
    mean_contribution = mean(contribution),
    sd_contribution = sd(contribution),
    .groups = 'drop'
  )
print(summary_stats)
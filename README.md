# SNP-Analysis-Head-Neck-Cancer
Comprehensive SNP variant analysis and visualization of whole-exome Sequencing (WES) data from head and neck cancer using R. Open-source project using dataset SRR32633603
# Make Folder 
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Load Libraries
library(VariantAnnotation)
library(ggplot2)
library(dplyr)
library(GenomicRanges)

# Load VCF File
vcf_file <- "filtered_snps.vcf"
vcf <- readVcf(vcf_file, genome = "hg38")
head(rowRanges(vcf))
summary(rowRanges(vcf))
length(vcf)

# Extract Variant Info 
variants <- rowRanges(vcf)
variant_data <- as.data.frame(mcols(variants))
variant_data$chr <- as.character(seqnames(variants))
variant_data$start <- start(variants)

# === Plot 1: QUAL Distribution (All Variants) ===
p1 <- ggplot(variant_data, aes(x = QUAL)) +
  geom_histogram(binwidth = 50, fill = "red", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Variant Quality Scores", x = "Quality Score (QUAL)", y = "Count")
ggsave("plots/01_qual_distribution_all_variants.png", plot = p1, width = 8, height = 5, dpi = 300)

# === Filter PASS Variants ===
pass_variants <- variant_data[variant_data$FILTER == "PASS", ]
write.csv(pass_variants, "filtered_pass_variants.csv", row.names = FALSE)

# === Plot 2: QUAL Distribution (PASS Only) ===
p2 <- ggplot(pass_variants, aes(x = QUAL)) +
  geom_histogram(binwidth = 50, fill = "green", color = "black") +
  theme_minimal() +
  labs(title = "Quality Distribution of PASS Variants", x = "Quality Score (QUAL)", y = "Count")
ggsave("plots/02_qual_distribution_pass_variants.png", plot = p2, width = 8, height = 5, dpi = 300)

# === Plot 3: Variant Count per Chromosome ===
variant_df <- as.data.frame(variant_data)
variant_df$chr <- as.character(variant_df$chr)
chr_counts <- variant_df %>%
  group_by(chr) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
p3 <- ggplot(chr_counts, aes(x = reorder(chr, -count), y = count)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_minimal() +
  labs(title = "Variant Count per Chromosome", x = "Chromosome", y = "Number of Variants") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/03_variant_count_per_chr.png", plot = p3, width = 8, height = 5, dpi = 300)

# === Plot 4: Top 30 Chromosomes by Variant Count ===
chr_counts_top <- chr_counts %>% top_n(30, count)
p4 <- ggplot(chr_counts_top, aes(x = reorder(chr, -count), y = count)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_minimal() +
  labs(title = "Top 30 Chromosomes by Variant Count", x = "Chromosome", y = "Number of Variants") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 10))
ggsave("plots/04_top30_chr_variant_count.png", plot = p4, width = 10, height = 6, dpi = 300)

# === Plot 5: Histogram of SNP QUAL Scores ===
qual_scores <- as.data.frame(mcols(vcf))$QUAL
p5 <- ggplot(data.frame(QUAL = qual_scores), aes(x = QUAL)) +
  geom_histogram(fill = "skyblue", color = "black", bins = 50) +
  theme_minimal() +
  labs(title = "Distribution of SNP Quality Scores", x = "QUAL Score", y = "Frequency")
ggsave("plots/05_snp_qual_hist.png", plot = p5, width = 8, height = 5, dpi = 300)

# === Plot 6: Transitions vs. Transversions ===
ref <- as.character(mcols(vcf)$REF)
alt <- sapply(mcols(vcf)$ALT, function(x) as.character(x[[1]]))
snp_df <- data.frame(REF = ref, ALT = alt)
snp_df$type <- ifelse(
  (snp_df$REF == "A" & snp_df$ALT == "G") | (snp_df$REF == "G" & snp_df$ALT == "A") |
  (snp_df$REF == "C" & snp_df$ALT == "T") | (snp_df$REF == "T" & snp_df$ALT == "C"),
  "Transition", "Transversion"
)
snp_type_counts <- snp_df %>% group_by(type) %>% summarise(count = n())
p6 <- ggplot(snp_type_counts, aes(x = type, y = count, fill = type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "SNP Type: Transitions vs. Transversions", x = "Type", y = "Count") +
  scale_fill_manual(values = c("Transition" = "blue", "Transversion" = "orange"))
ggsave("plots/06_snp_type_transition_vs_transversion.png", plot = p6, width = 6, height = 5, dpi = 300)

# === Plot 7: QUAL vs Position by Chromosome ===
positions <- start(rowRanges(vcf))
chromosomes <- as.character(seqnames(rowRanges(vcf)))
qual <- as.data.frame(mcols(vcf))$QUAL
plot_df <- data.frame(CHR = chromosomes, POS = positions, QUAL = qual)
p7 <- ggplot(plot_df, aes(x = POS, y = QUAL)) +
  geom_point(alpha = 0.6, color = "darkgreen", size = 1) +
  facet_wrap(~ CHR, scales = "free_x") +
  theme_minimal() +
  labs(title = "SNP Quality Scores Across Chromosomes", x = "Position", y = "QUAL")
ggsave("plots/07_qual_vs_position_faceted.png", plot = p7, width = 10, height = 8, dpi = 300)

# === Plot 8: Heatmap of SNP Density per 1Mb ===
df <- data.frame(chr = chromosomes, pos = positions)
df$bin <- floor(df$pos / 1e6)
bin_counts <- df %>%
  group_by(chr, bin) %>%
  summarise(SNP_count = n())
p8 <- ggplot(bin_counts, aes(x = bin, y = chr, fill = SNP_count)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap of SNP Density per 1Mb Region", x = "Genomic Bin (Mb)", y = "Chromosome") +
  theme_minimal()
ggsave("plots/08_snp_density_heatmap.png", plot = p8, width = 10, height = 6, dpi = 300)

# === Plot 9: Manhattan-like Plot ===
plot_df$CHR <- factor(plot_df$CHR, levels = unique(plot_df$CHR))
plot_df <- plot_df %>%
  group_by(CHR) %>%
  mutate(CHR_num = as.numeric(as.factor(CHR)),
         global_pos = POS + (CHR_num - 1) * 1e7)
p9 <- ggplot(plot_df, aes(x = global_pos, y = QUAL, color = CHR)) +
  geom_point(size = 1, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Manhattan-like Plot of SNP Quality", x = "Genomic Position (pseudo)", y = "QUAL Score") +
  theme(legend.position = "none")
ggsave("plots/09_manhattan_like_plot.png", plot = p9, width = 10, height = 6, dpi = 300)

# === Plot 10: Pie Chart of Mutation Types ===
snp_df <- snp_df[nchar(snp_df$REF) == 1 & nchar(snp_df$ALT) == 1, ]
snp_df$mutation <- paste0(snp_df$REF, ">", snp_df$ALT)
mutation_counts <- snp_df %>%
  group_by(mutation) %>%
  summarise(count = n())
p10 <- ggplot(mutation_counts, aes(x = "", y = count, fill = mutation)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "SNP Mutation Types (Pie Chart)") +
  scale_fill_brewer(palette = "Set3")
ggsave("plots/10_mutation_type_pie_chart.png", plot = p10

<img width="3000" height="1800" alt="09_manhattan_like_plot" src="https://github.com/user-attachments/assets/f42aa19c-024b-4174-93d2-aaaffbddb092" />
<img width="3000" height="1800" alt="08_snp_density_heatmap" src="https://github.com/user-attachments/assets/1f267d2c-5db2-447c-bfc4-8be24c5c11eb" />
<img width="3000" height="2400" alt="07_qual_vs_position_faceted" src="https://github.com/user-attachments/assets/e1c7b194-f2b0-4089-b1aa-b1f628d190c2" />
<img width="1800" height="1500" alt="06_snp_type_transition_vs_transversion" src="https://github.com/user-attachments/assets/0506b22c-5544-4e55-ac40-afacf92c4672" />
<img width="2400" height="1500" alt="05_snp_qual_hist" src="https://github.com/user-attachments/assets/4cf18499-76a3-47ca-ad10-7167857e32ba" />
<img width="3000" height="1800" alt="04_top30_chr_variant_count" src="https://github.com/user-attachments/assets/5744c828-1302-4f31-a617-0c594ab35468" />
<img width="2400" height="1500" alt="03_variant_count_per_chr" src="https://github.com/user-attachments/assets/348226fb-5a0b-44a5-a2d0-b0767c61010a" />
<img width="2400" height="1500" alt="02_qual_distribution_pass_variants" src="https://github.com/user-attachments/assets/c380b3ac-2f10-4659-bbb5-2969c01c48b4" />
<img width="2400" height="1500" alt="01_qual_distribution_all_variants" src="https://github.com/user-attachments/assets/2b716113-d5f4-497e-86fe-07a7b60d988d" />
<img width="704" height="354" alt="SNPs identified plot" src="https://github.com/user-attachments/assets/f2d5d271-1651-4108-897d-5e9abcbb61e8" />


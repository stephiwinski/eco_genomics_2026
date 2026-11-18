
# Genetic Diversity

library(tidyverse)

setwd("~/projects/eco_genomics_2026/population_genomics/mydata")

##we want to output results to myresults
OUTDIR="~/projects/eco_genomics_2026/population_genomics/myresults"

meta <- read_csv("sample_metadata.csv")

# one latitude per site, so we can order the sites north to south
sites <- meta %>%
  group_by(Population) %>%
  summarise(Latitude = mean(Latitude))

# Identify northernmost and southernmost populations
north_site <- sites %>%
  slice_max(Latitude, n = 1)

south_site <- sites %>%
  slice_min(Latitude, n = 1)

# Print the populations and their latitudes
north_site
south_site


## ASU ####
# Nucleotide diveristy


pi_ASU <- read_tsv("pi/ASU_pi_100kb.windowed.pi")


pi_ASU$Population <- "ASU"
pi_ASU$Population <- as.factor(pi_ASU$Population)


# summarise pi for each site
pi_summary_ASU <- pi_ASU %>%
  group_by(Population) %>%
  summarise(n_windows = n(),
            mean_pi   = mean(PI),
            median_pi = median(PI),
            .groups = "drop")
pi_summary_ASU


pi_ASU_plot <- pi_ASU %>%
  ggplot(aes(x = Population, y = PI, fill = Population)) +
  geom_violin(scale = "width", linewidth = 0.5, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3) +
  scale_fill_manual(values = c("ASU" = "#FFB6C1")) +
  scale_y_log10() +
  labs(x = "Population", y = "log10(Nucleotide diversity)") +
  theme_classic() +
  theme(legend.position = "none") + 
  coord_flip()

pi_ASU_plot


ggsave(paste0(OUTDIR, "/pi_ASU_violin.png"), plot = pi_ASU_plot, width = 6, height = 6)




# Helper: turn (scaffold, position) into a position along the whole genome by laying the scaffolds end to end, in the order in which they appear in the data.

add_genome_pos <- function(df) {
  scaffolds <- df %>% 
    group_by(CHROM) %>% ##group by CHROM
    summarise(length = max(POS)) %>% ## length = max position within a chromosome
    arrange(match(CHROM, unique(df$CHROM))) %>% ##make sure each is unique, arranged
    mutate(offset = cumsum(as.numeric(length)) - length, ###cumulatively number (helps with labels)
           shade  = factor(row_number() %% 2))    ##value for alternating colour 
  df %>%
    left_join(scaffolds, by = "CHROM") %>% ##join this data by given data with CHROM
    mutate(genome_pos = offset + POS) ###cumalative genome pos 
}

ASU_win <- pi_ASU %>%
  mutate(POS = (BIN_START + BIN_END) / 2) %>% 
  add_genome_pos() ## this is our new function


scan_plot_ASU <- ASU_win %>%
  ggplot(aes(genome_pos/ 1e6, PI, color=Population)) +
  geom_line(alpha=0.2) +
  scale_colour_manual(values="blue") +
  labs(x = paste("position (Mb)"), y = "Nucleotide diversity") +
  theme_classic()

scan_plot_ASU

ggsave(paste0(OUTDIR, "/pi/pi_along_scaffold_ASU.png"), plot = scan_plot_ASU, width = 16, height = 4)

## KRU ####
# Nucleotide diveristy

pi_KRU <- read_tsv("pi/KRU_pi_100kb.windowed.pi")

pi_KRU$Population <- "KRU"
pi_KRU$Population <- as.factor(pi_KRU$Population)


# summarize pi for each site
pi_summary_KRU <- pi_KRU %>%
  group_by(Population) %>%
  summarise(n_windows = n(),
            mean_pi   = mean(PI),
            median_pi = median(PI),
            .groups = "drop")
pi_summary_KRU


# violin box plot
pi_KRU_plot <- pi_KRU %>%
  ggplot(aes(x = Population, y = PI, fill = Population)) +
  geom_violin(scale = "width", linewidth = 0.5, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3) +
  scale_fill_manual(values = c("KRU" = "#00BFC4")) +
  scale_y_log10() +
  labs(x = "Population", y = "log10(Nucleotide diversity)") +
  theme_classic() +
  theme(legend.position = "none") + 
  coord_flip()

pi_KRU_plot

ggsave(paste0(OUTDIR, "/pi_KRU_violin.png"), plot = pi_KRU_plot, width = 6, height = 6)


# Helper: turn (scaffold, position) into a position along the whole genome by laying the scaffolds end to end, in the order in which they appear in the data.
add_genome_pos <- function(df) {
  scaffolds <- df %>% 
    group_by(CHROM) %>% ##group by CHROM
    summarise(length = max(POS)) %>% ## length = max position within a chromosome
    arrange(match(CHROM, unique(df$CHROM))) %>% ##make sure each is unique, arranged
    mutate(offset = cumsum(as.numeric(length)) - length, ###cumulatively number (helps with labels)
           shade  = factor(row_number() %% 2))    ##value for alternating colour 
  df %>%
    left_join(scaffolds, by = "CHROM") %>% ##join this data by given data with CHROM
    mutate(genome_pos = offset + POS) ###cumalative genome pos 
}

KRU_win <- pi_KRU %>%
  mutate(POS = (BIN_START + BIN_END) / 2) %>% 
  add_genome_pos() ## this is our new function


scan_plot_KRU <- KRU_win %>%
  ggplot(aes(genome_pos/ 1e6, PI, color=Population)) +
  geom_line(alpha=0.2) +
  scale_colour_manual(values="blue") +
  labs(x = paste("position (Mb)"), y = "Nucleotide diversity") +
  theme_classic()

scan_plot_KRU

ggsave(paste0(OUTDIR, "/pi/pi_along_scaffold_KRU.png"), plot = scan_plot_KRU, width = 16, height = 4)




# FST plot ####
pop1 <- "ASU"   #the same two populations as in the bash code
pop2 <- "KRU"

fst_win <- read_tsv(paste0("fst/", pop1, "_vs_", pop2, "_FST_100kb.windowed.weir.fst"))


# Helper: turn (scaffold, position) into a position along the whole genome by laying the scaffolds end to end, in the order in which they appear in the data.
add_genome_pos <- function(df) {
  scaffolds <- df %>% 
    group_by(CHROM) %>% ##group by CHROM
    summarise(length = max(POS)) %>% ## length = max position within a chromosome
    arrange(match(CHROM, unique(df$CHROM))) %>% ##make sure each is unique, arranged
    mutate(offset = cumsum(as.numeric(length)) - length, ###cumulatively number (helps with labels)
           shade  = factor(row_number() %% 2))    ##value for alternating colour 
  df %>%
    left_join(scaffolds, by = "CHROM") %>% ##join this data by given data with CHROM
    mutate(genome_pos = offset + POS) ###cumalative genome pos 
}

fst_win <- fst_win %>%
  mutate(POS = (BIN_START + BIN_END) / 2) %>% 
  add_genome_pos() ## this is our new function


# empirical outlier threshold: the top 1% of windows
fst_cut <- quantile(fst_win$WEIGHTED_FST, probs = 0.99)
fst_cut

fst_win <- mutate(fst_win, outlier = WEIGHTED_FST >= fst_cut)

ggplot(fst_win, aes(x = genome_pos / 1e6, y = WEIGHTED_FST)) +
  geom_point(aes(colour = shade), size = 0.8) +
  geom_point(data = filter(fst_win, outlier), colour = "red", size = 1.2) +
  geom_hline(yintercept = fst_cut, linetype = "dashed") +
  scale_colour_manual(values = c("grey35", "grey70"), guide = "none") +
  labs(x = "Position along the genome (Mb; scaffolds end to end)",
       y = "Weighted FST (10-kb windows)",
       title = paste(pop1, "vs", pop2)) +
  theme_classic()


ggsave(paste0(OUTDIR,"/fst_scan_", pop1, "_vs_", pop2, ".png"), width = 10, height = 4)


# saving outliers to a table 

fst_outliers <- fst_win %>%
  filter(outlier) %>%
  arrange(desc(WEIGHTED_FST)) %>%
  select(CHROM, BIN_START, BIN_END, N_VARIANTS, WEIGHTED_FST)

head(fst_outliers, 10)
write_tsv(fst_outliers, paste0("fst/", pop1, "_vs_", pop2, "_FST_100kb_outliers.txt"))


# Combined violin plot ####

pi_combined <- bind_rows(pi_KRU, pi_ASU) %>%
  mutate(Population = factor(Population, levels = c("ASU", "KRU")))

pi_combined_plot <- pi_combined %>%
  ggplot(aes(x = Population, y = PI, fill = Population)) +
  geom_violin(scale = "width", linewidth = 0.5, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3) +
  scale_fill_manual(values = c(
    "KRU" = "#00BFC4",
    "ASU" = "#FFB6C1"
  )) +
  scale_y_log10() +
  labs(x = "Population", y = "log10(Nucleotide diversity)") +
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()

pi_combined_plot

ggsave(paste0(OUTDIR, "/combined_violin_plot.png"), plot = pi_combined_plot, width = 6, height = 5)


# putting the graphs together ####

install.packages("patchwork")

library(patchwork)
library(ggplot2)

# Read FST data
pop1 <- "ASU"
pop2 <- "KRU"

fst_win <- read_tsv(
  paste0("fst/", pop1, "_vs_", pop2,"_FST_100kb.windowed.weir.fst"))

# Calculate midpoint of each 100kb window
ASU_win <- pi_ASU %>%
  mutate(POS = (BIN_START + BIN_END) / 2)

KRU_win <- pi_KRU %>%
  mutate(POS = (BIN_START + BIN_END) / 2)

fst_win <- fst_win %>%
  mutate(POS = (BIN_START + BIN_END) / 2)

# Combine positions to establish a shared scaffold order
all_pos <- bind_rows(
  ASU_win %>% select(CHROM, POS),
  KRU_win %>% select(CHROM, POS),
  fst_win %>% select(CHROM, POS)
)

scaffold_order <- unique(all_pos$CHROM)

# Calculate scaffold offsets
scaffold_offsets <- all_pos %>%
  group_by(CHROM) %>%
  summarise(scaffold_length = max(POS, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(match(CHROM, scaffold_order)) %>%
  mutate(offset = lag(cumsum(scaffold_length), default = 0))


# Add cumulative genome positions
add_genome_pos <- function(df) {
  df %>%
    left_join(scaffold_offsets %>% select(CHROM, offset), by = "CHROM") %>%
    mutate(genome_pos = offset + POS) %>%
    arrange(CHROM, POS)
}

ASU_win <- add_genome_pos(ASU_win)
KRU_win <- add_genome_pos(KRU_win)
fst_win <- add_genome_pos(fst_win)


# ASU log2 nucleotide diversity
pi_plot_ASU <- ASU_win %>%
  filter(!is.na(PI), PI > 0) %>%
  mutate(log2_PI = -log2(PI)) %>%
  ggplot(aes(x = genome_pos / 1e6, y = log2_PI)) +
  geom_line(aes(group = CHROM), color = "#FFB6C1", alpha = 0.8, linewidth = 0.3) +
  labs(title = "ASU", x = "Position along genome (Mb)", y = expression(log[2](pi))) +
  theme_classic()


# KRU log2 nucleotide diversity
pi_plot_KRU <- KRU_win %>%
  filter(!is.na(PI), PI > 0) %>%
  mutate(log2_PI = -log2(PI)) %>%
  ggplot(aes(x = genome_pos / 1e6, y = log2_PI)) +
  geom_line(aes(group = CHROM), color = "#00BFC4", alpha = 0.8, linewidth = 0.3) +
  labs(title = "KRU", x = "Position along genome (Mb)", y = expression(log[2](pi))) +
  theme_classic()


sum(ASU_win$PI == 0, na.rm = TRUE)
sum(KRU_win$PI == 0, na.rm = TRUE)


# Calculate 99th percentile FST threshold
fst_cut <- quantile(fst_win$WEIGHTED_FST, probs = 0.99)

# Identify outlier windows
fst_win <- fst_win %>%
  mutate(outlier = WEIGHTED_FST >= fst_cut)

# FST genome-wide line graph
fst_plot <- fst_win %>%
  ggplot(aes(x = genome_pos / 1e6, y = WEIGHTED_FST)) +
  geom_line(aes(group = CHROM), color = "grey35", linewidth = 0.3) +
  geom_point(data = filter(fst_win, outlier), color = "red",size = 1) +
  geom_hline(yintercept = fst_cut,linetype = "dashed", color = "red") +
  labs(title = "ASU vs KRU", x = "Position along genome (Mb)", y = expression(F[ST])) +
  theme_classic()


# Combine all three plots
combined_pi_fst <- (
  pi_plot_KRU /
    pi_plot_ASU /
    fst_plot) +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "A")


# Display combined figure
combined_pi_fst


# Save figure
ggsave(paste0(OUTDIR, "/combined_pi_fst_genome_scan.png"), plot = combined_pi_fst, width = 14, height = 10, dpi = 300)



# FST Range ####
fst_win %>%
  summarise(
    mean_FST = mean(WEIGHTED_FST, na.rm = TRUE),
    median_FST = median(WEIGHTED_FST, na.rm = TRUE),
    min_FST = min(WEIGHTED_FST, na.rm = TRUE),
    max_FST = max(WEIGHTED_FST, na.rm = TRUE),
    fst_cutoff = first(fst_cut)
  )



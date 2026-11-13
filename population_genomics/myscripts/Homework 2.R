
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


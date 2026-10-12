## Homework #1

setwd("~/projects/eco_genomics_2026/transcriptomics/mydata")

## Import the libraries

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)  

# Import the counts matrix
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)


# now filter out the OA/OW groups here and check if F11 subset works before DESeq
# script similar to the subset function below 

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~  treatment)

dim(dds)

dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 

# F11 ####

dds_AM_OWA_F11 <- subset(dds, select = treatment %in% c("AM", "OWA") & generation == "F11")
dds_AM_OWA_F11$treatment <- droplevels(dds_AM_OWA_F11$treatment)
dds_AM_OWA_F11$treatment <- relevel(dds_AM_OWA_F11$treatment, ref = "AM")

table(dds_AM_OWA_F11$treatment)   

dds_AM_OWA_F11 <- DESeq(dds_AM_OWA_F11)
resultsNames(dds_AM_OWA_F11)

res_OWAvsAMF11 <- results(dds_AM_OWA_F11, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAMF11 <- res_OWAvsAMF11[order(res_OWAvsAMF11$padj),]
head(res_OWAvsAMF11)  
summary(res_OWAvsAMF11)


#Volcano plot F11
volcano_df11 <- as.data.frame(res_OWAvsAMF11)

volcano_F11 <- volcano_df11 %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

p3 <- ggplot(volcano_F11,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Down" = "steelblue",
    "NS"   = "grey70",
    "Up"   = "firebrick"
  )) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    title = "F11",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = NULL
  )



# F0 ####

dds_AM_OWA_F0 <- subset(dds, select = treatment %in% c("AM", "OWA") & generation == "F0")
dds_AM_OWA_F0$treatment <- droplevels(dds_AM_OWA_F0$treatment)
dds_AM_OWA_F0$treatment <- relevel(dds_AM_OWA_F0$treatment, ref = "AM")

table(dds_AM_OWA_F0$treatment)   

dds_AM_OWA_F0 <- DESeq(dds_AM_OWA_F0)
resultsNames(dds_AM_OWA_F0)

res_OWAvsAMF0 <- results(dds_AM_OWA_F0, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAMF0 <- res_OWAvsAMF0[order(res_OWAvsAMF0$padj),]
head(res_OWAvsAMF0)  
summary(res_OWAvsAMF0)


#Volcano plot F0
volcano_df0 <- as.data.frame(res_OWAvsAMF0)

volcano_F0 <- volcano_df0 %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

p2 <- ggplot(volcano_F0,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Down" = "steelblue",
    "NS"   = "grey70",
    "Up"   = "firebrick"
  )) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    title = "F0",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = NULL
  )


# F4 ####

dds_AM_OWA_F4 <- subset(dds, select = treatment %in% c("AM", "OWA") & generation == "F4")
dds_AM_OWA_F4$treatment <- droplevels(dds_AM_OWA_F4$treatment)
dds_AM_OWA_F4$treatment <- relevel(dds_AM_OWA_F4$treatment, ref = "AM")

table(dds_AM_OWA_F4$treatment)   

dds_AM_OWA_F4 <- DESeq(dds_AM_OWA_F4)
resultsNames(dds_AM_OWA_F4)

res_OWAvsAMF4 <- results(dds_AM_OWA_F4, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAMF4 <- res_OWAvsAMF4[order(res_OWAvsAMF4$padj),]
head(res_OWAvsAMF4)  
summary(res_OWAvsAMF4)


#Volcano plot F4
volcano_df4 <- as.data.frame(res_OWAvsAMF4)

volcano_F4 <- volcano_df4 %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  ) 

p1 <- ggplot(volcano_F4,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Down" = "steelblue",
    "NS"   = "grey70",
    "Up"   = "firebrick"
  )) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    title = "F4",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = NULL
  )


# Combined plots ####

# make combined volcano plot figure

library(cowplot)
install.packages("patchwork")
library(patchwork)

combined_figure <- (p2 + p1 + p3) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')

print(combined_figure)

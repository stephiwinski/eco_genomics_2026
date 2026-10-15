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

# organizing data ####

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

p1 <- ggplot(volcano_F0,
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

p2 <- ggplot(volcano_F4,
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


# Combined plot volcano plot ####

# make combined volcano plot figure

library(cowplot)
install.packages("patchwork")
library(patchwork)

combined_figure <- (p1 + p2 + p3) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')&      
  scale_y_continuous(limits = c(0, 100))&
  scale_x_continuous(limits = c(-10, 10))


print(combined_figure)


# Euler plot comparing generations

# Euler plot ####


#### PLOT OVERLAPPING DEGS IN VENN EULER DIAGRAM


# For OWA vs AM F0
res_OWAvsAMF0 <- results(dds_AM_OWA_F0, name="treatment_OWA_vs_AM", alpha=0.05) # pull out the results for the contrast of interest
res_OWAvsAMF0 <- res_OWAvsAMF0[order(res_OWAvsAMF0$padj),] # order them by significance
res_OWAvsAMF0 <- res_OWAvsAMF0[!is.na(res_OWAvsAMF0$padj),] # get rid of any NAs
degs_OWAvsAMF0 <- row.names(res_OWAvsAMF0[res_OWAvsAMF0$padj < 0.05,]) # make a list of significant differentially expressed genes for this contrast

length(degs_OWAvsAMF0) # 3924

# For OWA vs AM F4
res_OWAvsAMF4 <- results(dds_AM_OWA_F4, name="treatment_OWA_vs_AM", alpha=0.05) # pull out the results for the contrast of interest
res_OWAvsAMF4 <- res_OWAvsAMF4[order(res_OWAvsAMF4$padj),] # order them by significance
res_OWAvsAMF4 <- res_OWAvsAMF4[!is.na(res_OWAvsAMF4$padj),] # get rid of any NAs
degs_OWAvsAMF4 <- row.names(res_OWAvsAMF4[res_OWAvsAMF4$padj < 0.05,]) # make a list of significant differentially expressed genes for this contrast

length(degs_OWAvsAMF4) # 91

# For OWA vs AM
res_OWAvsAMF11 <- results(dds_AM_OWA_F11, name="treatment_OWA_vs_AM", alpha=0.05) # pull out the results for the contrast of interest
res_OWAvsAMF11 <- res_OWAvsAMF11[order(res_OWAvsAMF11$padj),] # order them by significance
res_OWAvsAMF11 <- res_OWAvsAMF11[!is.na(res_OWAvsAMF11$padj),] # get rid of any NAs
degs_OWAvsAMF11 <- row.names(res_OWAvsAMF11[res_OWAvsAMF11$padj < 0.05,]) # make a list of significant differentially expressed genes for this contrast

length(degs_OWAvsAMF11) # 1403


library(eulerr)

# Total
length(degs_OWAvsAMF0)  # 3924
length(degs_OWAvsAMF4)  # 91
length(degs_OWAvsAMF11)  # 1403

# Intersections
length(intersect(degs_OWAvsAMF0,degs_OWAvsAMF4))  # 38
length(intersect(degs_OWAvsAMF0,degs_OWAvsAMF11))  # 232
length(intersect(degs_OWAvsAMF4,degs_OWAvsAMF11))  # 8

# Shared across all
intWA <- intersect(degs_OWAvsAMF0,degs_OWAvsAMF4)
length(intersect(degs_OWAvsAMF11,intWA)) # 2

# the two shared genes shared across all
intWA <- intersect(degs_OWAvsAMF0, degs_OWAvsAMF4)
shared_genes_all <- intersect(degs_OWAvsAMF11, intWA)
print(shared_genes_all)

# Number unique to each treatment

3924-38-232+2 # 3656 F0
91-38-8+2 # 47 F4 
1403-232-8+2 # 1165 F11

# Number shared in pairs of treatments

# Shared by exactly two treatments
38-2   # 36 F0 & F4 
232-2  # 230 F0 & F11
8-2    # 6 F4 & F11 


# assemble the results
# Note that the names are important and have to be specific to line up the diagram
fit1 <- euler(c("F0" = 3656, "F4" = 47, "F11" = 1165, "F0&F4" = 36, "F0&F11" = 230, "F4&F11" = 6, "F0&F4&F11" = 2))

# make the plot
plot(fit1,  lty = 1:3, quantities = TRUE)


plot(
  fit1,
  fills = c("#66C2A5", "#FFB6C1", "#8DA0CB"),
  edges = list(col = "black", lwd = 1.5),
  labels = list(col = "black", font = 2, cex = 1.2),
  quantities = list(col = "black", cex = 1), 
  )


#cross check with above lengths of DEGS: the four values, unique, shared with one other, shared with the second other, shared across all treatments, should sum to the length of DEGs for each treatment contrast to AM
3656+230+36+2 # 3924 total F0
47+36+2+6  # 91 total F4
1165+6+2+230    # 1403  total OA



# TopGO Analysis ####

# getting files

# F0

res_OWAvsAMF0.df <- as.data.frame(res_OWAvsAMF0)
res_OWAvsAMF0.df$fullID <- rownames(res_OWAvsAMF0.df)

parts <- strsplit(res_OWAvsAMF0.df$fullID, "::")

res_OWAvsAMF0.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWAvsAMF0.df,
  "myresults/F0_OWAvsAM_results.csv",
  row.names = FALSE
)


# F4

res_OWAvsAMF4.df <- as.data.frame(res_OWAvsAMF4)
res_OWAvsAMF4.df$fullID <- rownames(res_OWAvsAMF4.df)

parts <- strsplit(res_OWAvsAMF4.df$fullID, "::")

res_OWAvsAMF4.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWAvsAMF4.df,
  "myresults/F4_OWAvsAM_results.csv",
  row.names = FALSE
)

# F11

res_OWAvsAMF11.df <- as.data.frame(res_OWAvsAMF11)
res_OWAvsAMF11.df$fullID <- rownames(res_OWAvsAMF11.df)

parts <- strsplit(res_OWAvsAMF11.df$fullID, "::")

res_OWAvsAMF11.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWAvsAMF11.df,
  "myresults/F11_OWAvsAM_results.csv",
  row.names = FALSE
)


# prep to run the GO analyses:


library(topGO)

getwd()

mappingFile <- "trinotate_annotation_GOblastx_forTopGO.txt"

geneID2GO <- readMappings(
  file = mappingFile
)

cat("Genes in GO mapping:",
    length(geneID2GO),
    "\n")

# Filter GO terms by size

allGO <- table(unlist(geneID2GO))

keepTerms <- names(allGO)[
  allGO >= 5 &
    allGO <= 500
]

geneID2GO.filtered <- lapply(
  geneID2GO,
  function(x) intersect(x, keepTerms)
)

geneID2GO.filtered <- geneID2GO.filtered[
  lengths(geneID2GO.filtered) > 0
]


## Create a function to run a TopGO contrast ####

run_topGO_contrast <- function(
    infile,
    outfile,
    ontology = "BP",
    padj.cutoff = 0.05){
  
  deseq <- read.csv(
    infile,
    stringsAsFactors = FALSE
  )
  
  deseq <- subset(
    deseq,
    !is.na(padj)
  )
  
  geneList <- factor(
    as.integer(deseq$padj < padj.cutoff)
  )
  
  names(geneList) <- deseq$shortID
  
  cat("\nGenes tested:",
      length(geneList))
  
  cat("\nSignificant genes:",
      sum(geneList == 1),
      "\n")
  
  GOdata <- new(
    "topGOdata",
    ontology = ontology,
    allGenes = geneList,
    geneSelectionFun = function(x) x == 1,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO.filtered
  )
  
  resultWeight <- runTest(
    GOdata,
    algorithm = "weight01",
    statistic = "fisher"
  )
  
  GOresults <- GenTable(
    GOdata,
    weightFisher = resultWeight,
    orderBy = "weightFisher",
    topNodes = 100
  )
  
  write.csv(
    GOresults,
    outfile,
    row.names = FALSE
  )
  
  return(GOresults)
}



## Run our TopGO function for the three different contrasts

GO_F0 <- run_topGO_contrast(
  "myresults/F0_OWAvsAM_results.csv",
  "myresults/GO_F0_OWAvsAM_BP.csv"
)

GO_F4 <- run_topGO_contrast(
  "myresults/F4_OWAvsAM_results.csv",
  "myresults/GO_F4_OWAvsAM_BP.csv"
)

GO_F11 <- run_topGO_contrast(
  "myresults/F11_OWAvsAM_results.csv",
  "myresults/GO_F11_OWAvsAM_BP.csv"
)


## Make a bubble plot of the top 10 GO categories

# Read TopGO for OWA vs AM F0 ####
go <- read.csv(
  "myresults/GO_F0_OWAvsAM_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go$weightFisher <- gsub("^<\\s*", "", go$weightFisher)

go$weightFisher <- as.numeric(go$weightFisher)

# Convert counts to numeric
go$Significant <- as.numeric(go$Significant)

# Create -log10(p)
go$minusLogP <- -log10(go$weightFisher)

# Keep top 10 GO terms
go_top10 <- go %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_top10$Term <- factor(
  go_top10$Term,
  levels = rev(go_top10$Term)
)

# Bubble plot
ggplot(
  go_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top GO Terms: OWA vs AM F0",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Significant\nGenes"
  )



# Read TopGO for OWA vs AM F4 ####
go <- read.csv(
  "myresults/GO_F4_OWAvsAM_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go$weightFisher <- gsub("^<\\s*", "", go$weightFisher)

go$weightFisher <- as.numeric(go$weightFisher)

# Convert counts to numeric
go$Significant <- as.numeric(go$Significant)

# Create -log10(p)
go$minusLogP <- -log10(go$weightFisher)

# Keep top 10 GO terms
go_top10 <- go %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_top10$Term <- factor(
  go_top10$Term,
  levels = rev(go_top10$Term)
)

# Bubble plot
ggplot(
  go_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top GO Terms: OWA vs AM F4",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Significant\nGenes"
  )



# read TopGO for OWA vs AM F11 #### 
go <- read.csv(
  "myresults/GO_F11_OWAvsAM_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go$weightFisher <- gsub("^<\\s*", "", go$weightFisher)

go$weightFisher <- as.numeric(go$weightFisher)

# Convert counts to numeric
go$Significant <- as.numeric(go$Significant)

# Create -log10(p)
go$minusLogP <- -log10(go$weightFisher)

# Keep top 10 GO terms
go_top10 <- go %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_top10$Term <- factor(
  go_top10$Term,
  levels = rev(go_top10$Term)
)

# Bubble plot
ggplot(
  go_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top GO Terms: OWA vs AM F11",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Significant\nGenes"
  )


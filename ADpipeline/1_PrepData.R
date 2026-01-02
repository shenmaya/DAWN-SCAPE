# Load Libraries ----
library(ggplot2)
library(readxl)
library(readr)
library(locfdr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(biomaRt)
library(here)

# Set Working Directory ----
setwd(here("ADpipeline"))

sapply(paste0(here(), '/mycode/', list.files(paste0(here(), '/mycode/'))), source)

# Set Figure Save Path
dir.create('figs/')
figpath <- 'figs/1_PrepData/'
if (!dir.exists(figpath)) {
  dir.create(figpath)
}

#'#######################################################################
# Prepare Bulk Data for Network ----
#'#######################################################################

# Load bulk data
datExpr.reg = readRDS('data/haney/bulk_expression_regressed.rds')
datExpr = readRDS('data/haney/bulk_expression.rds')
datMeta = readRDS('data/haney/meta.rds')

# Slice the data
# Only want control brains 
control_datreg = datExpr.reg[,which(datMeta$Diagnosis=='CTL')]
control_dat = datExpr[,which(datMeta$Diagnosis=='CTL')]
# Include all brain regions: Frontal, Temporal, Parietal, BA17 (Occipital)
rawregions = as.vector(datMeta$region[which(datMeta$Diagnosis=='CTL')])
rawregions[which(rawregions%in%c('BA44-45','BA9','BA4-6','BA24'))] = 'Frontal'
rawregions[which(rawregions%in%c('BA38','BA20-37','BA41-42-22'))] = 'Temporal'
rawregions[which(rawregions%in%c('BA3-1-2-5','BA7','BA39-40'))] = 'Parietal'
regions = rawregions

dim(control_datreg)
dim(control_dat)
length(regions)

# Convert bulk gene expression from alias to ensembl 
# Not manually converting using mapIds because results in quite a few NAs... 
# 8037 NAs, 16799 non-NAs
# Using author-created table... no NAs
bulk_df <- read.csv('data/haney/DEG.csv', header=TRUE)
# Check that the gexp matrices genes (rows) are in the same order as the table
all(rownames(control_datreg) == bulk_df$alias) # TRUE
all(rownames(control_dat) == bulk_df$external_gene_name) # TRUE
rownames(control_datreg) <- bulk_df$ensembl_gene_id
rownames(control_dat) <- bulk_df$ensembl_gene_id

#'#######################################################################
# Prepare P-Value Sets ----
#'#######################################################################

# Use ensembl database
ensembl <- useMart("ensembl")
# Example: Select the human gene dataset
ensembl_human <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Load gene p-values ----
## DE (Random) P-Values ----
de_df <- as.data.frame(read_tsv('data/meta.anlz.ad_cntrl.tsv', show_col_types=F))
de_df <- de_df[c('ensembl_gene_id', 'pval.random', 'zval.random', 'fdr.random')]

# How many protein-coding genes?
de_annot <- getBM(attributes = c('ensembl_gene_id', 'gene_biotype', 'external_gene_name'),
                     filters = 'ensembl_gene_id',
                     values = de_df$ensembl_gene_id,
                     mart = ensembl_human)
de_protein_coding_genes <- de_annot[de_annot$gene_biotype == "protein_coding", ]
nrow(de_protein_coding_genes) #13103

## AP/Genetic Risk (Validation) P-Values ----
gen_df <- read.csv('data/TWAS validation-Table 1.csv')
# Remove .num extensions from ENSG IDs
gen_df$gene <- sub("\\..*", "", gen_df$gene)

# How many protein-coding genes?
gen_annot <- getBM(attributes = c('ensembl_gene_id', 'gene_biotype', 'external_gene_name'),
                  filters = 'ensembl_gene_id',
                  values = gen_df$gene,
                  mart = ensembl_human)
gen_protein_coding_genes <- gen_annot[gen_annot$gene_biotype == "protein_coding", ]
nrow(gen_protein_coding_genes) #14552

# Gene Names -----
gen_df$ensembl <- sub("\\..*", "", gen_df$gene)
gen_df$alias <- mapIds(org.Hs.eg.db, keys=gen_df$ensembl, keytype='ENSEMBL', column='SYMBOL')

de_df$ensembl <- de_df$ensembl_gene_id
de_df$alias <- mapIds(org.Hs.eg.db, keys=de_df$ensembl, keytype='ENSEMBL', column='SYMBOL')

# Setting AD Fixed Genes ----
fixedADGENgenes <- read.csv('data/fixedADGENgenes.csv')
# One gene has an extra space
fixedADGENgenes$Gene <- sapply(strsplit(fixedADGENgenes$Gene, ' '), "[[", 1)
# One gene is lowercase
fixedADGENgenes$Gene <- toupper(fixedADGENgenes$Gene)

fixedADGENgenes$Gene # 53 fixed genes
fixedADGENgenes$alias <- fixedADGENgenes$Gene
fixedADGENgenes$ensembl <- mapIds(org.Hs.eg.db, keys=fixedADGENgenes$Gene, keytype='SYMBOL', column='ENSEMBL')
fixedADGENgenes <- fixedADGENgenes[, c('alias', 'ensembl', 'Evidence')]

# Scenarios: As long as the gene exists in the bulk data, we can use the fixed gene...
# If missing AP p-value: set as 0
# If missing DE p-value: set as 0.5 (essentially no info)
length(intersect(rownames(control_dat), fixedADGENgenes$ensembl)) # 52 fixed AP genes in bulk exp matrix
# Missing fixed AP gene: MS4A2
fixedADGENgenes[fixedADGENgenes$ensembl == setdiff(fixedADGENgenes$ensembl, rownames(control_dat)),]

# Fill DE and AP p-values
de_df$fixed <- rep(F, nrow(de_df))
gen_df$fixed <- rep(F, nrow(gen_df))
de_df$orig_pval <- de_df$pval.random
gen_df$orig_pval <- gen_df$pvalue

gen_df <- gen_df[, c('ensembl', 'alias', 'fixed', 'pvalue', 'orig_pval')]
colnames(gen_df) <- c('ensembl', 'alias', 'fixed', 'pvals', 'orig_pval')
de_df <- de_df[, c('ensembl', 'alias', 'fixed', 'pval.random', 'orig_pval')]
colnames(de_df) <- c('ensembl', 'alias', 'fixed', 'pvals', 'orig_pval')

for (gene_ensembl in intersect(rownames(control_dat), fixedADGENgenes$ensembl)) {
  gene_alias <- fixedADGENgenes[fixedADGENgenes$ensembl == gene_ensembl,]$alias
  if (!(gene_ensembl %in% de_df$ensembl)) {
    de_df <- rbind(de_df, data.frame('ensembl'=gene_ensembl,
                                     'alias'=gene_alias,
                                     'fixed'=T,
                                     'pvals'=0.5, 
                                     'orig_pval'=NA)
    )
  }
  if (!(gene_ensembl %in% gen_df$ensembl)) {
    gen_df <- rbind(gen_df, data.frame('ensembl'=gene_ensembl,
                                     'alias'=gene_alias,
                                     'fixed'=T,
                                     'pvals'=0, 
                                     'orig_pval'=NA)
    )
  } else {
    gen_df[which(gen_df$ensembl == gene_ensembl), 'pvals'] <- 0
    gen_df[which(gen_df$ensembl == gene_ensembl), 'fixed'] <- T
  }
}
# Check that all AP fixed genes were correctly marked and corrected in gen_df
sum(gen_df$fixed) == length(intersect(rownames(control_dat), fixedADGENgenes$ensembl))
setdiff(gen_df[gen_df$fixed,]$ensembl, intersect(rownames(control_dat), fixedADGENgenes$ensembl))
setdiff(intersect(rownames(control_dat), fixedADGENgenes$ensembl), gen_df[gen_df$fixed,]$ensembl)
gen_df[gen_df$fixed,]$pvals

de_p <- de_df$pvals
gen_p <- gen_df$pvals # fixed p-values, i.e. some genes from literature set to 0
# gen_p <- gen_df$orig_pvals # original p-values

overlaid_hists(list(de_p, gen_p), 
               labels=c('DE', 'AP'), main_title='DE + AP P-Values')

dev.off()

#'#######################################################################
# Calibrate Z-Scores ----
#'#######################################################################

## Convert P-Values to Z-Scores ----
de_z <- qnorm(1-de_p)
sum(is.infinite(de_z)) # 1 Infinite
gen_z <- qnorm(1-gen_p) 
sum(is.infinite(gen_z)) # 55 Infinite with fixed, 9 infinite otherwise

overlaid_hists(list(de_z[is.finite(de_z)], gen_z[is.finite(gen_z)]), breaks=50,
               labels=c('DE', 'AP'), main_title='DE + AP Z-Scores',
               savepath=paste0(figpath, 'zscores_deap.png'))

## Correct P-Values using Efron's Correction ----
## DE (Random) P-Values ----
hist(de_z, breaks=100)
dev.off()
sort(de_z, decreasing=T)[1:20]
de_z_corr <- de_z
set.seed(42)
a <- floor(sort(unique(de_z_corr[is.finite(de_z_corr)])*2, decreasing=T)[2])/2 + rnorm(1, 0, 0.001) # 8 + noise
b <- ceiling(sort(unique(de_z_corr[is.finite(de_z_corr)])*2, decreasing=T)[1])/2 + rnorm(1, 0, 0.001) # 8.5 + noise
unif_seq <- seq(from=a, to=b, length.out=sum(is.infinite(de_z_corr))+2)
de_z_corr[is.infinite(de_z_corr)] <- unif_seq[2:(length(unif_seq)-1)]

sort(de_z_corr, decreasing=T)[1:20]
hist(de_z_corr, breaks=50)
de_z <- de_z_corr
de_df$zscores <- de_z
# locfdr
pdf(paste0(figpath, "dezscores_locfdr.pdf"), width = 8, height = 6)
de_locfdr <- locfdr(de_z, main='DE')
dev.off()
de_locfdr_mu <- de_locfdr$fp0['mlest','delta'] # 0.8920404
de_locfdr_sigma <- de_locfdr$fp0['mlest','sigma'] # 1.579717
de_z_transf <- (de_z - de_locfdr_mu)/de_locfdr_sigma
de_transf_locfdr <- locfdr(de_z_transf)
de_transf_locfdr_mu <- de_transf_locfdr$fp0['mlest','delta'] 
de_transf_locfdr_sigma <- de_transf_locfdr$fp0['mlest','sigma']
str(round(c(de_transf_locfdr_mu, de_transf_locfdr_sigma), 4)) # 0, 1
de_df$zscores_cal <- de_z_transf
de_df$pvals_cal <- 1-pnorm(de_z_transf) 
dev.off()

overlaid_hists(list(de_df$zscores, de_df$zscores_cal), breaks=50,
               labels=c('zscores', 'cal zscores'), main_title='DE Z-Scores',
               savepath=paste0(figpath, 'zscores_de.png'))

overlaid_hists(list(de_df$pvals, de_df$pvals_cal), breaks=50,
               labels=c('pvals', 'cal pvals'), main_title='DE P-Values',
               savepath=paste0(figpath, 'pvals_de.png'))

before_label <- paste0('before: mean = ', round(de_locfdr_mu, 2), ', sd = ', round(de_locfdr_sigma, 2))
after_label <- paste0('after: mean = ', round(de_transf_locfdr_mu, 2), ', sd = ', round(de_transf_locfdr_sigma, 2))
ggplot(data.frame(zscore=c(de_z, de_z_transf),
                  Calibration = factor(c(rep('Before', length(de_z)),
                                         rep('After', length(de_z_transf))),
                                       levels=c('Before', 'After'))),
       aes(x=zscore, fill = Calibration)) +
  geom_histogram(alpha = 0.5, position="identity", bins=60)+
  ggtitle('DE')+
  xlab('Z-Score')+
  ylab('Count')+
  # geom_vline(xintercept=0.05, color = "black", linewidth=0.7)+
  annotate("text", x = 2.2, y = 650, label = before_label,
           hjust = 0, vjust = 0, size = 2.9, color = "#F8766D") +
  annotate("text", x = 1.5, y = 850, label = after_label,
           hjust = 0, vjust = 0, size = 2.9, color = "#00BFC4") +
  coord_cartesian(xlim=c(-5,10))
ggsave(paste0(figpath, 'dezscores_cal.pdf'),width=5,height=3,dpi=700)

## AP (Validation) P-Values ----
hist(gen_z, breaks=100)
sort(gen_z, decreasing=T)[1:70]
gen_z_corr <- gen_z
set.seed(42)
a <- floor(sort(unique(gen_z_corr[is.finite(gen_z_corr)])*2, decreasing=T)[2])/2 + rnorm(1, 0, 0.001) # 5.5 + noise
b <- ceiling(sort(unique(gen_z_corr[is.finite(gen_z_corr)])*2, decreasing=T)[1])/2 + rnorm(1, 0, 0.001) # 8.5 + noise
unif_seq <- seq(from=a, to=b, length.out=sum(is.infinite(gen_z_corr))+2)
gen_z_corr[is.infinite(gen_z_corr)] <- sample(unif_seq[2:(length(unif_seq)-1)], replace=F)
hist(gen_z_corr, breaks=100)
sort(gen_z_corr, decreasing=T)[1:70]
sort(gen_z, decreasing=T)[1:70]
gen_z_transf <- gen_z_corr[gen_z_corr > 4]
gen_z_transf <- (gen_z_transf-min(gen_z_transf))/(max(gen_z_transf) - min(gen_z_transf))
gen_z_transf <- 4+(gen_z_transf*(6-4))
gen_z_corr[gen_z_corr > 4] <- gen_z_transf
hist(gen_z_corr, breaks=100)

gen_z <- gen_z_corr
gen_df$zscores <- gen_z

# locfdr
pdf(paste0(figpath, "apzscores_locfdr.pdf"), width = 8, height = 6)
gen_locfdr <- locfdr(gen_z, main='AP')
dev.off()
gen_locfdr_mu <- gen_locfdr$fp0['mlest','delta'] # 0.08886073
gen_locfdr_sigma <- gen_locfdr$fp0['mlest','sigma'] # 1.039575
gen_z_transf <- (gen_z - gen_locfdr_mu)/gen_locfdr_sigma
gen_df$zscores_cal <- gen_z_transf
gen_df$pvals_cal <- 1-pnorm(gen_z_transf) 

# No calibration/correction for AP z-scores
overlaid_hists(list(gen_df$zscores, gen_df$zscores_cal), breaks=50,
               labels=c('zscores', 'cal zscores'), main_title='AP Z-Scores',
               savepath=paste0(figpath, 'zscores_ap.png'))

overlaid_hists(list(gen_df$pvals, gen_df$pvals_cal), breaks=50,
               labels=c('pvals', 'cal pvals'), main_title='AP P-Values',
               savepath=paste0(figpath, 'pvals_ap.png'))

gen_label <- paste0('mean = ', round(gen_locfdr_mu, 2), ', sd = ', round(gen_locfdr_sigma, 2))
ggplot(data.frame(zscore=c(gen_z)),
       aes(x=zscore)) +
  geom_histogram(alpha = 0.5, position="identity", bins=60)+
  ggtitle('AP')+
  xlab('Z-Score')+
  ylab('Count')+
  geom_vline(xintercept=0.05, color = "black", linewidth=0.7)+
  annotate("text", x = 2.5, y = 800, label = gen_label,
           hjust = 0, vjust = 0, size = 2.9, color = "black") +
  coord_cartesian(xlim=c(-5,10))
ggsave(paste0(figpath, 'apzscores.pdf'),width=3.97,height=3,dpi=700)

# Compare DE + AP Distributions ----
overlaid_hists(list(de_df$pvals_cal, gen_df$pvals_cal), breaks=50,
               labels=c('DE', 'AP'), main_title='DE + AP P-Values: Calibrated',
               savepath=paste0(figpath, 'pvals_deap_cal.png'))

overlaid_hists(list(de_df$zscores_cal, gen_df$zscores_cal), breaks=50,
               labels=c('DE', 'AP'), main_title='DE + AP Z-Scores: Calibrated',
               savepath=paste0(figpath, 'zscores_deap_cal.png'))

# Generate Histograms ----
overlaid_hists(list(de_p, gen_p), 
               labels=c('DE', 'AP'), main_title='DE + AP P-Values')

overlaid_hists(list(de_df$zscores, de_df$zscores_cal), breaks=50,
               labels=c('zscores', 'cal zscores'), main_title='DE Z-Scores',
               savepath=paste0(figpath, 'zscores_de.png'))

overlaid_hists(list(de_df$pvals, de_df$pvals_cal), breaks=50,
               labels=c('pvals', 'cal pvals'), main_title='DE P-Values',
               savepath=paste0(figpath, 'pvals_de.png'))

overlaid_hists(list(gen_df$zscores, gen_df$zscores_cal), breaks=50,
               labels=c('zscores', 'cal zscores'), main_title='AP Z-Scores',
               savepath=paste0(figpath, 'zscores_ap.png'))

overlaid_hists(list(gen_df$pvals, gen_df$pvals_cal), breaks=50,
               labels=c('pvals', 'cal pvals'), main_title='AP P-Values',
               savepath=paste0(figpath, 'pvals_ap.png'))

#'#######################################################################
# Select Genes ----
#'#######################################################################

## Get intersection of genes ----
# DE genes and bulk data set are from same study
selgenes <- intersect(intersect(rownames(control_datreg), gen_df$ensembl), de_df$ensembl)
str(selgenes) # 11772

gen_df <- gen_df[match(selgenes, gen_df$ensembl),]
de_df <- de_df[match(selgenes, de_df$ensembl),]

control_datreg <- control_datreg[selgenes, ]
control_dat <- control_dat[selgenes, ]

plot(gen_df$zscores, de_df$zscores_cal, pch=16, cex=0.3)

#'#######################################################################
# Save Files ----
#'#######################################################################

write.csv(gen_df, 'data/gen_pvals.csv')
write.csv(de_df, 'data/de_pvals.csv')
saveRDS(control_datreg, 'data/bulk_reg.rds')
saveRDS(control_dat, 'data/bulk.rds')

overlaid_hists(list(de_df$pvals, gen_df$pvals), breaks=50,
               labels=c('DE', 'AP'), main_title='Common Genes: DE + AP P-Values',
               xlab_title='P-Value',
               savepath=paste0(figpath, 'commongenes_pvals_degen.png'))
overlaid_hists(list(de_df$zscores, gen_df$zscores), breaks=50,
               labels=c('DE', 'AP'), main_title='Common Genes: DE + AP Z-Scores',
               xlab_title='Z-Score',
               savepath=paste0(figpath, 'commongenes_zscores_degen.png'))
overlaid_hists(list(de_df$pvals, de_df$pvals_cal), breaks=50,
               labels=c('P-Value', 'Calib P-Value'), main_title='DE Common Genes: P-Values vs Calibrated',
               xlab_title='P-Value',
               savepath=paste0(figpath, 'commongenes_pvals_de.png'))
overlaid_hists(list(de_df$zscores, de_df$zscores_cal), breaks=50,
               labels=c('Z-Score', 'Calib Z-Score'), main_title='DE Common Genes: Z-Scores vs Calibrated',
               xlab_title='P-Value',
               savepath=paste0(figpath, 'commongenes_zscores_de.png'))

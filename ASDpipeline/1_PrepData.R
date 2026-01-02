# Load Libraries ----
library(ggplot2)
library(readxl)
library(locfdr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Set Working Directory ----
setwd(here("ASDpipeline"))

sapply(paste0(here(), '/mycode/', list.files(paste0(here(), '/mycode/'))), source)

# Set Figure Save Path
figpath <- 'figs/1_PrepData/'
if (!dir.exists(figpath)) {
  dir.create(figpath)
}

#'#######################################################################
# Prepare Bulk Data for Network ----
#'#######################################################################

# Load bulk data
datExpr.reg = readRDS('./data/haney/bulk_expression_regressed.rds')
datExpr = readRDS('./data/haney/bulk_expression.rds')
datMeta = readRDS('./data/haney/meta.rds')

# Slice the data
# Only want control brains 
control_datreg = datExpr.reg[,which(datMeta$Diagnosis=='CTL')]
control_dat = datExpr[,which(datMeta$Diagnosis=='CTL')]
# Only want certain brain regions: Frontal, Temporal, BA17 (Occipital)
rawregions = as.vector(datMeta$region[which(datMeta$Diagnosis=='CTL')])
rawregions[which(rawregions%in%c('BA44-45','BA9','BA4-6','BA24'))] = 'Frontal'
rawregions[which(rawregions%in%c('BA38','BA20-37','BA41-42-22'))] = 'Temporal'
rawregions[which(rawregions%in%c('BA3-1-2-5','BA7','BA39-40'))] = 'Parietal'
control_datreg = control_datreg[,which(rawregions!='Parietal')]
control_dat = control_dat[,which(rawregions!='Parietal')]
regions = rawregions[which(rawregions!='Parietal')]

#'#######################################################################
# Prepare P-Value Sets ----
#'#######################################################################

## DE P-Values ----
DE <- read.csv('./data/haney/DEG.csv', header=TRUE)
# Calculate P-Values (Non-FDR-Corrected P-Values)
DE$pvals <- DE$Whole_Cortex_FDR_P * rank(DE$Whole_Cortex_FDR_P) / (max(DE$Whole_Cortex_FDR_P) * length(DE$Whole_Cortex_FDR_P))
de_p <- DE$pvals
str(de_p) # 24836
hist(de_p)

## AP/Genetic Risk P-Values ----
TADA <- read_excel('./data/supplementary_table_11.xlsx')
TADA$pvals <- TADA$p_TADA_ASD
gen_p <- TADA$pvals
str(gen_p) # 18128
hist(gen_p)

overlaid_hists(list(de_p, gen_p), 
               labels=c('DE', 'AP'), main_title='DE + AP P-Values')

#'#######################################################################
# Calibrate Z-Scores ----
#'#######################################################################

## Convert P-Values to Z-Scores ----
de_z <- qnorm(1-de_p)
sum(is.infinite(de_z)) # 0 Infinite
gen_z <- qnorm(1-gen_p) 
sum(is.infinite(gen_z)) # 12 Infinite

hist(de_z, breaks=50)
hist(gen_z, breaks=50)
overlaid_hists(list(de_z, gen_z[is.finite(gen_z)]), breaks=50,
               labels=c('DE', 'AP'), main_title='DE + AP Z-Scores',
               savepath=paste0(figpath, 'zscores_degen.png'))

## Correct P-Values using Efron's Correction ----
# DE
# No infinite z-scores to deal with 
DE$zscores <- de_z
pdf(paste0(figpath, "dezscores_locfdr.pdf"), width = 8, height = 6)
de_locfdr <- locfdr(de_z, main='DE')
dev.off()
de_locfdr_mu <- de_locfdr$fp0['mlest','delta'] 
de_locfdr_sigma <- de_locfdr$fp0['mlest','sigma']
str(round(c(de_locfdr_mu, de_locfdr_sigma), 3)) # 0.537, 1.477
de_z_transf <- (de_z - de_locfdr_mu)/de_locfdr_sigma
de_transf_locfdr <- locfdr(de_z_transf)
de_transf_locfdr_mu <- de_transf_locfdr$fp0['mlest','delta'] 
de_transf_locfdr_sigma <- de_transf_locfdr$fp0['mlest','sigma']
str(round(c(de_transf_locfdr_mu, de_transf_locfdr_sigma), 4)) # 0, 1
DE$zscores_cal <- de_z_transf
DE$pvals_cal <- 1-pnorm(de_z_transf) 

overlaid_hists(list(DE$zscores, DE$zscores_cal), breaks=50,
               labels=c('zscores', 'cal zscores'), main_title='DE Z-Scores',
               savepath=paste0(figpath, 'zscores_de.png'))

overlaid_hists(list(DE$pvals, DE$pvals_cal), breaks=50,
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
  annotate("text", x = 2.2, y = 1000, label = before_label,
           hjust = 0, vjust = 0, size = 2.9, color = "#F8766D") +
  annotate("text", x = 1, y = 1500, label = after_label,
           hjust = 0, vjust = 0, size = 2.9, color = "#00BFC4") +
  # geom_text(aes(x = 2.2, y = 1000, label = before_label), 
  #           hjust = 0, vjust = 0, size = 3, color = "#F8766D")+
  # geom_text(aes(x = 1, y = 1500, label = after_label), 
  #           hjust = 0, vjust = 0, size = 3, color = "#00BFC4")+
  coord_cartesian(xlim=c(-5,10))
ggsave(paste0(figpath, 'dezscores_cal.pdf'),width=5,height=3,dpi=700)

# AP 
# 12 infinite z-scores to deal with 
hist(gen_z, breaks=100)
sort(gen_z, decreasing=T)[1:70]
gen_z_corr <- gen_z
set.seed(42)
a <- floor(sort(unique(gen_z_corr[is.finite(gen_z_corr)])*2, decreasing=T)[2])/2 + rnorm(1, 0, 0.001) # 7.5 + noise
b <- ceiling(sort(unique(gen_z_corr[is.finite(gen_z_corr)])*2, decreasing=T)[1])/2 + rnorm(1, 0, 0.001) # 8.5 + noise 
unif_seq <- seq(from=a, to=b, length.out=sum(is.infinite(gen_z_corr))+2)
gen_z_corr[is.infinite(gen_z_corr)] <- sample(unif_seq[2:(length(unif_seq)-1)], replace=F)
hist(gen_z_corr, breaks=100)

gen_z <- gen_z_corr
TADA$zscores <- gen_z
pdf(paste0(figpath, "apzscores_locfdr.pdf"), width = 8, height = 6)
gen_locfdr <- locfdr(gen_z, main='AP')
dev.off()
gen_locfdr_mu <- gen_locfdr$fp0['mlest','delta'] 
gen_locfdr_sigma <- gen_locfdr$fp0['mlest','sigma']
str(round(c(gen_locfdr_mu, gen_locfdr_sigma), 3)) # 0.018, 1.018
# No calibration/correction for AP z-scores

gen_label <- paste0('mean = ', round(gen_locfdr_mu, 2), ', sd = ', round(gen_locfdr_sigma, 2))
ggplot(data.frame(zscore=c(gen_z)),
       aes(x=zscore)) +
  geom_histogram(alpha = 0.5, position="identity", bins=60)+
  ggtitle('AP')+
  xlab('Z-Score')+
  ylab('Count')+
  geom_vline(xintercept=0.05, color = "black", linewidth=0.7)+
  annotate("text", x = 2.5, y = 1000, label = gen_label,
           hjust = 0, vjust = 0, size = 2.9, color = "black") +
  coord_cartesian(xlim=c(-5,10))
ggsave(paste0(figpath, 'apzscores.pdf'),width=3.97,height=3,dpi=700)

# Recreate overlaid histograms with corrected AP z-sscores
overlaid_hists(list(de_z, gen_z), breaks=50,
               labels=c('DE', 'AP'), main_title='DE + AP Z-Scores (Corrected)',
               savepath=paste0(figpath, 'zscores_degen_corr.png'))

#'#######################################################################
# Select Genes ----
#'#######################################################################

# Convert bulk gene expression from alias to ensembl 
all(rownames(control_datreg) == DE$external_gene_name)
all(rownames(control_dat) == DE$external_gene_name)
rownames(control_datreg) <- DE$ensembl_gene_id
rownames(control_dat) <- DE$ensembl_gene_id

## Get intersection of genes ----
# DE genes and bulk data set are from same study
selgenes <- intersect(rownames(control_datreg), TADA$gene_id)
str(selgenes) # 15633

gen_df <- TADA[match(selgenes, TADA$gene_id),][c('gene_id', 'gene', 'gene_gencodeV33', 'gene_gnomad', 
                                                 'pvals', 'zscores')]
colnames(gen_df) <- c('ensembl', 'gene', 'gene_gencodeV33', 'gene_gnomad', 'pvals', 'zscores')

de_df <- DE[match(selgenes, DE$ensembl_gene_id),][c('ensembl_gene_id', 'external_gene_name', 'hgnc_symbol',
                                                    'pvals', 'zscores', 'zscores_cal', 'pvals_cal')]
colnames(de_df) <- c('ensembl', 'gene', 'gene_hgnc', 'pvals', 'zscores', 'zscores_cal', 'pvals_cal')

control_datreg <- control_datreg[selgenes, ]
control_dat <- control_dat[selgenes, ]

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
               labels=c('P-Value', 'Calib P-Value'), main_title='Common Genes: P-Values vs Calibrated',
               xlab_title='P-Value',
               savepath=paste0(figpath, 'commongenes_pvals_de.png'))
overlaid_hists(list(de_df$zscores, de_df$zscores_cal), breaks=50,
               labels=c('Z-Score', 'Calib Z-Score'), main_title='Common Genes: Z-Scores vs Calibrated',
               xlab_title='P-Value',
               savepath=paste0(figpath, 'commongenes_zscores_de.png'))

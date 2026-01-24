library(doBy)
library(igraph)
library(here)

setwd(here("tutorial"))

sapply(paste0(here(), '/mycode/', list.files(paste0(here(), '/mycode/'))), source)

# Load files ----
bulk <- readRDS('data/bulk.rds')
gen_df <- readRDS('data/gen_df.rds')
rownames(gen_df) <- gen_df$ensembl
de_df <- readRDS('data/de_df.rds')
rownames(de_df) <- de_df$ensembl

gen_pvals <- gen_df$pvals
names(gen_pvals) <- gen_df$ensembl
de_pvals <- de_df$pvals
names(de_pvals) <- de_df$ensembl

# Construct Network ----
screen_genes <- unique(c(names(gen_pvals[gen_pvals < quantile(gen_pvals, 0.5)]), 
                         names(de_pvals[de_pvals < quantile(de_pvals, 0.5)])))

best.lambda <- 0.25
pnsgraph <- construct_graph(t(bulk), gen_pvals, pvalthresh = NULL, corthresh = 0.35, pvalscreen = screen_genes, 
                            lambda = best.lambda, savefile=paste0('outputs/pns_lambda',best.lambda, '.rds'), standardize=T, verbose=T)
cleangraph <- clean_graph(pnsgraph$graphfinal, degcut=1, mincompsize=10, degdivide=5, savefile=paste0('outputs/pns_lambda',best.lambda, '.rds'), verbose=T)

set.seed(42)
layout <- layout_with_fr(cleangraph$net)

pdf(paste0("outputs/graph_fr.pdf"), height = 5, width = 5)
par(mar=c(0,0,0,0)+.1)
plot.igraph(cleangraph$net, vertex.size=V(cleangraph$net)$size*2,
            layout=layout,
            vertex.label=NA,
            vertex.color='orange',
            vertex.frame.color='black',
            edge.width = 0.75,
            edge.color = 'gray')
dev.off()

# Run Joint-HMRF ----
genes <- rownames(cleangraph$ne)
dezscores <- de_df[genes,]$zscores
genzscores <- gen_df[genes,]$zscores
graph <- cleangraph$ne
deinit <- geninit <- 1.25
detrim <- gentrim <- 4
numiter <- 100

Istart <- 1*((dezscores > deinit) & (genzscores > geninit)) + 2*((dezscores > deinit) & (genzscores <= geninit)) + 3*((dezscores <= deinit) & (genzscores > geninit)) + 4*((dezscores <= deinit) & (genzscores <= geninit))
table(Istart)

# Run Joint-HMRF on one initialization set
# For running Joint-HMRF on multiple initializations, see run_hmrf_parallel and analysis files 3_RunJointHMRF.R
hmrf_res <- run_hmrf(genes, dezscores, genzscores, graph, deinit, geninit, detrim, gentrim,
                     numiter, b0_lb=-20, b1_ub=10, savefile='outputs/hmrf_res.rds', verbose=T)
table(hmrf_res$Iupdate)

# Obtain final states
t <- ifelse(hmrf_res$logppost_hist[numiter] >= hmrf_res$logppost_hist[numiter-1], numiter, numiter-1)

states <- final_states(genes, graph, hmrf_res$state_hist[t,], 
                                  dezscores, genzscores,
                           hmrf_res$param_hist[t, 1], hmrf_res$param_hist[t, 2], hmrf_res$param_hist[t, 3], 
                           hmrf_res$param_hist[t, 4], hmrf_res$param_hist[t, 5], 
                           hmrf_res$param_hist[t, 6], hmrf_res$param_hist[t, 7], hmrf_res$param_hist[t, 8], hmrf_res$param_hist[t, 9],
                           hmrf_res$param_hist[t, 10], hmrf_res$param_hist[t, 11])$Iupdate
names(states) <- genes
hmrf_res$Ifinal <- states
saveRDS(hmrf_res, 'outputs/hmrf_res.rds')

# State Color Palette

# Color Scheme
moduleTypecolors = c('Etiological'='#E66100',
                     'Emergent'='#1A85FF',
                     'Other'='#C5B6E0')
state_colorref <- c(moduleTypecolors[['Etiological']], moduleTypecolors[['Emergent']], moduleTypecolors[['Other']], moduleTypecolors[['Other']])
# Get color for each gene based on state
state_colors <- state_colorref[states]
names(state_colors) <- names(states)
# Make sure genes are in same order as graph nodes
state_colors <- state_colors[V(cleangraph$net)$name]

pdf(paste0("outputs/graph_fr_state.pdf"), height = 5, width = 5)
par(mar=c(0,0,0,0)+.1)
plot.igraph(cleangraph$net, vertex.size=V(cleangraph$net)$size*2,
            layout=layout,
            vertex.label=NA,
            vertex.color=state_colors,
            vertex.frame.color='black',
            edge.width = 0.75,
            edge.color = 'gray')
dev.off()

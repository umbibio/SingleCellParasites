library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(RColorBrewer)
library(sctransform)

## Reading the data
## These files were prepared using the python script `get_h5ad_data_kz.ipynb`. 


input.dir <- '~/work/ToxoplasmaGondii/scRNAseqToxoAtlas/kz/'

rh384.expr.file <- 'rh384_expression_filtered.csv'
rh96.expr.file  <- 'rh96_expression_filtered.csv'
pru.expr.file   <- 'pru_expression_filtered.csv'
me49.expr.file  <- 'me49_expression_filtered.csv'
pru.falc.expr.file   <- 'pru_falc_expression_filtered.csv'

rh384.obs.file <- 'rh384_obs_filtered.csv'
rh96.obs.file  <- 'rh96_obs_filtered.csv'
pru.obs.file   <- 'pru_obs_filtered.csv'
me49.obs.file  <- 'me49_obs_filtered.csv'
pru.obs.expr.file   <- 'pru_falc_obs_filtered.csv'

rh384.var.file <- 'rh384_var_filtered.csv'
rh96.var.file  <- 'rh96_var_filtered.csv'
pru.var.file   <- 'pru_var_filtered.csv'
me49.var.file  <- 'me49_var_filtered.csv'
pru.var.expr.file   <- 'pru_falc_var_filtered.csv'

### Plasmodium
input.dir.plasmodium <- "/Users/kouroshz/work/ToxoPlasmaGondii/scRNAseqPlasmodAtlas/Expression_Matrices/Smartseq2/"
plasmodium.count.file <- 'SS2_counts.csv'
plasmodium.pheno.file <- 'SS2_pheno.csv'



## Reciprocal orthologs
GT1.Pberghei <- read.xlsx('~/work/ToxoPlasmaGondii/orthologs/rec_GT1.vs.P.bergei.xlsx')
GT1.ME49 <- read.xlsx('~/work/ToxoPlasmaGondii/orthologs/rec_GT1.vs.ME49.xlsx')

GT1.ME49.Pb <- inner_join(GT1.ME49, GT1.Pberghei, by='query_id')
GT1.ME49.Pb <- GT1.ME49.Pb %>% dplyr::select('query_id', contains("subject_id"))

colnames(GT1.ME49.Pb) <- c('GT1', 'ME49', 'PB')

## Reading in the Toxo data
rh384.expr <- read.csv(paste(input.dir, rh384.expr.file, sep = ''))
rh384.var <- read.csv(paste(input.dir, rh384.var.file, sep = ''))
rh384.obs <- read.csv(paste(input.dir, rh384.obs.file, sep = ''))

rh96.expr <- read.csv(paste(input.dir, rh96.expr.file, sep = ''))
rh96.var <- read.csv(paste(input.dir, rh96.var.file, sep = ''))
rh96.obs <- read.csv(paste(input.dir, rh96.obs.file, sep = ''))

pru.expr <- read.csv(paste(input.dir, pru.expr.file, sep = ''))
pru.var <- read.csv(paste(input.dir, pru.var.file, sep = ''))
pru.obs <- read.csv(paste(input.dir, pru.obs.file, sep = ''))


me49.expr <- read.csv(paste(input.dir, me49.expr.file, sep = ''))
me49.var <- read.csv(paste(input.dir, me49.var.file, sep = ''))
me49.obs <- read.csv(paste(input.dir, me49.obs.file, sep = ''))

## Reading Plasmodium data
plasmodium.count <- read.csv(paste(input.dir.plasmodium, plasmodium.count.file, sep = ''))
plasmodium.pheno <- read.csv(paste(input.dir.plasmodium, plasmodium.pheno.file, sep = ''))


## Processing the RH data
processCounts <- function(expr){
  cols <- expr[,1]
  expr <- t(expr[,-1])
  colnames(expr) <- cols
  return(expr)
}

rh384.expr <- processCounts(rh384.expr)
rh96.expr <- processCounts(rh96.expr)
pru.expr <- processCounts(pru.expr)
me49.expr <- processCounts(me49.expr)


genes <- plasmodium.count$X
plasmodium.count <- plasmodium.count[,-1]
rownames(plasmodium.count) <- genes


## Get common genes
pb.genes    <- data.frame(PB = rownames(plasmodium.count), stringsAsFactors = F)
rh348.genes <- data.frame(GT1 = rownames(rh384.expr), stringsAsFactors = F)
rh96.genes  <- data.frame(GT1 = rownames(rh96.expr), stringsAsFactors = F)
pru.genes   <- data.frame(ME49 = rownames(pru.expr), stringsAsFactors = F)
me49.genes  <- data.frame(ME49 = rownames(me49.expr), stringsAsFactors = F)

## ortholog genes that exist in all datasets
GT1.ME49.Pb <- GT1.ME49.Pb %>% 
  dplyr::filter(GT1 %in% rh348.genes$GT1 & GT1 %in% rh96.genes$GT1 &
                  ME49 %in% pru.genes$ME49 & ME49 %in% me49.genes$ME49 &
                  PB %in% pb.genes$PB) 

## Identify genes with orthologs
pb.genes <- inner_join(pb.genes, GT1.ME49.Pb, by = 'PB')
rh348.genes <- inner_join(rh348.genes, GT1.ME49.Pb, by = 'GT1')
rh96.genes <- inner_join(rh96.genes, GT1.ME49.Pb, by = 'GT1')
pru.genes <- inner_join(pru.genes, GT1.ME49.Pb, by = 'ME49')
me49.genes <- inner_join(me49.genes, GT1.ME49.Pb, by = 'ME49')

## Take the subset of expression datasets
ortho.ind <- rownames(rh384.expr) %in% rh348.genes$GT1
rh384.expr.sub <- rh384.expr[ortho.ind, ]

ortho.ind <- rownames(rh96.expr) %in% rh96.genes$GT1
rh96.expr.sub <- rh96.expr[ortho.ind, ]

ortho.ind <- rownames(pru.expr) %in% pru.genes$ME49
pru.expr.sub <- pru.expr[ortho.ind, ]

ortho.ind <- rownames(me49.expr) %in% me49.genes$ME49
me49.expr.sub <- me49.expr[ortho.ind, ]

ortho.ind <- rownames(plasmodium.count) %in% pb.genes$PB
pb.expr.sub <- plasmodium.count[ortho.ind, ]

## Write all the genes using GT1 id
rownames(pru.expr.sub) <- pru.genes$GT1[match(rownames(pru.expr.sub), pru.genes$ME49)]
rownames(me49.expr.sub) <- me49.genes$GT1[match(rownames(me49.expr.sub), me49.genes$ME49)]
rownames(pb.expr.sub) <- pb.genes$GT1[match(rownames(pb.expr.sub), pb.genes$PB)]


## sort out species names and conditions
pru.obs$spp <- paste('PRU', gsub(' ', '', as.character(pru.obs$dpi)), sep = '')
pru.obs$NAME <- paste(pru.obs$spp, pru.obs$index, sep = '_')

me49.obs$spp <- paste('ME49', gsub(' ', '', as.character(me49.obs$dpi)), sep = '')
me49.obs$NAME <- paste(me49.obs$spp, me49.obs$index, sep = '_')

rh384.obs$spp <- 'RH384'
rh384.obs$NAME <- paste(rh384.obs$spp, rh384.obs$index, sep = '_')

rh96.obs$spp <- 'RH96'
rh96.obs$NAME <- paste(rh96.obs$spp, rh96.obs$index, sep = '_')

plasmodium.pheno$spp <- plasmodium.pheno$ShortenedLifeStage2
plasmodium.pheno$NAME <- paste(plasmodium.pheno$spp, plasmodium.pheno$sample_id, sep = '_')


rh384.ind <- match(colnames(rh384.expr.sub), rh384.obs$index)
rh96.ind <- match(colnames(rh96.expr.sub), rh96.obs$index)
pru.ind <- match(colnames(pru.expr.sub), pru.obs$index)
me49.ind <- match(colnames(me49.expr.sub), me49.obs$index)
pb.ind <- match(colnames(pb.expr.sub), plasmodium.pheno$sample_id)


colnames(rh384.expr.sub) <- rh384.obs$NAME[rh384.ind]
colnames(rh96.expr.sub) <- rh96.obs$NAME[rh96.ind]
colnames(pru.expr.sub) <- pru.obs$NAME[pru.ind]
colnames(me49.expr.sub) <- me49.obs$NAME[me49.ind]
colnames(pb.expr.sub) <- plasmodium.pheno$NAME[pb.ind]

all.inds <- c(rh384.ind, rh96.ind, pru.ind, me49.ind, pb.ind)
all.samples <- data.frame(Sample = rep(NA, length(all.inds)), 
                          spp = rep(NA, length(all.inds)), 
                          cell.cycle = rep(NA, length(all.inds)), stringsAsFactors = F)

all.samples$Sample <- c(as.character(rh384.obs$index[rh384.ind]),
                        as.character(rh96.obs$index[rh96.ind]),
                        as.character(pru.obs$index[pru.ind]),
                        as.character(me49.obs$index[me49.ind]),
                        as.character(plasmodium.pheno$sample_id[pb.ind]))

all.samples$spp <- c(as.character(rh384.obs$spp[rh384.ind]),
                     as.character(rh96.obs$spp[rh96.ind]),
                     as.character(pru.obs$spp[pru.ind]),
                     as.character(me49.obs$spp[me49.ind]),
                     as.character(plasmodium.pheno$spp[pb.ind]))

all.samples$cell.cycle <- c(gsub("\"", "", as.character(rh384.obs$cell_cycle[rh384.ind])),
                            rep('NA', length(rh96.ind)),
                            gsub("\"", "", as.character(pru.obs$cell_cycle[pru.ind])),
                            gsub("\"", "", as.character(me49.obs$cell_cycle[me49.ind])),
                            rep('NA', length(pb.ind)))




## functions
##
prep_S.O <- function(S.O){
  S.O <- NormalizeData(S.O, normalization.method = "LogNormalize", scale.factor = 10000)
  S.O <- FindVariableFeatures(S.O, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(S.O)
  S.O <- ScaleData(S.O, features = all.genes)
  S.O <- RunPCA(S.O, features = VariableFeatures(object = S.O))
  S.O <- FindNeighbors(S.O, dims = 1:13)
  S.O <- FindClusters(S.O, resolution = 0.1)
  S.O <- RunUMAP(S.O, dims = 1:13)
  return(S.O)
}

getNormExpr <- function(S.O){
  expr.norm <- as.matrix(S.O[["RNA"]]@data)
  rownames(expr.norm) <- gsub('-', '_',rownames(expr.norm))
  expr.norm <- as.data.frame(expr.norm) 
  expr.norm <- expr.norm %>% dplyr::mutate(GeneID = rownames(expr.norm))
  expr.norm <- expr.norm %>% gather(key = Sample, value = expr, -GeneID)
  expr.norm <- expr.norm  %>% group_by(GeneID) %>% 
    dplyr::mutate(quantile = ifelse(expr <= quantile(expr)[2], 1, 
                                    ifelse(expr <= quantile(expr)[3], 2, 
                                           ifelse(expr <= quantile(expr)[4], 3,4))))
  return(expr.norm)
}




S.O.RH348 <- CreateSeuratObject(counts = rh384.expr.sub, min.cells = 3, min.features = 200)
S.O.RH96 <- CreateSeuratObject(counts = rh96.expr.sub, min.cells = 3, min.features = 200)
S.O.PRU <- CreateSeuratObject(counts = pru.expr.sub, min.cells = 3, min.features = 200)
S.O.ME49 <- CreateSeuratObject(counts = me49.expr.sub, min.cells = 3, min.features = 200)
S.O.PB <- CreateSeuratObject(counts = pb.expr.sub, min.cells = 3, min.features = 200)



### Combine the datasets
S.O.combined <- merge(x = S.O.RH348, y = S.O.RH96)
S.O.combined <- merge(x = S.O.combined, y = S.O.PRU)
S.O.combined <- merge(x = S.O.combined, y = S.O.ME49)
S.O.combined <- merge(x = S.O.combined, y = S.O.PB)

RenameCells(S.O.combined, add.cell.id = c('RH348','RH96','PRU', 'ME49', 'PM'))

VlnPlot(S.O.combined, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

## unique(S.O.combined@meta.data$orig.ident)
cells.sub <- c('RH384', "Ring", "Schizont", "Trophozoite")
S.O.combined <- subset(S.O.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 )
S.O.combined <- subset(S.O.combined, subset = orig.ident %in% cells.sub)


S.O.combined <- prep_S.O(S.O.combined)
combined.expr <- getNormExpr(S.O.combined)

## colors
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

umap <- S.O.combined[['umap']]@cell.embeddings
umap <- as.data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) 
umap$cell.type <- unlist(lapply(strsplit(umap$Sample, split = '_'), `[[`, 1))
umap <- umap %>% ungroup() %>% dplyr::mutate_at(vars(-contains('UMAp')), funs(factor))
umap$SampleName <- unlist(lapply(as.character(umap$Sample), function(x) gsub('^([^_]+)_', '', x)))
umap$cell.cycle <- all.samples$cell.cycle[match(umap$SampleName, all.samples$Sample)]
#umap$cell.cycle[umap$cell.cycle == 'G1 a'] = 'G1'
#umap$cell.cycle[umap$cell.cycle == 'G1 b'] = 'G1'
p1 <- ggplot(umap, aes(x=UMAP_1,y=UMAP_2)) + 
  geom_point(aes(
    fill = cell.type, color = cell.cycle), shape=21, size = 2, stroke = 1.2) + 
  scale_fill_manual(breaks = unique(umap$cell.type), values = col_vector) + 
  theme_bw()

plot(p1)




pc <- S.O.combined[['pca']]@cell.embeddings
pc <- as.data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) 
pc$cell.type <- unlist(lapply(strsplit(pc$Sample, split = '_'), `[[`, 1))
pc <- pc %>% ungroup() %>% dplyr::mutate_at(vars(-contains('PC')), funs(factor))


pc$SampleName <- unlist(lapply(as.character(pc$Sample), function(x) gsub('^([^_]+)_', '', x)))

pc$cell.cycle <- all.samples$cell.cycle[match(pc$SampleName, all.samples$Sample)]

pc$cell.cycle[pc$cell.cycle == 'G1 a'] = 'G1'
pc$cell.cycle[pc$cell.cycle == 'G1 b'] = 'G1'


p2 <- ggplot(pc, aes(x=PC_2,y=PC_3)) + 
  geom_point(aes(
    fill = cell.type, color = cell.cycle), shape=21, size = 2, stroke = 1.2) + 
  scale_fill_manual(breaks = unique(pc$cell.type), values = col_vector) + 
  theme_bw()

plot(p2)

## Get highly variable genes
S.O.combined <- FindVariableFeatures(S.O.combined, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(S.O.combined), 10)
plot1 <- VariableFeaturePlot(S.O.combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


S.O.combined <- FindNeighbors(S.O.combined, dims = 1:10)
S.O.combined <- FindClusters(S.O.combined, resolution = 0.5)
DimPlot(S.O.combined, reduction = "umap")

sub.cells <- c('PRUDay0', 'Ring', 'Schizont', 'Trophozoite')
#sub.cells <- c('RH384', 'Schizont')
#sub.cells <- c('RH384', 'Trophozoite')

pc.filter <- pc %>% 
  dplyr::filter(pc$cell.type %in% sub.cells)


for(item in unique(as.character(pc.filter$cell.type))[-1]){
  pc.filter <- pc.filter %>% 
    mutate(!!item := ifelse(pc.filter$cell.type %in% 
                              c(unique(as.character(pc.filter$cell.type))[1], item), 1, 0))
}

pc.filter <- pc.filter %>% pivot_longer(cols = all_of(sub.cells[-1]), names_to = "comp", 
                                        values_to = "values")

pc.filter <- pc.filter %>% dplyr::filter(values == 1)

p1 <- ggplot(pc.filter, aes(x=PC_1,y=PC_2)) + 
  geom_point(aes(
    fill = cell.cycle), shape=21, size = 2) + 
  scale_fill_manual(breaks = unique(pc$cell.cycle), values = col_vector) + 
  facet_wrap(comp~.) + 
  theme_bw()

plot(p1)


p2 <- ggplot(pc.filter, aes(x=PC_2,y=PC_3)) + 
  geom_point(aes(
    fill = cell.cycle, color = cell.type), shape=21, size = 2, stroke = 1.1) + 
  scale_fill_manual(breaks = unique(pc$cell.cycle), values = col_vector) + 
  facet_grid(comp~.) + 
  theme_bw()

plot(p2)

## perform a PCA on genes

combined.expr$SampleName <- unlist(lapply(as.character(combined.expr$Sample), function(x) gsub('^([^_]+)_', '', x)))
combined.expr$spp <- unlist(lapply(as.character(combined.expr$Sample), function(x) gsub('_.*', '', x)))

combined.expr$cell.cycle <- all.samples$cell.cycle[match(combined.expr$SampleName, all.samples$Sample)]

highly.var <- combined.expr %>% group_by(GeneID) %>% summarise(var = sd(expr))
highly.var <- highly.var %>% 
  mutate(var.quantile = ifelse(var <= quantile(var)[2], 1,
                               ifelse(var <= quantile(var)[3], 2,
                                      ifelse(var <= quantile(var)[4], 3,4))))
combined.expr <- left_join(combined.expr, highly.var, by = 'GeneID')

combined.expr.filt <- combined.expr %>% dplyr::filter(var.quantile >=2)
combined.expr.filt.mat <- combined.expr.filt %>% 
  pivot_wider(-c(quantile:var.quantile), names_from = 'Sample', values_from = 'expr')
genes <- combined.expr.filt.mat[,1]

pr.genes <- prcomp(combined.expr.filt.mat[,-1], center = T, scale. = F)



## Cluster the data
XX <- combined.expr.filt.mat[,-1]
XX.scale <- scale(XX)
km.res <- kmeans(XX, centers = 5, iter.max = 50, nstart = 25)
length(km.res$cluster)
combined.expr.filt.mat$clusters <- km.res$cluster

library(umap)
gene.umap = umap(combined.expr.filt.mat[,-1])
gene.umap2d <- data.frame(gene.umap$layout)
colnames(gene.umap2d) <- c('UMAP1', 'UMAP2')
gene.umap2d$clusters <- as.factor(km.res$cluster)

p1 <- ggplot(gene.umap2d, aes(x=UMAP1,y=UMAP2)) + 
  geom_point(aes(
    fill = clusters), shape=21, size = 2) + 
  scale_fill_manual(breaks = unique(gene.umap2d$clusters), values = col_vector) + 
  theme_bw()

plot(p1)

pr.df <- data.frame(pr.genes$x)
pr.df$clusters <- as.factor(km.res$cluster)

p1 <- ggplot(pr.df, aes(x=PC2,y=PC3)) + 
  geom_point(aes(
    fill = clusters), shape=21, size = 2) + 
  scale_fill_manual(breaks = unique(pr.df$clusters), values = col_vector) + 
  theme_bw()

plot(p1)


prep_geneset.190 <- function(GeneSet.190){
  ## Make less lenghty names
  colnames(GeneSet.190) <- 
    gsub("development", "devel", 
         gsub('tgo.*[[:digit:]]\\.', '', 
              gsub('Tachyzoite', 'Tachy', 
                   gsub('Bradyzoite', 'Brady', 
                        gsub('tachyzoite', 'Tachy',  
                             gsub('bradyzoite', 'Brady',
                                  gsub('\\,', '', 
                                       gsub('\\(.*', '', 
                                            gsub('_', '.', 
                                                 gsub('-', '.', colnames(GeneSet.190)))))))))))
  
  GeneSet.190 <- GeneSet.190 %>% gather(key = GeneSet, value = GeneName) %>% na.omit() %>% 
    group_by(GeneSet) %>% summarise(genes = list(as.character(GeneName)), total = n())
  return(GeneSet.190)
}

fisherEnrichment <- function(GeneSet, GeneSet.list){
  ## This is a Hack to avoid nested for loops
  GeneSet$all <- 2932
  GeneSet.list$all <- 2932
  
  XX <- full_join(GeneSet, GeneSet.list, by = 'all')
  XX <- XX %>% rowwise() %>% mutate(overlap = length(intersect(c(genes.x), c(genes.y))),
                                    overlap.genes = list(intersect(c(genes.x), c(genes.y))))
  XX <- XX %>% rowwise() %>% 
    mutate(pvalue = fisher.test(matrix(c(overlap, total.x - overlap, total.y - overlap, 
                                         all - (total.x + total.y - overlap) ), byrow = T, ncol = 2, nrow = 2),
                                alternative = "greater")$p.value)
  
  return(XX)
}

GeneSet.190 <- read.xlsx('./Input/GeneSets/GSEA sytenic list plus extras 9.25.19.xlsx')
GeneSet.190 <- prep_geneset.190(GeneSet.190)
#GeneSet.190 <- GeneSet.190 %>% unnest(genes) %>% dplyr::select(-c('total'))
#colnames(GeneSet.190) <- c('GeneSet', 'GeneName')

## What are these genes enriched in
Candidates <- combined.expr.filt.mat %>% dplyr::select(GeneID, clusters)
GeneSet.list <- Candidates %>% group_by(clusters) %>% 
  summarise(cluster = clusters[1], genes = list(as.character(GeneID)), total = n())
  
FE <- fisherEnrichment(GeneSet.190, GeneSet.list)  

FE <- FE %>% dplyr::select(GeneSet, cluster, pvalue) 
  
FE.filt <- FE %>% dplyr::filter(pvalue < 0.01) %>% arrange(cluster, pvalue)

FE.filt.top2 <- FE.filt %>% group_by(cluster) %>% 
  summarise(top.GO = paste(GeneSet[1:min(n(), 2)], collapse = '/'))

exp.clust.df <- combined.expr.filt.mat %>% ungroup() %>% 
  dplyr::select(-GeneID) %>% group_by(clusters) %>%
  summarise_all(mean)

exp.clust.df <- left_join(exp.clust.df, FE.filt.top2, by = c('clusters' = 'cluster'))

exp.clust.df <- exp.clust.df %>% 
  pivot_longer(-c(clusters,top.GO), names_to = 'cells', values_to = 'ave_expr')
exp.clust.df$spp <-unlist(lapply(as.character(exp.clust.df$cells), function(x) gsub('_.*', '', x)))

exp.clust.df <- exp.clust.df %>% arrange(desc(clusters), spp, desc(ave_expr))
exp.clust.df.h <- exp.clust.df %>%
  mutate(clusters = factor(clusters, levels = unique(clusters)),
         cells = factor(cells, levels = unique(cells)),
         spp = factor(spp, levels = unique(spp)),
         top_GO = factor(top.GO, levels = unique(top.GO)),
         ave_expr = ave_expr)

hm.palette <- colorRampPalette(brewer.pal(9, 'Blues'), space='Lab')

p <- ggplot(exp.clust.df.h, aes(x = cells, y = as.numeric(clusters), fill = ave_expr)) + 
  geom_tile() + 
  scale_y_continuous(breaks = 1:length(levels(exp.clust.df.h$clusters)),
                     labels = levels(exp.clust.df.h$clusters),
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(levels(exp.clust.df.h$clusters)),
                                         labels = levels(exp.clust.df.h$top_GO))) +
  scale_x_discrete(expand=c(0,0)) +
  ylab("clusters") + xlab("cells") + 
  #scale_fill_gradientn(colours = hm.palette(100)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none") +
  facet_grid(.~spp, scales = "free_x")
plot(p)  
  

library(metacell)
library(Seurat)
library(viridis)
library(tidyverse)
#library(clustree)
library(factoextra)
#library(cluster)


## R code used for MetaCell Clustering

setwd("/SET/WORKING/DIRECTORY")

# Set up MetaCell folders
if(!dir.exists("./metacell_full_db/")) dir.create("./metacell_full_db/") 
scdb_init("./metacell_full_db/", force_reinit=T) 

if(!dir.exists("./metacell_full_figs")) dir.create("./metacell_full_figs") 
scfigs_init("./metacell_full_figs") 

# Import seurat object
seurat_data <- readRDS("./Data/CD8set_sct_int_noTCR(2).Rds")


########### PREP DATA ###################

# Check where cutoffs should be

table(seurat_data@meta.data$hash.ID)
HTOHeatmap(seurat_data, assay = "HTO", ncells = 5000)

to.check <- as_tibble(seurat_data@meta.data, rownames = "cell.code")

to.check %>%
  ggplot(aes(x = mito, fill = HTO_classification.global))+
  geom_density()+
  geom_vline(xintercept = 0.075)+
  facet_wrap(~HTO_classification.global, ncol = 1)

to.check %>%
  ggplot(aes(x = rp, fill = HTO_classification.global))+
  geom_density()+
  geom_vline(xintercept = 0.075)+
  facet_wrap(~HTO_classification.global, ncol = 1)


to.check %>%
  ggplot(aes(x = (nCount_RNA), fill = HTO_classification.global))+
  geom_density()+
  scale_x_log10()+
  geom_vline(xintercept = 2000)+
  facet_wrap(~HTO_classification.global, ncol = 1)

# Set Cutoffs for each patient

## NSCLC 36
to.check %>% 
  filter(HTO_classification.global == "Singlet") %>% 
  filter(mito < 0.075) %>% 
  filter(nCount_RNA > 2000) %>% 
  filter(orig.ident == "6671") %>% 
  pull(cell.code) -> N36

to.check %>% 
  filter(HTO_classification.global == "Singlet") %>% 
  filter(mito < 0.075) %>% 
  filter(nCount_RNA > 2000) %>% 
  filter(orig.ident == "6789_S1") %>% 
  filter(HTO_classification %in% c("HTO1", "HTO2", "HTO4")) %>% 
  pull(cell.code) %>% 
  c(N36) -> N36


## NSCLC 37
to.check %>% 
  filter(HTO_classification.global == "Singlet") %>% 
  filter(mito < 0.075) %>% 
  filter(nCount_RNA > 2000) %>% 
  filter(orig.ident == "6789_S1") %>% 
  filter(HTO_classification %in% c("HTO6", "HTO7", "HTO8")) %>% 
  pull(cell.code) -> N37

## NSCLC 20
to.check %>% 
  filter(HTO_classification.global == "Singlet") %>% 
  filter(mito < 0.075) %>% 
  filter(nCount_RNA > 2000) %>% 
  filter(orig.ident == "6817_S1") %>% 
  filter(HTO_classification %in% c("HTO1", "HTO2", "HTO4")) %>% 
  pull(cell.code) -> N20


## NSCLC 31
to.check %>% 
  filter(HTO_classification.global == "Singlet") %>% 
  filter(mito < 0.075) %>% 
  filter(nCount_RNA > 2000) %>% 
  filter(orig.ident == "6817_S1") %>% 
  filter(HTO_classification %in% c("HTO6", "HTO7", "HTO8")) %>% 
  pull(cell.code) -> N31

## NSCLC 41
to.check %>% 
  filter(HTO_classification.global == "Singlet") %>% 
  filter(mito < 0.075) %>% 
  filter(nCount_RNA > 2000) %>% 
  filter(orig.ident == "7109_S2") %>% 
  filter(HTO_classification %in% c("HTO1", "HTO2", "HTO4")) %>% 
  pull(cell.code) -> N41


## NSCLC 52
to.check %>% 
  filter(HTO_classification.global == "Singlet") %>% 
  filter(mito < 0.075) %>% 
  filter(nCount_RNA > 2000) %>% 
  filter(orig.ident == "7109_S2") %>% 
  filter(HTO_classification %in% c("HTO6", "HTO7", "HTO8")) %>% 
  pull(cell.code) -> N52

## NSCLC 50
to.check %>% 
  filter(HTO_classification.global == "Singlet") %>% 
  filter(mito < 0.075) %>% 
  filter(nCount_RNA > 2000) %>% 
  filter(orig.ident == "7114_S1") %>% 
  filter(HTO_classification %in% c("HTO1", "HTO2", "HTO4")) %>% 
  pull(cell.code) -> N50


## NSCLC 57
to.check %>% 
  filter(HTO_classification.global == "Singlet") %>% 
  filter(mito < 0.075) %>% 
  filter(nCount_RNA > 2000) %>% 
  filter(orig.ident == "7114_S1") %>% 
  filter(HTO_classification %in% c("HTO6", "HTO7", "HTO8")) %>% 
  pull(cell.code) -> N57

# Subset data

to.check %>% 
  select(cell.code) %>%
  mutate(patient = case_when(cell.code %in% N20 ~ "NSCLC_20",
                             cell.code %in% N31 ~ "NSCLC_31",
                             cell.code %in% N36 ~ "NSCLC_36",
                             cell.code %in% N37 ~ "NSCLC_37",
                             cell.code %in% N41 ~ "NSCLC_41",
                             cell.code %in% N50 ~ "NSCLC_50",
                             cell.code %in% N52 ~ "NSCLC_52",
                             cell.code %in% N57 ~ "NSCLC_57",
                             TRUE ~ "other")) %>% 
  deframe() -> patient.meta

seurat_data <- AddMetaData(object = seurat_data, metadata = patient.meta, col.name = "Patient")

seurat_data_filt <- subset(seurat_data, subset = Patient != "other")


## Generate metacell mat object
sce <- as.SingleCellExperiment(seurat_data_filt, assay = "RNA")
mat <- scm_import_sce_to_mat(sce)
scdb_add_mat(id = "NSCLC_data",  mat = mat)

## Clean-up
remove(sce)
remove(seurat_data)

write_rds(seurat_data_filt, "Data/CD8_NSCLC_filt.rds")

########### START METACELL PREP ###################

## Plot UMI distribution
mcell_plot_umis_per_cell("NSCLC_data", min_umis_cutoff = 2000)


## Get Bad genes
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
ig_genes = c(grep("^IGJ", nms, v=T), 
             grep("^IGH",nms,v=T),
             grep("^IGK", nms, v=T), 
             grep("^IGL", nms, v=T))

tcr_genes = c(grep("^TRAV", nms, v=T), 
             grep("^TRBV",nms,v=T),
             grep("^TRAJ", nms, v=T), 
             grep("^TCBJ", nms, v=T))

bad_genes = unique(c(tcr_genes, grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes))


## Ignore bad genes
mcell_mat_ignore_genes(new_mat_id="NSCLC_filt", mat_id="NSCLC_data", bad_genes, reverse=F) 


# Make gstat object
mcell_add_gene_stat(mat_id = "NSCLC_filt", gstat_id = "NSCLC_gs", force = T)
gstat <- scdb_gstat("NSCLC_gs")
dim(gstat)

# generate feats_gset and plot stats
mcell_gset_filter_varmean(gstat_id = "NSCLC_gs", gset_id = "NSCLC_feats", T_vm=0.12, force_new=T)
mcell_gset_filter_cov(gstat_id = "NSCLC_gs", gset_id = "NSCLC_feats", T_tot=100, T_top3=2)

# Check length feature genes and plot statistics
feats_gset <- scdb_gset("NSCLC_feats")
length(names(feats_gset@gene_set))
mcell_plot_gstats(gstat_id = "NSCLC_gs", "NSCLC_feats")


# calculate gene-gene corrections from feature gene list
gene.anchors = c('PCNA', 'TOP2A', 'TXN', 'HSP90AB1', 'FOS')
mcell_mat_rpt_cor_anchors(mat_id = "NSCLC_filt", gene_anchors = gene.anchors, cor_thresh = 0.1, gene_anti = c(),
                          tab_fn = "./Output/g2g_correlations_lat.txt", sz_cor_thresh = 0.2)


### Read correlation matrix and generate  gene-set to use for gene-clustering
gcor.mat <- read.table("./Output/g2g_correlations_lat.txt", header = T)
foc.genes <- apply(gcor.mat[,-1], 1, which.max)
gset <- gset_new_gset(sets = foc.genes, desc = "Cell cycle and stress correlated genes")
scdb_add_gset("All_corr_lateral_genes", gset)

# Check amount of clusters with elbow method
fviz_nbclust(t(gcor.mat), kmeans, method = "wss", k.max = 6) + theme_minimal() + ggtitle("the Elbow Method")
ggsave(filename = "./metacell_full_figs/elbow_lat.pdf", width = 5, height = 4)


### generate mat object containing only correlated genes
mcell_mat_ignore_genes(new_mat_id = "NSCLC_filt_lateral", mat_id = "NSCLC_filt",
                       ig_genes = names(foc.genes), reverse = T)

### cluster  genes and identify modules | Set here to 10 clusters
mcell_gset_split_by_dsmat(gset_id = "All_corr_lateral_genes" , mat_id = "NSCLC_filt_lateral", K = 20)

feats_gset <- scdb_gset("All_corr_lateral_genes")
mat = scdb_mat("NSCLC_filt_lateral")

# Plot heatmaps of gene-gene correlations per module
mcell_plot_gset_cor_mats(gset_id = "All_corr_lateral_genes", scmat_id = "NSCLC_filt_lateral")

# Manually inspect and annotate

## 1 = Loosely correlated
## 2 = Transcription stuff  loose correlated   
## 3 = Transcription, nucleus, EZH2
## 4 = Cell cycle, transcription
## 5 =  Potential cell cycle
## 6 = Chaperones, stress
## 7 = Large cell cycle
## 8 = Cell cycle
## 9 = Ribosomal
## 10 = Random
## 11 = DNA replication stuff
## 12 = T cell activation
## 13 = HEat shock proteins
## 14 = Immune related
## 15 = Histones
## 16 = Immune related
## 17 = Heatshock proteins
## 18 = Immune
## 19 = Immune
## 20 = Immune

# write out general annotations
# feats_gset <- scdb_gset("All_corr_diff_genes")
# enframe(feats_gset@gene_set, "gene", "cluster") %>% 
#   mutate(annotation = case_when( cluster %in% c(2) ~ "ribosome", 
#                                  cluster %in% c(5) ~ "translation",
#                                  cluster %in% c(7, 20) ~ "Cell cycle",
#                                  cluster %in% c(10, 19, 26,29,30,33,39,40,47,48) ~ "Immune genes",
#                                  cluster %in% c(12) ~ "respiration",
#                                  cluster %in% c(17,21,24) ~ "Fibroblast genes",
#                                  cluster %in% c(28) ~ "IFN response",
#                                  cluster %in% c(44) ~ "melanocyte",
#                                  cluster %in% c(49) ~ "Hemoglobin",
#                                  cluster %in% c(50) ~ "mitochondrium",
#                                  TRUE ~ "ambigious")) %>% 
#   write_tsv("./misc_data/gene_modules.tsv")

# set gene to remove for clustering to lateral
to.lateral <- names(feats_gset@gene_set[feats_gset@gene_set %in% c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17)])
to.lateral <- names(feats_gset@gene_set)
test <- rep(1, length(to.lateral))
names(test) <- to.lateral
scdb_add_gset("lateral",gset_new_gset(test, "lateral"))


# filter feats gset for lateral genes
lateral_gset = scdb_gset("lateral")
feats_gset <- scdb_gset("NSCLC_feats")

feats_gset = gset_new_restrict_gset(feats_gset, lateral_gset, inverse = T, desc = "cgraph gene markers w/o lateral genes")
scdb_add_gset(id = "NSCLC_feats_filt", gset = feats_gset)
feats_gset <- scdb_gset("NSCLC_feats_filt")


# Run MetaCell clustering
mcell_add_cgraph_from_mat_bknn(mat_id = "NSCLC_filt", gset_id = "NSCLC_feats_filt", graph_id = "NSCLC_graph", K=100, dsamp=T)
mcell_coclust_from_graph_resamp(coc_id = "NSCLC_coc", graph_id = "NSCLC_graph", min_mc_size=30, p_resamp=0.75, n_resamp=500)
mcell_mc_from_coclust_balanced(mc_id = "NSCLC_MC", coc_id =  "NSCLC_coc", mat_id = "NSCLC_filt", K=30, min_mc_size=60, alpha=2)

# Import clustering and add colors
mc <- scdb_mc("NSCLC_MC")
length(names(mc@mc))
length(names(mc@annots))
mc@colors <- viridis(length(mc@annots))
scdb_add_mc("NSCLC_MC",mc)

# Plot MetaCell 2D projection
mcell_mc2d_force_knn(mc2d_id = "NSCLC_mc2d" ,mc_id =  "NSCLC_MC", "NSCLC_graph", ignore_mismatch = T)
tgconfig::set_param("mcell_mc2d_height",800, "metacell")
tgconfig::set_param("mcell_mc2d_width",800, "metacell")
mcell_mc2d_plot("NSCLC_mc2d")


lfp = log2(mc@mc_fp)

lfp %>% 
  as.data.frame() %>% 
  as_tibble(rownames = "gene") -> lfp


mcell_mc_reorder_hc(mc_id = "NSCLC_MC", gene_left =  "SELL")

lfp %>% 
  filter(gene == "ENTPD1") %>% 
  pivot_longer(cols = !gene) %>% 
  arrange(value) %>% 
  pull(name) -> ordering.mc

mc.reord <- mc_reorder(mc, as.numeric(ordering.mc))
scdb_add_mc("NSCLC_MC_reordered", mc.reord)

# Plot MetaCell 2D projection
mcell_mc2d_force_knn(mc2d_id = "NSCLC_mc2d_r" ,mc_id =  "NSCLC_MC_reordered", "NSCLC_graph", ignore_mismatch = T)
tgconfig::set_param("mcell_mc2d_height",800, "metacell")
tgconfig::set_param("mcell_mc2d_width",800, "metacell")
mcell_mc2d_plot("NSCLC_mc2d_r")


###### HIERARCHICAL CLUSTERING #######

mc_hc = mcell_mc_hclust_confu(mc_id = "NSCLC_MC_reordered", graph_id = "NSCLC_graph")

mc_sup = mcell_mc_hierarchy(mc_id = "NSCLC_MC_reordered", mc_hc = mc_hc, T_gap = 0.04)

write_rds(mc_sup, "Data/CD8_NSCLC_Sup_MC.rds")

mcell_mc_plot_hierarchy(mc_id = "NSCLC_MC_reordered", graph_id = "NSCLC_graph", 
                        mc_order = mc_hc$order, 
                        sup_mcs = mc_sup, 
                        width = 3600, height = 7200, 
                        min_nmc=2, show_mc_ids = T)



###### quick plots #######

lfp = log2(mc.reord@mc_fp)

plt = function(gene1, gene2, lfp, colors) 
{
  plot(lfp[gene1, ], lfp[gene2, ], pch=21, cex=3, bg=colors, xlab=gene1, ylab=gene2)
  text(lfp[gene1, ], lfp[gene2, ], colnames(lfp))
  
}

genes1 = c('ENTPD1', 'ENTPD1', 'ENTPD1', 'ENTPD1')
genes2 = c('ITGAE', 'CXCL13', 'TCF7', 'PDCD1')
par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
for (i in seq_along(genes1)) {
  plt(gene1 = genes1[i], gene2 = genes2[i], lfp = lfp, colors = mc@colors)
}



genes1 = c('TCF7', 'IL2RA', 'GZMB', 'CD27')
genes2 = c('TOX', 'CXCL13', 'TBX21', 'GZMK')
par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
for (i in seq_along(genes1)) {
  plt(gene1 = genes1[i], gene2 = genes2[i], lfp = lfp, colors = mc@colors)
}



genes1 = c('ENTPD1', 'ENTPD1', 'ENTPD1', 'ENTPD1')
genes2 = c('ACACA', 'LDHB', 'LDHA', 'MAT2A')
par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
for (i in seq_along(genes1)) {
  plt(gene1 = genes1[i], gene2 = genes2[i], lfp = lfp, colors = mc@colors)
}




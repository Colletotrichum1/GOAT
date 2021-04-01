### load CAZYme, effector, SM backBone genes and plot heatmap for each of them
strain = "gd10"
sm.geneID <- read.table(paste0("/Users/alexwang/0data/0mango/inf_associated_gene/SM/",strain,"_SM.geneClass"), fill = T, header = T)[,1]
sm.type <- read.table(paste0("/Users/alexwang/0data/0mango/inf_associated_gene/SM/",strain,"_SM.geneClass"), fill = T, header = T)[,2]
effector.geneID <- read.table(paste0("/Users/alexwang/0data/0mango/inf_associated_gene/secretome/",strain,"_effector.geneID"))[,1]
cazy.geneID <- read.table(paste0("/Users/alexwang/0data/0mango/inf_associated_gene/cazy/",strain,"_cazy.geneID"),fill=T)[,1]

accessory.gene <- read.table(paste0("/Users/alexwang/0data/0mango/accessory/",strain,"_accessory.geneID"))[,1]
###### plot tpm Heatmap
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

plot.Tpm.heat <- function(strain, type, h=3) {
  tpm.df <- read.table(paste0("/Users/alexwang/0data/0mango/transcriptome/tpm/",strain,"_tpm.txt"), header = T)
  geneName <- rownames(tpm.df)
  if(strain %in% c("gd10", "yn56")){
    colnames(tpm.df) <- c("coni1", "coni2", "coni3", "hypha1", "hypha2", "hypha3",
                          "inf3d1", "inf3d2", "inf3d3", "inf5d1", "inf5d2", "inf5d3")
    tpm.df <- tpm.df %>% mutate(coni=(coni1+coni2+coni3)/3, hypha=(hypha1+hypha2+hypha3)/3,
                                inf3d=(inf3d1+inf3d2+inf3d3)/3, inf5d=(inf5d1+inf5d2+inf5d3)/3) %>% .[,13:16]
  }else if(strain == "gz15"){
    colnames(tpm.df) <- c("hypha1", "hypha2","inf5d1", "inf5d2", "inf5d3")
    tpm.df <- tpm.df %>% mutate(hypha=(hypha1+hypha2)/2,inf5d=(inf5d1+inf5d2+inf5d3)/3) %>% .[,1:5]
  }else{
    colnames(tpm.df) <- c("hypha1", "hypha2", "hypha3","inf5d1", "inf5d2", "inf5d3")
    tpm.df <- tpm.df %>% mutate(hypha=(hypha1+hypha2+hypha3)/3,inf5d=(inf5d1+inf5d2+inf5d3)/3) %>% .[,1:6]
  }
  rownames(tpm.df) <- geneName
  if(type=="effector"){
    mat <- as.matrix(tpm.df[rownames(tpm.df) %in% effector.geneID, ])
  }else if (type=="sm"){
    mat <- as.matrix(tpm.df[rownames(tpm.df) %in% sm.geneID, ])
  }else if(type=="cazy"){
    mat <- as.matrix(tpm.df[rownames(tpm.df) %in% cazy.geneID, ])
  }
  # scale
  mat <- t(apply(mat, 1, scale))
  mat.scaled <- na.omit(mat)
  idx <- which(rownames(mat) %in% rownames(mat.scaled))
  if (strain %in% c("gd10", "yn56")) {
    colnames(mat.scaled) <- c("Conidia", "Hypha", "3 dpi", "5 dpi")
  }else if(strain == "gz15"){
    colnames(mat.scaled) <- c("Hypha1", "Hypha2", "5 dpi1", "5 dpi2", "5 dpi3")
  }else{
    colnames(mat.scaled) <- c("Hypha1", "Hypha2", "Hypha3", "5 dpi1", "5 dpi2", "5 dpi3")
  }
  # pdf(file = paste0("~/0data/0mango/Figures/heamap/",strain,"_",type,".pdf"), width = 7, height = 10)
  ## group by accessory and core chromosomes
  g <- ifelse(rownames(mat.scaled) %in% accessory.gene, "Accessory", "Core")
  cc <- c(rgb(84,178,84,maxColorValue = 255),rgb(14,14,14,maxColorValue = 255), rgb(217,37,42, maxColorValue = 255))
  ### some global parameters can be set with ht_opt() method
  # ht_opt("heatmap_row_names_gp" = gpar(fontsize = 5))
  if(length(which(g == "Accessory")) == 1 || type=="sm"){
    if(type=="sm"){
      ra <- rowAnnotation(foo=anno_mark(at=1:nrow(mat.scaled), labels = sm.type[idx]), gp=gpar(fontsize=0.1))
      Heatmap(mat.scaled, cluster_rows = T, cluster_columns = F, column_title_side = "top", 
              show_heatmap_legend = T,column_names_rot = T,col=colorRamp2(c(-1.5, 0, 1.5), colors = cc), column_names_centered = T,
              border = T, row_names_gp = gpar(fontsize=4), right_annotation = ra,
              heatmap_legend_param = list(title="", at=c(-2,2), labels=c("low", "high"),legend_height = unit(2, "cm")))
    }else{
      Heatmap(mat.scaled, cluster_rows = T, cluster_columns = F, column_title_side = "top", 
              show_heatmap_legend = T,column_names_rot = T,col=colorRamp2(c(-1.5, 0, 1.5), colors = cc), column_names_centered = T,
              border = T, row_names_gp = gpar(fontsize=2),
              heatmap_legend_param = list(title="", at=c(-2,2), labels=c("low", "high"),legend_height = unit(2, "cm")))
    }
  }else{
    ### split accessory and core chromosomes
    h1 <- Heatmap(mat.scaled[which(g=="Core"),], cluster_rows = T, cluster_columns = F, column_title_side = "top", row_title = "Core",
                  show_heatmap_legend = T,column_names_rot = T,col=colorRamp2(c(min(mat.scaled),(min(mat.scaled) + max(mat.scaled))/2, max(mat.scaled)), colors = cc), column_names_centered = T,
                  border = T, row_names_gp = gpar(fontsize=2),
                  heatmap_legend_param = list(title="", at=c(-2,2), labels=c("low", "high"),legend_height = unit(2, "cm")))
    h2 <- Heatmap(mat.scaled[which(g=="Accessory"),], cluster_rows = T, cluster_columns = F,row_title = "Accessory",
                  show_heatmap_legend = F,column_names_rot = T,col=colorRamp2(c(min(mat.scaled),(min(mat.scaled) + max(mat.scaled))/2, max(mat.scaled)), colors = cc), column_names_centered = T,
                  border = T, height = unit(2, "cm"),row_names_gp = gpar(fontsize=4))
    heat.list <- h1 %v% h2
    draw(heat.list)
  }
  # dev.off()
}
# effector, sm, cazy
plot.Tpm.heat(strain = strain, type = "sm", h = 3)


######################################################################## plot orthogroups gene count
geneCount <- read.delim("/Users/alexwang/0data/0mango/transcriptome/ortholog_heat/cazy/cazy_count2.tsv", header = F, row.names = "V1")
# only keep orthogroups with different gene counts in each strains, remove conserved orthogroups
geneCount <- apply(geneCount, 1, function(x){
  if (length(unique(x))!=1) {
    return(as.vector(x))
  }
})
diff.og <- t(as.data.frame(Filter(Negate(is.null), geneCount)))
colnames(diff.og) <- c("fj11", "gd10", "gz15", "hn47", "qz3", "yn55", "yn56")

mat.scaled <- t(apply(diff.og, 1, scale))
colnames(mat.scaled) <- c("fj11", "gd10", "gz15", "hn47", "qz3", "yn55", "yn56")

cc <- c(rgb(8,186,255,maxColorValue = 255), "white", rgb(255,128,2, maxColorValue = 255))
# Heatmap(mat.scaled, cluster_rows = T, cluster_columns = F, column_title_side = "top", 
#         show_heatmap_legend = T,column_names_rot = T,col=colorRamp2(c(-1.5, 0, 1.5), colors = cc), column_names_centered = T,
#         border = T, row_names_gp = gpar(fontsize=2),
#         heatmap_legend_param = list(title="", at=c(-2,2), labels=c("low", "high"),legend_height = unit(2, "cm")))

### do not scale
Heatmap(diff.og, cluster_rows = T, cluster_columns = F, column_title_side = "top", 
        show_heatmap_legend = T,column_names_rot = T,col=colorRamp2(c(0,2,6), colors = cc), column_names_centered = T,
        border = T, row_names_gp = gpar(fontsize=2),
        heatmap_legend_param = list(title="Gene number", at=c(0,1,2,3,4,5,6), labels=c(0,1,2,3,4,5,6),legend_height = unit(3, "cm")))


######################################################### plot gene heatmap on accessory chromosomes

plot.heat.geneID <- function(strain, geneID){
  tpm.df <- read.table(paste0("/Users/alexwang/0data/0mango/transcriptome/tpm/",strain,"_tpm.txt"), header = T)
  geneID <- read.table(geneID)[,1]
  if(strain %in% c("gd10", "yn56")){
    ## get 5days dpi
    tpm.df <- tpm.df[rownames(tpm.df) %in% geneID, c(4,5,6,
                                                     10,11,12)]
  }else if(strain == "gz15") {
    tpm.df <- tpm.df[rownames(tpm.df) %in% geneID, c(1,2,3,
                                                     4,5)]
  }else{
    tpm.df <- tpm.df[rownames(tpm.df) %in% geneID, c(1,2,3,
                                                     4,5,6)]
  }
  mat <- as.matrix(tpm.df)
  mat.scaled <- t(apply(mat, 1, scale))
  if(strain == "gz15"){
    colnames(mat.scaled) <- c("Hypha-1", "Hypha-2", "5 dpi-1", "5 dpi-2", "5 dpi-3")
    
  }else{
    colnames(mat.scaled) <- c("Hypha-1", "Hypha-2", "Hypha-3", "5 dpi-1", "5 dpi-2", "5 dpi-3")
    
  }
  cc <- c(rgb(84,178,84,maxColorValue = 255),rgb(14,14,14,maxColorValue = 255), rgb(217,37,42, maxColorValue = 255))
  Heatmap(mat.scaled, cluster_rows = T, cluster_columns = F, column_title_side = "top", 
          show_heatmap_legend = T,column_names_rot = T,col=colorRamp2(c(-1.5, 0, 1.5), colors = cc), column_names_centered = T,
          border = T, row_names_gp = gpar(fontsize=10),
          heatmap_legend_param = list(title="", at=c(-2,2), labels=c("low", "high"),legend_height = unit(2, "cm")))
}

## plot all heat map
for(i in c("fj11", "yn55", "gd10", "gz15", "hn47", "yn56")){
  # pdf(paste0(i, "_mini_sp.pdf"))
  ppi=720
  png(paste0(i, "_mini_sp.png"),width = 6*ppi, height = 10*ppi, res = ppi)
  strain = i
  geneID = paste0("/Users/alexwang/0data/0mango/transcriptome/DEG/附属染色体上的差异表达基因/",strain,"_inf5d_Vs_hypha.geneID.accessory.specific")
  print(plot.heat.geneID(strain, geneID))
  dev.off()
}

########################################################################## plot YN56, GD10 all stage heatmap inf 3d up-regulated

inf3.up.yn56.coni <- read.table("/Users/alexwang/0data/0mango/transcriptome/DEG/yn56_inf3d_Vs_coni.up.geneID")[,1]
inf3.up.yn56.hypha <- read.table("/Users/alexwang/0data/0mango/transcriptome/DEG/yn56_inf3d_Vs_hypha.up.geneID")[,1]

inf3.up.yn56 <- intersect(inf3.up.yn56.coni, inf3.up.yn56.hypha)

strain = "yn56"
tpm.df <- read.table(paste0("/Users/alexwang/0data/0mango/transcriptome/tpm/",strain,"_tpm.txt"))
accessory.gene <- read.table(paste0("/Users/alexwang/0data/0mango/accessory/",strain,"_accessory.geneID"))[,1]
mat <- as.matrix(tpm.df[rownames(tpm.df) %in% inf3.up.yn56, ])
mat <- t(apply(mat, 1, scale))
mat.scaled <- na.omit(mat)
colnames(mat.scaled) <- c("Coni1", "Coni2","Coni3","Hypha1", "Hypha2", "Hypha3", "3 dpi1", "3 dpi2", "3 dpi3", "5 dpi1", "5 dpi2", "5 dpi3")
### plot
g <- ifelse(rownames(mat.scaled) %in% accessory.gene, "*", "")
# accessory annotation
ha1 <- rowAnnotation(foo = anno_text(g, gp = gpar(fontsize = 4)))
# 5d up regulated annotation
inf5.up <- read.table(paste0("/Users/alexwang/0data/0mango/transcriptome/DEG/",strain,"_inf5d_Vs_inf3d.up.geneID"))[,1]
inf5.down <- read.table(paste0("/Users/alexwang/0data/0mango/transcriptome/DEG/",strain,"_inf5d_Vs_inf3d.down.geneID"))[,1]
idx.up <- which(rownames(mat.scaled) %in% inf5.up)
idx.down <- which(rownames(mat.scaled) %in% inf5.down)

ha2 <- rowAnnotation(foo = anno_mark(at = c(idx.up, idx.down), 
                                     labels = c(rep("up", length(idx.up)), 
                                                rep("down", length(idx.down)))),
                                      gp = gpar(fontsize = 1, lwd=0.2))

cc <- c(rgb(84,178,84,maxColorValue = 255),rgb(14,14,14,maxColorValue = 255), rgb(217,37,42, maxColorValue = 255))
# rownames(mat.scaled)[277]
ppi=720
png(paste0(strain, "_complexHeatmap.png"),width = 12*ppi, height = 35*ppi, res = ppi)
Heatmap(mat.scaled, cluster_rows = T, cluster_columns = F, column_title_side = "top",column_names_rot = T,
            col=colorRamp2(c(min(mat.scaled),(min(mat.scaled) + max(mat.scaled))/2, max(mat.scaled)), colors = cc), column_names_centered = T,
              border = T, row_names_gp = gpar(fontsize=2),
              # annotation
              left_annotation = ha1,
              right_annotation = ha2,
              # legend
              show_heatmap_legend = T, heatmap_legend_param = list(title="", at=c(-2,2), labels=c("low", "high"),legend_height = unit(1, "cm")))
dev.off()
###################################### all genes on mini chromosomes
library(dplyr)
strain = "gd10"
tpm.df <- read.table(paste0("/Users/alexwang/0data/0mango/transcriptome/tpm/",strain,"_tpm.txt"))
colnames(tpm.df) <- c("coni1", "coni2", "coni3", "hypha1", "hypha2", "hypha3",
                      "inf3d1", "inf3d2", "inf3d3", "inf5d1", "inf5d2", "inf5d3")
geneName <- rownames(tpm.df)
tpm.df <- tpm.df %>% mutate(coni=(coni1+coni2+coni3)/3, hypha=(hypha1+hypha2+hypha3)/3,
                            inf3d=(inf3d1+inf3d2+inf3d3)/3, inf5d=(inf5d1+inf5d2+inf5d3)/3) %>% .[,13:16]
rownames(tpm.df) <- geneName
# select accessory genes
accessory.gene <- read.table(paste0("/Users/alexwang/0data/0mango/accessory/",strain,"_accessory.geneID"))[,1]
mat <- as.matrix(tpm.df[rownames(tpm.df) %in% accessory.gene, ])
mat <- t(apply(mat, 1, scale))
mat.scaled <- na.omit(mat)
# colnames(mat.scaled) <- c("Coni1", "Coni2","Coni3","Hypha1", "Hypha2", "Hypha3", "3 dpi1", "3 dpi2", "3 dpi3", "5 dpi1", "5 dpi2", "5 dpi3")
colnames(mat.scaled) <- c("Coni", "Hypha", "3 dpi", "5 dpi")
## read up-regulated geneID   3d vs coni
inf3.up <- read.table(paste0("/Users/alexwang/0data/0mango/transcriptome/DEG/",strain,"_inf3d_Vs_coni.up.geneID"))[,1]
idx1 <- which(rownames(mat.scaled) %in% inf3.up)
# 5d vs coni
inf5.up <- read.table(paste0("/Users/alexwang/0data/0mango/transcriptome/DEG/",strain,"_inf5d_Vs_coni.up.geneID"))[,1]
idx2 <- which(rownames(mat.scaled) %in% inf5.up)

ha1 <- rowAnnotation(foo = anno_mark(at = c(idx1, idx2), 
                                     labels = c(rep("inf3d up", length(idx1)), 
                                                rep("inf5d up", length(idx2)))),
                     gp = gpar(fontsize = 2, lwd=0.5))

cc <- c(rgb(84,178,84,maxColorValue = 255),rgb(14,14,14,maxColorValue = 255), rgb(217,37,42, maxColorValue = 255))
ppi=720
png(paste0(strain, "_complexHeatmap_mini_mergeReplicates.png"),width = 12*ppi, height = 20*ppi, res = ppi)
Heatmap(mat.scaled, cluster_rows = T, cluster_columns = F, column_title_side = "top",column_names_rot = T,
        col=colorRamp2(c(min(mat.scaled),(min(mat.scaled) + max(mat.scaled))/2, max(mat.scaled)), colors = cc), column_names_centered = T,
        border = T, row_names_gp = gpar(fontsize=3), 
        # annotation
        right_annotation = ha1,
        # legend
        show_heatmap_legend = T, heatmap_legend_param = list(title="", at=c(-2,2), labels=c("low", "high"),legend_height = unit(1, "cm")))
dev.off()
### only display up regulated genes
mat.up <- rbind(mat.scaled[rownames(mat.scaled) %in% inf3.up,, drop=FALSE], mat.scaled[rownames(mat.scaled) %in% inf5.up, ])
mat.up <- mat.up[!duplicated(mat.up), ]
ppi=720
png(paste0(strain, "_complexHeatmap_mini_mergeReplicates_upregulated.png"),width = 12*ppi, height = 20*ppi, res = ppi)
Heatmap(mat.up, cluster_rows = T, cluster_columns = F, column_title_side = "top",column_names_rot = T,
        col=colorRamp2(c(min(mat.up),(min(mat.up) + max(mat.up))/2, max(mat.up)), colors = cc), column_names_centered = T,
        border = T, row_names_gp = gpar(fontsize=6), 
        # legend
        show_heatmap_legend = T, heatmap_legend_param = list(title="", at=c(-2,2), labels=c("low", "high"),legend_height = unit(1, "cm")))
dev.off()
#################################################################
## FJ11 YN56 comparasion
tpm.df <- read.table("/Users/alexwang/0data/0mango/accessory/yn55_fj11_anchors/yn55_fj11.anchors.tpm.txt")
colnames(tpm.df) <- c("yn55", "yn55_Rep1", "yn55_Rep2", "yn55_Rep3", "fj11", "fj11_rep1", "fj11_rep2", "fj11_rep3")
rownames(tpm.df) <- paste0(tpm.df$yn55, "-",tpm.df$fj11)
tpm.df <- tpm.df[, c(2,3,4,6,7,8)]
mat <- t(apply(tpm.df, 1, scale))
mat.scaled <- na.omit(mat)
csep <- read.table("/Users/alexwang/0data/0mango/inf_associated_gene/CSEP/fj11.CSEP.geneID")[,1]
mat.csep.id <- unlist(lapply(stringi::stri_split(rownames(mat.scaled), regex = "-"), function(x) x[2]))
idx<- which( mat.csep.id %in% csep)

ha1 <- rowAnnotation(foo = anno_mark(at = idx, 
                                     labels = rep("CSEP", length(idx))),
                     gp = gpar(fontsize = 2, lwd=0.2))

cc <- c(rgb(84,178,84,maxColorValue = 255),rgb(14,14,14,maxColorValue = 255), rgb(217,37,42, maxColorValue = 255))
cc = RColorBrewer::brewer.pal(3,"PRGn")
pdf("fj11-yn55_comparasion.pdf",width = 16, height = 12)
Heatmap(mat.scaled, cluster_rows = F, cluster_columns = F, column_names_side = "bottom",column_names_rot = T,row_names_side = "left",
        col=colorRamp2(c(min(mat.scaled), (min(mat.scaled)+max(mat.scaled))/2, max(mat.scaled)), colors = cc), column_names_centered = T,
        border = T, row_names_gp = gpar(fontsize=10),
        # annotation
        right_annotation = ha1,
        # legend
        show_heatmap_legend = T, heatmap_legend_param = list(title="", at=c(-2,2), labels=c("low", "high"),legend_height = unit(1, "cm")))
dev.off()



\


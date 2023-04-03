# load packages
rm(list=ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(pheatmap)
library(plyr)
library(data.table)
library(ggplot2)
library(AnnotationHub)
library(org.Bt.eg.db)
library(clusterProfiler)
library(GenomicFeatures)
library(DESeq2)
library(ggpubr)
library(ComplexHeatmap)

#=========================================================================================================
#=================================== load all bovine scRNA-seq data ========================================
setwd("/home/bovine/merge/blastoid")
exp1<-Read10X(data.dir="filtered_feature_bc_matrix/") # load cellranger files
# load blastoid 10x scRNA-seq data in our research 
cow1<-CreateSeuratObject(counts = exp1, min.cells = 3, min.features = 300,project = "Blastoid")
cow1@meta.data$tech<-"10xseq"
cow1@meta.data$celltype<-"Blastoid"
cow1@meta.data$data<-"data0"
cow1[["percent.mt"]] <- PercentageFeatureSet(cow1, pattern = "^MT-")
cow1@meta.data$sample<-as.matrix(cow1@active.ident)                
cow1 <- subset(x = cow1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
cow1 <- subset(x = cow1, subset = nCount_RNA > 5000 & nCount_RNA < 30000)
# load data1
exp1<-read.table("data1/data1_TPM_out.txt",header = T)
meta1<-read.table("data1/mata-data1.txt",header = T,sep = "\t")
rownames(meta1)<-meta1$cellname
meta1<-meta1[,-1]
en2sy<-read.table("jie/cow/Ensembl2symble_bovine.txt")
colnames(en2sy)<-c("Ensemblid","symbol")
en2sy<-as.data.frame(en2sy)
fid=as.character(rownames(exp1))
geneLists=data.frame(Ensemblid=fid)
results_up=merge(geneLists,en2sy,by='Ensemblid',all.x=T)
results_up$symbol<-as.matrix(results_up$symbol)
results_up$Ensemblid<-as.matrix(results_up$Ensemblid)
results_up$symbol[is.na(results_up$symbol)] = results_up$Ensemblid[is.na(results_up$symbol)]
results_up<-results_up[order(results_up$Ensemblid),]
exp1<-exp1[order(rownames(exp1)),]
exp1$symbol<-results_up$symbol
genelist<-unique(exp1$symbol[duplicated(exp1$symbol)])
for (i in 1:length(genelist)) {
  iindex=which(exp1$symbol==genelist[i])
  exp1[iindex[1],1:169]<-colSums(exp1[iindex,1:169])
  exp1<-exp1[-iindex[-1],]
  print(i)
}
rownames(exp1)<-exp1$symbol
exp1<-exp1[,-170]
meta1<-meta1[rownames(meta1)%in%colnames(exp1),]
data1<-CreateSeuratObject(counts = exp1, meta.data =meta1, project = "data1_29511234")
data1[["percent.mt"]] <- PercentageFeatureSet(data1, pattern = "^MT-")
data1@meta.data$sample<-as.matrix(data1@active.ident)                
data1 <- subset(x = data1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
# load data2
exp2<-read.table("data2/data2_TPM_out.txt",header = T,sep = "\t")
meta2<-read.table("data2/meta-data2.txt",header = T,sep = "\t")
rownames(meta2)<-meta2$cellname
meta2<-meta2[,-1]
fid=as.character(rownames(exp2))
geneLists=data.frame(Ensemblid=fid)
results_up=merge(geneLists,en2sy,by='Ensemblid',all.x=T)
results_up$symbol<-as.matrix(results_up$symbol)
results_up$Ensemblid<-as.matrix(results_up$Ensemblid)
results_up$symbol[is.na(results_up$symbol)] = results_up$Ensemblid[is.na(results_up$symbol)]
results_up<-results_up[order(results_up$Ensemblid),]
exp2<-exp2[order(rownames(exp2)),]
exp2$symbol<-results_up$symbol
genelist<-unique(exp2$symbol[duplicated(exp2$symbol)])
for (i in 1:length(genelist)) {
  iindex=which(exp2$symbol==genelist[i])
  exp2[iindex[1],1:(ncol(exp2)-1)]<-colSums(exp2[iindex,1:(ncol(exp2)-1)])
  exp2<-exp2[-iindex[-1],]
  print(i)
}
rownames(exp2)<-exp2$symbol
exp2<-exp2[,-ncol(exp2)]
meta2<-meta2[rownames(meta2)%in%colnames(exp2),]
data2<-CreateSeuratObject(counts = exp2, meta.data =meta2,project = "data2_Jiang")
data2[["percent.mt"]] <- PercentageFeatureSet(data2, pattern = "^MT-")
data2@meta.data$sample<-as.matrix(data2@active.ident)                
data2 <- subset(x = data2, subset =  nFeature_RNA > 200 & percent.mt < 20)
# load data3
countlist<-list.files("countfile")
for (i in 1:length(countlist)) {
  counti<-read.table(file = paste("countfile/",countlist[i],sep = ""),header = F)
  counti<-head(counti,-5)
  colnames(counti)<-c("geneid",strsplit(countlist[i],"_")[[1]][1])
  if (i==1) {
    data3 = counti
  }else{
    data3=merge(x = data3, y = counti, by = "geneid")
  }
  print(countlist[i])
}
write.table(x=data3,file = "data3_Count_out.txt",sep = "\t",row.names = F)
txdb <- makeTxDbFromGFF("../Bos_taurus.UMD3.1.94.gtf", format = "gtf", circ_seqs = character())
ebg <- exonsBy(txdb, by="gene")
exon_gene_sizes <- sum(width(reduce(ebg)))
rownames(data3)<-data3$geneid
data3<-data3[,-1]
data3<-data3[rowSums(data3)>0,]
genelen=exon_gene_sizes[names(exon_gene_sizes)%in%rownames(data3)] 
if (identical(names(genelen),rownames(data3))) {
  tmp <- data3 / genelen
  TPM <- 1e6 * t(t(tmp) / colSums(tmp))
}
colSums(tpms)
write.table(as.data.frame(TPM), file="data3_TPM_out.txt",sep="\t",row.names = F)
exp3<-read.table("data3/data3_TPM_out.txt",header = T,sep = "\t")
meta3<-read.table("data3/meta-data3.txt",header = T,sep = "\t")
rownames(meta3)<-meta3$cellname
meta3<-meta3[,-1]
en2sy<-as.data.frame(en2sy)
fid=as.character(rownames(exp3))
geneLists=data.frame(Ensemblid=fid)
results_up=merge(geneLists,en2sy,by='Ensemblid',all.x=T)
results_up$symbol<-as.matrix(results_up$symbol)
results_up$Ensemblid<-as.matrix(results_up$Ensemblid)
results_up$symbol[is.na(results_up$symbol)] = results_up$Ensemblid[is.na(results_up$symbol)]
results_up<-results_up[order(results_up$Ensemblid),]
exp3<-exp3[order(rownames(exp3)),]
rownames(exp3)<-results_up$symbol
exp3$symbol<-results_up$symbol
genelist<-unique(exp3$symbol[duplicated(exp3$symbol)])
for (i in 1:length(genelist)) {
  iindex=which(exp3$symbol==genelist[i])
  exp3[iindex[1],1:98]<-colSums(exp3[iindex,1:98])
  exp3<-exp3[-iindex[-1],]
  print(i)
}
rownames(exp3)<-exp3$symbol
exp3<-exp3[,-99]
meta3<-meta3[rownames(meta3)%in%colnames(exp3),]
colnames(meta3)<-c("data","celltype")
data3<-CreateSeuratObject(counts = exp3, meta.data =meta3, project = "data3_35971640")
data3[["percent.mt"]] <- PercentageFeatureSet(data3, pattern = "^MT-")
data3@meta.data$sample<-as.matrix(data3@active.ident)                
data3 <- subset(x = data3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

#=========================================================================================================
#================================= merge all data and quanlity control ===================================
cow_all<-merge(x = cow1,y = c(data1,data2,data3))
cowall.list <- SplitObject(cow_all, split.by = "data")
for (i in 1:length(x = cowall.list)) {
  cowall.list[[i]] <- NormalizeData(object = cowall.list[[i]], verbose = FALSE)
  cowall.list[[i]] <- FindVariableFeatures(object = cowall.list[[i]],
                                                  selection.method = "vst", verbose = FALSE)
  print(i)
}
cow_all.anchors <- FindIntegrationAnchors(object.list = cowall.list, dims = 1:80,k.anchor = 20,k.filter = 40, reduction = "rpca") 
cow_all <- IntegrateData(anchorset = cow_all.anchors, dims = 1:80)          
cow_all <- ScaleData(cow_all)                                               
cow_all <- RunPCA(object = cow_all,npcs=100)
pdf(file = "cow_all_pca_use.pdf")
ElbowPlot(cow_all, ndims = 100)
dev.off()
# tsne choose best dims and pca
pcause=90
orig.ident<-cow_all@active.ident
cow_all <- FindNeighbors(object = cow_all, dims = 1:pcause)                         
cow_all <- FindClusters(object = cow_all, resolution = 0.6)                         
cow_all <- RunUMAP(object = cow_all, dims = 1:pcause)                               
pdf(file = paste("cow_all_tsne_cluster",pcause,".pdf",sep = ""),width = 8,height = 6)
DimPlot(object = cow_all,reduction = "umap",label = T,pt.size = 1.3)
dev.off()
cow_all@meta.data$cluster<-cow_all@active.ident
pdf(file = paste("cow_all_tsne_sample",pcause,".pdf",sep = ""),width = 12,height = 8)
DimPlot(object = cow_all,reduction = "umap",group.by = "data",pt.size = 1.3)
dev.off()
cow_all@meta.data$celltypeless=cow_all$celltype
pdf(file = paste("cow_all_tsne_celltype",pcause,".pdf",sep = ""),width = 12,height = 8)
DimPlot(object = cow_all,reduction = "umap",group.by = "celltype",pt.size = 1.3)
dev.off() 
# data3 IVF_cells annotation
cow_all@meta.data$celltype2 = cow_all$celltype
cow_all$celltype2[cow_all$data=="data3_35971640"&cow_all$celltype=="IVF_blastocyst"] = "data3_IVF_blastocyst"
# Figure 2A
pdf(file = paste("Figure_2A",pcause,".pdf",sep = ""),width = 11,height = 8)
DimPlot(object = cow_all,reduction = "umap",group.by = "celltype2",pt.size = 1.0,
           cols = color_bar)
dev.off()
cow_all@meta.data$ori_celltype = cow_all@active.ident

#=========================================================================================================
#=============================== cell type annotation through marker genes  ==============================
markergenes<-read.csv("Three germ layers markers.csv",header = T)
markers<-as.vector(as.matrix(markergenes))
markers<-unique(markers)
markers<-markers[!markers==""]
data.usage <- DotPlot(cow_all,features = markers, assay = 'RNA')$data
data.anno<-cbind(markers,"celltype")
for (i in 1:ncol(markergenes)) {
  data.anno[markers%in%markergenes[,i],2]=colnames(markergenes)[i]
}
colnames(data.anno)<-c("features.plot","label")
data.anno<-as.data.frame(data.anno)
df.plot <- plyr::join(data.usage,data.anno)
df.plot<-na.omit(df.plot)
#cluster
data.cluster<-"a"
for (i in levels(df.plot$id)){
  i_list<-c(data.usage$avg.exp.scaled[data.usage$id == i],data.usage$pct.exp[data.usage$id == i])
  if (data.cluster=="a") {
    data.cluster=i_list
  }else{
    data.cluster=rbind(data.cluster,i_list)
  }
}
row.names(data.cluster)<-levels(df.plot$id)
d <- dist(data.cluster, method = "euclidean")
fit2 <- hclust(d, method="ward.D")
pdf(file = "cow_all_markergene_cluster.pdf",width = 8,height =6 )
plot(fit2, col = "black", col.main = "#45ADA8", col.lab = "#7C8071",
     col.axis = "#F38630", lwd = 3, lty = 1, 
     axes = F, hang = -1)
dev.off()
cluster_order<-row.names(data.cluster)[fit2$order]
df.plot$id <- factor(df.plot$id,levels = cluster_order)
p <- ggplot(df.plot,aes(x=features.plot,y = as.numeric(id),size = (pct.exp^2)/100, color = avg.exp.scaled))+
  geom_point() + scale_size("% detected", range = c(0,6)) +
  scale_color_gradientn(colours = c("white","red"),
                        guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                        name = "Average\nexpression") + cowplot::theme_cowplot() +
  ylab("") + xlab("Markers") + theme_bw() + 
  scale_y_continuous(breaks = 1:length(levels(df.plot$id)),labels = levels(df.plot$id),sec.axis = dup_axis())+
  facet_grid(~label, scales="free_x",space = "free")+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf(file = "cow_all_marker_dotplot.pdf",width = 14,height = 6)
print(p)
dev.off()

#=========================================================================================================
#============================= Count DEGs and cell numbers in each cell type  ============================
cow_all.markers <- FindAllMarkers(object = cow_all, only.pos = F)
cow_all.markers<-cow_all.markers[cow_all.markers$p_val_adj<0.05,]
cow_all.markers %>% group_by(cluster) %>% top_n(25, avg_logFC) -> top25_up
cow_all.markers %>% group_by(cluster) %>% top_n(25, -avg_logFC) -> top25_down
top25_up<-top25_up[top25_up$avg_logFC>0,]
top25_down<-top25_down[top25_down$avg_logFC<0,]
write.table(x=top25_up,file = "Bovine_merged_top25_UP_DEGs.txt",sep = "\t",row.names = F)
write.table(x=top25_down,file = "Bovine_merged_top25_Down_DEGs.txt",sep = "\t",row.names = F)
# show marker genes expression in each cell types
markers_select = c("SOX2","SOX17","GATA2")
markers = markers_select
markers<-unique(markers)
markers<-markers[!markers==""]
pdf(file = "cow_all_original_marker_annotation_selectmarkers.pdf",width = 12,height = 6)
DotPlot(cow_all, features = markers, assay = 'RNA', cols = c('white', 'red')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
dev.off()
data.usage <- DotPlot(cow_all,features = markers, assay = 'RNA')$data
data.anno<-cbind(markers,"celltype")
for (i in 1:ncol(markergenes)) {
  data.anno[markers%in%markergenes[,i],2]=colnames(markergenes)[i]
}
colnames(data.anno)<-c("features.plot","label")
data.anno<-as.data.frame(data.anno)
df.plot <- plyr::join(data.usage,data.anno)
df.plot<-na.omit(df.plot)
data.cluster<-"a"
for (i in levels(df.plot$id)){
  i_list<-c(data.usage$avg.exp.scaled[data.usage$id == i],data.usage$pct.exp[data.usage$id == i])
  if (data.cluster=="a") {
    data.cluster=i_list
  }else{
    data.cluster=rbind(data.cluster,i_list)
  }
}
row.names(data.cluster)<-levels(df.plot$id)
d <- dist(data.cluster, method = "euclidean")
fit2 <- hclust(d, method="ward.D")
# Figure 2B
pdf(file = "Figure_2B.pdf",width = 8,height =6 )
plot(fit2, col = "black", col.main = "#45ADA8", col.lab = "#7C8071",
     col.axis = "#F38630", lwd = 3, lty = 1, 
     axes = F, hang = -1)
dev.off()
cluster_order<-row.names(data.cluster)[fit2$order]
df.plot$id <- factor(df.plot$id,levels = cluster_order)
p <- ggplot(df.plot,aes(x=features.plot,y = as.numeric(id),size = (pct.exp^2)/100, color = avg.exp.scaled))+
  geom_point() + scale_size("% detected", range = c(0,6)) + 
  scale_color_gradientn(colours = c("white","red"),
                        guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                        name = "Average\nexpression") + cowplot::theme_cowplot() +
  ylab("") + xlab("Markers") + theme_bw() + 
  scale_y_continuous(breaks = 1:length(levels(df.plot$id)),labels = levels(df.plot$id),sec.axis = dup_axis())+
  facet_grid(~label, scales="free_x",space = "free")+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf(file = "cow_all_marker_dotplot_selectmarkers.pdf",width = 14,height = 6)
print(p)
dev.off()
# merged cells percent in each cell types
# Figure 2F
cow_all@meta.data$celltype_pre<-cow_all$celltype
cow_all$celltype_pre<-plyr::mapvalues(x = cow_all$celltype_pre, from = c("IVF_2 cell","IVF_16 cell","IVF_8 cell","IVF_zygote"), to = rep("Pre_lineage",4))
cow_merge_bar_all<-cbind(as.matrix(cow_all@active.ident),as.matrix(cow_all$celltype_pre),as.matrix(cow_all$data))
cow_merge_bar1<-cow_merge_bar_all[cow_merge_bar_all[,1]%in%c(2,3,6),]
cow_merge_bar2<-cow_merge_bar_all[cow_merge_bar_all[,2]%in%c("data2_IVF_blastocyst","data3_IVF_blastocyst"),]
celltypelist<-names(table(cow_merge_bar1[,2]))
for (i in levels(cow_all@active.ident)) {
  sign<-table(cow_merge_bar_all[cow_merge_bar_all[,1]==i,2])
  singletype_table<-cbind(rep(i,length(sign)),names(sign),as.matrix(sign))
  rownames(singletype_table)<-NULL
  if (i==0) {
    table_type<-singletype_table
  }else{
    table_type<-rbind(table_type,singletype_table)
  }
}
colnames(table_type)<-c("Celltype","Cancertype","Cellnumber")
table_sample_type<-as.data.frame(table_type)
table_sample_type$Cancertype<-as.factor(table_sample_type$Cancertype)
table_sample_type$Cellnumber<-as.numeric(as.matrix(table_sample_type$Cellnumber))
ce<-ddply(table_sample_type,"Celltype",transform,percent_weight=Cellnumber/sum(Cellnumber)*100)
ce<-ddply(ce,"Cancertype",transform,percent_weight2=Cellnumber/sum(Cellnumber))
ce<-ddply(ce,"Celltype",transform,percent_weight3=percent_weight2/sum(percent_weight2)*100)
ce$y_position<-as.factor(as.numeric(ce$Cancertype)*2)
p_bar2<-ggplot(ce,aes(x=Celltype,y=percent_weight3,fill=Cancertype))+geom_bar(stat = "identity")+scale_x_discrete(limits=c("5","4","3","2","1","0"))+
  theme_bw()+geom_text(aes(label=paste(round(percent_weight3),"%",sep = ""), y=percent_weight3), vjust=0)+coord_flip()
ggsave("Figure_2F.pdf", p_bar2, width=10 ,height=8)
for (i in levels(cow_all@active.ident)) {                 
  sign<-table(cow_merge_bar_all[cow_merge_bar_all[,1]==i,3]) #2 for celltype 3 for data type
  singletype_table<-cbind(rep(i,length(sign)),names(sign),as.matrix(sign))
  rownames(singletype_table)<-NULL
  if (i==0) {
    table_type<-singletype_table
  }else{
    table_type<-rbind(table_type,singletype_table)
  }   
}
colnames(table_type)<-c("Celltype","Cancertype","Cellnumber")
table_sample_type<-as.data.frame(table_type)
table_sample_type$Cancertype<-as.factor(table_sample_type$Cancertype)
table_sample_type$Cellnumber<-as.numeric(as.matrix(table_sample_type$Cellnumber))
p_bar<-ggplot(table_sample_type,aes(x=Celltype,y=Cellnumber,fill=Cancertype))+geom_bar(stat = "identity")+theme_bw()
ggsave("cowall_datatype_count.pdf", p_bar, width=10 ,height=12)
ce<-ddply(table_sample_type,"Celltype",transform,percent_weight=Cellnumber/sum(Cellnumber)*100)
p_bar2<-ggplot(ce,aes(x=Celltype,y=percent_weight,fill=Cancertype))+geom_bar(stat = "identity")+theme_bw()
ggsave("cowall_datatype_percent.pdf", p_bar2, width=10 ,height=8)

#=========================================================================================================
#=============================== Cell trajectory analysis through monocle3  ==============================
library(monocle3)
library(tidyverse)
library(patchwork)
library(spdep)
library(cowplot)
library(loomR)
cow_all@meta.data<-cow_all@meta.data[,!colnames(cow_all@meta.data)=="tech"]
file.remove("bovine_all_S4toh5ad.loom")
sdata.loom <- as.loom(x = cow_all, filename = "bovine_all_S4toh5ad.loom", verbose = FALSE)
sdata.loom$close_all()
adatas = sc.read_loom("bovine_all_S4toh5ad.loom", sparse=True, cleanup=False, dtype='float32')
adatas.write('./bovine_all_seurat2scanpy.h5ad')
cell_metadata <- cow_all@meta.data
data<-GetAssayData(cow_all,assay = "RNA")
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds) 
cds <- reduce_dimension(cds, preprocess_method="PCA")
cds <- cluster_cells(cds) 
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(cow_all, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed 
cds@int_colData$reducedDims$UMAP <- cds.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") 
cds <- learn_graph(cds,verbose = T,use_partition = F,close_loop = F)
# Figure 2I
pdf(file = "pseudotime_celltype.pdf",height = 8,width = 9)
plot_cells(cds,
           color_cells_by = "celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE, 
           label_branch_points=T,
           group_label_size=4,
           label_cell_groups = F,
           cell_size=1.5)
dev.off()
cds_backup = cds
cds = order_cells(cds)
pdf(file = "pseudotime_color.pdf",height = 8,width = 10)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)
dev.off()
pdf(file = "pseudotime_color2.pdf",height = 8,width = 10)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=F,
           label_branch_points=T,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)
dev.off()

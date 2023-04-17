# I read it in on the cluster using the anndata package
#w <- anndata::read_h5ad("ad_worm_aging.h5ad")
library(tidyverse)
library(Seurat)
library(SeuratDisk)
#library(SeuratData)


# to have a quick look without actually loading the dataset
#hfile <- Connect("D:/Harry_Jones/worms.h5Seurat")
# to actually load the Seurat object
worms <- LoadH5Seurat("D:/Harry_Jones/worms.h5Seurat")

# Using the convert function from SeuratDisk. This worked fine on the cluster
# Convert("ad_worm_aging.h5ad", dest = "h5seurat")
#worms <- LoadH5Seurat("D:/Harry_Jones/ad_worm_aging.h5seurat") # this does not seem to work
#worms_filt <- subset(worms, nFeature_RNA > 100)


#==========
# Some QC
#==========
# for human and mouse we can identify mitochondrial genes  (prefixed with "MT").
# C elegans MT genes are named differently 
# From WormBase parasite
# https://parasite.wormbase.org/biomart/martview/a635bf7189ad83f47d6778648d67fd2e
mt_genes <- c("ctb-1","nduo-6","ndfl-4","nduo-1","atp-6","nduo-2","ctc-3","nduo-4",
							"ctc-1","ctc-2","nduo-3","nduo-5","MTCE.1","MTCE.2","MTCE.5","MTCE.6",
							"MTCE.7","MTCE.8","MTCE.9","MTCE.10","MTCE.13","MTCE.14","MTCE.15",
							"MTCE.17","MTCE.18","MTCE.19","MTCE.20","MTCE.22","MTCE.24","MTCE.27",
							"MTCE.28","MTCE.29","MTCE.30","MTCE.32","MTCE.33","MTCE.36")

mt_genes <- mt_genes[mt_genes %in% rownames(worms@assays$RNA@counts)]
worms$percent.MT <- PercentageFeatureSet(worms, features = mt_genes)

# rrn seems to be ribosomal rRNA
grep("^rrn",rownames(worms@assays$RNA@counts),value = TRUE) # there may be more than these 3 but we'll go with this for now
worms$percent.Ribosomal <- PercentageFeatureSet(worms, pattern="^rrn")

# https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html#QC
worms$largest_count <- apply(
	worms@assays$RNA@counts,
	2,
	max
)

worms$largest_index <- apply(
	worms@assays$RNA@counts,
	2,
	which.max
)

worms$largest_gene <- rownames(worms)[worms$largest_index] 
worms$percent.Largest.Gene <- 100 * worms$largest_count/worms$nCount_RNA




VlnPlot(
	worms, 
	features = c("nCount_RNA", "nFeature_RNA", "percent.MT", "percent.Ribosomal", "percent.Largest.Gene")
)

VlnPlot(
	worms, 
	features = c("nCount_RNA", "nFeature_RNA", "percent.MT", "percent.Ribosomal", "percent.Largest.Gene")
) +
scale_y_log10()


qc.metrics <- as_tibble(
	worms[[]],
	rownames="Cell.Barcode"
) 

head(qc.metrics)

qc.metrics %>%
	arrange(percent.MT) %>%
	ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.MT)) + 
	geom_point() + 
	scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
	ggtitle("Example of plotting QC metrics") +
	geom_hline(yintercept = 750) +
	geom_hline(yintercept = 2000) +
	scale_x_log10() + scale_y_log10()

qc.metrics  <- qc.metrics %>%
	mutate(complexity=log10(nFeature_RNA) / log10(nCount_RNA))  

complexity.lm <- lm(log10(qc.metrics$nFeature_RNA)~log10(qc.metrics$nCount_RNA))

complexity.lm

qc.metrics <- qc.metrics %>%
	mutate(
		complexity_diff = log10(nFeature_RNA) - ((log10(qc.metrics$nCount_RNA)*complexity.lm$coefficients[2])+complexity.lm$coefficients[1])
	)

qc.metrics %>%
	ggplot(aes(x=complexity_diff)) +
	geom_density(fill="seagreen")


complexity_scale <- min(c(max(qc.metrics$complexity_diff),0-min(qc.metrics$complexity_diff))) 

qc.metrics %>%
	mutate(complexity_diff=replace(complexity_diff,complexity_diff< -0.1,-0.1)) %>%
	ggplot(aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=complexity_diff)) +
	geom_point(size=0.5) +
	geom_abline(slope=complexity.lm$coefficients[2], intercept = complexity.lm$coefficients[1]) +
	scale_colour_gradient2(low="blue2",mid="grey",high="red2")

qc.metrics %>%
	ggplot(aes(x=complexity_diff, y=percent.Largest.Gene)) +
	geom_point() 


largest_gene_list <- qc.metrics %>%
	group_by(largest_gene) %>%
	count() %>%
	arrange(desc(n)) 

largest_gene_list

largest_gene_list %>%
	filter(n>814) %>%
	pull(largest_gene) -> largest_genes_to_plot

qc.metrics %>%
	filter(largest_gene %in% largest_genes_to_plot) %>%
	mutate(largest_gene=factor(largest_gene, levels=largest_genes_to_plot)) %>%
	arrange(largest_gene) %>%
	ggplot(aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=largest_gene)) +
	geom_point(size=1) +
	scale_colour_manual(values=c("black",RColorBrewer::brewer.pal(9,"Set1")))

qc.metrics %>%
	filter(largest_gene %in% largest_genes_to_plot) %>%
	mutate(largest_gene=factor(largest_gene, levels=largest_genes_to_plot)) %>%
	arrange(largest_gene) %>%
	ggplot(aes(x=complexity_diff, y=percent.Largest.Gene, colour=largest_gene)) +
	geom_point()+
	scale_colour_manual(values=c("black",RColorBrewer::brewer.pal(9,"Set1")))

qc.metrics %>%
	arrange(percent.MT) %>%
	ggplot(aes(x=complexity_diff, y=percent.Largest.Gene, colour=percent.MT)) +
	geom_point() +
	scale_colour_gradient(low="grey", high="red2")

# doesn't look like MT genes are a big problem

qc.metrics %>%
	arrange(percent.Ribosomal) %>%
	ggplot(aes(x=complexity_diff, y=percent.Largest.Gene, colour=percent.Ribosomal)) +
	geom_point() +
	scale_colour_gradient(low="grey", high="red2")

qc.metrics %>%
	ggplot(aes(percent.MT)) + 
	geom_histogram(binwidth = 0.3, fill="seagreen", colour="black") +
	ggtitle("Distribution of Percentage Mitochondrion") +
	geom_vline(xintercept = 8)

qc.metrics %>%
	ggplot(aes(percent.Largest.Gene)) + 
	geom_histogram(binwidth = 0.7, fill="seagreen", colour="black") +
	ggtitle("Distribution of Percentage Largest Gene") +
	geom_vline(xintercept = 30)

qc.metrics %>%
	ggplot(aes(nFeature_RNA)) + 
	geom_histogram(fill="seagreen", binwidth = 50, colour="black") +
	ggtitle("Distribution of nFeature_RNA") +
	geom_vline(xintercept = 2000) +
	geom_vline(xintercept = 180)

qc.metrics %>%
	ggplot(aes(x = timepoint, y=nFeature_RNA)) + 
	geom_violin(fill="seagreen", colour="black") +
	ggtitle("Distribution of nFeature_RNA") 

qc.metrics %>%
	ggplot(aes(x = timepoint, y=nFeature_RNA)) + 
	geom_violin(fill="seagreen", colour="black") +
	ggtitle("Distribution of nFeature_RNA") +
	scale_y_log10() +
	geom_hline(yintercept = c(180, 2000))

# Filtering
worms <- subset(
	worms,
	nFeature_RNA>180 & 
		nFeature_RNA < 2000 & 
		percent.MT < 8 & 
		percent.Largest.Gene < 30
) 

#worms <- NormalizeData(worms, normalization.method = "LogNormalize")

gene.expression <- apply(worms@assays$RNA@data,1,mean)
gene.expression <- sort(gene.expression, decreasing = TRUE)
head(gene.expression, n=50)
# we probably want to remove those first couple - one is rRNA

# housekeeping genes
#tba-1, Y45F10D.4 and pmp-3
ggplot(mapping = aes(worms@assays$RNA@data["tba-1",])) + 
	geom_histogram(binwidth = 0.05, fill="seagreen", colour="black") + 
	ggtitle("")
# that's a lot of 0 values 
ggplot(mapping = aes(worms@assays$RNA@data["Y45F10D.4",])) + 
	geom_histogram(binwidth = 0.05, fill="seagreen", colour="black") + 
	ggtitle("")
ggplot(mapping = aes(worms@assays$RNA@data["pmp-3",])) + 
	geom_histogram(binwidth = 0.05, fill="seagreen", colour="black") + 
	ggtitle("")


as_tibble(
	worms@assays$RNA@data[,1:100]
) %>%
	pivot_longer(
		cols=everything(),
		names_to="cell",
		values_to="expression"
	) %>%
	ggplot(aes(x=expression, group=cell)) +
	geom_density() +
	coord_cartesian(ylim=c(0,0.6), xlim=c(0,3))

# way too many 0 values to use the simple normalisation
# Normalise again, this time using a centered log ratio transformation - 
# more similar to the sort of size factor based normalisation which is used for many RNA-Seq experiments. The margin=2 option means that it normalises per cell instead of per gene
worms <- NormalizeData(worms, normalization.method = "CLR", margin = 2)

save(worms, file = "D:/Harry_Jones/clr_norm_worms.rds")


as_tibble(
	worms@assays$RNA@data[,1:100]
) %>%
	pivot_longer(
		cols=everything(),
		names_to="cell",
		values_to="expression"
	) %>%
	ggplot(aes(x=expression, group=cell)) +
	geom_density() +
	coord_cartesian(ylim=c(0,0.6), xlim=c(0,3))
# That looks better

normalisation.qc <- tibble(
	pc95 = apply(worms[["RNA"]]@data,2,quantile,0.95),
	measured = apply(worms[["RNA"]]@data,2,function(x)(100*sum(x!=0))/length(x))
) 

normalisation.qc %>% 
	ggplot(aes(x=measured,y=pc95))+
	geom_point()+
	ggtitle("Normalisation of data")

# We should then do some cell cycle correction but I don't have a set of marker genes 
# so we could come back to this.
# http://www.wormbook.org/chapters/www_celldivision/celldivision.html

worms <- FindVariableFeatures(
	worms, 
	selection.method = "vst", 
	nfeatures=500
)

variance.data <- as_tibble(HVFInfo(worms),rownames = "Gene") 

variance.data <- variance.data %>% 
	mutate(hypervariable=Gene %in% VariableFeatures(worms)
	)

head(variance.data, n=10)

variance.data %>% 
	ggplot(aes(log(mean),log(variance),color=hypervariable)) + 
	geom_point() + 
	scale_color_manual(values=c("black","red"))


worms <- ScaleData(worms,features=rownames(data))
save(worms, file = "D:/Harry_Jones/clr_norm_scaled_worms.rds")


RunPCA(worms,features=VariableFeatures(worms)) -> worms
DimPlot(worms,reduction="pca")

DimPlot(worms,reduction="pca", group.by = "largest_gene", label = TRUE, label.size = 3) + NoLegend()
DimPlot(worms,reduction="pca", dims = c(3,4))
# these look odd

ElbowPlot(worms)

DimHeatmap(worms,dims=1:20, cells=500)
#DimHeatmap(worms,dims=21:30, cells=500)

# worms <- RunTSNE(
# 	worms,
# 	dims=1:20,
# 	perplexity=10,
# 	check_duplicates = FALSE
# ) # this took a long time

#DimPlot(worms, reduction = "tsne", pt.size = 1) + ggtitle("tSNE with Perplexity 10")

# worms <- RunTSNE(
# 	worms,
# 	dims=1:9,
# 	perplexity=200,
# 	check_duplicates = FALSE
# ) # as did this

#DimPlot(worms, reduction = "tsne", pt.size = 1) + ggtitle("tSNE with Perplexity 200")

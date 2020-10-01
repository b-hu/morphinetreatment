setwd("/local/Bo/Remap_allhearts/clustering_allhearts/")

load("/local/Bo/Remap_allhearts/clustering_allhearts/all.hearts.merged.150co.noRFP.40pc.Robj")
all.hearts <- SetIdent(all.hearts,value = "orig.ident")
morphcomp <- SubsetData(all.hearts,ident.use = c("Hr8","Hr11","Hr15","Hr16","Hr18","Hr22"))
rm(all.hearts)
morphcomp <- SCTransform(morphcomp,variable.features.n = 3500,vars.to.regress = c("percent.mito") )
morphcomp <- NormalizeData(object = morphcomp)
morphcomp <- FindVariableFeatures(object = morphcomp, 
                              selection.method = "vst", nfeatures = 2000)
morphcomp <- ScaleData(object = morphcomp,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mito"))


morphcomp <- RunPCA(object = morphcomp, features = VariableFeatures(object = morphcomp),npcs = 100)

ElbowPlot(object = morphcomp,ndims = 100)

morphcomp <- RunUMAP(morphcomp, dims = 1:7)
morphcomp <- FindNeighbors(morphcomp, dims = 1:7)
morphcomp <- FindClusters(morphcomp, verbose = FALSE)

morphcomp <- SetIdent(morphcomp,value = "basic.Cell.type")
morphcomp <- SubsetData(morphcomp,ident.remove = "Dead cells")
morphcomp <- SCTransform(morphcomp,variable.features.n = 3500,vars.to.regress = c("percent.mito") )
morphcomp <- NormalizeData(object = morphcomp)
morphcomp <- FindVariableFeatures(object = morphcomp, 
                                  selection.method = "vst", nfeatures = 1000)
morphcomp <- ScaleData(object = morphcomp,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mito"))

morphcomp <- RunPCA(object = morphcomp, features = VariableFeatures(object = morphcomp),npcs = 100)
ElbowPlot(object = morphcomp,ndims = 100)
morphcomp <- RunUMAP(morphcomp, dims = 1:10)
morphcomp <- FindNeighbors(morphcomp, dims = 1:10)
#morphcomp <- FindClusters(morphcomp, verbose = FALSE)
#morphcomp <- SetIdent(morphcomp,value = "basic.Cell.type")
#morphcomp <- SetIdent(morphcomp,cells = WhichCells(morphcomp,expression = basic.Cell.type == "Dead cells"),value = "Endocardium")

#p1 <- FeaturePlot(morphcomp,features = c("myl7","hbba1","acta2","cd74a"),cols = c("lightgrey","red"),reduction = "umap") + NoLegend()
#p2 <- DimPlot(object = morphcomp, reduction = "umap",label = TRUE) + NoLegend()
#p3 <- DimPlot(object = morphcomp, reduction = "umap",label = F, group.by = "basic.Cell.type")
p4 <- DimPlot(object = morphcomp, reduction = "umap",label = F, group.by = "morphine")
plot_grid(DimPlot(morphcomp,group.by = "basic.Cell.type",label = F),p4)

plot_grid(p3,p4,p1,p2)

save(morphcomp,file = "../../morph.comp/morphcomp2.Robj")
pdf("../../morph.comp/morphcomp.umaps2.pdf",width = 13,height = 4)
plot_grid(DimPlot(morphcomp,group.by = "basic.Cell.type",label = F),p4)
dev.off()

#FeaturePlot(morphcomp,c("epdl1", "id2a", "c1qa", "grn1", "grn2", "lgmn"),cols = c("grey","red"))

#function summarySE####
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#calc ####

morph.counts <- data.frame(table(morphcomp@meta.data$orig.ident, morphcomp@meta.data$basic.Cell.type))
colnames(morph.counts) <- c("Library","Cell.type","Freq")
morph.counts <- merge(x=morph.counts, y= unique(data.frame(Library=morphcomp@meta.data$orig.ident,
                                                            morph=morphcomp@meta.data$morphine,
                                                           time=morphcomp@meta.data$time)),
                       by.x = "Library", by.y = "Library",all.y = F)

for (i in unique(morph.counts$Library)) {
  morph.counts$ratio[morph.counts$Library == i] <-  morph.counts$Freq[morph.counts$Library == i] / sum(morph.counts$Freq[morph.counts$Library == i])
}

ord <- names(sort(summary(as.factor(morphcomp@meta.data$basic.Cell.type)),decreasing = T))

morph.plot <- list()
for (j in c("3dpi","7dpi")) {
  temp <- morph.counts[morph.counts$time == j,]
  for (i in ord) {
    morph.plot[[paste(j,i)]] <- 
      ggplot(temp[temp$Cell.type==i,], aes(x=morph,y=ratio,fill="grey")) +
      geom_bar(position=position_dodge(), stat="identity",show.legend = FALSE) +
      ylab("% of all cells") +
      xlab("")+
      ggtitle(paste(j,i))
  }
}
pdf(file = "../morph.comp/testsplit.pdf",width = 18,height = 22)
plot_grid(plotlist = morph.plot)
dev.off()

morph.plot <- list()
  for (i in ord) {
    morph.plot[[i]] <- ggplot(temp[temp$Cell.type==i,], aes(x=morph,y=ratio,fill="grey")) +
      geom_bar(position=position_dodge(), stat="identity",show.legend = FALSE) +
      ylab("% of all cells") +
      xlab("")+
      ggtitle(paste(j,i))
  }
pdf(file = "../morph.comp/split2.pdf",width = 18,height = 22)
plot_grid(plotlist = morph.plot)
dev.off()






morph.calc <- summarySE(data = morph.counts,measurevar = "ratio" ,groupvars = c("morph","Cell.type"))
morph.calc$ratio <- morph.calc$ratio *100
morph.calc$se <- morph.calc$se * 100
morph.plot <- list()
for (i in ord) {
  morph.plot[[i]] <- 
    ggplot(morph.calc[morph.calc$Cell.type == i,], aes(x=morph,y=ratio,fill="grey")) +
    geom_bar(position=position_dodge(), stat="identity",show.legend = FALSE) +
    geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
    #coord_cartesian(ylim = c(0, max(deep.nppc[deep.nppc$Cell.type == i,]$norm.ratio)+0.05         ),) +  
    ylab("% of all cells") +
    ggtitle(i)
}


p.values <- list()
for (i in ord) {
test <- t.test(
morph.counts[morph.counts$Cell.type == i,][morph.counts[morph.counts$Cell.type == i,]$morph == "Morph", ]$ratio,
morph.counts[morph.counts$Cell.type == i,][morph.counts[morph.counts$Cell.type == i,]$morph == "Ctrl", ]$ratio)
p.values[[i]] <- test$p.value
}
p.values <- as.data.frame(unlist(p.values))
colnames(p.values) <- "p.values"
p.values$"Significant?" <- "No" 
library(gridExtra)


plot_grid( plot_grid(plotlist = morph.plot),tableGrob(p.values))

DotPlot(all.hearts,features = c("cx40.8"),group.by = "fibro.CM.subtypes")

- aldh1a2
- postnb
- gata4
- cxcl12a
- mpeg1.1
- tgfb1b
- pdgfb
- cxcr4b
- col1a1
- fh1a
load("morph.comp/morphcomp.Robj")

mor.genes <- c("aldh1a2","postnb","gata4","cxcl12a",
               "mpeg1.1","tgfb1b","pdgfba","pdgfbb","cxcr4b","col1a1a","fn1a","actb1","actb2","pcna"
)
mor.comp.genes <- DotPlot(morphcomp,features =mor.genes,group.by = "morphine",split.by = "orig.ident",cols = rep("red",11) )
mor.comp.genes <- mor.comp.genes$data
mor.comp.genes$treatment <- 
  sapply(mor.comp.genes$id,  function(x) unlist(strsplit(as.character(x),"_"))[1]  )
mor.comp.genes$library <- 
  sapply(mor.comp.genes$id,  function(x) unlist(strsplit(as.character(x),"_"))[2]  )
mor.comp.genes <- summarySE(data = mor.comp.genes,measurevar = "avg.exp" ,groupvars = c("treatment","features.plot"))
mor.comp.genes[,c(1,3,4,7)]
write.csv(mor.comp.genes[,c(1,3,7,8)],file = "/local/Bo/Remap_allhearts/morph.comp/morphcomp.raw.csv")

mor.plot <- list()
for (i in mor.genes) {
  mor.plot[[i]] <- 
ggplot(mor.comp.genes[mor.comp.genes$features.plot == i,], aes(x=treatment,y=avg.exp,fill = treatment)) +
  geom_bar(stat = "identity",position = "dodge") +
  scale_fill_manual(values = c("red","blue") ) +
  geom_errorbar(aes(ymin=avg.exp-se, ymax=avg.exp+se),width=.2,position=position_dodge(.9)) +
  #geom_text(aes(label=time),position=position_dodge(width=0.9), vjust=0.5,hjust = -0.25, angle = 90) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,face = "bold")) +
  #scale_y_continuous(limits = c(0,0.3)) +
  ylab(paste0(i," expression")) +
  ggtitle(i)
  }
combine_plots(plotlist = mor.plot)

mor.plot <- list()
for (i in mor.genes) {
  mor.plot[[i]] <- 
    ggplot(mor.comp.genes[mor.comp.genes$features.plot == i,], aes(x=treatment,y=avg.exp,fill = treatment)) +
    geom_bar(stat = "identity",position = "dodge") +
    scale_fill_manual(values = c("red","blue") ) +
    geom_errorbar(aes(ymin=avg.exp-se, ymax=avg.exp+se),width=.2,position=position_dodge(.9)) +
    #geom_text(aes(label=time),position=position_dodge(width=0.9), vjust=0.5,hjust = -0.25, angle = 90) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,face = "bold")) +
    #scale_y_continuous(limits = c(0,0.3)) +
    ylab(paste0(i," expression")) +
    ggtitle(i)
  
}
combine_plots(plotlist = mor.plot)

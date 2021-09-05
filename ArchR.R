#ArchR 
#ssh fat02
#cd /share2/pub/zhenggw/zhenggw/HemeFragments/
conda activate R410
R
library(ArchR)
library(parallel)
# mclapply

#make arrow files
addArchRGenome("hg19") # hg38, mm9, mm10
addArchRThreads(threads = 16)

pathFragments <- "/share2/pub/zhenggw/zhenggw/HemeFragments/"

inputFiles <- list.files(pathFragments, pattern = ".gz",
        full.names = TRUE)
    names(inputFiles) <- gsub(".fragments.tsv.gz", "", list.files(pathFragments,
        pattern = ".gz"))
    inputFiles <- inputFiles[!grepl(".tbi", inputFiles)]
    inputFiles

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later 信噪比，根据TSS富集分数进行计算 
  filterFrags = 1000, #ArchR默认会过滤TSS富集得分低于4或唯一比对数小于1000（也就是保留TSS富集得分大于4且唯一比对数大于1000的细胞）#细胞质控
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# > ArrowFiles
# [1] "scATAC_BMMC_R1.arrow"      "scATAC_CD34_BMMC_R1.arrow"
# [3] "scATAC_PBMC_R1.arrow"

#添加doublets信息，后续去除
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

#创建ArchRProject 允许多个Arrow文件整理到单个项目之中
projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)



#一些操作
# paste0("Memory Size = ", round(object.size(projHeme1) / 10^6, 3), " MB")
# # "Memory Size = 37.477 MB"

# 我们还可以检查当前的ArchRProject中存放了哪些矩阵数据，这些数据一旦增加后就可以在下游分析中使用。
# getAvailableMatrices(projHeme1)
# # "GeneScoreMatrix" "TileMatrix"
# head(projHeme1$cellNames)
# head(projHeme1$Sample)
# quantile(projHeme1$TSSEnrichment)
#projHeme1$[TAB]

# 从ArchRProject中提取部分细胞
# 以前学习的R语言向量/矩阵/数据框提取数据的方法可以无缝的迁移到ArchRProject上，但区别在于ArchR并不会直接提取数据，而是构建一个新的ArchRProject对象。

# 最简单的方法是根据位置和细胞名
# # 根据位置
# projHeme1[1:100, ]
# # 根据细胞名
# projHeme1[projHeme1$cellNames[1:100], ]

# 复杂一些就是先根据逻辑判断提取细胞名，接着根据细胞名提取数据。例如挑选scATAC_BMMC_R1的细胞，或是TSSEnrichment大于8的细胞。

# # sample name
# idxSample <- BiocGenerics::which(projHeme1$Sample %in% "scATAC_BMMC_R1")
# cellsSample <- projHeme1$cellNames[idxSample]
# projHeme1[cellsSample, ]
# # TSS enrichment score
# idxPass <- which(projHeme1$TSSEnrichment >= 8)
# cellsPass <- projHeme1$cellNames[idxPass]
# projHeme1[cellsPass, ]

# 用getCellColData提取我们需要的两列，其中nFrages需要进行log10运算
df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))

p <- ggPoint(
    x = df[,1],
    y = df[,2],
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

p

# png("TSS-vs-Frags.png") 
# plot(p)
# dev.off()
# getwd()

library(ggplot2)
ggsave("TSS-vs-Frags.pdf")


# ArchR提供了小提琴图(violin plot)和山脊图(ridge plot)用来展示不同组之间的信息。这两种类型的图可以用一个函数plotGroups进行绘制。除了根据样本进行分组绘图外，还可以使用下游分析得到的分组信息（例如聚类）。
# 根据TSS富集得分为每个样本绘制山脊图。设置plotAs = "ridges"
p1 <- plotGroups(
    ArchRProj = projHeme1,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "ridges"
   )
p1

#绘制小提琴图 plotAs = "violin"
p2 <- plotGroups(
    ArchRProj = projHeme1,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p2

p3 <- plotGroups(
    ArchRProj = projHeme1,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "ridges"
   )
p3

p4 <- plotGroups(
    ArchRProj = projHeme1,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p4

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 4, height = 4)

# 绘制样本的TSS富集谱和Fragment大小分布
#Fragments大小分布
p1 <- plotFragmentSizes(ArchRProj = projHeme1)
p1

#TSS富集谱
p2 <- plotTSSEnrichment(ArchRProj = projHeme1)
p2

plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)

# #这两个图有问题，报错如下，不知道是不是hdf5的原因
# Error in .safelapply(seq_along(uniqGroups), function(x) { :
# Error Found Iteration 1 :
#         [1] "Error in H5Fopen(file) : HDF5. File accessibility. Unable to open file.\n"
#         <simpleError in H5Fopen(file): HDF5. File accessibility. Unable to open file.>
# Error Found Iteration 2 :
#         [1] "Error in H5Fopen(file) : HDF5. File accessibility. Unable to open file.\n"
#         <simpleError in H5Fopen(file): HDF5. File accessibility. Unable to open file.>


# --------------------------------------------------------
# 目前的问题
# Cairo
# parallel的mclapply容易出问题                                       !!!!!解决办法：library(parallel)
# HDF5. File accessibility. Unable to open file.                    !!!!!解决办法：chmod 777 *.arrow




#-------------------------------------------------------------------------------------------------------------

# 保存和加载ArchRProject
saveArchRProject(ArchRProj = projHeme1, outputDirectory = "Save-ProjHeme1", load = FALSE)
# #load等于False不会改变当前环境的projHeme 如果想要将当前的项目复制到新的目录，可以设置load=TRUE

projHeme1 <- readRDS("./Save-ProjHeme1/Save-ArchR-Project.rds")

#------------------------------------------------------------------------------------------------
#从ArchRProject过滤doublets
projHeme2 <- filterDoublets(projHeme1)  #projHemeTmp <- filterDoublets(projHeme1, filterRatio = 1.5) 提高filterRatio会过滤更多的细胞


#------------------------------------------------------------------------------------------------
#ArchR降维分析 隐语义(Latent Semantic Indexing)迭代 （增加LSI的迭代次数，也可以用来处理批次效应）
projHeme2 <- addIterativeLSI(
    ArchRProj = projHeme2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

#使用Harmony矫正批次效应
projHeme2 <- addHarmony(
    ArchRProj = projHeme2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

#使用ArchR聚类
projHeme2 <- addClusters(
    input = projHeme2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

head(projHeme2$Clusters)
table(projHeme2$Clusters)
#   C1  C10  C11  C12   C2   C3   C4   C5   C6   C7   C8   C9
# 1575  720 1221 1384 1079  307  388  432  353 1276  816  699

library(pheatmap)

cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)

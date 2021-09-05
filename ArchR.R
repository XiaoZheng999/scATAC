#ArchR 
#ssh fat02
#cd /share2/pub/zhenggw/zhenggw/HemeFragments/
conda activate R410
R
library(ArchR)
addArchRThreads(threads = 16)
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

#创建ArchRProject ArchRProject对象允许我们将多个Arrow文件整理到单个项目之中
projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
library(lattice)
#source() #to source function in a seperate file
prepareData <- function(luckyGenes, prDes){
  geneList <- prDat[luckyGenes,]
  #geneVec <- as.vector (geneList) #change each element of list to a vector
  geneMatrix <- as.matrix(geneList)
  geneVec <- as.vector (t(geneMatrix))
  
  geneName <- rep(luckyGenes , each = 39)
  
  jDat<- data.frame(prDes, gExp = geneVec, geneName)
  return(jDat)
}

makeStripplot <- function(jDat, pch, cex){
  stripplot(gExp ~ devStage | geneName, jDat, group = gType, jitter.data = TRUE, auto.key = TRUE, type = c('p', 'a'), grid = TRUE, pch = pch, cex = cex)
}


prDat <- read.table("/Users/fdorri/Documents/UBC/courses/STAT540/workspace/stat540_2014/examples/photoRec/data/GSE4051_data.tsv")
str(prDat, max.level = 0)
prDes <- readRDS("/Users/fdorri/Documents/UBC/courses/STAT540/workspace/stat540_2014/examples/photoRec/data/GSE4051_design.rds")
str(prDes)

  
(luckyGenes <- c("1419655_at","1438815_at"))
jDat <- prepareData(luckyGenes, prDes)
str(jDat)

head(jDat)
tail(jDat)
makeStripplot(jDat, pch = 20, cex = 1.3)

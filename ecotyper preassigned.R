#Run ecotyper analysis based off what ecotyper could discover
#READ in cutpoints and format (OS cutpoints already filtered to significance)
OScutpoints <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/optimalcutpointsOS.csv", header = TRUE, row.names = 1)
OScutpoints$ID <- NULL
CSScutpoints <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/optimalcutpointsCSS.csv", header = TRUE, row.names = 1)
CSScutpoints$OS.months <- NULL
CSScutpoints$Bioinformatic.CSS_status <- NULL
DFScutpoints <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/optimalcutpointsDFS.csv", header = TRUE, row.names = 1)
DFScutpoints$Bioinformatic.DFS_months <- NULL
DFScutpoints$Bioinformatic.DFS_status <- NULL

#read in ecotyper assignments
Ecotyperassign <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/ecotyper_output ALL/Carcinoma_Ecotypes/Ecotype_Assignment.csv")
Ecotyperassign$ID <-strtrim(Ecotyperassign$ID, 12)
rownames(Ecotyperassign) <- Ecotyperassign$ID

EcotyperassignOS <- Ecotyperassign[rownames(Ecotyperassign) %in% rownames(OScutpoints),]
EcotyperassignCSS <- Ecotyperassign[rownames(Ecotyperassign) %in% rownames(CSScutpoints),]
EcotyperassignDFS <- Ecotyperassign[rownames(Ecotyperassign) %in% rownames(DFScutpoints),]

#merge ecotypes with OS data
#Merge the ecotyper TIME type data with the cutpoints
EcotyperassignOS$ID <- rownames(EcotyperassignOS)
OScutpoints$ID <- row.names(OScutpoints)
mergeECO_OS <- merge.data.frame(OScutpoints, EcotyperassignOS, by = "ID")
EcotyperassignCSS$ID <- rownames(EcotyperassignCSS)
CSScutpoints$ID <- row.names(CSScutpoints)
mergeECO_CSS <- merge.data.frame(CSScutpoints, EcotyperassignCSS, by = "ID")
EcotyperassignDFS$ID <- rownames(EcotyperassignDFS)
DFScutpoints$ID <- row.names(DFScutpoints)
mergeECO_DFS <- merge.data.frame(DFScutpoints, EcotyperassignDFS, by = "ID")

#Remove ID column
rownames(mergeECO_OS) <- mergeECO_OS$ID
rownames(mergeECO_CSS) <- mergeECO_CSS$ID
rownames(mergeECO_DFS)<- mergeECO_DFS$ID
mergeECO_OS$ID <- NULL
mergeECO_CSS$ID <- NULL
mergeECO_DFS$ID <- NULL


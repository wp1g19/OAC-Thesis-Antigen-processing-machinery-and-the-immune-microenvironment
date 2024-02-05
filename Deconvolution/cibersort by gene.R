#read in cibsort data
cibersort <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/CIBERSORTx_Final_Simple.csv", header = TRUE, row.names = 1)
rownames(cibersort) <- gsub("/",".", rownames(cibersort))
rownames(cibersort) <- gsub("-",".", rownames(cibersort))
#READ in cutpoints and format (OS cutpoints already filtered to significance)
OScutpoints <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/optimalcutpointsOS.csv", header = TRUE, row.names = 1)
OScutpoints$ID <- NULL
CSScutpoints <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/optimalcutpointsCSS.csv", header = TRUE, row.names = 1)
CSScutpoints$OS.months <- NULL
CSScutpoints$Bioinformatic.CSS_status <- NULL
DFScutpoints <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/optimalcutpointsDFS.csv", header = TRUE, row.names = 1)
DFScutpoints$Bioinformatic.DFS_months <- NULL
DFScutpoints$Bioinformatic.DFS_status <- NULL

# merge cibersort and OS/CSS/DFS data
rownames(cibersort) <- strtrim(rownames(cibersort), 12)

#filter ecotyper data to OS,CSS,DFS data 
cibersortOS <- cibersort[rownames(cibersort) %in% rownames(OScutpoints),]
cibersortCSS <- cibersort[rownames(cibersort) %in% rownames(CSScutpoints),]
cibersortDFS <- cibersort[rownames(cibersort) %in% rownames(DFScutpoints),]

cibersortOS$ID <- rownames(cibersortOS)
OScutpoints$ID <- row.names(OScutpoints)
mergeciber_OS <- merge.data.frame(OScutpoints, cibersortOS, by = "ID")
cibersortCSS$ID <- rownames(cibersortCSS)
CSScutpoints$ID <- row.names(CSScutpoints)
mergeciber_CSS <- merge.data.frame(CSScutpoints, cibersortCSS, by = "ID")
cibersortDFS$ID <- rownames(cibersortDFS)
DFScutpoints$ID <- row.names(cibersortDFS)
mergeciber_DFS <- merge.data.frame(DFScutpoints, cibersortDFS, by = "ID")

#Remove ID column
rownames(mergeciber_OS) <- mergeciber_OS$ID
rownames(mergeciber_CSS) <- mergeciber_CSS$ID
rownames(mergeciber_DFS)<- mergeciber_DFS$ID
mergeciber_OS$ID <- NULL
mergeciber_CSS$ID <- NULL
mergeciber_DFS$ID <- NULL

#make sure numeric columns are numeric
i <- c(28:37)   
mergeciber_OS[ , i] <- apply(mergeciber_OS[ , i], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))
i <- c(45:54) 
mergeciber_CSS[ , i] <- apply(mergeciber_CSS[ , i], 2,            # Specify own function within apply
                             function(x) as.numeric(as.character(x)))
mergeciber_DFS[ , i] <- apply(mergeciber_DFS[ , i], 2,            # Specify own function within apply
                              function(x) as.numeric(as.character(x)))


library(zoo)
#split dataframe by gene survival association
#setup colnames for downstream analysis in graphpad
colnamesciber <- c("High:B cells",
                   "High:Plasma cells",
                   "High:T cells CD8",
                   "High:T cells CD4",
                   "High:T cells regulatory (Tregs)",
                   "High:NK cells",
                   "High:Macrophages",
                   "High:Dendritic cells",
                   "High:Granulocytes",
                   "High:Monocytes",
                   "Low:B cells",
                   "Low:Plasma cells",
                   "Low:T cells CD8",
                   "Low:T cells CD4",
                   "Low:T cells regulatory (Tregs)",
                   "Low:NK cells",
                   "Low:Macrophages",
                   "Low:Dendritic cells",
                   "Low:Granulocytes",
                   "Low:Monocytes")
colnamesciberopposite <- c("Low:B cells",
                           "Low:Plasma cells",
                           "Low:T cells CD8",
                           "Low:T cells CD4",
                           "Low:T cells regulatory (Tregs)",
                           "Low:NK cells",
                           "Low:Macrophages",
                           "Low:Dendritic cells",
                           "Low:Granulocytes",
                           "Low:Monocytes",
                           "High:B cells",
                   "High:Plasma cells",
                   "High:T cells CD8",
                   "High:T cells CD4",
                   "High:T cells regulatory (Tregs)",
                   "High:NK cells",
                   "High:Macrophages",
                   "High:Dendritic cells",
                   "High:Granulocytes",
                   "High:Monocytes"
                   )


HLA_A_OS <-split(mergeciber_OS, mergeciber_OS$HLA_A)
HLA_A_High_OS <- HLA_A_OS[["high"]]
HLA_A_Low_OS <- HLA_A_OS[["low"]]
HLA_A_High_OS <- HLA_A_High_OS[,-1:-27]
HLA_A_Low_OS <- HLA_A_Low_OS[,-1:-27]
HLA_A_OS_split <- data.frame(HLA_A_High_OS, cbind(zoo(, 1:nrow(HLA_A_High_OS)), as.zoo(HLA_A_Low_OS)))
colnames(HLA_A_OS_split) <- colnamesciber

HLA_B_OS <-split(mergeciber_OS, mergeciber_OS$HLA_B)
HLA_B_High_OS <- HLA_B_OS[["high"]]
HLA_B_Low_OS <- HLA_B_OS[["low"]]
HLA_B_High_OS <- HLA_B_High_OS[,-1:-27]
HLA_B_Low_OS <- HLA_B_Low_OS[,-1:-27]
HLA_B_OS_split <- data.frame(HLA_B_High_OS, cbind(zoo(, 1:nrow(HLA_B_High_OS)), as.zoo(HLA_B_Low_OS)))
colnames(HLA_B_OS_split) <- colnamesciber

HLA_E_OS <-split(mergeciber_OS, mergeciber_OS$HLA_E)
HLA_E_High_OS <- HLA_E_OS[["high"]]
HLA_E_Low_OS <- HLA_E_OS[["low"]]
HLA_E_High_OS <- HLA_E_High_OS[,-1:-27]
HLA_E_Low_OS <- HLA_E_Low_OS[,-1:-27]
HLA_E_OS_split <- data.frame(HLA_E_Low_OS, cbind(zoo(, 1:nrow(HLA_E_Low_OS)), as.zoo(HLA_E_High_OS)))
colnames(HLA_E_OS_split) <- colnamesciberopposite


HLA_G_OS <-split(mergeciber_OS, mergeciber_OS$HLA_G)
HLA_G_High_OS <- HLA_G_OS[["high"]]
HLA_G_Low_OS <- HLA_G_OS[["low"]]
HLA_G_High_OS <- HLA_G_High_OS[,-1:-27]
HLA_G_Low_OS <- HLA_G_Low_OS[,-1:-27]
HLA_G_OS_split <- data.frame(HLA_G_Low_OS, cbind(zoo(, 1:nrow(HLA_G_Low_OS)), as.zoo(HLA_G_High_OS)))
colnames(HLA_G_OS_split) <- colnamesciberopposite


HLA_DPA1_OS <-split(mergeciber_OS, mergeciber_OS$HLA_DPA1)
HLA_DPA1_High_OS <- HLA_DPA1_OS[["high"]]
HLA_DPA1_Low_OS <- HLA_DPA1_OS[["low"]]
HLA_DPA1_High_OS <- HLA_DPA1_High_OS[,-1:-27]
HLA_DPA1_Low_OS <- HLA_DPA1_Low_OS[,-1:-27]
HLA_DPA1_OS_split <- data.frame(HLA_DPA1_High_OS, cbind(zoo(, 1:nrow(HLA_DPA1_High_OS)), as.zoo(HLA_DPA1_Low_OS)))
colnames(HLA_DPA1_OS_split) <- colnamesciber


HLA_DQA1_OS <-split(mergeciber_OS, mergeciber_OS$HLA_DQA1)
HLA_DQA1_High_OS <- HLA_DQA1_OS[["high"]]
HLA_DQA1_Low_OS <- HLA_DQA1_OS[["low"]]
HLA_DQA1_High_OS <- HLA_DQA1_High_OS[,-1:-27]
HLA_DQA1_Low_OS <- HLA_DQA1_Low_OS[,-1:-27]
HLA_DQA1_OS_split <- data.frame(HLA_DQA1_High_OS, cbind(zoo(, 1:nrow(HLA_DQA1_High_OS)), as.zoo(HLA_DQA1_Low_OS)))
colnames(HLA_DQA1_OS_split) <- colnamesciber


HLA_DRA_OS <-split(mergeciber_OS, mergeciber_OS$HLA_DRA)
HLA_DRA_High_OS <- HLA_DRA_OS[["high"]]
HLA_DRA_Low_OS <- HLA_DRA_OS[["low"]]
HLA_DRA_High_OS <- HLA_DRA_High_OS[,-1:-27]
HLA_DRA_Low_OS <- HLA_DRA_Low_OS[,-1:-27]
HLA_DRA_OS_split <- data.frame(HLA_DRA_Low_OS, cbind(zoo(, 1:nrow(HLA_DRA_Low_OS)), as.zoo(HLA_DRA_High_OS)))
colnames(HLA_DRA_OS_split) <- colnamesciberopposite

HLA_DRB1_OS <-split(mergeciber_OS, mergeciber_OS$HLA_DRB1)
HLA_DRB1_High_OS <- HLA_DRB1_OS[["high"]]
HLA_DRB1_Low_OS <- HLA_DRB1_OS[["low"]]
HLA_DRB1_High_OS <- HLA_DRB1_High_OS[,-1:-27]
HLA_DRB1_Low_OS <- HLA_DRB1_Low_OS[,-1:-27]
HLA_DRB1_OS_split <- data.frame(HLA_DRB1_High_OS, cbind(zoo(, 1:nrow(HLA_DRB1_High_OS)), as.zoo(HLA_DRB1_Low_OS)))
colnames(HLA_DRB1_OS_split) <- colnamesciber

HLA_DRB5_OS <-split(mergeciber_OS, mergeciber_OS$HLA_DRB5)
HLA_DRB5_High_OS <- HLA_DRB5_OS[["high"]]
HLA_DRB5_Low_OS <- HLA_DRB5_OS[["low"]]
HLA_DRB5_High_OS <- HLA_DRB5_High_OS[,-1:-27]
HLA_DRB5_Low_OS <- HLA_DRB5_Low_OS[,-1:-27]
HLA_DRB5_OS_split <- data.frame(HLA_DRB5_High_OS, cbind(zoo(, 1:nrow(HLA_DRB5_High_OS)), as.zoo(HLA_DRB5_Low_OS)))
colnames(HLA_DRB5_OS_split) <- colnamesciber

LGMN_OS <-split(mergeciber_OS, mergeciber_OS$LGMN)
LGMN_High_OS <- LGMN_OS[["high"]]
LGMN_Low_OS <- LGMN_OS[["low"]]
LGMN_High_OS <- LGMN_High_OS[,-1:-27]
LGMN_Low_OS <- LGMN_Low_OS[,-1:-27]
LGMN_OS_split <- data.frame(LGMN_High_OS, cbind(zoo(, 1:nrow(LGMN_High_OS)), as.zoo(LGMN_Low_OS)))
colnames(LGMN_OS_split) <- colnamesciber

MR1_OS <-split(mergeciber_OS, mergeciber_OS$MR1)
MR1_High_OS <- MR1_OS[["high"]]
MR1_Low_OS <- MR1_OS[["low"]]
MR1_High_OS <- MR1_High_OS[,-1:-27]
MR1_Low_OS <- MR1_Low_OS[,-1:-27]
MR1_OS_split <- data.frame(MR1_Low_OS, cbind(zoo(, 1:nrow(MR1_Low_OS)), as.zoo(MR1_High_OS)))
colnames(MR1_OS_split) <- colnamesciberopposite

PSMB10_OS <-split(mergeciber_OS, mergeciber_OS$PSMB10)
PSMB10_High_OS <- PSMB10_OS[["high"]]
PSMB10_Low_OS <- PSMB10_OS[["low"]]
PSMB10_High_OS <- PSMB10_High_OS[,-1:-27]
PSMB10_Low_OS <- PSMB10_Low_OS[,-1:-27]
PSMB10_OS_split <- data.frame(PSMB10_Low_OS, cbind(zoo(, 1:nrow(PSMB10_Low_OS)), as.zoo(PSMB10_High_OS)))
colnames(PSMB10_OS_split) <- colnamesciberopposite

PSMB8_OS <-split(mergeciber_OS, mergeciber_OS$PSMB8)
PSMB8_High_OS <- PSMB8_OS[["high"]]
PSMB8_Low_OS <- PSMB8_OS[["low"]]
PSMB8_High_OS <- PSMB8_High_OS[,-1:-27]
PSMB8_Low_OS <- PSMB8_Low_OS[,-1:-27]
PSMB8_OS_split <- data.frame(PSMB8_Low_OS, cbind(zoo(, 1:nrow(PSMB8_Low_OS)), as.zoo(PSMB8_High_OS)))
colnames(PSMB8_OS_split) <- colnamesciber

PSMB9_OS <-split(mergeciber_OS, mergeciber_OS$PSMB9)
PSMB9_High_OS <- PSMB9_OS[["high"]]
PSMB9_Low_OS <- PSMB9_OS[["low"]]
PSMB9_High_OS <- PSMB9_High_OS[,-1:-27]
PSMB9_Low_OS <- PSMB9_Low_OS[,-1:-27]
PSMB9_OS_split <- data.frame(PSMB9_Low_OS, cbind(zoo(, 1:nrow(PSMB9_Low_OS)), as.zoo(PSMB9_High_OS)))
colnames(PSMB9_OS_split) <- colnamesciberopposite

RFX5_OS <-split(mergeciber_OS, mergeciber_OS$RFX5)
RFX5_High_OS <- RFX5_OS[["high"]]
RFX5_Low_OS <- RFX5_OS[["low"]]
RFX5_High_OS <- RFX5_High_OS[,-1:-27]
RFX5_Low_OS <- RFX5_Low_OS[,-1:-27]
RFX5_OS_split <- data.frame(RFX5_Low_OS, cbind(zoo(, 1:nrow(RFX5_Low_OS)), as.zoo(RFX5_High_OS)))
colnames(RFX5_OS_split) <- colnamesciberopposite

SPPL2A_OS <-split(mergeciber_OS, mergeciber_OS$SPPL2A)
SPPL2A_High_OS <- SPPL2A_OS[["high"]]
SPPL2A_Low_OS <- SPPL2A_OS[["low"]]
SPPL2A_High_OS <- SPPL2A_High_OS[,-1:-27]
SPPL2A_Low_OS <- SPPL2A_Low_OS[,-1:-27]
SPPL2A_OS_split <- data.frame(SPPL2A_High_OS, cbind(zoo(, 1:nrow(SPPL2A_High_OS)), as.zoo(SPPL2A_Low_OS)))
colnames(SPPL2A_OS_split) <- colnamesciber

TAPBPL_OS <-split(mergeciber_OS, mergeciber_OS$TAPBPL)
TAPBPL_High_OS <- TAPBPL_OS[["high"]]
TAPBPL_Low_OS <- TAPBPL_OS[["low"]]
TAPBPL_High_OS <- TAPBPL_High_OS[,-1:-27]
TAPBPL_Low_OS <- TAPBPL_Low_OS[,-1:-27]
TAPBPL_OS_split <- data.frame(TAPBPL_High_OS, cbind(zoo(, 1:nrow(TAPBPL_High_OS)), as.zoo(TAPBPL_Low_OS)))
colnames(TAPBPL_OS_split) <- colnamesciber

TAPBP_OS <-split(mergeciber_OS, mergeciber_OS$TAPBP)
TAPBP_High_OS <- TAPBP_OS[["high"]]
TAPBP_Low_OS <- TAPBP_OS[["low"]]
TAPBP_High_OS <- TAPBP_High_OS[,-1:-27]
TAPBP_Low_OS <- TAPBP_Low_OS[,-1:-27]
TAPBP_OS_split <- data.frame(TAPBP_High_OS, cbind(zoo(, 1:nrow(TAPBP_High_OS)), as.zoo(TAPBP_Low_OS)))
colnames(TAPBP_OS_split) <- colnamesciber

ERAP2_OS <-split(mergeciber_OS, mergeciber_OS$ERAP2)
ERAP2_High_OS <- ERAP2_OS[["high"]]
ERAP2_Low_OS <- ERAP2_OS[["low"]]
ERAP2_High_OS <- ERAP2_High_OS[,-1:-27]
ERAP2_Low_OS <- ERAP2_Low_OS[,-1:-27]
ERAP2_OS_split <- data.frame(ERAP2_High_OS, cbind(zoo(, 1:nrow(ERAP2_High_OS)), as.zoo(ERAP2_Low_OS)))
colnames(ERAP2_OS_split) <- colnamesciber

ERAP2_OS <-split(mergeciber_OS, mergeciber_OS$ERAP2)
ERAP2_High_OS <- ERAP2_OS[["high"]]
ERAP2_Low_OS <- ERAP2_OS[["low"]]
ERAP2_High_OS <- ERAP2_High_OS[,-1:-27]
ERAP2_Low_OS <- ERAP2_Low_OS[,-1:-27]
ERAP2_OS_split <- data.frame(ERAP2_High_OS, cbind(zoo(, 1:nrow(ERAP2_High_OS)), as.zoo(ERAP2_Low_OS)))
colnames(ERAP2_OS_split) <- colnamesciber

CTSS_OS <-split(mergeciber_OS, mergeciber_OS$CTSS)
CTSS_High_OS <- CTSS_OS[["high"]]
CTSS_Low_OS <- CTSS_OS[["low"]]
CTSS_High_OS <- CTSS_High_OS[,-1:-27]
CTSS_Low_OS <- CTSS_Low_OS[,-1:-27]
CTSS_OS_split <- data.frame(CTSS_High_OS, cbind(zoo(, 1:nrow(CTSS_High_OS)), as.zoo(CTSS_Low_OS)))
colnames(CTSS_OS_split) <- colnamesciber

CTSL_OS <-split(mergeciber_OS, mergeciber_OS$CTSL)
CTSL_High_OS <- CTSL_OS[["high"]]
CTSL_Low_OS <- CTSL_OS[["low"]]
CTSL_High_OS <- CTSL_High_OS[,-1:-27]
CTSL_Low_OS <- CTSL_Low_OS[,-1:-27]
CTSL_OS_split <- data.frame(CTSL_Low_OS, cbind(zoo(, 1:nrow(CTSL_Low_OS)), as.zoo(CTSL_High_OS)))
colnames(CTSL_OS_split) <- colnamesciberopposite

CSDE1_OS <-split(mergeciber_OS, mergeciber_OS$CSDE1)
CSDE1_High_OS <- CSDE1_OS[["high"]]
CSDE1_Low_OS <- CSDE1_OS[["low"]]
CSDE1_High_OS <- CSDE1_High_OS[,-1:-27]
CSDE1_Low_OS <- CSDE1_Low_OS[,-1:-27]
CSDE1_OS_split <- data.frame(CSDE1_Low_OS, cbind(zoo(, 1:nrow(CSDE1_Low_OS)), as.zoo(CSDE1_High_OS)))
colnames(CSDE1_OS_split) <- colnamesciberopposite

CIITA_OS <-split(mergeciber_OS, mergeciber_OS$CIITA)
CIITA_High_OS <- CIITA_OS[["high"]]
CIITA_Low_OS <- CIITA_OS[["low"]]
CIITA_High_OS <- CIITA_High_OS[,-1:-27]
CIITA_Low_OS <- CIITA_Low_OS[,-1:-27]
CIITA_OS_split <- data.frame(CIITA_Low_OS, cbind(zoo(, 1:nrow(CIITA_Low_OS)), as.zoo(CIITA_High_OS)))
colnames(CIITA_OS_split) <- colnamesciberopposite

CD74_OS <-split(mergeciber_OS, mergeciber_OS$CD74)
CD74_High_OS <- CD74_OS[["high"]]
CD74_Low_OS <- CD74_OS[["low"]]
CD74_High_OS <- CD74_High_OS[,-1:-27]
CD74_Low_OS <- CD74_Low_OS[,-1:-27]
CD74_OS_split <- data.frame(CD74_Low_OS, cbind(zoo(, 1:nrow(CD74_Low_OS)), as.zoo(CD74_High_OS)))
colnames(CD74_OS_split) <- colnamesciberopposite

CD1D_OS <-split(mergeciber_OS, mergeciber_OS$CD1D)
CD1D_High_OS <- CD1D_OS[["high"]]
CD1D_Low_OS <- CD1D_OS[["low"]]
CD1D_High_OS <- CD1D_High_OS[,-1:-27]
CD1D_Low_OS <- CD1D_Low_OS[,-1:-27]
CD1D_OS_split <- data.frame(CD1D_Low_OS, cbind(zoo(, 1:nrow(CD1D_Low_OS)), as.zoo(CD1D_High_OS)))
colnames(CD1D_OS_split) <- colnamesciberopposite

CALR_OS <-split(mergeciber_OS, mergeciber_OS$CALR)
CALR_High_OS <- CALR_OS[["high"]]
CALR_Low_OS <- CALR_OS[["low"]]
CALR_High_OS <- CALR_High_OS[,-1:-27]
CALR_Low_OS <- CALR_Low_OS[,-1:-27]
CALR_OS_split <- data.frame(CALR_Low_OS, cbind(zoo(, 1:nrow(CALR_Low_OS)), as.zoo(CALR_High_OS)))
colnames(CALR_OS_split) <- colnamesciberopposite

library(tidyverse)
dfs <- list(HLA_A_OS_split = HLA_A_OS_split,
            HLA_B_OS_split = HLA_B_OS_split,
            HLA_E_OS_split = HLA_E_OS_split,
            HLA_G_OS_split = HLA_G_OS_split,
            HLA_DPA1_OS_split  = HLA_DPA1_OS_split,
            HLA_DQA1_OS_split = HLA_DQA1_OS_split,
            HLA_DRA_OS_split = HLA_DRA_OS_split,
            HLA_DRB1_OS_split = HLA_DRB1_OS_split,
            HLA_DRB5_OS_split = HLA_DRB5_OS_split,
            LGMN_OS_split = LGMN_OS_split,
            MR1_OS_split = MR1_OS_split,
            PSMB10_OS_split = PSMB10_OS_split,
            PSMB8_OS_split = PSMB8_OS_split,
            PSMB9_OS_split = PSMB9_OS_split,
            RFX5_OS_split = RFX5_OS_split,
            SPPL2A_OS_split = SPPL2A_OS_split,
            TAPBPL_OS_split = TAPBPL_OS_split,
            TAPBP_OS_split = TAPBP_OS_split,
            ERAP2_OS_split = ERAP2_OS_split,
            ERAP2_OS_split = ERAP2_OS_split,
            CTSS_OS_split = CTSS_OS_split,
            CTSL_OS_split = CTSL_OS_split,
            CSDE1_OS_split = CSDE1_OS_split,
            CIITA_OS_split = CIITA_OS_split,
            CD74_OS_split = CD74_OS_split,
            CD1D_OS_split = CD1D_OS_split,
            CALR_OS_split = CALR_OS_split)
setwd("C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/OSCibermerge/OS")
walk2(dfs, paste0(names(dfs), ".csv"), write_csv)

#Quantile analysis for gene approaching significance

#read in quantile data
UpperquantileCSDE1 <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileCSDE1.csv", header = TRUE, row.names =1)
LowerquantileCSDE1 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileCSDE1.csv", header = TRUE, row.names =1)
UpperquantileHLA.E <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileHLA.E.csv", header = TRUE, row.names =1)
LowerquantileHLA.E <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileHLA.E.csv", header = TRUE, row.names =1)
UpperquantileTAPBP <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileTAPBP.csv", header = TRUE, row.names =1)
LowerquantileTAPBP <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileTAPBP.csv", header = TRUE, row.names =1)
UpperquantilePSMB9 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantilePSMB9.csv", header = TRUE, row.names =1)
LowerquantilePSMB9 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantilePSMB9.csv", header = TRUE, row.names =1)
UpperquantileSPPL2A <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileSPPL2A.csv", header = TRUE, row.names =1)
LowerquantileSPPL2A <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileSPPL2A.csv", header = TRUE, row.names =1)
UpperquantileHLA.DRB1 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileHLA.DRB1.csv", header = TRUE, row.names =1)
LowerquantileHLA.DRB1 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileHLA.DRB1.csv", header = TRUE, row.names =1)
UpperquantileHLA.DQA1 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileHLA.DQA1.csv", header = TRUE, row.names =1)
LowerquantileHLA.DQA1 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileHLA.DQA1.csv", header = TRUE, row.names =1)
UpperquantileCIITA <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileCIITA.csv", header = TRUE, row.names =1)
LowerquantileCIITA <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileCIITA.csv", header = TRUE, row.names =1)
UpperquantileMR1 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileMR1.csv", header = TRUE, row.names =1)
LowerquantileMR1 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileMR1.csv", header = TRUE, row.names =1)
UpperquantileCD1D <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileCD1D.csv", header = TRUE, row.names =1)
LowerquantileCD1D <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileCD1D.csv", header = TRUE, row.names =1)
UpperquantileCALR <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileCALR.csv", header = TRUE, row.names =1)
LowerquantileCALR <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileCALR.csv", header = TRUE, row.names =1)
UpperquantilePSMB8 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantilePSMB8.csv", header = TRUE, row.names =1)
LowerquantilePSMB8 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantilePSMB8.csv", header = TRUE, row.names =1)
UpperquantilePSMB9 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantilePSMB9.csv", header = TRUE, row.names =1)
LowerquantilePSMB9 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantilePSMB9.csv", header = TRUE, row.names =1)
UpperquantileHLA.DRA <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileHLA.DRA.csv", header = TRUE, row.names =1)
LowerquantileHLA.DRA <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileHLA.DRA.csv", header = TRUE, row.names =1)
UpperquantileHLA.DPA1 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileHLA.DPA1.csv", header = TRUE, row.names =1)
LowerquantileHLA.DPA1 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileHLA.DPA1.csv", header = TRUE, row.names =1)
UpperquantileLGMN <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileLGMN.csv", header = TRUE, row.names =1)
LowerquantileLGMN <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileLGMN.csv", header = TRUE, row.names =1)
UpperquantileCTSL <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileCTSL.csv", header = TRUE, row.names =1)
LowerquantileCTSL <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileCTSL.csv", header = TRUE, row.names =1)
UpperquantilePSMB10 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantilePSMB10.csv", header = TRUE, row.names =1)
LowerquantilePSMB10 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantilePSMB10.csv", header = TRUE, row.names =1)
UpperquantileCD74 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileCD74.csv", header = TRUE, row.names =1)
LowerquantileCD74 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileCD74.csv", header = TRUE, row.names =1)
UpperquantileERAP2 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileERAP2.csv", header = TRUE, row.names =1)
LowerquantileERAP2 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileERAP2.csv", header = TRUE, row.names =1)
UpperquantileERAP1 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileERAP1.csv", header = TRUE, row.names =1)
LowerquantileERAP1 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileERAP1.csv", header = TRUE, row.names =1)
UpperquantileTAPBPL <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileTAPBPL.csv", header = TRUE, row.names =1)
LowerquantileTAPBPL <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileTAPBPL.csv", header = TRUE, row.names =1)
UpperquantileCTSS <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileCTSS.csv", header = TRUE, row.names =1)
LowerquantileCTSS <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileCTSS.csv", header = TRUE, row.names =1)
UpperquantileHLA.DRB5 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileHLA.DRB5.csv", header = TRUE, row.names =1)
LowerquantileHLA.DRB5 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileHLA.DRB5.csv", header = TRUE, row.names =1)
UpperquantileHLA.A <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileHLA.A.csv", header = TRUE, row.names =1)
LowerquantileHLA.A <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileHLA.A.csv", header = TRUE, row.names =1)
UpperquantileHLA.B <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileHLA.B.csv", header = TRUE, row.names =1)
LowerquantileHLA.B <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileHLA.B.csv", header = TRUE, row.names =1)
UpperquantileHLA.G <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileHLA.G.csv", header = TRUE, row.names =1)
LowerquantileHLA.G <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileHLA.G.csv", header = TRUE, row.names =1)
UpperquantileRFX5 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/UpperquantileRFX5.csv", header = TRUE, row.names =1)
LowerquantileRFX5 <- read.csv( file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/LowerquantileRFX5.csv", header = TRUE, row.names =1)



#make up quantile defining dataframes
DUpperquantileCSDE1 <- as.data.frame(rownames(UpperquantileCSDE1))
DUpperquantileCSDE1$Quantile <- "Upper_Quantile"
colnames(DUpperquantileCSDE1) <- c("ID", "Quantile")
DLowerquantileCSDE1 <- as.data.frame(rownames(LowerquantileCSDE1))
DLowerquantileCSDE1$Quantile <- "Lower_Quantile"
colnames(DLowerquantileCSDE1) <- c("ID", "Quantile")
CSDE1quantiles <- rbind(DUpperquantileCSDE1,DLowerquantileCSDE1)

DUpperquantileHLA.E <- as.data.frame(rownames(UpperquantileHLA.E))
DUpperquantileHLA.E$Quantile <- "Upper_Quantile"
colnames(DUpperquantileHLA.E) <- c("ID", "Quantile")
DLowerquantileHLA.E <- as.data.frame(rownames(LowerquantileHLA.E))
DLowerquantileHLA.E$Quantile <- "Lower_Quantile"
colnames(DLowerquantileHLA.E) <- c("ID", "Quantile")
HLA.Equantiles <- rbind(DUpperquantileHLA.E,DLowerquantileHLA.E)

DUpperquantileTAPBP <- as.data.frame(rownames(UpperquantileTAPBP))
DUpperquantileTAPBP$Quantile <- "Upper_Quantile"
colnames(DUpperquantileTAPBP) <- c("ID", "Quantile")
DLowerquantileTAPBP <- as.data.frame(rownames(LowerquantileTAPBP))
DLowerquantileTAPBP$Quantile <- "Lower_Quantile"
colnames(DLowerquantileTAPBP) <- c("ID", "Quantile")
TAPBPquantiles <- rbind(DUpperquantileTAPBP,DLowerquantileTAPBP)

DUpperquantilePSMB9 <- as.data.frame(rownames(UpperquantilePSMB9))
DUpperquantilePSMB9$Quantile <- "Upper_Quantile"
colnames(DUpperquantilePSMB9) <- c("ID", "Quantile")
DLowerquantilePSMB9 <- as.data.frame(rownames(LowerquantilePSMB9))
DLowerquantilePSMB9$Quantile <- "Lower_Quantile"
colnames(DLowerquantilePSMB9) <- c("ID", "Quantile")
PSMB9quantiles <- rbind(DUpperquantilePSMB9,DLowerquantilePSMB9)

DUpperquantileSPPL2A <- as.data.frame(rownames(UpperquantileSPPL2A))
DUpperquantileSPPL2A$Quantile <- "Upper_Quantile"
colnames(DUpperquantileSPPL2A) <- c("ID", "Quantile")
DLowerquantileSPPL2A <- as.data.frame(rownames(LowerquantileSPPL2A))
DLowerquantileSPPL2A$Quantile <- "Lower_Quantile"
colnames(DLowerquantileSPPL2A) <- c("ID", "Quantile")
SPPL2Aquantiles <- rbind(DUpperquantileSPPL2A,DLowerquantileSPPL2A)

DUpperquantileHLA.DRB1 <- as.data.frame(rownames(UpperquantileHLA.DRB1))
DUpperquantileHLA.DRB1$Quantile <- "Upper_Quantile"
colnames(DUpperquantileHLA.DRB1) <- c("ID", "Quantile")
DLowerquantileHLA.DRB1 <- as.data.frame(rownames(LowerquantileHLA.DRB1))
DLowerquantileHLA.DRB1$Quantile <- "Lower_Quantile"
colnames(DLowerquantileHLA.DRB1) <- c("ID", "Quantile")
HLA.DRB1quantiles <- rbind(DUpperquantileHLA.DRB1,DLowerquantileHLA.DRB1)

DUpperquantileHLA.DQA1 <- as.data.frame(rownames(UpperquantileHLA.DQA1))
DUpperquantileHLA.DQA1$Quantile <- "Upper_Quantile"
colnames(DUpperquantileHLA.DQA1) <- c("ID", "Quantile")
DLowerquantileHLA.DQA1 <- as.data.frame(rownames(LowerquantileHLA.DQA1))
DLowerquantileHLA.DQA1$Quantile <- "Lower_Quantile"
colnames(DLowerquantileHLA.DQA1) <- c("ID", "Quantile")
HLA.DQA1quantiles <- rbind(DUpperquantileHLA.DQA1,DLowerquantileHLA.DQA1)

DUpperquantileCIITA <- as.data.frame(rownames(UpperquantileCIITA))
DUpperquantileCIITA$Quantile <- "Upper_Quantile"
colnames(DUpperquantileCIITA) <- c("ID", "Quantile")
DLowerquantileCIITA <- as.data.frame(rownames(LowerquantileCIITA))
DLowerquantileCIITA$Quantile <- "Lower_Quantile"
colnames(DLowerquantileCIITA) <- c("ID", "Quantile")
CIITAquantiles <- rbind(DUpperquantileCIITA,DLowerquantileCIITA)

DUpperquantileMR1 <- as.data.frame(rownames(UpperquantileMR1))
DUpperquantileMR1$Quantile <- "Upper_Quantile"
colnames(DUpperquantileMR1) <- c("ID", "Quantile")
DLowerquantileMR1 <- as.data.frame(rownames(LowerquantileMR1))
DLowerquantileMR1$Quantile <- "Lower_Quantile"
colnames(DLowerquantileMR1) <- c("ID", "Quantile")
MR1quantiles <- rbind(DUpperquantileMR1,DLowerquantileMR1)


DUpperquantileCD1D <- as.data.frame(rownames(UpperquantileCD1D))
DUpperquantileCD1D$Quantile <- "Upper_Quantile"
colnames(DUpperquantileCD1D) <- c("ID", "Quantile")
DLowerquantileCD1D <- as.data.frame(rownames(LowerquantileCD1D))
DLowerquantileCD1D$Quantile <- "Lower_Quantile"
colnames(DLowerquantileCD1D) <- c("ID", "Quantile")
CD1Dquantiles <- rbind(DUpperquantileCD1D,DLowerquantileCD1D)

DUpperquantileERAP2 <- as.data.frame(rownames(UpperquantileERAP2))
DUpperquantileERAP2$Quantile <- "Upper_Quantile"
colnames(DUpperquantileERAP2) <- c("ID", "Quantile")
DLowerquantileERAP2 <- as.data.frame(rownames(LowerquantileERAP2))
DLowerquantileERAP2$Quantile <- "Lower_Quantile"
colnames(DLowerquantileERAP2) <- c("ID", "Quantile")
ERAP2quantiles <- rbind(DUpperquantileERAP2,DLowerquantileERAP2)

DUpperquantileCALR <- as.data.frame(rownames(UpperquantileCALR))
DUpperquantileCALR$Quantile <- "Upper_Quantile"
colnames(DUpperquantileCALR) <- c("ID", "Quantile")
DLowerquantileCALR <- as.data.frame(rownames(LowerquantileCALR))
DLowerquantileCALR$Quantile <- "Lower_Quantile"
colnames(DLowerquantileCALR) <- c("ID", "Quantile")
CALRquantiles <- rbind(DUpperquantileCALR,DLowerquantileCALR)

DUpperquantilePSMB9 <- as.data.frame(rownames(UpperquantilePSMB9))
DUpperquantilePSMB9$Quantile <- "Upper_Quantile"
colnames(DUpperquantilePSMB9) <- c("ID", "Quantile")
DLowerquantilePSMB9 <- as.data.frame(rownames(LowerquantilePSMB9))
DLowerquantilePSMB9$Quantile <- "Lower_Quantile"
colnames(DLowerquantilePSMB9) <- c("ID", "Quantile")
PSMB9quantiles <- rbind(DUpperquantilePSMB9,DLowerquantilePSMB9)

DUpperquantileHLA.DRA <- as.data.frame(rownames(UpperquantileHLA.DRA))
DUpperquantileHLA.DRA$Quantile <- "Upper_Quantile"
colnames(DUpperquantileHLA.DRA) <- c("ID", "Quantile")
DLowerquantileHLA.DRA <- as.data.frame(rownames(LowerquantileHLA.DRA))
DLowerquantileHLA.DRA$Quantile <- "Lower_Quantile"
colnames(DLowerquantileHLA.DRA) <- c("ID", "Quantile")
HLA.DRAquantiles <- rbind(DUpperquantileHLA.DRA,DLowerquantileHLA.DRA)

DUpperquantileHLA.DPA1 <- as.data.frame(rownames(UpperquantileHLA.DPA1))
DUpperquantileHLA.DPA1$Quantile <- "Upper_Quantile"
colnames(DUpperquantileHLA.DPA1) <- c("ID", "Quantile")
DLowerquantileHLA.DPA1 <- as.data.frame(rownames(LowerquantileHLA.DPA1))
DLowerquantileHLA.DPA1$Quantile <- "Lower_Quantile"
colnames(DLowerquantileHLA.DPA1) <- c("ID", "Quantile")
HLA.DPA1quantiles <- rbind(DUpperquantileHLA.DPA1,DLowerquantileHLA.DPA1)

DUpperquantileLGMN <- as.data.frame(rownames(UpperquantileLGMN))
DUpperquantileLGMN$Quantile <- "Upper_Quantile"
colnames(DUpperquantileLGMN) <- c("ID", "Quantile")
DLowerquantileLGMN <- as.data.frame(rownames(LowerquantileLGMN))
DLowerquantileLGMN$Quantile <- "Lower_Quantile"
colnames(DLowerquantileLGMN) <- c("ID", "Quantile")
LGMNquantiles <- rbind(DUpperquantileLGMN,DLowerquantileLGMN)

DUpperquantileCTSL <- as.data.frame(rownames(UpperquantileCTSL))
DUpperquantileCTSL$Quantile <- "Upper_Quantile"
colnames(DUpperquantileCTSL) <- c("ID", "Quantile")
DLowerquantileCTSL <- as.data.frame(rownames(LowerquantileCTSL))
DLowerquantileCTSL$Quantile <- "Lower_Quantile"
colnames(DLowerquantileCTSL) <- c("ID", "Quantile")
CTSLquantiles <- rbind(DUpperquantileCTSL,DLowerquantileCTSL)

DUpperquantilePSMB10 <- as.data.frame(rownames(UpperquantilePSMB10))
DUpperquantilePSMB10$Quantile <- "Upper_Quantile"
colnames(DUpperquantilePSMB10) <- c("ID", "Quantile")
DLowerquantilePSMB10 <- as.data.frame(rownames(LowerquantilePSMB10))
DLowerquantilePSMB10$Quantile <- "Lower_Quantile"
colnames(DLowerquantilePSMB10) <- c("ID", "Quantile")
PSMB10quantiles <- rbind(DUpperquantilePSMB10,DLowerquantilePSMB10)

DUpperquantileCD74 <- as.data.frame(rownames(UpperquantileCD74))
DUpperquantileCD74$Quantile <- "Upper_Quantile"
colnames(DUpperquantileCD74) <- c("ID", "Quantile")
DLowerquantileCD74 <- as.data.frame(rownames(LowerquantileCD74))
DLowerquantileCD74$Quantile <- "Lower_Quantile"
colnames(DLowerquantileCD74) <- c("ID", "Quantile")
CD74quantiles <- rbind(DUpperquantileCD74,DLowerquantileCD74)

DUpperquantilePSMB8 <- as.data.frame(rownames(UpperquantilePSMB8))
DUpperquantilePSMB8$Quantile <- "Upper_Quantile"
colnames(DUpperquantilePSMB8) <- c("ID", "Quantile")
DLowerquantilePSMB8 <- as.data.frame(rownames(LowerquantilePSMB8))
DLowerquantilePSMB8$Quantile <- "Lower_Quantile"
colnames(DLowerquantilePSMB8) <- c("ID", "Quantile")
PSMB8quantiles <- rbind(DUpperquantilePSMB8,DLowerquantilePSMB8)

DUpperquantileERAP1 <- as.data.frame(rownames(UpperquantileERAP1))
DUpperquantileERAP1$Quantile <- "Upper_Quantile"
colnames(DUpperquantileERAP1) <- c("ID", "Quantile")
DLowerquantileERAP1 <- as.data.frame(rownames(LowerquantileERAP1))
DLowerquantileERAP1$Quantile <- "Lower_Quantile"
colnames(DLowerquantileERAP1) <- c("ID", "Quantile")
ERAP1quantiles <- rbind(DUpperquantileERAP1,DLowerquantileERAP1)

DUpperquantileTAPBPL <- as.data.frame(rownames(UpperquantileTAPBPL))
DUpperquantileTAPBPL$Quantile <- "Upper_Quantile"
colnames(DUpperquantileTAPBPL) <- c("ID", "Quantile")
DLowerquantileTAPBPL <- as.data.frame(rownames(LowerquantileTAPBPL))
DLowerquantileTAPBPL$Quantile <- "Lower_Quantile"
colnames(DLowerquantileTAPBPL) <- c("ID", "Quantile")
TAPBPLquantiles <- rbind(DUpperquantileTAPBPL,DLowerquantileTAPBPL)

DUpperquantileCTSS <- as.data.frame(rownames(UpperquantileCTSS))
DUpperquantileCTSS$Quantile <- "Upper_Quantile"
colnames(DUpperquantileCTSS) <- c("ID", "Quantile")
DLowerquantileCTSS <- as.data.frame(rownames(LowerquantileCTSS))
DLowerquantileCTSS$Quantile <- "Lower_Quantile"
colnames(DLowerquantileCTSS) <- c("ID", "Quantile")
CTSSquantiles <- rbind(DUpperquantileCTSS,DLowerquantileCTSS)

DUpperquantileHLA.DRB5 <- as.data.frame(rownames(UpperquantileHLA.DRB5))
DUpperquantileHLA.DRB5$Quantile <- "Upper_Quantile"
colnames(DUpperquantileHLA.DRB5) <- c("ID", "Quantile")
DLowerquantileHLA.DRB5 <- as.data.frame(rownames(LowerquantileHLA.DRB5))
DLowerquantileHLA.DRB5$Quantile <- "Lower_Quantile"
colnames(DLowerquantileHLA.DRB5) <- c("ID", "Quantile")
HLA.DRB5quantiles <- rbind(DUpperquantileHLA.DRB5,DLowerquantileHLA.DRB5)

DUpperquantileHLA.A <- as.data.frame(rownames(UpperquantileHLA.A))
DUpperquantileHLA.A$Quantile <- "Upper_Quantile"
colnames(DUpperquantileHLA.A) <- c("ID", "Quantile")
DLowerquantileHLA.A <- as.data.frame(rownames(LowerquantileHLA.A))
DLowerquantileHLA.A$Quantile <- "Lower_Quantile"
colnames(DLowerquantileHLA.A) <- c("ID", "Quantile")
HLA.Aquantiles <- rbind(DUpperquantileHLA.A,DLowerquantileHLA.A)

DUpperquantileHLA.B <- as.data.frame(rownames(UpperquantileHLA.B))
DUpperquantileHLA.B$Quantile <- "Upper_Quantile"
colnames(DUpperquantileHLA.B) <- c("ID", "Quantile")
DLowerquantileHLA.B <- as.data.frame(rownames(LowerquantileHLA.B))
DLowerquantileHLA.B$Quantile <- "Lower_Quantile"
colnames(DLowerquantileHLA.B) <- c("ID", "Quantile")
HLA.Bquantiles <- rbind(DUpperquantileHLA.B,DLowerquantileHLA.B)

DUpperquantileHLA.G <- as.data.frame(rownames(UpperquantileHLA.G))
DUpperquantileHLA.G$Quantile <- "Upper_Quantile"
colnames(DUpperquantileHLA.G) <- c("ID", "Quantile")
DLowerquantileHLA.G <- as.data.frame(rownames(LowerquantileHLA.G))
DLowerquantileHLA.G$Quantile <- "Lower_Quantile"
colnames(DLowerquantileHLA.G) <- c("ID", "Quantile")
HLA.Gquantiles <- rbind(DUpperquantileHLA.G,DLowerquantileHLA.G)

DUpperquantileRFX5 <- as.data.frame(rownames(UpperquantileRFX5))
DUpperquantileRFX5$Quantile <- "Upper_Quantile"
colnames(DUpperquantileRFX5) <- c("ID", "Quantile")
DLowerquantileRFX5 <- as.data.frame(rownames(LowerquantileRFX5))
DLowerquantileRFX5$Quantile <- "Lower_Quantile"
colnames(DLowerquantileRFX5) <- c("ID", "Quantile")
RFX5quantiles <- rbind(DUpperquantileRFX5,DLowerquantileRFX5)


#left join cibersort data to quantile data
cibersortleft <- cibersort
cibersortleft$ID <- rownames(cibersortleft)
CSDE1quantilecibersort <- left_join(CSDE1quantiles, cibersortleft)
HLA.Equantilecibersort <- left_join(HLA.Equantiles, cibersortleft)
TAPBPquantilecibersort <- left_join(TAPBPquantiles, cibersortleft)
SPPL2Aquantilecibersort <- left_join(SPPL2Aquantiles, cibersortleft)
HLA.DRB1quantilecibersort <- left_join(HLA.DRB1quantiles, cibersortleft)
HLA.DQA1quantilecibersort <- left_join(HLA.DQA1quantiles, cibersortleft)
CITTAquantilecibersort <- left_join(CIITAquantiles, cibersortleft)
MR1quantilecibersort <- left_join(MR1quantiles, cibersortleft)
CD1Dquantilecibersort <- left_join(CD1Dquantiles, cibersortleft)
ERAP2Dquantilecibersort <- left_join(ERAP2quantiles, cibersortleft)
CALRDquantilecibersort <- left_join(CALRquantiles, cibersortleft)
PSMB8Dquantilecibersort <- left_join(PSMB8quantiles, cibersortleft)
PSMB9Dquantilecibersort <- left_join(PSMB9quantiles, cibersortleft)
PSMB10Dquantilecibersort <- left_join(PSMB10quantiles, cibersortleft)
HLA.DRADquantilecibersort <- left_join(HLA.DRAquantiles, cibersortleft)
HLA.DPA1Dquantilecibersort <- left_join(HLA.DPA1quantiles, cibersortleft)
LGMNDquantilecibersort <- left_join(LGMNquantiles, cibersortleft)
CTSLDquantilecibersort <- left_join(CTSLquantiles, cibersortleft)
CD74Dquantilecibersort <- left_join(CD74quantiles, cibersortleft)
CTSLDquantilecibersort <- left_join(CTSLquantiles, cibersortleft)
ERAP1Dquantilecibersort <- left_join(ERAP1quantiles, cibersortleft)
TAPBPLDquantilecibersort <- left_join(TAPBPLquantiles, cibersortleft)
CTSSDquantilecibersort <- left_join(CTSSquantiles, cibersortleft)
HLA.DRB5Dquantilecibersort <- left_join(HLA.DRB5quantiles, cibersortleft)
HLA.ADquantilecibersort <- left_join(HLA.Aquantiles, cibersortleft)
HLA.BDquantilecibersort <- left_join(HLA.Bquantiles, cibersortleft)
HLA.GDquantilecibersort <- left_join(HLA.Gquantiles, cibersortleft)
RFX5Dquantilecibersort <- left_join(RFX5quantiles, cibersortleft)
#merge with cluster data
clusterHEATMAP <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/Column_Annotation_Heatmap_Cibersort.csv", header = TRUE)
clusterHEATMAP <- clusterHEATMAP[c("X", "cluster")]

#MAtch up cluster data
CSDE1quantilecibersort <- merge(CSDE1quantilecibersort, clusterHEATMAP, by = 1)
HLA.Equantilecibersort <- merge(HLA.Equantilecibersort, clusterHEATMAP, by = 1)
TAPBPquantilecibersort <- merge(TAPBPquantilecibersort, clusterHEATMAP, by = 1)
SPPL2Aquantilecibersort <- merge(SPPL2Aquantilecibersort, clusterHEATMAP, by = 1)
HLA.DRB1quantilecibersort <- merge(HLA.DRB1quantilecibersort, clusterHEATMAP, by = 1)
HLA.DQA1quantilecibersort <- merge(HLA.DQA1quantilecibersort, clusterHEATMAP, by = 1)
CITTAquantilecibersort<- merge(CITTAquantilecibersort, clusterHEATMAP, by = 1)
MR1quantilecibersort <- merge(MR1quantilecibersort, clusterHEATMAP, by = 1)
CD1Dquantilecibersort <- merge(CD1Dquantilecibersort, clusterHEATMAP, by = 1)
ERAP2Dquantilecibersort <- merge(ERAP2Dquantilecibersort, clusterHEATMAP, by = 1)
CALRDquantilecibersort <- merge(CALRDquantilecibersort, clusterHEATMAP, by = 1)
PSMB8Dquantilecibersort <- merge(PSMB8Dquantilecibersort, clusterHEATMAP, by = 1)
PSMB9Dquantilecibersort <- merge(PSMB9Dquantilecibersort, clusterHEATMAP, by = 1)
PSMB10Dquantilecibersort <- merge(PSMB10Dquantilecibersort, clusterHEATMAP, by = 1)
HLA.DRADquantilecibersort <- merge(HLA.DRADquantilecibersort, clusterHEATMAP, by = 1)
LGMNDquantilecibersort <- merge(LGMNDquantilecibersort, clusterHEATMAP, by = 1)
CTSLDquantilecibersort <- merge(CTSLDquantilecibersort, clusterHEATMAP, by = 1)
CD74Dquantilecibersort <- merge(CD74Dquantilecibersort, clusterHEATMAP, by = 1)
CTSLDquantilecibersort <- merge(CTSLDquantilecibersort, clusterHEATMAP, by = 1)
HLA.DPA1Dquantilecibersort <- merge(HLA.DPA1Dquantilecibersort, clusterHEATMAP, by = 1)
ERAP1Dquantilecibersort <- merge(ERAP1Dquantilecibersort, clusterHEATMAP, by = 1)
TAPBPLDquantilecibersort <- merge(TAPBPLDquantilecibersort, clusterHEATMAP, by = 1)
CTSSDquantilecibersort <- merge(CTSSDquantilecibersort, clusterHEATMAP, by = 1)
HLA.DRB5Dquantilecibersort <- merge(HLA.DRB5Dquantilecibersort, clusterHEATMAP, by = 1)
HLA.GDquantilecibersort <- merge(HLA.GDquantilecibersort, clusterHEATMAP, by = 1)
HLA.BDquantilecibersort <- merge(HLA.BDquantilecibersort, clusterHEATMAP, by = 1)
HLA.ADquantilecibersort <- merge(HLA.ADquantilecibersort, clusterHEATMAP, by = 1)
RFX5Dquantilecibersort <- merge(RFX5Dquantilecibersort, clusterHEATMAP, by = 1)

dfs <- list(CSDE1quantilecibersort = CSDE1quantilecibersort,
            HLA.Equantilecibersort = HLA.Equantilecibersort,
            TAPBPquantilecibersort = TAPBPquantilecibersort,
            SPPL2Aquantilecibersort = SPPL2Aquantilecibersort,
            HLA.DRB1quantilecibersort = HLA.DRB1quantilecibersort,
            HLA.DQA1quantilecibersort = HLA.DQA1quantilecibersort,
            CITTAquantilecibersort = CITTAquantilecibersort,
            MR1quantilecibersort = MR1quantilecibersort,
            CD1Dquantilecibersort = CD1Dquantilecibersort,
            ERAP2Dquantilecibersort = ERAP2Dquantilecibersort,
            CALRDquantilecibersort = CALRDquantilecibersort,
            PSMB8Dquantilecibersort = PSMB8Dquantilecibersort,
            PSMB9Dquantilecibersort = PSMB9Dquantilecibersort,
            PSMB10Dquantilecibersort = PSMB10Dquantilecibersort,
            HLA.DRADquantilecibersort = HLA.DRADquantilecibersort,
            LGMNDquantilecibersort = LGMNDquantilecibersort,
            CTSLDquantilecibersort = CTSLDquantilecibersort,
            CD74Dquantilecibersort = CD74Dquantilecibersort,
            CTSLDquantilecibersort = CTSLDquantilecibersort,
            HLA.DPA1Dquantilecibersort = HLA.DPA1Dquantilecibersort,
            ERAP1Dquantilecibersort = ERAP1Dquantilecibersort,
            TAPBPLDquantilecibersort = TAPBPLDquantilecibersort,
            CTSSDquantilecibersort = CTSSDquantilecibersort,
            HLA.DRB5Dquantilecibersort = HLA.DRB5Dquantilecibersort,
            HLA.GDquantilecibersort = HLA.GDquantilecibersort,
            HLA.BDquantilecibersort = HLA.BDquantilecibersort,
            HLA.ADquantilecibersort = HLA.ADquantilecibersort,
            CALRDquantilecibersort = CALRDquantilecibersort,
            RFX5Dquantilecibersort = RFX5Dquantilecibersort)
setwd("C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/OSCibermerge/Quantile")
walk2(dfs, paste0(names(dfs), ".csv"), write_csv)
  
CALRDquantilecibersort2 <- CALRDquantilecibersort[,2:11]

modelList<-list()
for(i in 2:11){
  fmla <- formula(paste(names(CALRDquantilecibersort2)[i], " ~ Quantile"))
  modelList[[i]]<-wilcox.test(fmla, data = CALRDquantilecibersort2, paired = FALSE)
}


modelList<-list()
for(i in 2:11){
  fmla <- formula(paste(names(CALRDquantilecibersort2)[i], " ~ Quantile"))
  modelList[[i]]<-wilcox.test(fmla, data = CALRDquantilecibersort2)
}
  
RFX5Dquantilecibersort2 <- RFX5Dquantilecibersort[,2:12]

cols <- names(RFX5Dquantilecibersort2)[2:ncol(RFX5Dquantilecibersort2)]

all_test <- lapply(cols, function(x) 
  kruskal.test(reformulate("Quantile", x), data = RFX5Dquantilecibersort2))

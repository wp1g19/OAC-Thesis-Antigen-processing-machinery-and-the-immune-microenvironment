  #Combined OCCAMS TCGA survival analysis
  
  #Load library
  library("tidyverse")
  library(tidyr)
  library(survival)
  library(survminer)
  library(RTCGA)
  library(forestmodel)
  library(finalfit)
  library(survivalAnalysis)
  #Now I can combine the TMM rna seq data with the combined occams TCGa clinical data to perform survival analysis
  
  #read in the combined TMM ran seq data
  CombinedTMMdata <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGATMMMergedBatchCorrectedAPMCSDE1OutliersCut.csv", header = TRUE,row.names = 1)
  CombinedTMMdata <- as.data.frame(t(CombinedTMMdata))
  rownames(CombinedTMMdata) <- strtrim(rownames(CombinedTMMdata), 12)
  #read in the combined clinical data
  
  clin <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGAMergedclinfinal.csv", header = TRUE, row.names = 1)
  rownames(clin) <- gsub("/", ".", rownames(clin))
  rownames(clin) <- gsub("-", ".", rownames(clin))
  rownames(clin) <- strtrim(rownames(clin), 12)
  #Merge the clinical and TMM data
  TMMclinMerge <- as.data.frame(merge(clin, CombinedTMMdata, by = 0))
  rownames(TMMclinMerge) <- TMMclinMerge$Row.names
  TMMclinMerge$Row.names <- NULL
  
  colnames(TMMclinMerge) <- gsub("-", "_", colnames(TMMclinMerge))
  #optimal cutpoint overall survival analysis
  setwd("C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/Visualisation")
  
  
  #make the optimal cut
  resOS.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = c("RFX5",
                  "CTSS",
                  "CD1D",
                  "MR1",
                  "RFXANK",
                  "SPPL2A",
                  "RFXAP",
                  "CD74",
                  "CIITA",
                  "PSMB10",
                  "LGMN",
                  "CTSL",
                  "TAPBPL",
                  "CALR",
                  "ERAP1",
                  "ERAP2",
                  "CANX",
                  "PDIA3",
                  "B2M",
                  "HLA_E",
                  "HLA_DMB",
                  "HLA_DPA1",
                  "HLA_DPB1",
                  "HLA_DOA",
                  "HLA_DMA",
                  "PSMB9",
                  "HLA_DQA2",
                  "HLA_B",
                  "HLA_DQB2",
                  "HLA_C",
                  "HLA_DRA",
                  "TAP2",
                  "PSMB8",
                  "HLA_DRB5",
                  "HLA_DQA1",
                  "TAP1",
                  "HLA_DRB1",
                  "TAPBP",
                  "HLA_G",
                  "HLA_A",
                  "CSDE1"), minprop = 0.15
  )
  
  
  #make the optimal cuts for CoxPH
  HLA_A.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "HLA_A", minprop = 0.2
  )
  
  HLA_A.surv_count.cat <- surv_categorize(HLA_A.surv_count.cut) 
  
  HLA_B.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "HLA_B", minprop = 0.2
  )
  HLA_B.surv_count.cat <- surv_categorize(HLA_B.surv_count.cut) 
  
  HLA_C.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "HLA_C", minprop = 0.2
  )
  
  HLA_C.surv_count.cat <- surv_categorize(HLA_C.surv_count.cut) 
  
  HLA_E.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "HLA_E", minprop = 0.2
  )
  
  HLA_E.surv_count.cat <- surv_categorize(HLA_E.surv_count.cut) 
  
  RFX5.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "RFX5", minprop = 0.2
  )
  
  RFX5.surv_count.cat <- surv_categorize(RFX5.surv_count.cut) 
  
  CTSS.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "CTSS", minprop = 0.2
  )
  
  CTSS.surv_count.cat <- surv_categorize(CTSS.surv_count.cut) 
  
  CD1D.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "CD1D", minprop = 0.2
  )
  
  CD1D.surv_count.cat <- surv_categorize(CD1D.surv_count.cut) 
  
  MR1.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "MR1", minprop = 0.2
  )
  
  MR1.surv_count.cat <- surv_categorize(MR1.surv_count.cut) 
  
  RFXANK.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "RFXANK", minprop = 0.2
  )
  
  RFXANK.surv_count.cat <- surv_categorize(RFXANK.surv_count.cut) 
  
  SPPL2A.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "SPPL2A", minprop = 0.2
  )
  
  SPPL2A.surv_count.cat <- surv_categorize(SPPL2A.surv_count.cut) 
  
  RFXAP.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "RFXAP", minprop = 0.2
  )
  
  RFXAP.surv_count.cat <- surv_categorize(RFXAP.surv_count.cut) 
  
  CD74.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "CD74", minprop = 0.2
  )
  
  CD74.surv_count.cat <- surv_categorize(CD74.surv_count.cut) 
  
  CIITA.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "CIITA", minprop = 0.2
  )
  
  CIITA.surv_count.cat <- surv_categorize(CIITA.surv_count.cut) 
  
  PSMB10.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "PSMB10", minprop = 0.2
  )
  
  PSMB10.surv_count.cat <- surv_categorize(PSMB10.surv_count.cut) 
  
  LGMN.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "LGMN", minprop = 0.2
  )
  
  LGMN.surv_count.cat <- surv_categorize(LGMN.surv_count.cut) 
  
  CTSL.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "CTSL", minprop = 0.2
  )
  
  CTSL.surv_count.cat <- surv_categorize(CTSL.surv_count.cut) 
  
  TAPBPL.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "TAPBPL", minprop = 0.2
  )
  
  TAPBPL.surv_count.cat <- surv_categorize(TAPBPL.surv_count.cut) 
  
  CALR.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "CALR", minprop = 0.2
  )
  
  CALR.surv_count.cat <- surv_categorize(CALR.surv_count.cut) 
  
  ERAP1.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "ERAP1", minprop = 0.2
  )
  
  ERAP1.surv_count.cat <- surv_categorize(ERAP1.surv_count.cut) 
  
  ERAP2.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "ERAP2", minprop = 0.2
  )
  
  ERAP2.surv_count.cat <- surv_categorize(ERAP2.surv_count.cut) 
  
  CANX.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "CANX", minprop = 0.2
  )
  
  CANX.surv_count.cat <- surv_categorize(CANX.surv_count.cut) 
  
  PDIA3.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "PDIA3", minprop = 0.2
  )
  
  PDIA3.surv_count.cat <- surv_categorize(PDIA3.surv_count.cut) 
  
  B2M.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "B2M", minprop = 0.2
  )
  
  B2M.surv_count.cat <- surv_categorize(B2M.surv_count.cut) 
  
  HLA_DMB.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "HLA_DMB", minprop = 0.2
  )
  
  HLA_DMB.surv_count.cat <- surv_categorize(HLA_DMB.surv_count.cut) 
  
  HLA_DPA1.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "HLA_DPA1", minprop = 0.2
  )
  
  HLA_DPA1.surv_count.cat <- surv_categorize(HLA_DPA1.surv_count.cut) 
  
  
  HLA_DOA.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "HLA_DOA", minprop = 0.2
  )
  
  HLA_DOA.surv_count.cat <- surv_categorize(HLA_DOA.surv_count.cut) 
  
  HLA_DMA.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "HLA_DMA", minprop = 0.2
  )
  
  HLA_DMA.surv_count.cat <- surv_categorize(HLA_DMA.surv_count.cut) 
  
  PSMB9.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "PSMB9", minprop = 0.2
  )
  
  PSMB9.surv_count.cat <- surv_categorize(PSMB9.surv_count.cut) 
  
  
  
  HLA_DQA2.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "HLA_DQA2", minprop = 0.2
  )
  
  HLA_DQA2.surv_count.cat <- surv_categorize(HLA_DQA2.surv_count.cut) 
  
  HLA_DQB2.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "HLA_DQB2", minprop = 0.2
  )
  
  HLA_DQB2.surv_count.cat <- surv_categorize(HLA_DQB2.surv_count.cut) 
  
  HLA_DRA.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "HLA_DRA", minprop = 0.2
  )
  
  HLA_DRA.surv_count.cat <- surv_categorize(HLA_DRA.surv_count.cut) 
  
  TAP2.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "TAP2", minprop = 0.2
  )
  
  TAP2.surv_count.cat <- surv_categorize(TAP2.surv_count.cut) 
  
  PSMB8.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "PSMB8", minprop = 0.2
  )
  
  PSMB8.surv_count.cat <- surv_categorize(PSMB8.surv_count.cut) 
  
  HLA_DRB5.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "HLA_DRB5", minprop = 0.2
  )
  
  HLA_DRB5.surv_count.cat <- surv_categorize(HLA_DRB5.surv_count.cut) 
  
  HLA_DQA1.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "HLA_DQA1", minprop = 0.2
  )
  
  HLA_DQA1.surv_count.cat <- surv_categorize(HLA_DQA1.surv_count.cut) 
  
  TAP1.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "TAP1", minprop = 0.2
  )
  
  TAP1.surv_count.cat <- surv_categorize(TAP1.surv_count.cut) 
  
  
  TAPBP.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "TAPBP", minprop = 0.2
  )
  
  TAPBP.surv_count.cat <- surv_categorize(TAPBP.surv_count.cut) 
  
  HLA_DRB1.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "HLA_DRB1", minprop = 0.2
  )
  
  HLA_DRB1.surv_count.cat <- surv_categorize(HLA_DRB1.surv_count.cut) 
  
  HLA_G.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "HLA_G", minprop = 0.2
  )
  
  HLA_G.surv_count.cat <- surv_categorize(HLA_G.surv_count.cut) 
  
  CSDE1.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "CSDE1", minprop = 0.2
  )
  
  CSDE1.surv_count.cat <- surv_categorize(CSDE1.surv_count.cut) 
  
  PTPN2.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "PTPN2", minprop = 0.2
  )
  
  PTPN2.surv_count.cat <- surv_categorize(PTPN2.surv_count.cut) 
  
  SMYD3.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "SMYD3", minprop = 0.2
  )
  
  SMYD3.surv_count.cat <- surv_categorize(SMYD3.surv_count.cut) 
  
  NLRC5.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "NLRC5", minprop = 0.1
  )
  
  NLRC5.surv_count.cat <- surv_categorize(NLRC5.surv_count.cut) 
  
  
  IRF1.surv_count.cut <- surv_cutpoint(
    TMMclinMerge,
    time = "OS.months",
    event = "Vital_status",
    variables = "IRF1", minprop = 0.1
  )
  
  IRF1.surv_count.cat <- surv_categorize(IRF1.surv_count.cut) 
  
  OptimalOSoptimalcut <- do.call("cbind", list(HLA_A.surv_count.cat,
                        HLA_B.surv_count.cat,
                        HLA_C.surv_count.cat,
                        HLA_C.surv_count.cat,
                        HLA_E.surv_count.cat,
                        HLA_G.surv_count.cat,
                        HLA_DMA.surv_count.cat,
                        HLA_DMB.surv_count.cat,
                        HLA_DOA.surv_count.cat,
                        HLA_DPA1.surv_count.cat,
                        HLA_DQA1.surv_count.cat,
                        HLA_DQA2.surv_count.cat,
                        HLA_DQA2.surv_count.cat,
                        HLA_DQB2.surv_count.cat,
                        HLA_DRA.surv_count.cat,
                        HLA_DRB1.surv_count.cat,
                        HLA_DRB5.surv_count.cat,
                        LGMN.surv_count.cat,
                        MR1.surv_count.cat,
                        PDIA3.surv_count.cat,
                        PSMB10.surv_count.cat,
                        PSMB8.surv_count.cat,
                        PSMB9.surv_count.cat,
                        RFX5.surv_count.cat,
                        RFXANK.surv_count.cat,
                        RFXAP.surv_count.cat,
                        SPPL2A.surv_count.cat,
                        TAP1.surv_count.cat,
                        TAP2.surv_count.cat,
                        TAPBPL.surv_count.cat,
                        TAPBP.surv_count.cat,
                        ERAP1.surv_count.cat,
                        ERAP2.surv_count.cat,
                        CTSL.surv_count.cat,
                        CTSS.surv_count.cat,
                        CSDE1.surv_count.cat,
                        CIITA.surv_count.cat,
                        CD74.surv_count.cat,
                        CD1D.surv_count.cat,
                        CD74.surv_count.cat,
                        CANX.surv_count.cat,
                        CALR.surv_count.cat,
                        B2M.surv_count.cat,
                        PTPN2.surv_count.cat,
                        SMYD3.surv_count.cat,
                        NLRC5.surv_count.cat,
                        IRF1.surv_count.cat))
  
  
  OptimalOSoptimalcut2 <- OptimalOSoptimalcut[, !duplicated(colnames(OptimalOSoptimalcut), fromLast = FALSE)] 
  
  write.csv(OptimalOSoptimalcut2, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/CombinedOSoptimalcut.csv")
  OptimalOSoptimalcut2 <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/CombinedOSoptimalcut.csv", header = TRUE, row.names = 1)
  OptimalOSoptimalcut2$Vital_status
  covariates <- c(   "RFX5",
                     "CTSS",
                     "CD1D",
                     "MR1",
                     "RFXANK",
                     "SPPL2A",
                     "RFXAP",
                     "CD74",
                     "CIITA",
                     "PSMB10",
                     "LGMN",
                     "CTSL",
                     "TAPBPL",
                     "CALR",
                     "ERAP1",
                     "ERAP2",
                     "CANX",
                     "PDIA3",
                     "B2M",
                     "HLA_E",
                     "HLA_DMB",
                     "HLA_DPA1",
                     "HLA_DOA",
                     "HLA_DMA",
                     "PSMB9",
                     "HLA_DQA2",
                     "HLA_B",
                     "HLA_DQB2",
                     "HLA_C",
                     "HLA_DRA",
                     "TAP2",
                     "PSMB8",
                     "HLA_DRB5",
                     "HLA_DQA1",
                     "TAP1",
                     "HLA_DRB1",
                     "TAPBP",
                     "HLA_G",
                     "HLA_A",
                     "CSDE1",
                     "PTPN2",
                     "SMYD3",
                     "NLRC5",
                     "IRF1")
  
  #Univariate clinical model
  univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(OS.months, Vital_status)~', x)))
  
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = OptimalOSoptimalcut2, na.action = na.omit)})
  
  forest_model(model_list = univ_models,
               covariates = covariates,
               merge_models = T,
               recalculate_height = 10,
               recalculate_width = 3)
  
  covariates <- c("HLA_A",
                "CSDE1")
  OptimalOSoptimalcut2$CSDE1 <- as.factor(OptimalOSoptimalcut2$CSDE1)
  OptimalOSoptimalcut2$CSDE1 = relevel(OptimalOSoptimalcut2$CSDE1, ref = "low")
  OptimalOSoptimalcut2$SPPL2A <- as.factor(OptimalOSoptimalcut2$SPPL2A)
  OptimalOSoptimalcut2$SPPL2A = relevel(OptimalOSoptimalcut2$SPPL2A, ref = "low")
  OptimalOSoptimalcut2$RFXAP <- as.factor(OptimalOSoptimalcut2$RFXAP)
  OptimalOSoptimalcut2$RFXAP = relevel(OptimalOSoptimalcut2$RFXAP, ref = "low")
  OptimalOSoptimalcut2$PSMB10 <- as.factor(OptimalOSoptimalcut2$PSMB10)
  OptimalOSoptimalcut2$PSMB10 = relevel(OptimalOSoptimalcut2$PSMB10, ref = "low")
  
  #Univariate clinical model
  univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(OS.months, Vital_status)~', x)))
  
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = OptimalOSoptimalcut2, na.action = na.omit)})
  
  
  
  forest_model(model_list = univ_models,
               covariates = covariates,
               merge_models = T,
               recalculate_height = 10,
               recalculate_width = 6, )
  #select the genes with univariate survival significance
  siggenes <- c("RFX5",
                "CTSS",
                "CD1D",
                "MR1",
                "SPPL2A",
                "CD74",
                "CIITA",
                "PSMB10",
                "LGMN",
                "CTSL",
                "TAPBPL",
                "CALR",
                "ERAP1",
                "ERAP2",
                "HLA_E",
                "HLA_DPA1",
                "PSMB9",
                "HLA_B",
                "HLA_DRA",
                "PSMB8",
                "HLA_DRB5",
                "HLA_DQA1",
                "HLA_DRB1",
                "TAPBP",
                "HLA_G",
                "HLA_A",
                "CSDE1")
  
  #select univariate significant genes for a multivariate analysis
  OptimalOSoptimalcutsig <- OptimalOSoptimalcut2[,colnames(OptimalOSoptimalcut2) %in% siggenes ]
  
  #merge significant genes with clinical features
  optsiggene <- OptimalOSoptimalcutsig
  optsiggene$ID <- rownames(optsiggene)
  write.csv(optsiggene, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/optimalcutpointsOS.csv")
  clin2 <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGAMergedclinfinal.csv", header = TRUE, row.names = 1)
  rownames(clin2) <- gsub("/",".", rownames(clin2))
  clin2$ID <- rownames(clin2)
  
  OptimalOSoptimalcutsigclin <- left_join(clin2, optsiggene,"ID")
  colnames(OptimalOSoptimalcutsigclin) <- gsub(".x","", colnames(OptimalOSoptimalcutsigclin))
  colnames(OptimalOSoptimalcutsigclin) <- gsub(".y","", colnames(OptimalOSoptimalcutsigclin))
  OptimalOSoptimalcutsigclin2 <- OptimalOSoptimalcutsigclin[, !duplicated(colnames(OptimalOSoptimalcutsigclin), fromLast = FALSE)] 
  OptimalOSoptimalcutsigclin <- OptimalOSoptimalcutsigclin2
  colnames(OptimalOSoptimalcutsigclin)[colnames(OptimalOSoptimalcutsigclin) == 'S'] <- 'Sex'
  library(finalfit)
  dependent_os  <- "Surv(OS.months, Vital_status)"
  OptimalOSoptimalcutsigclin$OS.months
  explanatory   <- c(     "Age", "Sex","pT",
                          "pN",
                          "pM",
                          "HLA_A",
                          "CSDE1")
  
  OptimalOSoptimalcutsigclin %>% 
    finalfit(dependent_os, explanatory)
  
  OptimalOSoptimalcutsigclin %>% 
    finalfit(dependent_os, explanatory, add_dependent_label = FALSE, metrics = TRUE)
  OptimalOSoptimalcutsigclin %>% 
  hr_plot(dependent_os, explanatory)


csde1mult <- coxph(Surv(OS.months, Vital_status)~ Age + Sex + pT + pN + pM + CSDE1, data = OptimalOSoptimalcutsigclin, na.action = na.omit)




#make the optimal cut
HLA_A.surv_count.cut <- surv_cutpoint(
  TMMclinMerge,
  time = "OS.months",
  event = "Vital_status",
  variables = "HLA_A", minprop = 0.2
)

HLA_A.surv_count.cut.sum <- as.data.frame(summary(HLA_A.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_A.surv_count.cut, "HLA_A", palette = "npg", bins = 176, xlim=c(250, 1000))
dev.copy(png,filename="HLA_A.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_A.surv_count.cat <- surv_categorize(HLA_A.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Vital_status) ~ HLA_A,
               data = HLA_A.surv_count.cat)
# extract p value
HLA_Apvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,40), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 10,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE, surv.median.line = "hv" # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_A_OSKm.png");
while (dev.cur()>1) dev.off()

#make the optimal cut
HLA_A.surv_count.cut <- surv_cutpoint(
  TMMclinMerge,
  time = "OS.months",
  event = "Vital_status",
  variables = "HLA_A", minprop = 0.2
)

HLA_A.surv_count.cut.sum <- as.data.frame(summary(HLA_A.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_A.surv_count.cut, "HLA_A", palette = "npg", bins = 176, xlim=c(250, 1000))
dev.copy(png,filename="HLA_A.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_A.surv_count.cat <- surv_categorize(HLA_A.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Vital_status) ~ HLA_A,
               data = HLA_A.surv_count.cat)
# extract p value
HLA_Apvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,40), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 10,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE, surv.median.line = "hv" # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_A_OSKm.png");
while (dev.cur()>1) dev.off()

#make the optimal cut
HLA_B.surv_count.cut <- surv_cutpoint(
  TMMclinMerge,
  time = "OS.months",
  event = "Vital_status",
  variables = "HLA_B", minprop = 0.15
)

HLA_B.surv_count.cut.sum <- as.data.frame(summary(HLA_B.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_B.surv_count.cut, "HLA_B", palette = "npg", bins = 176, xlim=c(250, 1000))
dev.copy(png,filename="HLA_B.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_B.surv_count.cat <- surv_categorize(HLA_B.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Vital_status) ~ HLA_B,
               data = HLA_B.surv_count.cat)
# extract p value
HLA_Bpvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,80), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_B_OSKm.png");
while (dev.cur()>1) dev.off()

#make the optimal cut
HLA_C.surv_count.cut <- surv_cutpoint(
  TMMclinMerge,
  time = "OS.months",
  event = "Vital_status",
  variables = "HLA_C", minprop = 0.2
)

HLA_C.surv_count.cut.sum <- as.data.frame(summary(HLA_C.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_C.surv_count.cut, "HLA_C", palette = "npg", bins = 176, xlim=c(250, 1000))
dev.copy(png,filename="HLA_C.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_C.surv_count.cat <- surv_categorize(HLA_C.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Vital_status) ~ HLA_C,
               data = HLA_C.surv_count.cat)
# extract p value
HLA_Cpvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,80), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_C_OSKm.png");
while (dev.cur()>1) dev.off()

#make the optimal cut
HLA_E.surv_count.cut <- surv_cutpoint(
  TMMclinMerge,
  time = "OS.months",
  event = "Vital_status",
  variables = "HLA_E"
)

HLA_E.surv_count.cut.sum <- as.data.frame(summary(HLA_E.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_E.surv_count.cut, "HLA_E", palette = "npg", bins = 176, xlim=c(250, 1000))
dev.copy(png,filename="HLA_E.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_E.surv_count.cat <- surv_categorize(HLA_E.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Vital_status) ~ HLA_E,
               data = HLA_E.surv_count.cat)
# extract p value
HLA_Epvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,80), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_E_OSKm.png");
while (dev.cur()>1) dev.off()

#make the optimal cut
HLA_E.surv_count.cut <- surv_cutpoint(
  TMMclinMerge,
  time = "OS.months",
  event = "Vital_status",
  variables = "HLA_E"
)

HLA_E.surv_count.cut.sum <- as.data.frame(summary(HLA_E.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(HLA_E.surv_count.cut, "HLA_E", palette = "npg", bins = 176, xlim=c(250, 1000))
dev.copy(png,filename="HLA_E.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
HLA_E.surv_count.cat <- surv_categorize(HLA_E.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Vital_status) ~ HLA_E,
               data = HLA_E.surv_count.cat)
# extract p value
HLA_Epvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,80), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="HLA_E_OSKm.png");
while (dev.cur()>1) dev.off()

#make the optimal cut
CSDE1.surv_count.cut <- surv_cutpoint(
  TMMclinMerge,
  time = "OS.months",
  event = "Vital_status",
  variables = "CSDE1", minprop = 0.15
)

CSDE1.surv_count.cut.sum <- as.data.frame(summary(CSDE1.surv_count.cut))

#Plot the cutpoint 

# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(CSDE1.surv_count.cut, "CSDE1", palette = "npg", bins = 176, xlim=c(250, 1000))
dev.copy(png,filename="CSDE1.OS.surv_count.cut.png");
while (dev.cur()>1) dev.off()

#categorise the cutpoint 
CSDE1.surv_count.cat <- surv_categorize(CSDE1.surv_count.cut) 

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Vital_status) ~ CSDE1,
               data = CSDE1.surv_count.cat)
# extract p value
CSDE1pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,60), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(),  font.x = 25, font.y = 25, font.legend = 25, font.caption = 25,  font.tickslab = 20, pval.size = 10, # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE, surv.median.line = "hv" # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="CSDE1_OSKm.png");
while (dev.cur()>1) dev.off()


#Cancer specific survival

TMMclinMergeCSS <- subset(TMMclinMerge, Bioinformatic.CSS_status!="2" & Bioinformatic.CSS_status!="#N/A")

TMMclinMergeCSS <-TMMclinMergeCSS %>% drop_na(Bioinformatic.CSS_status)
TMMclinMergeCSS$Bioinformatic.CSS_status <- as.numeric(TMMclinMergeCSS$Bioinformatic.CSS_status)
resCSS.surv_count.cut <- surv_cutpoint(
  TMMclinMergeCSS,
  time = "OS.months",
  event = "Bioinformatic.CSS_status",
  variables = c("RFX5",
                "CTSS",
                "CD1D",
                "MR1",
                "RFXANK",
                "SPPL2A",
                "RFXAP",
                "CD74",
                "CIITA",
                "PSMB10",
                "LGMN",
                "CTSL",
                "TAPBPL",
                "CALR",
                "ERAP1",
                "ERAP2",
                "CANX",
                "PDIA3",
                "B2M",
                "HLA_E",
                "HLA_DMB",
                "HLA_DPA1",
                "HLA_DOA",
                "HLA_DMA",
                "PSMB9",
                "HLA_DQA2",
                "HLA_B",
                "HLA_DQB2",
                "HLA_C",
                "HLA_DRA",
                "TAP2",
                "PSMB8",
                "HLA_DRB5",
                "HLA_DQA1",
                "TAP1",
                "HLA_DRB1",
                "TAPBP",
                "HLA_G",
                "HLA_A",
                "CSDE1",
                "PTPN2",
                "SMYD3",
                "NLRC5",
                "IRF1"), minprop = 0.10
)

res.cat <- surv_categorize(resCSS.surv_count.cut)
head(res.cat)
write.csv(res.cat, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/optimalcutpointsCSS.csv")
cssclin <- TMMclinMergeCSS[,c("Age", "Sex", "pN", "pT", "pM")]
res.cat <- merge(res.cat, cssclin, by = 0)
rownames(res.cat) <- res.cat$Row.names
res.cat$Row.names <- NULL
res.cat$Bioinformatic.DFS_status <- as.numeric(res.cat$Bioinformatic.CSS_status)

library(survival)
library(RTCGA)
fit <- survfit(Surv(OS.months, Bioinformatic.CSS_status) ~ HLA_E,
               data = res.cat)
# extract p value
CSDE1pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,60), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE, font.x = 17, font.y = 17, font.legend = 17, font.caption = 12,  font.tickslab = 12,
  surv.median.line = "hv",  # show bars instead of names in text annotations
  # in legend of risk table
)  


#save the image file
dev.copy(png,filename="CSDE1_CSSKm.png");
while (dev.cur()>1) dev.off()



covariates <- c(   "RFX5",
                   "CTSS",
                   "CD1D",
                   "MR1",
                   "RFXANK",
                   "SPPL2A",
                   "RFXAP",
                   "CD74",
                   "CIITA",
                   "PSMB10",
                   "LGMN",
                   "CTSL",
                   "TAPBPL",
                   "CALR",
                   "ERAP1",
                   "ERAP2",
                   "CANX",
                   "PDIA3",
                   "B2M",
                   "HLA_E",
                   "HLA_DMB",
                   "HLA_DPA1",
                   "HLA_DOA",
                   "HLA_DMA",
                   "PSMB9",
                   "HLA_DQA2",
                   "HLA_B",
                   "HLA_DQB2",
                   "HLA_C",
                   "HLA_DRA",
                   "TAP2",
                   "PSMB8",
                   "HLA_DRB5",
                   "HLA_DQA1",
                   "TAP1",
                   "HLA_DRB1",
                   "TAPBP",
                   "HLA_G",
                   "HLA_A",
                   "CSDE1",
                   "PTPN2",
                   "SMYD3",
                   "NLRC5",
                   "IRF1")

#Univariate clinical model
univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(OS.months, Bioinformatic.CSS_status)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = res.cat, na.action = na.omit)})

forest_model(model_list = univ_models,
             covariates = covariates,
             merge_models = T,
             recalculate_height = 10,
             recalculate_width = 3)


covariates <- c(   "RFX5",
                   "RFXAP",
                   "CTSS",
                   "MR1",
                   "SPPL2A",
                   "CD74",
                   "CIITA",
                   "PSMB10",
                   "LGMN",
                   "CTSL",
                   "TAPBPL",
                   "CALR",
                   "ERAP1",
                   "ERAP2",
                   "HLA_E",
                   "HLA_DPA1",
                   "PSMB9",
                   "HLA_B",
                   "HLA_DRA",
                   "PSMB8",
                   "HLA_DRB5",
                   "HLA_DQA1",
                   "HLA_DRB1",
                   "TAPBP",
                   "HLA_G",
                   "HLA_A",
                   "CSDE1")
res.cat$RFXAP <- as.factor(res.cat$RFXAP)
res.cat$RFXAP = relevel(res.cat$RFXAP, ref = "low")

#Univariate clinical model
univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(OS.months, Bioinformatic.CSS_status)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = res.cat, na.action = na.omit)})

forest_model(model_list = univ_models,
             covariates = covariates,
             merge_models = T,
             recalculate_height = 10,
             recalculate_width = 3)
#multivariate 
dependent_os  <- "Surv(OS.months, Bioinformatic.CSS_status)"

explanatory   <- c(    "CSDE1",
                       "HLA-E" ,
                       "SPPL2A",
                       "HLA-DPA1", 
                       "PSMB9",
                       "HLA-B" ,
                       "HLA-DRA", 
                       "PSMB8" ,
                       "HLA-DRB5",
                       "HLA-DQA1",
                       "HLA-DRB1",
                       "TAPBP" ,
                       "HLA-G",
                       "HLA-A",
                       "ERAP2",
                       "ERAP1",
                       "CALR"  ,    
                       "TAPBPL" ,     
                       "CTSL"    ,  
                       "LGMN"     , 
                       "PSMB10"    ,  
                       "CIITA"  ,
                       "CD74"  ,
                       "MR1" ,
                       "CTSS" ,     
                       "RFX5"   )

res.cat %>% 
  finalfit(dependent_os, explanatory)

res.cat %>% 
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE, metrics = TRUE)




TMMclinMergeDFS <- subset(TMMclinMerge, Bioinformatic.DFS_status!="#N/A")

TMMclinMergeDFS <-TMMclinMergeDFS %>% drop_na(Bioinformatic.DFS_status)
TMMclinMergeCSS$Bioinformatic.DFS_status <- as.numeric(TMMclinMergeCSS$Bioinformatic.DFS_status)
resDFS.surv_count.cut <- surv_cutpoint(
  TMMclinMergeCSS,
  time = "Bioinformatic.DFS_months",
  event = "Bioinformatic.DFS_status",
  variables = c(   "RFX5",
                   "CTSS",
                   "CD1D",
                   "MR1",
                   "RFXANK",
                   "SPPL2A",
                   "RFXAP",
                   "CD74",
                   "CIITA",
                   "PSMB10",
                   "LGMN",
                   "CTSL",
                   "TAPBPL",
                   "CALR",
                   "ERAP1",
                   "ERAP2",
                   "CANX",
                   "PDIA3",
                   "B2M",
                   "HLA_E",
                   "HLA_DMB",
                   "HLA_DPA1",
                   "HLA_DOA",
                   "HLA_DMA",
                   "PSMB9",
                   "HLA_DQA2",
                   "HLA_B",
                   "HLA_DQB2",
                   "HLA_C",
                   "HLA_DRA",
                   "TAP2",
                   "PSMB8",
                   "HLA_DRB5",
                   "HLA_DQA1",
                   "TAP1",
                   "HLA_DRB1",
                   "TAPBP",
                   "HLA_G",
                   "HLA_A",
                   "CSDE1",
                   "PTPN2",
                   "SMYD3",
                   "NLRC5",
                   "IRF1"), minprop = 0.10
)

resDFS.cat <- surv_categorize(resDFS.surv_count.cut)
head(resDFS.cat)
write.csv(resDFS.cat, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/optimalcutpointsDFS.csv")
resDFS.cat <- merge(resDFS.cat, cssclin, by = 0)
library(survival)
library(RTCGA)
fit <- survfit(Surv(Bioinformatic.DFS_months, Bioinformatic.DFS_status) ~ CSDE1,
               data = resDFS.cat)



covariates <- c(   "RFX5",
                   "CTSS",
                   "CD1D",
                   "MR1",
                   "RFXANK",
                   "SPPL2A",
                   "RFXAP",
                   "CD74",
                   "CIITA",
                   "PSMB10",
                   "LGMN",
                   "CTSL",
                   "TAPBPL",
                   "CALR",
                   "ERAP1",
                   "ERAP2",
                   "CANX",
                   "PDIA3",
                   "B2M",
                   "HLA_E",
                   "HLA_DMB",
                   "HLA_DPA1",
                   "HLA_DOA",
                   "HLA_DMA",
                   "PSMB9",
                   "HLA_DQA2",
                   "HLA_B",
                   "HLA_DQB2",
                   "HLA_C",
                   "HLA_DRA",
                   "TAP2",
                   "PSMB8",
                   "HLA_DRB5",
                   "HLA_DQA1",
                   "TAP1",
                   "HLA_DRB1",
                   "TAPBP",
                   "HLA_G",
                   "HLA_A",
                   "CSDE1",
                   "PTPN2",
                   "SMYD3",
                   "NLRC5",
                   "IRF1")

#Univariate clinical model
univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(Bioinformatic.DFS_months, Bioinformatic.DFS_status)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = resDFS.cat, na.action = na.omit)})

forest_model(model_list = univ_models,
             covariates = covariates,
             merge_models = T,
             recalculate_height = 10,
             recalculate_width = 3)


covariates <- c(  )
res.cat$CD1D <- as.factor(res.cat$CD1D)
resDFS.cat$PSMB10 = relevel(resDFS.cat$PSMB10, ref = "low")
#Univariate clinical model
univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(Bioinformatic.DFS_months, Bioinformatic.DFS_status)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = resDFS.cat, na.action = na.omit)})

forest_model(model_list = univ_models,
             covariates = covariates,
             merge_models = T,
             recalculate_height = 10,
             recalculate_width = 3)



#multivariate 
dependent_os  <- "Surv(Bioinformatic.DFS_months, Bioinformatic.DFS_status)"


explanatory   <- c(    "RFX5",
                       "SPPL2A",
                       "pT",
                       "pN",
                       "pM",
                       "CSDE1")
 


resDFS.cat %>% 
  finalfit(dependent_os, explanatory)

resDFS.cat %>% 
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE, metrics = TRUE)









fit <- survfit(Surv(OS.months, Bioinformatic.CSS_status) ~ CSDE1,
               data = res.cat)
# extract p value
CSDE1pvalue <- surv_pvalue(fit)$pval.txt
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = FALSE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,60), 
  xlab = ("Months"),# present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 20,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE,  font.x = 25, font.y = 25, font.legend = 25, font.caption = 25,  font.tickslab = 20, pval.size = 10,
  surv.median.line = "hv", 
  # show bars instead of names in text annotations
  # in legend of risk table
)  
#save the image file
dev.copy(png,filename="CSDE1_DFSKm.png");
while (dev.cur()>1) dev.off()









# Here I analyse the correlation data of the combined batch corrected TCGA and OCCAMS datasets for clinical outcomes

#read in the TMM combined batch corrected data 
CombinedTMMdata <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGATMMMergedBatchCorrectedAPMCSDE1OutliersCut.csv", header = TRUE,row.names = 1)

#Now create a clinical dataframe of the OCCAMS and TCGA data

library(TCGAbiolinks)
library(dplyr)

ESCAclin <- GDCquery_clinic(project = "TCGA-ESCA", save.csv = TRUE, type = "clinical")
ESCAclin$submitter_id <- gsub("\\-", ".", ESCAclin$submitter_id)
#Filter the Clinical data to only our IDs of interest
sampleids <- colnames(CombinedTMMdata)

ESCAclinfilt <- as.data.frame(ESCAclin[ESCAclin$submitter_id %in% sampleids, ])
#query clinical supplement for more data
ESCAclinquery <- GDCquery(
  project = "TCGA-ESCA", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab"
)
GDCdownload(ESCAclinquery)
clinical.BCRtab.all <- GDCprepare(ESCAclinquery)
names(clinical.BCRtab.all)
nte <- clinical.BCRtab.all$clinical_nte_esca
patientsup <- clinical.BCRtab.all$clinical_patient_esca
patientfollowupnte <- clinical.BCRtab.all$clinical_follow_up_v4.0_nte_esca
patientfollowup <- clinical.BCRtab.all$clinical_follow_up_v4.0_esca
ESCAclinsupplement <- Reduce(function(x, y) merge(x, y, all=TRUE), clinical.BCRtab.all)
cols <- as.data.frame(colnames(ESCAclinsupplement))



#Select the clinical features we want to assess in analysis

Clinfeatures <- c("submitter_id",
                  "days_to_diagnosis",
                  "days_to_last_follow_up",
                  "days_to_death",
                  "days_to_birth",
                  "days_to_recurrence",
                  "gender",
                  "vital_status",
                  "ajcc_pathologic_m",
                  "ajcc_pathologic_n",
                  "ajcc_pathologic_t"
                  )

ESCAclinfilt2 <-  as.data.frame(ESCAclinfilt[,colnames(ESCAclinfilt) %in% Clinfeatures])
ESCAclinfilt3 <- as.data.frame(ESCAclinfilt2)
ESCAclinfilt3 <- ESCAclinfilt3 %>% 
  mutate(OS.Days = coalesce(days_to_death, days_to_last_follow_up))
ESCAclinfilt3 <- ESCAclinfilt3 %>% 
  mutate(OS.months = round(ESCAclinfilt3$OS.Days/30.417, digit=0))

ESCAclinfilt3$vital_status <- ifelse(ESCAclinfilt3$vital_status == "Dead", yes = 1, no = 0)

ESCAclinfilt3$days_to_birth <- as.numeric(gsub("\\-", "", ESCAclinfilt3$days_to_birth))

ESCAclinfilt3 <- ESCAclinfilt3 %>% 
  mutate(Age = round(ESCAclinfilt3$days_to_birth/365.242199, digit=0))

#read in the TCGA paper data for DSS and PFS
ESCApaperclin <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/ESCA tcga Paper clin.csv", header = TRUE)
rownames(ESCApaperclin) <- ESCApaperclin$submitter_id
ESCApaperclinfilt <- as.data.frame(ESCApaperclin[rownames(ESCApaperclin) %in% sampleids, ])

ESCApaperMergeGDCclinfilt <- left_join(x = ESCApaperclinfilt,y =  ESCAclinfilt3)
#Read in the OCCAMS data

OCCAMSclin <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMS clinical.csv", header = TRUE)
OCCAMSclin2 <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMS.clinical.csv", header = TRUE)
#filter occams clinical to samples of interest
row.names(OCCAMSclin) <- OCCAMSclin$ID
rownames(OCCAMSclin) <- gsub("/", ".", rownames(OCCAMSclin))
OCCAMSclinfilt <- as.data.frame(OCCAMSclin[rownames(OCCAMSclin) %in% sampleids, ])
row.names(OCCAMSclin2) <- OCCAMSclin2$ID
rownames(OCCAMSclin2) <- gsub("/", ".", rownames(OCCAMSclin2))
sampleids <- rownames(OCCAMSclinfilt)
OCCAMSclinfilt2 <- as.data.frame(OCCAMSclin2[rownames(OCCAMSclin2) %in% sampleids, ])
OCCAMSclinfilt3 <- as.data.frame(OCCAMSclinfilt2[rownames(OCCAMSclinfilt2) %in% rownames(OCCAMSclinfilt), ])
OCCAMSclinfilt4 <- as.data.frame(OCCAMSclinfilt)
OCCAMSclinfilt4$DI.ageAtDiagnosis <- OCCAMSclinfilt3$DI.ageAtDiagnosis
OCCAMSclinfilt4$RP.MStage.DistantMetastasis

Clinfeatures2 <- c("ID",
                  "OS.days",
                  "DI.PatientGender",
                  "FE.Patient.Died",
                  "RP.Nstage.RP.TNM7",
                  "RP.TStage.PrimaryTumour",
                  "RP.MStage.DistantMetastasis",
                  "Bioinformatic.DFS_status",
                  "Bioinformatic.CSS_status",
                  "Bioinformatic.DFS_days",
                  "DI.ageAtDiagnosis"
)

OCCAMSclinfilt5 <-  as.data.frame(OCCAMSclinfilt4[,colnames(OCCAMSclinfilt4) %in% Clinfeatures2])


OCCAMSclinfilt5 <- OCCAMSclinfilt5 %>% 
  mutate(OS.months = round(OCCAMSclinfilt5$OS.days/30.417, digit=0))


# Now I must remold the order and names of the clinical parameters in each dataframe to produce a single rbinded dataframe.

colnames(ESCApaperMergeGDCclinfilt)
colnames(OCCAMSclinfilt5)
colnames(ESCApaperMergeGDCclinfilt) <- c("ID",
                                         "type",
                                         "Bioinformatic.DFS_status",
                                         "Bioinformatic.DFS_days",
                                         "Bioinformatic.CSS_status",
                                         "Bioinformatic.CSS_days",
                                         "days_to_diagnosis",
                                         "days_to_last_follow_up",
                                         "pT",
                                         "days_to_recurrence",
                                         "pN",
                                         "pM",
                                         "Sex",
                                         "Vital_status",
                                         "days_to_birth",
                                         "days_to_death",
                                         "OS.Days",
                                         "OS.months",
                                         "Age")
colnames(OCCAMSclinfilt5) <- c("ID",
                               "Vital_status",
                               "Age",
                               "Sex",
                               "pN",
                               "pM",
                               "pT",
                               "OS.days",
                               "Bioinformatic.CSS_status",
                               "Bioinformatic.DFS_status",
                               "Bioinformatic.DFS_days",
                               "OS.months")


OCCAMSclinfilt5 <- OCCAMSclinfilt5 %>% 
  mutate(Bioinformatic.DFS_months = round(OCCAMSclinfilt5$Bioinformatic.DFS_days/30.417, digit=0))

ESCApaperMergeGDCclinfilt <- ESCApaperMergeGDCclinfilt %>% 
  mutate(Bioinformatic.DFS_months = round(ESCApaperMergeGDCclinfilt$Bioinformatic.DFS_days/30.417, digit=0))


Clinfeatures3 <- c("ID",
                   "Vital_status",
                   "Age",
                   "Sex",
                   "pN",
                   "pM",
                   "pT",
                   "Bioinformatic.CSS_status",
                   "Bioinformatic.DFS_status",
                   "Bioinformatic.DFS_months",
                   "OS.months")
OCCAMSclinfiltfinal <-  as.data.frame(OCCAMSclinfilt5[,colnames(OCCAMSclinfilt5) %in% Clinfeatures3])
ESCAclinfinal <-  as.data.frame(ESCApaperMergeGDCclinfilt[,colnames(ESCApaperMergeGDCclinfilt) %in% Clinfeatures3])

#write out the final clinical dataframes and correct in excel before reading in and merging

write.csv(OCCAMSclinfiltfinal, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSclinfiltfinal.csv")
write.csv(ESCAclinfinal, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/ESCAclinfinal.csv")

OCCAMSclinfiltfinal <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSclinfiltfinal.csv", header = TRUE, row.names = 1)
ESCAclinfinal <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/ESCAclinfinal.csv", header = TRUE, row.names = 1)

sampleids <- colnames(CombinedTMMdata)

OCCAMSTCGAClinMergedFinal <- rbind(OCCAMSclinfiltfinal, ESCAclinfinal)

#write out the final merged clinical dataframe
write.csv(OCCAMSTCGAClinMergedFinal, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGAMergedclinfinal.csv")

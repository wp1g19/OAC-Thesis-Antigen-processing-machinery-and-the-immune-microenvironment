library(TCGAbiolinks)
library(dplyr)
library(maftools)
library(readr)

#Download the TCGA mutation data and place into a maf
tcgaCNquery <- GDCquery(project = "TCGA-ESCA", 
                         legacy = FALSE,
                         data.category = "Copy Number Variation",
                         data.type = "Copy Number Segment",
                         barcode =  c("TCGA-L7-A6VZ", "TCGA-IG-A4QS", "TCGA-L5-A8NV", "TCGA-R6-A6XG", "TCGA-R6-A6DQ", "TCGA-L5-A8NE", "TCGA-L5-A8NT", "TCGA-L5-A4OQ", "TCGA-2H-A9GH", "TCGA-R6-A8W8", "TCGA-L5-A8NJ", "TCGA-V5-AASX", "TCGA-L5-A4OP", "TCGA-L5-A8NF", "TCGA-L5-A8NM", "TCGA-L5-A8NR", "TCGA-2H-A9GN", "TCGA-L5-A8NS", "TCGA-L5-A4OI", "TCGA-2H-A9GO", "TCGA-R6-A8W5", "TCGA-L5-A88V", "TCGA-L5-A8NI", "TCGA-L5-A4OW", "TCGA-L5-A4OJ", "TCGA-L5-A8NW", "TCGA-R6-A6L6", "TCGA-L5-A891", "TCGA-RE-A7BO", "TCGA-IC-A6RE", "TCGA-JY-A6F8", "TCGA-JY-A93E", "TCGA-VR-AA4D", "TCGA-R6-A6KZ", "TCGA-S8-A6BV", "TCGA-X8-AAAR", "TCGA-L5-A4OO", "TCGA-JY-A6FB", "TCGA-JY-A93C", "TCGA-2H-A9GK", "TCGA-2H-A9GF", "TCGA-2H-A9GI", "TCGA-2H-A9GL", "TCGA-2H-A9GM", "TCGA-2H-A9GR","TCGA-JY-A6FH", "TCGA-JY-A93D", "TCGA-JY-A938", "TCGA-L5-A4OE", "TCGA-L5-A4OF", "TCGA-L5-A4OG","TCGA-L5-A4OH", "TCGA-L5-A4ON", "TCGA-L5-A4OR", "TCGA-L5-A4OS", "TCGA-L5-A4OT", "TCGA-L5-A4OU", "TCGA-L5-A4OX", "TCGA-L5-A8NG", "TCGA-L5-A8NH", "TCGA-L5-A8NL", "TCGA-L5-A8NN", "TCGA-L5-A8NU", "TCGA-L5-A43C", "TCGA-L5-A43I", "TCGA-L5-A88Y", "TCGA-M9-A5M8", "TCGA-Q9-A6FW", "TCGA-R6-A6DN", "TCGA-R6-A6L4", "TCGA-R6-A6XQ", "TCGA-R6-A6Y2", "TCGA-R6-A8WC", "TCGA-R6-A8WG", "TCGA-V5-A7RB", "TCGA-V5-A7RE", "TCGA-V5-AASW", "TCGA-VR-A8EQ", "TCGA-ZR-A9CJ"))
GDCdownload(tcgaCNquery)

tcgaCN <- GDCprepare(tcgaCNquery, summarizedExperiment = FALSE)

genes <- c("TP53",
           "ARID1A",
           "NLRC5",
           "CSDE1",
           "CIITA",
           "ERBB2",
           "CDKN2A",
           "CANX",
           "ERAP1",
           "CD1C",
           "HLA-DRB1",
           "ERAP2",
           "CD1A",
           "SPPL2A",
           "TAPBP",
           "CD1B",
           "PMS2",
           "HLA-DQA1",
           "PDIA3",
           "LGMN",
           "CTSS",
           "CD1D",
           "PSMB9",
           "MR1",
           "TAPBPL",
           "TAP2",
           "HLA-A",
           "HLA-DQA2",
           "HLA-B",
           "HLA-DPA1",
           "TAP1",
           "HLA-DQB2",
           "HLA-DMA",
           "HLA-DRB5",
           "CALR",
           "HLA-DPB1",
           "B2M",
           "CD74",
           "RFXANK",
           "RFX5",
           "HLA-DRA",
           "HLA-C",
           "HLA-G",
           "PSMB8",
           "HLA-DOA",
           "IRF1",
           "CTSL",
           "HLA-DMB",
           "HLA-E",
           "PSMB10",
           "RFXAP")


#Write out copy number segment for GISTIC analysis

write.csv(tcgaCN, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Copy number/TCGA copy number segment.csv")


#Download the TCGA mutation data and place into a maf
tcgaCNquery <- GDCquery(project = "TCGA-ESCA", 
                        legacy = FALSE,
                        data.category = "Copy Number Variation",
                        data.type = "Gene Level Copy Number",
                        workflow.type = "ASCAT2",
                        barcode =  c("TCGA-L7-A6VZ", "TCGA-IG-A4QS", "TCGA-L5-A8NV", "TCGA-R6-A6XG", "TCGA-R6-A6DQ", "TCGA-L5-A8NE", "TCGA-L5-A8NT", "TCGA-L5-A4OQ", "TCGA-2H-A9GH", "TCGA-R6-A8W8", "TCGA-L5-A8NJ", "TCGA-V5-AASX", "TCGA-L5-A4OP", "TCGA-L5-A8NF", "TCGA-L5-A8NM", "TCGA-L5-A8NR", "TCGA-2H-A9GN", "TCGA-L5-A8NS", "TCGA-L5-A4OI", "TCGA-2H-A9GO", "TCGA-R6-A8W5", "TCGA-L5-A88V", "TCGA-L5-A8NI", "TCGA-L5-A4OW", "TCGA-L5-A4OJ", "TCGA-L5-A8NW", "TCGA-R6-A6L6", "TCGA-L5-A891", "TCGA-RE-A7BO", "TCGA-IC-A6RE", "TCGA-JY-A6F8", "TCGA-JY-A93E", "TCGA-VR-AA4D", "TCGA-R6-A6KZ", "TCGA-S8-A6BV", "TCGA-X8-AAAR", "TCGA-L5-A4OO", "TCGA-JY-A6FB", "TCGA-JY-A93C", "TCGA-2H-A9GK", "TCGA-2H-A9GF", "TCGA-2H-A9GI", "TCGA-2H-A9GL", "TCGA-2H-A9GM", "TCGA-2H-A9GR","TCGA-JY-A6FH", "TCGA-JY-A93D", "TCGA-JY-A938", "TCGA-L5-A4OE", "TCGA-L5-A4OF", "TCGA-L5-A4OG","TCGA-L5-A4OH", "TCGA-L5-A4ON", "TCGA-L5-A4OR", "TCGA-L5-A4OS", "TCGA-L5-A4OT", "TCGA-L5-A4OU", "TCGA-L5-A4OX", "TCGA-L5-A8NG", "TCGA-L5-A8NH", "TCGA-L5-A8NL", "TCGA-L5-A8NN", "TCGA-L5-A8NU", "TCGA-L5-A43C", "TCGA-L5-A43I", "TCGA-L5-A88Y", "TCGA-M9-A5M8", "TCGA-Q9-A6FW", "TCGA-R6-A6DN", "TCGA-R6-A6L4", "TCGA-R6-A6XQ", "TCGA-R6-A6Y2", "TCGA-R6-A8WC", "TCGA-R6-A8WG", "TCGA-V5-A7RB", "TCGA-V5-A7RE", "TCGA-V5-AASW", "TCGA-VR-A8EQ", "TCGA-ZR-A9CJ"))
GDCdownload(tcgaCNquery)

tcgaCN <- GDCprepare(tcgaCNquery, summarizedExperiment = FALSE)

esca.subtype <- TCGAquery_subtype(tumor = "esca")

ploidytcga <- esca.subtype[c("patient", "Absolute extract ploidy")]

tcgaCNfilt <- tcgaCN[ tcgaCN$gene_name %in% genes ,]

#remove the min max chromosome after finding no indication of LOH

tcgaCNfilt <- select(tcgaCNfilt, -contains("min"))
tcgaCNfilt <- select(tcgaCNfilt, -contains("max"))

# Remove other genomic identifiers other than gene name
tcgaCNfilt <- select(tcgaCNfilt, -contains(c("chromosome", "gene_id", "start", "end" )))

#merge the ploidy with the copy number
TtcgaCNfilt <- as.data.frame(t(tcgaCNfilt))
rownames(TtcgaCNfilt) <- strtrim(rownames(TtcgaCNfilt), 12)
colnames(TtcgaCNfilt) <- TtcgaCNfilt[1,]
TtcgaCNfilt <- TtcgaCNfilt[-1,]

rownames(ploidytcga) <- ploidytcga$patient

ploidyCNtcgafilt <- ploidyCNtcga[ ploidyCNtcga$patient %in% rownames(TtcgaCNfilt) ,]
ploidyCNtcgafilt <- select(ploidyCNtcgafilt, -contains(c("V", "Row.names" )))
rownames(ploidyCNtcgafilt) <- ploidyCNtcgafilt$patient

ploidyCNtcgafinal <- merge(ploidyCNtcgafilt, TtcgaCNfilt, by = 0, all = TRUE)

ploidyCNtcgafinal <- select(ploidyCNtcgafinal, -contains(c("Row.names" )))

# Write the ploidy data to csv and produce gain loss neutral to the specification set out here https://cancer.sanger.ac.uk/cosmic/help/cnv/overview

write.csv(ploidyCNtcgafinal, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Copy number/ploidyCNtcgafinal.csv")

ploidyCNtcgafinal <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Copy number/ploidyCNtcgafinal.csv", header = TRUE, row.names = 1)

ploidyCNtcgafinal <- ploidyCNtcgafinal[,-1]

icgcSimpleMutationToMAF(
  "C:/Users/wp1g19/Documents/Revisit bioinformatics/Copy number/ICGC/copy_number_somatic_mutation.ESAD-UK.tsv",
  basename = NA,
  MAFobj = FALSE,
  clinicalData = NULL,
  removeDuplicatedVariants = TRUE,
  addHugoSymbol = TRUE
)

clusterfractional1 <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/CIBERSORTx_Fractional.csv", header = TRUE, row.names = 1)
rownames(clusterfractional1) <- gsub("/",".", rownames(clusterfractional1))
rownames(clusterfractional1) <- gsub("-",".", rownames(clusterfractional1))

wr <- merge(clustersFract, clusterfractional1, by = 0)
 write.csv(x = wr, file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/CIBERSORTx_Fractional_withclusters.csv")
 
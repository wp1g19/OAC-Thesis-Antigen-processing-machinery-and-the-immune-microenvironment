#Correlation plot of Cibersort results
library(ggplot2)
library(gridExtra)
library(ggpubr)
#read in data
CIBERSORT_COMP <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Digital cytometry/TPM rerun/Methyl CIBERSORT comparisons.csv", header = TRUE, row.names = 1)
CIBERSORT_COMP$T.cells.CD8.CIBERSORT.B.mode
CIBERSORT_COMP.lm <- lm(CD8..MethylCIBERSORT. ~ T.cells.CD8.CIBERSORT.B.mode, data = CIBERSORT_COMP)
CIBERSORT_COMP.lm2 <- lm(CD8..MethylCIBERSORT. ~ 0 + T.cells.CD8.CIBERSORT.B.mode, data = CIBERSORT_COMP) # Adding the 0 term tells the lm() to fit the line through the origin
summary(CIBERSORT_COMP.lm)
summary(CIBERSORT_COMP.lm2)


lm.scatter1 <- ggplot(CIBERSORT_COMP, aes(x=CD8..MethylCIBERSORT., y=T.cells.CD8.CIBERSORT.B.mode)) + 
  geom_point(color='#2980B9', size = 4) + xlim(c(0, 0.3)) + 
  geom_abline(intercept=0, slope=0.9, color='#2C3E50', size=1.1) + 
  labs(title='A')+
  stat_cor(label.y = 0.7)+ 
  stat_regline_equation(label.y = 0.6) 


CIBERSORT_COMP.lm3 <- lm(CD8..MethylCIBERSORT. ~ 0 + T.cells.CD8.CIBERSORTx.B.mode.Combined, data = CIBERSORT_COMP)
summary(CIBERSORT_COMP.lm3)

lm.scatter2 <- ggplot(CIBERSORT_COMP, aes(x=CD8..MethylCIBERSORT., y=T.cells.CD8.CIBERSORTx.B.mode.Combined)) + 
  geom_point(color='#2980B9', size = 4) + xlim(c(0, 0.3)) + 
  geom_abline(intercept=0, slope=0.3858, color='#2C3E50', size=1.1) + 
  labs(title='A')+
  stat_cor(label.y = 0.7)+ 
  stat_regline_equation(label.y = 0.6) 

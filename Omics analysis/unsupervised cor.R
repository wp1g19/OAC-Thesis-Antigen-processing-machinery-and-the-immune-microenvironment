
#unsupervised correlations 
TCGAOCCAMSmergedComBatBatchC <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGATMMMergedBatchCorrectedAPMCSDE1OutliersCut.csv", header = TRUE, row.names = 1)
mat <- t(TCGAOCCAMSmergedComBatBatchC)
library(parallel)
library(doParallel)
cores <- makeCluster(detectCores(), type='PSOCK') # grabs max available

cores <- 6  # explicitly choose number of cores

cl <- makeCluster(getOption('cl.cores', cores))
registerDoParallel(cl)
registerDoSEQ()
on.exit(stopCluster(cl))

res <- foreach(i = seq_len(ncol(mat)),
               .combine = rbind,
               .multicombine = TRUE,
               .inorder = FALSE,
               .packages = c('data.table', 'doParallel')) %dopar% {
                 cor(mat[,i], mat, method = 'pearson')
               }

rownames(res) <- colnames(res)


genes <- c(   "RFX5",
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
                   "HLA-E",
                   "HLA-DMB",
                   "HLA-DPA1",
                   "HLA-DOA",
                   "HLA-DMA",
                   "PSMB9",
                   "HLA-DQA2",
                   "HLA-B",
                   "HLA-DQB2",
                   "HLA-C",
                   "HLA-DRA",
                   "TAP2",
                   "PSMB8",
                   "HLA-DRB5",
                   "HLA-DQA1",
                   "TAP1",
                   "HLA-DRB1",
                   "TAPBP",
                   "HLA-G",
                   "HLA-A",
                   "CSDE1",
                   "PTPN2",
                   "SMYD3",
                   "NLRC5",
                   "IRF1",
                   "PDCD1",
              "CD274",
              "PAF1",
              "CTLA4",
              "LAG3",
              "HAVCR2")
resAPM <- res[,colnames(res) %in% genes ]
resAPM2 <- res[row.names(res) %in% genes,colnames(res) %in% genes ]
resAPM2 <- res[row.names(res) %in% genes,colnames(res) %in% genes]
rescsde1 <- as.data.frame(t((TCGAOCCAMSmergedComBatBatchC)))
res <- cor.test(rescsde1$CSDE1, rescsde1$`HLA-A`, 
                method = "pearson")
res

#apm genes
rescsde1 <- rescsde1[,colnames(rescsde1) %in% genes]
library(psych)
cor_test_mat <- corr.test(rescsde1)$p    # Apply corr.test function
cor_test_mat                         # Print matrix of p-values
cor_mat <- cor(rescsde1)                # Correlation matrix of example data
cor_mat     
library("corrplot") 
corrplot(cor_mat,                    # Draw corrplot with p-values
         p.mat = cor_test_mat,
         insig = "p-value")
library("ggcorrplot")                # Load ggcorrplot package

ggcorrplot(cor_mat)                  # Draw ggcorrplot



#+++++++++++++++++++++++
# Helper Functions
#+++++++++++++++++++++++

# Get lower triangle of the correlation matrix
.get_lower_tri <- function(cormat, show.diag = FALSE) {
  if (is.null(cormat)) {
    return(cormat)
  }
  cormat[upper.tri(cormat)] <- NA
  if (!show.diag) {
    diag(cormat) <- NA
  }
  return(cormat)
}

# Get upper triangle of the correlation matrix
.get_upper_tri <- function(cormat, show.diag = FALSE) {
  if (is.null(cormat)) {
    return(cormat)
  }
  cormat[lower.tri(cormat)] <- NA
  if (!show.diag) {
    diag(cormat) <- NA
  }
  return(cormat)
}

.remove_diag <- function(cormat) {
  if (is.null(cormat)) {
    return(cormat)
  }
  diag(cormat) <- NA
  cormat
}
# hc.order correlation matrix
.hc_cormat_order <- function(cormat, hc.method = "complete") {
  dd <- stats::as.dist((1 - cormat) / 2)
  hc <- stats::hclust(dd, method = hc.method)
  hc$order
}

.no_panel <- function() {
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  )
}


# Convert a tbl to matrix
.tibble_to_matrix <- function(x) {
  x <- as.data.frame(x)
  rownames(x) <- x[, 1]
  x <- x[, -1]
  as.matrix(x)
}

#APM genes from outliers cut

outliercutAPM <- read.csv(file = "C:/Users/wp1g19/Documents/Revisit bioinformatics/Combined/OCCAMSTCGATMMMergedBatchCorrectedAPMCSDE1OutliersCut.csv", header = TRUE, row.names = 1)

outliercutAPM2 <- outliercutAPM[ genes,]

#quick rewrite of ggcorplot to label the significant values not the insignificant values

# function body
#quick rewrite of ggcorplot to label the significant values not the insignificant values

# function body
ggcorrplot <- function(corr,
                       method = c("square", "circle"),
                       type = c("full", "lower", "upper"),
                       ggtheme = ggplot2::theme_minimal,
                       title = "",
                       show.legend = TRUE,
                       legend.title = "Corr",
                       show.diag = NULL,
                       colors = c("blue", "white", "red"),
                       outline.color = "gray",
                       hc.order = FALSE,
                       hc.method = "complete",
                       lab = FALSE,
                       lab_col = "black",
                       lab_size = 4,
                       p.mat = NULL,
                       sig.level = 0.05,
                       insig = c("pch", "blank"),
                       pch = 4,
                       pch.col = "black",
                       pch.cex = 5,
                       tl.cex = 12,
                       tl.col = "black",
                       tl.srt = 45,
                       digits = 2,
                       as.is = FALSE) {
  type <- match.arg(type)
  method <- match.arg(method)
  insig <- match.arg(insig)
  if (is.null(show.diag)) {
    if (type == "full") {
      show.diag <- TRUE
    } else {
      show.diag <- FALSE
    }
  }
  
  if (inherits(corr, "cor_mat")) {
    # cor_mat object from rstatix
    cor.mat <- corr
    corr <- .tibble_to_matrix(cor.mat)
    p.mat <- .tibble_to_matrix(attr(cor.mat, "pvalue"))
  }
  
  if (!is.matrix(corr) & !is.data.frame(corr)) {
    stop("Need a matrix or data frame!")
  }
  corr <- as.matrix(corr)
  
  corr <- base::round(x = corr, digits = digits)
  
  if (hc.order) {
    ord <- .hc_cormat_order(corr, hc.method = hc.method)
    corr <- corr[ord, ord]
    if (!is.null(p.mat)) {
      p.mat <- p.mat[ord, ord]
      p.mat <- base::round(x = p.mat, digits = digits)
    }
  }
  
  if (!show.diag) {
    corr <- .remove_diag(corr)
    p.mat <- .remove_diag(p.mat)
  }
  
  # Get lower or upper triangle
  if (type == "lower") {
    corr <- .get_lower_tri(corr, show.diag)
    p.mat <- .get_lower_tri(p.mat, show.diag)
  } else if (type == "upper") {
    corr <- .get_upper_tri(corr, show.diag)
    p.mat <- .get_upper_tri(p.mat, show.diag)
  }
  
  # Melt corr and pmat
  corr <- reshape2::melt(corr, na.rm = TRUE, as.is = as.is)
  colnames(corr) <- c("Var1", "Var2", "value")
  corr$pvalue <- rep(NA, nrow(corr))
  corr$signif <- rep(NA, nrow(corr))
  
  if (!is.null(p.mat)) {
    p.mat <- reshape2::melt(p.mat, na.rm = TRUE)
    corr$coef <- corr$value
    corr$pvalue <- p.mat$value
    corr$signif <- as.numeric(p.mat$value <= sig.level)
    p.mat <- subset(p.mat, p.mat$value < sig.level)
    if (insig == "blank") {
      corr$value <- corr$value * corr$signif
    }
  }
  
  
  corr$abs_corr <- abs(corr$value) * 10
  
  # heatmap
  p <-
    ggplot2::ggplot(
      data = corr,
      mapping = ggplot2::aes_string(x = "Var1", y = "Var2", fill = "value")
    )
  
  # modification based on method
  if (method == "square") {
    p <- p +
      ggplot2::geom_tile(color = outline.color)
  } else if (method == "circle") {
    p <- p +
      ggplot2::geom_point(
        color = outline.color,
        shape = 21,
        ggplot2::aes_string(size = "abs_corr")
      ) +
      ggplot2::scale_size(range = c(4, 10)) +
      ggplot2::guides(size = "none")
  }
  
  # adding colors
  p <- p + ggplot2::scale_fill_gradient2(
    low = colors[1],
    high = colors[3],
    mid = colors[2],
    midpoint = 0,
    limit = c(-1, 1),
    space = "Lab",
    name = legend.title
  )
  
  # depending on the class of the object, add the specified theme
  if (class(ggtheme)[[1]] == "function") {
    p <- p + ggtheme()
  } else if (class(ggtheme)[[1]] == "theme") {
    p <- p + ggtheme
  }
  
  
  p <- p +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = tl.srt,
        vjust = 1,
        size = tl.cex,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(size = tl.cex)
    ) +
    ggplot2::coord_fixed()
  
  label <- round(x = corr[, "value"], digits = digits)
  if (!is.null(p.mat) & insig == "blank") {
    ns <- corr$pvalue > sig.level
    if (sum(ns) > 0) label[ns] <- " "
  }
  
  # matrix cell labels
  if (lab) {
    p <- p +
      ggplot2::geom_text(
        mapping = ggplot2::aes_string(x = "Var1", y = "Var2"),
        label = label,
        color = lab_col,
        size = lab_size
      )
  }
  
  # matrix cell glyphs
  if (!is.null(p.mat) & insig == "pch") {
    p <- p + ggplot2::geom_point(
      data = p.mat,
      mapping = ggplot2::aes_string(x = "Var1", y = "Var2"),
      shape = pch,
      size = pch.cex,
      color = pch.col
    )
  }
  
  # add titles
  if (title != "") {
    p <- p +
      ggplot2::ggtitle(title)
  }
  
  # removing legend
  if (!show.legend) {
    p <- p +
      ggplot2::theme(legend.position = "none")
  }
  
  # removing panel
  p <- p +
    .no_panel()
  p
}



#' Compute the matrix of correlation p-values
#'
#' @param x numeric matrix or data frame
#' @param ... other arguments to be passed to the function cor.test.
#' @rdname ggcorrplot
#' @export

cor_pmat <- function(x, ...) {
  
  # initializing values
  mat <- as.matrix(x)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  
  # creating the p-value matrix
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- stats::cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  
  # name rows and columns of the p-value matrix
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  
  # return the final matrix
  p.mat
}

library(psych)
cor_test_mat <- corr.test(t(outliercutAPM2))$p    # Apply corr.test function
cor_test_mat                         # Print matrix of p-values
cor_mat <- cor(t(outliercutAPM2))              # Correlation matrix of example data
cor_mat     
library("corrplot") 
corrplot(cor_mat,                    # Draw corrplot with p-values
         p.mat = cor_test_mat,
         insig = "p-value")
ggcorrplot(cor_mat) 

p1 <- ggcorrplot(cor_mat,                  # Draw ggcorrplot with p-values
                 p.mat = cor_test_mat, type ="lower",show.legend = TRUE, lab = TRUE, sig.level = 0.05, pch = 1,pch.cex = 0,pch.col = "red", lab_size = 1 )
p1 + ggtitle("Correlation of APM genes and CSDE1 n = 176")

#TILS correlation
library(ggplot2)
library(dplyr)
cor <- Hmisc::rcorr(t(outliercutAPM2) %>% as.matrix())
nm = rownames(cor$r)
m = t(combn(nm, 2))
d = cbind(data.frame(m), R = cor$r[m], P = cor$P[m])
d$label = round(d$R, 2)
d$label[d$P < 0.05] = paste0(d$label[d$P < 0.05], "\n*")
d$X1 = factor(d$X1, nm)
d$X2 = factor(d$X2, rev(nm))

graphics.off()
ggplot(d, aes(X2, X1, fill = R, label = label)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = .02
  ) +
  geom_text(color = ifelse(d$R > 0.35, "black", "black"),size = 4, lineheight = .8) + 
  theme_bw() +
  coord_equal()+ theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + ggtitle("Correlation of MHC genes and markers of T cell exhaustion") + theme(text = element_text(size = 13)) 


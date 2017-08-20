#' Convert TPM to FPKM
#'
#' Convert TPM to FPKM based on the given gene length. This is
#' useful when gene expression in neurons are represented as TPM, 
#' e.g. when quantifying with RSEM.
#'
#' @param tpm A matrix of gene expression represented as TPM. Each row represent one gene and each column represent one sample
#' @param geneLength A vector with the same length as the row number of tpm matrix, representing length of genes in bp
#' @return A matrix of FPKM
#' @export
tpmToFpkm <- function(tpm, geneLength){
	count <- apply(tpm, 2, function(e){ e/1000*geneLength })
	fpkm <- apply(count, 2, function(e){ e*1000*1000000/sum(e)/geneLength })
	return(fpkm)
}

#' Data standardization
#'
#' Standardize FPKM input by log10-transformation followed
#' by central scaled according to gene expression average 
#' and sd in the training data. Missing expression values 
#' are randomly imputated based on the training data average
#' and sd.
#'
#' @param fpkm A matrix of gene expression represented as FPKM. Each row represent one gene (named by EnsemblID) and each column represent one sample
#' @return A matrix of standardized FPKM
#' @export
dataStandardize <- function(fpkm, standardize=TRUE){
	num.detected <- sum(rownames(fpkm) %in% names(mean_markers))
	if(num.detected < length(mean_markers)/1.5){
		stop(paste0("Too few marker genes. Only ", num.detected, " marker genes are available while at least ", ceiling(length(mean_markers)/1.5), " is required."))
	}
	if(standardize){
		fpkm <- log10(fpkm + 1)
	}
	
	fpkm.filled <- t(sapply(names(mean_markers), function(gene){
		if(sum(rownames(fpkm)==gene)>0){
			return(as.numeric(fpkm[gene,]))
		} else{
			return(rnorm(ncol(fpkm), mean_markers[gene], sd_markers[gene]))
		}
	}))
	if(standardize){
		fpkm.std <- apply(fpkm.filled, 2, function(e){ (e - mean_markers[rownames(fpkm.filled)]) / sd_markers[rownames(fpkm.filled)] })
	} else{
		fpkm.std <- fpkm.filled
	}
	colnames(fpkm.std) <- colnames(fpkm)
	return(fpkm.std)
}

#' Neuron Maturity Index (NMI) estimation
#'
#' Taking gene expressions of neuron samples (in FPKM), this
#' function predicts the Neuron Maturity Index (NMI) which
#' is an estimate of degree of neuron maturation. The Index
#' ranges between 0 and 1. The more it is close to one, the
#' more mature the neuron sample is. The models were trained
#' using the mature and immature neuron single cell RNA-seq
#' data in Darmanis, et al. PNAS, 2015.
#'
#' @param fpkm A matrix of gene expression represented by FPKM. Each row represent one gene (named by EnsemblID) and each column represent one sample
#' @param standardize To determine whether standardization should be applied to fpkm matrix first
#' @param returnFpkmStd To determine whether the standardlized FPKM of samples should be returned.
#' @return A list with at least two components. The first component (modularNMI) is a matrix of modular NMI, with 109 columns each of which shows NMI of one module. The second component (overallNMI) is a matrix of overall NMI, with three columns representing NMI-all, NMI-discriminable and Neuron Activity Index (NAI). The third element is a matrix of standardized FPKM if returnFpkmStd is TRUE.
#' @export
predictNMI <- function(fpkm, standardize = TRUE, returnFpkmStd = FALSE){
	require("Matrix")
	fpkm <- dataStandardize(fpkm, standardize = standardize)
	fpkm.std <- rbind(rep(1, ncol(fpkm)), fpkm)
	pred <- t(fpkm.std) %*% coefs
	pred <- exp(pred) / (exp(pred) + 1)
	
	pred.overall <- as.numeric(pred %*% as.matrix(modules_info$weight)) / sum(modules_info$weight)
	pred.discriminable <- as.numeric(pred[,modules_info$discriminable] %*% as.matrix(modules_info$weight[modules_info$discriminable])) / sum(modules_info$weight[modules_info$discriminable])
	pred.activity <- as.numeric(pred[,modules_info$discriminable & modules_info$matureHigh] %*% as.matrix(modules_info$weight[modules_info$discriminable & modules_info$matureHigh])) / sum(modules_info$weight[modules_info$discriminable & modules_info$matureHigh])
	pred.integrated <- data.frame(overall = pred.overall, discriminable = pred.discriminable, activity = pred.activity)
	rownames(pred.integrated) <- colnames(fpkm)
	
	if(returnFpkmStd){
		res <- list(modularNMI = pred, overallNMI = pred.integrated, fpkmStd = fpkm.std)
	} else{
		res <- list(modularNMI = pred, overallNMI = pred.integrated)
	}
	return(res)
}

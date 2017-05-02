#' estimateSex Function
#'
#' This function estimates the sex of DNA methylation microarray data frmo the distribution of beta-value and detection p-value of methylation sites from X chromosome and Y chromosome.
#' @param beta-value	beta-value matrix 
#' @param detecP	detection p-value matrix
#' @param beta.intervals.X	Intervals of beta-value of X chromosome probes within each of which the proportion is calculated. Default is 'seq(0,1,0.1)', which is 0 to 0.1, 0.1 to 0.2, ..., 0.9 to 1.0. 
#' @param beta.intervals.Y	Intervals of beta-value of Y chromosome probes within each of which the proportion is calculated. Default is the same as that of 'beta.intervals.X'
#' @param p.intervals.X	Intervals of detection p-value of X chromosome probes within each of which the proportion is calculated. Default is -Inf to 1e-18m 1e-18 to 1e-5, 1e-5 to 1e-2, 1e-2 to 0.
#' @param p.intervals.Y	Intervals of detection p-value of Y chromosome probes within each of which the proportion is calculated. Default is the same as that of 'p.intervals.X'
#' @param return_with_reference	Logical value to indicate whether or not to return the data derived from references data (Default to FALSE).
#' @keywords estimateSex
#' @export
#' @examples



estimateSex <- function(beta.value=NULL, detecP=NULL, beta.intervals.X=seq(0,1,0.1), beta.intervals.Y=seq(0,1,0.1), p.intervals.X=c(-18,-5,-2,0), p.intervals.Y=c(-18,-5,-2,0), return_with_reference=FALSE){

	# re-compile proportion table for reference samples
	tab.prop.reference <- .get.proportion_table.base_prop(beta.intervals.X, beta.intervals.Y, p.intervals.X, p.intervals.Y)

	DNAm_columns <- colnames(tab.prop.reference)[grep(colnames(tab.prop.reference), pattern="^p.", invert=TRUE)]
	DNAm_columns.X <- colnames(tab.prop.reference)[grep(colnames(tab.prop.reference), pattern="^X:", invert=FALSE)]
	DNAm_columns.Y <- colnames(tab.prop.reference)[grep(colnames(tab.prop.reference), pattern="^Y:", invert=FALSE)]
	pval_columns <- colnames(tab.prop.reference)[grep(colnames(tab.prop.reference), pattern="^p.")]
	pval_columns.X <- colnames(tab.prop.reference)[grep(colnames(tab.prop.reference), pattern="^p.X:")]
	pval_columns.Y <- colnames(tab.prop.reference)[grep(colnames(tab.prop.reference), pattern="^p.Y:")]
	columns.X <- c(DNAm_columns.X, pval_columns.X)
	columns.Y <- c(DNAm_columns.Y, pval_columns.Y)

	# obtain proportion table for test samples
	tab.prop.test <- get.proportion_table(beta.value, detecP, beta.intervals.X, beta.intervals.Y, p.intervals.X, p.intervals.Y)
	samples.test <- rownames(tab.prop.test)

	# merge reference set and test set.
	tab.prop.all <- rbind(tab.prop.reference, tab.prop.test)
	samples.all <- rownames(tab.prop.all)

	# start prediction by kmeans from PCA results


	pca.all.X <- prcomp(tab.prop.all[, columns.X], scale.=FALSE)
	pca.all.Y <- prcomp(tab.prop.all[, columns.Y], scale.=FALSE)
	kmeans.X <- kmeans(pca.all.X$x, 2)
	kmeans.Y <- kmeans(pca.all.Y$x, 2)

	cluster_name.male.X <- names(table(kmeans.X$cluster[reference.male]))[order(table(kmeans.X$cluster[reference.male]), decreasing=TRUE)[1]]
	cluster_name.male.Y <- names(table(kmeans.Y$cluster[reference.male]))[order(table(kmeans.Y$cluster[reference.male]), decreasing=TRUE)[1]]
	cluster_name.female.X <- names(table(kmeans.X$cluster[reference.female]))[order(table(kmeans.X$cluster[reference.female]), decreasing=TRUE)[1]]
	cluster_name.female.Y <- names(table(kmeans.Y$cluster[reference.female]))[order(table(kmeans.Y$cluster[reference.female]), decreasing=TRUE)[1]]

	test.male.X <- names(kmeans.X$cluster)[kmeans.X$cluster == cluster_name.male.X & (!names(kmeans.X$cluster) %in% reference.all)]
	test.male.Y <- names(kmeans.Y$cluster)[kmeans.Y$cluster == cluster_name.male.Y & (!names(kmeans.X$cluster) %in% reference.all)]
	test.female.X <- names(kmeans.X$cluster)[kmeans.X$cluster == cluster_name.female.X & (!names(kmeans.X$cluster) %in% reference.all)]
	test.female.Y <- names(kmeans.Y$cluster)[kmeans.Y$cluster == cluster_name.female.Y & (!names(kmeans.X$cluster) %in% reference.all)]

	tab.pred <- data.frame(row.names=samples.all, X.PC1=pca.all.X$x[samples.all, "PC1"], Y.PC1=pca.all.Y$x[samples.all, "PC1"], predicted.X=NA, predicted.Y=NA, predicted=NA)
	tab.pred[, "predicted.X"] <- ifelse(kmeans.X$cluster[samples.all] == cluster_name.male.X, "M", "F")
	tab.pred[, "predicted.Y"] <- ifelse(kmeans.Y$cluster[samples.all] == cluster_name.male.Y, "M", "F")
	tab.pred[, "predicted"] <- ifelse(tab.pred[samples.all, "predicted.X"]==tab.pred[samples.all, "predicted.Y"], tab.pred[samples.all, "predicted.X"], "N")
	tab.pred[samples.all %in% reference.male, "predicted"] <- "M"
	tab.pred[samples.all %in% reference.female, "predicted"] <- "F"

	#tab.pred$group <- "reference"
	#tab.pred[samples.test, "group"] <- "test"

	sex_estimation <- list(test=tab.pred[samples.test, c("X.PC1", "Y.PC1", "predicted.X", "predicted.Y", "predicted")] )

	if (return_with_reference) {
		sex_estimation$reference <- tab.pred[reference.all, c("X.PC1", "Y.PC1", "predicted.X", "predicted.Y", "predicted")]
	}

	return (sex_estimation)
}




.get.proportion_table.base_prop <- function(beta.intervals.X=seq(0,1,0.1), beta.intervals.Y=seq(0,1,0.1), p.intervals.X=c(-18,-5,-2,0), p.intervals.Y=c(-18,-5,-2,0)){
	tab.prop.beta.X <- .merge_intervals(tab.base_prop=tab.base_prop.beta.X, intervals=beta.intervals.X)
	tab.prop.beta.Y <- .merge_intervals(tab.base_prop=tab.base_prop.beta.Y, intervals=beta.intervals.Y)
	tab.prop.p.X <- .merge_intervals(tab.base_prop=tab.base_prop.p.X, intervals=p.intervals.X)
	tab.prop.p.Y <- .merge_intervals(tab.base_prop=tab.base_prop.p.Y, intervals=p.intervals.Y)

	colnames(tab.prop.beta.X) <- paste("X", colnames(tab.prop.beta.X), sep=":")
	colnames(tab.prop.beta.Y) <- paste("Y", colnames(tab.prop.beta.Y), sep=":")
	colnames(tab.prop.p.X) <- paste("p.X", colnames(tab.prop.p.X), sep=":")
	colnames(tab.prop.p.Y) <- paste("p.Y", colnames(tab.prop.p.Y), sep=":")

	tab.prop <- cbind(tab.prop.beta.X, tab.prop.beta.Y, tab.prop.p.X, tab.prop.p.Y)
	return (tab.prop)
}



.merge_intervals <- function(tab.base_prop = NULL, intervals=NULL){
	base_prop.columns <- as.numeric(colnames(tab.base_prop))
	#tab.prop <- data.frame(row.names=rownames(tab.base_prop), matrix(nrow=nrow(tab.base_prop), ncol=(length(intervals)-1)))
	tab.prop <- data.frame(row.names=rownames(tab.base_prop))

	prev <- as.numeric(intervals[1])
	for (i in 2:length(intervals)) {
		curr <- as.numeric(intervals[i])
		cols <- as.character(base_prop.columns[base_prop.columns>prev & base_prop.columns <= curr])
		interval_name <- sprintf("(%g,%g]", prev, curr)
		tab.prop[, interval_name] <- rowSums(tab.base_prop[, cols])
		prev <- curr
	}
	return(tab.prop)
}


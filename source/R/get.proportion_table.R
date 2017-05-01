#' get.proportion_table Function
#'
#' This function returns the beta-value and detection p-value proportion table within each intervals specified by parameters.
#' @param beta-value	beta-value matrix 
#' @param detecP	detection p-value matrix
#' @param beta.intervals.X	Intervals of beta-value of X chromosome probes within each of which the proportion is calculated. Default is 'seq(0,1,0.1)', which is 0 to 0.1, 0.1 to 0.2, ..., 0.9 to 1.0. 
#' @param beta.intervals.Y	Intervals of beta-value of Y chromosome probes within each of which the proportion is calculated. Default is the same as that of 'beta.intervals.X'
#' @param p.intervals.X	Intervals of detection p-value of X chromosome probes within each of which the proportion is calculated. Default is -Inf to 1e-18m 1e-18 to 1e-5, 1e-5 to 1e-2, 1e-2 to 0.
#' @param p.intervals.Y	Intervals of detection p-value of Y chromosome probes within each of which the proportion is calculated. Default is the same as that of 'p.intervals.X'
#' @param samples	list of samples to calculate beta/detection p-value proportions. If NULL (default) all samples in the beta-value and detection p-value matrices are used.
#' @keywords	proportion_table
#' @export
#' @examples
#' proportion_table <- get.proportion_table(beta.value=beta, detecP=detecP)

get.proportion_table <- function(beta.value=NULL, detecP=NULL, beta.intervals.X=seq(0,1,0.1), beta.intervals.Y=seq(0,1,0.1), p.intervals.X=c(-18,-5,-2,0), p.intervals.Y=c(-18,-5,-2,0), samples=NULL){
	if (is.null(dim(beta.value))) {
		beta.value <- as.matrix(beta.value)
		colnames(beta.value) <- ifelse(is.null(samples), "1", samples)
	}
	if (is.null(dim(detecP))) {
		detecP <- as.matrix(detecP)
		colnames(detecP) <- "1"
		colnames(detecP) <- ifelse(is.null(samples), "1", samples)
	}

	if (is.null(samples)) {
		samples <- colnames(beta.value)
	} else {
		beta.value <- beta.value[, samples]
		detecP <- detecP[, samples]
	}

	if (is.null(rownames(beta.value)) | is.null(rownames(detecP))) {
		cat("cannot figure out probe names.\n")
		return (NULL)
	}

	probes.X <- probes.noSNP.X[probes.noSNP.X %in% rownames(beta.value)]
	probes.Y <- probes.noSNP.Y[probes.noSNP.Y %in% rownames(beta.value)]


	tab.prop.test.beta.X <- do.call(rbind, lapply(as.list(samples), function(x) table(cut(beta.value[probes.X, x], breaks=beta.intervals.X, include.lowest=FALSE))/length(probes.noSNP.X)))
	tab.prop.test.beta.Y <- do.call(rbind, lapply(as.list(samples), function(x) table(cut(beta.value[probes.Y, x], breaks=beta.intervals.Y, include.lowest=FALSE))/length(probes.noSNP.Y)))
	tab.prop.test.p.X <- do.call(rbind, lapply(as.list(samples), function(x) table(cut(log10(detecP[probes.X, x]+(1e-17)), breaks=p.intervals.X))/length(probes.noSNP.X)))
	tab.prop.test.p.Y <- do.call(rbind, lapply(as.list(samples), function(x) table(cut(log10(detecP[probes.Y, x]+(1e-17)), breaks=p.intervals.Y))/length(probes.noSNP.Y)))

	colnames(tab.prop.test.beta.X) <- paste("X", colnames(tab.prop.test.beta.X), sep=":")
	colnames(tab.prop.test.beta.Y) <- paste("Y", colnames(tab.prop.test.beta.Y), sep=":")
	colnames(tab.prop.test.p.X) <- paste("p.X", colnames(tab.prop.test.p.X), sep=":")
	colnames(tab.prop.test.p.Y) <- paste("p.Y", colnames(tab.prop.test.p.Y), sep=":")

	tab.prop.test <- cbind(tab.prop.test.beta.X, tab.prop.test.beta.Y, tab.prop.test.p.X, tab.prop.test.p.Y)
	rownames(tab.prop.test) <- samples
	return (tab.prop.test)
}


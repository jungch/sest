#' A plotSexEstimation Function
#'
#' This function generates a scatter plot using 'ggplot' showing the first principal component of two principal component analysis results (pincipal component analysis using the data from X chromosome and the data from Y chromosome)
#' @param sex_estimation        result of estimateSex function 
#' @param samples       the list of samples to be plotted. If NULL (default) all samples in the estimateSex result are plotted.
#' @param main  character string to be used as the plot title. Default is an empty string.
#' @param filename      When specified, a PDF is generated.
#' @param include_reference     Logical value to specify whether to plot the reference data  (Default is FALSE).
#' @keywords plotSexEstimation
#' @export
#' @examples
#' plotSexEstimation(sex_estimation = sexEst)

plotSexEstimation <- function(sex_estimation = NULL, samples=NULL, main="", filename=NULL, include_reference=FALSE){
	if (is.null(samples)) {
		samples <- rownames(sex_estimation$test)
	}

	reference_color <- "#404040"
	reference_shape <- c(M=17, F=19)

	tab.test <- sex_estimation$test[samples,]

	tab.reference <- NULL
	if ("reference" %in% names(sex_estimation)) {
		tab.reference <- sex_estimation$reference
	}

	if(include_reference) {
		if (! "reference" %in% names(sex_estimation)) {
			cat("warining: no information on 'reference' set\n")
		}
		ref.alpha=0.5
	} else {
		# when include_reference==FALSE but there is PCA information for reference samples, plot them invisibly.
		ref.alpha=0
	}

	tmp.tab <- rbind(tab.test, tab.reference)

	if (is.character(filename)) {
		pdf(filename)
	}

	def.par <- par(no.readonly = TRUE) # save default, for resetting...

	layout(matrix(c(1,2), ncol=1), heights=c(9,1))
	plot(0,0, type="n", xlim=range(tmp.tab$X.PC1), ylim=range(tmp.tab$Y.PC1), xlab="X.PC1", ylab="Y.PC1", main=main)

	if(include_reference) {
		points(tab.reference[, "X.PC1"], tab.reference[, "Y.PC1"], col=alpha(reference_color, ref.alpha), pch=reference_shape[tab.reference$predicted])
	}

	test.male <- rownames(tab.test[tab.test$predicted=="M",])
	test.female <- rownames(tab.test[tab.test$predicted=="F",])
	test.N <- rownames(tab.test[tab.test$predicted=="N",])
	
	points(tab.test[, "X.PC1"], tab.test[, "Y.PC1"], col=color.annots$sex.colors[tab.test$predicted], pch=color.annots$sex.shapes[tab.test$predicted])

	legend_text <- c("Male", "Female", "N")
	legend_pch <- color.annots$sex.shapes[legend_text]
	legend_col <- color.annots$sex.colors[legend_text]
	if(include_reference) {
		legend_text <- c(legend_text, "", "Male reference", "Female reference")
		legend_pch <- c(legend_pch, NA, reference_shape[c("M", "F")])
		legend_col <- c(legend_col, NA, alpha(rep(reference_color, 2), ref.alpha))
	}
	par(mai=c(0,0,0,0))
	plot.new()
	legend("top", legend=legend_text, pch=legend_pch, col=legend_col, ncol=3, bty="n")

	par(def.par)
	if (is.character(filename)) {
		dev.off()
	}
}


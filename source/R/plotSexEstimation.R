#' A plotSexEstimation Function
#'
#' This function generates a scatter plot using 'ggplot' showing the first principal component of two principal component analysis results (pincipal component analysis using the data from X chromosome and the data from Y chromosome)
#' @param sex_estimation	result of estimateSex function 
#' @param samples	the list of samples to be plotted. If NULL (default) all samples in the estimateSex result are plotted.
#' @param main	character string to be used as the plot title. Default is an empty string.
#' @param filename	When specified, a PDF is generated.
#' @param include_reference	Logical value to specify whether to plot the reference data  (Default is FALSE).
#' @keywords plotSexEstimation
#' @export
#' @examples
#' plotSexEstimation(sex_estimation = sexEst)

plotSexEstimation <- function(sex_estimation = NULL, samples=NULL, main="", filename=NULL, include_reference=FALSE){
	reference_color <- "#404040"
	reference_shape <- c(M=0, F=2)
	g_point <- ggplot() + theme_bw()
	if(include_reference) {
		if ("reference" %in% names(sex_estimation)) {
			tab.reference <- sex_estimation$reference
			tab.reference$sex <- tab.reference$predicted

			g_point <- g_point + geom_point(data=tab.reference, aes(x=X.PC1, y=Y.PC1, shape=predicted), col=reference_color, alpha=0.8, cex=2) + labs(shape="reference set") + scale_shape_manual(values=reference_shape)
		} else {
			cat("warining: no information on 'reference' set\n")
		}
	} else {
		if ("reference" %in% names(sex_estimation)) {
			tab.reference <- sex_estimation$reference
			tab.reference$sex <- tab.reference$predicted

			g_point <- g_point + geom_point(data=tab.reference, aes(x=X.PC1, y=Y.PC1, shape=predicted), col=reference_color, alpha=0, cex=2) +guides(shape=FALSE)
 		}
	}

	if (is.null(samples)) {
		samples <- rownames(sex_estimation$test)
	}

	tab.test <- sex_estimation$test[samples,]
	test.male <- rownames(tab.test[tab.test$predicted=="M",])
	test.female <- rownames(tab.test[tab.test$predicted=="F",])
	test.N <- rownames(tab.test[tab.test$predicted=="N",])

	g_point <- g_point + geom_point(data=tab.test, aes(x=X.PC1, y=Y.PC1, col=predicted), cex=2) + geom_point(cex=2) + scale_color_manual(values=color.annots$sex.colors)

	if (is.character(filename)) {
		pdf(filename)
	}

	print(g_point + ggtitle(main))

	if (is.character(filename)) {
		dev.off()
	}
}



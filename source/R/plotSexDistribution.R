#' A plotSexDistribution Function
#'
#' This function generates a line chart using 'ggplot' showing the distribution of beta-values of chrX and chrY and detection p-values of chrY.
#' @param beta-value    beta-value matrix 
#' @param detecP        detection p-value matrix
#' @param beta.intervals.X      Intervals of beta-value of X chromosome probes within each of which the proportion is calculated. Default is 'seq(0,1,0.1)', which is 0 to 0.1, 0.1 to 0.2, ..., 0.9 to 1.0. 
#' @param beta.intervals.Y      Intervals of beta-value of Y chromosome probes within each of which the proportion is calculated. Default is the same as that of 'beta.intervals.X'
#' @param p.intervals.X Intervals of detection p-value of X chromosome probes within each of which the proportion is calculated. Default is -Inf to 1e-18m 1e-18 to 1e-5, 1e-5 to 1e-2, 1e-2 to 0.
#' @param p.intervals.Y Intervals of detection p-value of Y chromosome probes within each of which the proportion is calculated. Default is the same as that of 'p.intervals.X'
#' @param samples	the list of samples to be plotted. If NULL (default) all samples in the beta.value matrix are plotted.
#' @param main	character string to be used as the plot title. Default is an empty string.
#' @param filename	When specified, a PDF is generated.
#' @param include_reference	Logical value to specify whether to plot the reference data  (Default is FALSE).
#' @param color	Logical value to specify whether to plot the reference data  (Default is '#984ea3', one of the shades of purple).
#' @keywords plotSexDistribution
#' @export
#' @examples
#' plotSexDistribution(beta.value=beta, detecP=pval)


plotSexDistribution <- function(beta.value=NULL, detecP=NULL, beta.intervals.X=seq(0,1,0.1), beta.intervals.Y=seq(0,1,0.1), p.intervals.X=c(-18,-5,-2,0), p.intervals.Y=c(-18,-5,-2,0), samples=NULL, filename=NULL, include_reference=TRUE, color=NULL){

	# prepare the plot for reference samples
	tab.prop.reference <- .get.proportion_table.base_prop(beta.intervals.X, beta.intervals.Y, p.intervals.X, p.intervals.Y)

	DNAm_columns <- colnames(tab.prop.reference)[grep(colnames(tab.prop.reference), pattern="^p.", invert=TRUE)]
	DNAm_columns.X <- colnames(tab.prop.reference)[grep(colnames(tab.prop.reference), pattern="^X:", invert=FALSE)]
	DNAm_columns.Y <- colnames(tab.prop.reference)[grep(colnames(tab.prop.reference), pattern="^Y:", invert=FALSE)]
	pval_columns <- colnames(tab.prop.reference)[grep(colnames(tab.prop.reference), pattern="^p.")]
	pval_columns.X <- colnames(tab.prop.reference)[grep(colnames(tab.prop.reference), pattern="^p.X:")]
	pval_columns.Y <- colnames(tab.prop.reference)[grep(colnames(tab.prop.reference), pattern="^p.Y:")]
	columns.X <- c(DNAm_columns.X, pval_columns.X)
	columns.Y <- c(DNAm_columns.Y, pval_columns.Y)

	chr_annotate <- data.frame(row.names=c("chrX", "chrY"), xpos = rep(Inf, 2), ypos = rep(Inf, 2), annotateText = c("chrX", "chrY"), hjustvar = rep(1, 2), vjustvar = rep(1, 2))

	g.ref.X.beta <- ggplot() + ylim(0,0.45) + theme_bw()
	g.ref.Y.beta <- ggplot() + ylim(0,0.45) + theme_bw()
	g.ref.Y.p <- ggplot() + ylim(0,1.0) + theme_bw()
	if(include_reference) {
		tab.prop.reference.melted <- melt(cbind(sample=rownames(tab.prop.reference), tab.prop.reference))

		colnames(tab.prop.reference.melted) <- c("sample", "range", "proportion")
		tab.prop.reference.melted$sample <- as.character(tab.prop.reference.melted$sample)
		tab.prop.reference.melted$range <- factor(tab.prop.reference.melted$range, levels=colnames(tab.prop.reference))
		tab.prop.reference.melted$sex <- factor(ifelse(tab.prop.reference.melted$sample %in% reference.male, "Male", "Female"), levels=c("Male", "Female"))

		g.ref.X.beta <- g.ref.X.beta + geom_line(data=tab.prop.reference.melted[tab.prop.reference.melted$range %in% DNAm_columns.X, ], aes(x=range, y=proportion, group=sample, col=sex), alpha=0.01, lty=2) + geom_boxplot(data=tab.prop.reference.melted[tab.prop.reference.melted$range %in% DNAm_columns.X, ], aes(x=range, y=proportion, fill=sex), alpha=0.2)

		g.ref.Y.beta <- g.ref.Y.beta + ylim(0,0.45) + geom_line(data=tab.prop.reference.melted[tab.prop.reference.melted$range %in% DNAm_columns.Y, ], aes(x=range, y=proportion, group=sample, col=sex), alpha=0.01, lty=2) + geom_boxplot(data=tab.prop.reference.melted[tab.prop.reference.melted$range %in% DNAm_columns.Y, ], aes(x=range, y=proportion, fill=sex), alpha=0.2)

		g.ref.Y.p <- g.ref.Y.p + ylim(0,1.00) + geom_line(data=tab.prop.reference.melted[tab.prop.reference.melted$range %in% pval_columns.Y, ], aes(x=range, y=proportion, group=sample, col=sex), alpha=0.01, lty=2) + geom_boxplot(data=tab.prop.reference.melted[tab.prop.reference.melted$range %in% pval_columns.Y, ], aes(x=range, y=proportion, fill=sex), alpha=0.2)
	}


	# obtain proportion table for test samples
	if (is.null(samples)) {
		samples <- rownames(beta.value)
	}
	if (length(samples) == 1) {
		tmp.beta <- data.frame(row.names=rownames(beta.value))
		tmp.beta[, samples] <- beta.value[, samples]
		beta.value <- tmp.beta
		tmp.pval <- data.frame(row.names=rownames(detecP))
		tmp.pval[, samples] <- detecP[, samples]
		detecP <- tmp.pval
	} else {
		beta.value <- beta.value[, samples]
		detecP <- detecP[, samples]
	}	


	tab.prop.test <- get.proportion_table(beta.value, detecP, beta.intervals.X, beta.intervals.Y, p.intervals.X, p.intervals.Y)
	tab.prop.test.melted <- melt(tab.prop.test)
	colnames(tab.prop.test.melted) <- c("sample", "range", "proportion")
	tab.prop.test.melted$sample <- as.character(tab.prop.test.melted$sample)
	tab.prop.test.melted$range <- factor(as.character(tab.prop.test.melted$range), levels=c(DNAm_columns.X, DNAm_columns.Y, pval_columns.X, pval_columns.Y))

	if (is.null(color)) {
		color <- "#984ea3"
	}

	g.ref.X.beta <- g.ref.X.beta + geom_text(data=chr_annotate["chrX",], aes(x=xpos, y=ypos, hjust=hjustvar, vjust=vjustvar,label=annotateText)) + scale_color_manual(values=color.annots$sex.colors) + scale_fill_manual(values=color.annots$sex.colors) + scale_x_discrete(labels=beta_range) + xlab("beta-value range") + ylab("proportion") + guides(fill=FALSE, col=FALSE) + theme(axis.text.x = element_text(angle=30, hjust=1))
	g.ref.Y.beta <- g.ref.Y.beta + geom_text(data=chr_annotate["chrY",], aes(x=xpos, y=ypos, hjust=hjustvar, vjust=vjustvar,label=annotateText)) + scale_color_manual(values=color.annots$sex.colors) + scale_fill_manual(values=color.annots$sex.colors) + scale_x_discrete(labels=beta_range) + xlab("beta-value range") + ylab("proportion") + guides(fill=FALSE, col=FALSE) + theme(axis.text.x = element_text(angle=30, hjust=1))
	g.ref.Y.p <- g.ref.Y.p + geom_text(data=chr_annotate["chrY",], aes(x=xpos, y=ypos, hjust=hjustvar, vjust=vjustvar,label=annotateText)) + scale_color_manual(values=color.annots$sex.colors) + scale_fill_manual(values=color.annots$sex.colors) + scale_x_discrete(labels=pval_range) + xlab("detection p-value range (log10)") + ylab("proportion") + theme(legend.position="bottom") + guides(col=FALSE) + labs(fill="reference set")


	#g_box <- g_box + geom_point(data=tab.prop.test.melted, aes(x=range, y=proportion), pch=4, col=color, alpha=0.8)
	d.X.beta <- geom_line(data=tab.prop.test.melted[tab.prop.test.melted$range %in% DNAm_columns.X,], aes(x=range, y=proportion, group=sample), col=color)
	d.Y.beta <- geom_line(data=tab.prop.test.melted[tab.prop.test.melted$range %in% DNAm_columns.Y,], aes(x=range, y=proportion, group=sample), col=color)
	#d.X.p <- geom_line(data=tab.prop.test.melted[tab.prop.test.melted$range %in% pval_columns.X,], aes(x=range, y=proportion, group=sample), col=color, alpha=0.5)
	d.Y.p <- geom_line(data=tab.prop.test.melted[tab.prop.test.melted$range %in% pval_columns.Y,], aes(x=range, y=proportion, group=sample), col=color)

	g.X.beta <- g.ref.X.beta + d.X.beta
	g.Y.beta <- g.ref.Y.beta + d.Y.beta
	g.Y.p <- g.ref.Y.p + d.Y.p

	if (is.character(filename)) {
		pdf(filename)
	}

	print(grid.arrange(g.X.beta, g.Y.beta, g.Y.p, ncol=1))
		
	if (is.character(filename)) {
		dev.off()
	}
}


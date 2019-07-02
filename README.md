# sEst
An R-package for estimating the sex from DNA methylation microarray data

## Installation
### 1) Click on the package source file ('sest.tar').

### 2) Click on the 'Download' button.

### 3) Install the package from source:
a) In command line:
```
R CMD INSTALL '/path/to/package source'
```
b) Or, in R environment:
```{r}
install.packages(pkgs="/path/to/package source", repos=NULL, type="source")
```
### Alternatively:
```{r}
install.packages(pkgs="https://github.com/jungch/sest/raw/master/sest.tar", repos=NULL, type="source")
```

## Example codes
[Here](example/sEst_example.md) is an example scenario of estimating the sex of HM450k data samples.
This example uses a publicly available dataset from GEO (GSE51032).


## Functions:
### 'estimateSex'
##### usage:
```{r}
estimateSex(beta.value = NULL, detecP = NULL,
  beta.intervals.X = seq(0, 1, 0.1),
  beta.intervals.Y = seq(0, 1, 0.1),
  p.intervals.X = c(-18, -5, -2, 0),
  p.intervals.Y = c(-18, -5, -2, 0),
  return_with_reference = FALSE)
```
##### Arguments:
name | description
-|-
beta-value | beta-value matrix
detecP | detection p-value matrix
beta.intervals.X | Intervals of beta-value of X chromosome probes within each of which the proportion is calculated. Default is 'seq(0,1,0.1)', which is 0 to 0.1, 0.1 to 0.2, ..., 0.9 to 1.0.
beta.intervals.Y | Intervals of beta-value of Y chromosome probes within each of which the proportion is calculated. Default is the same as that of 'beta.intervals.X'
p.intervals.X | Intervals of detection p-value of X chromosome probes within each of which the proportion is calculated. Default is -Inf to 1e-18, 1e-18 to 1e-5, 1e-5 to 1e-2, 1e-2 to 0.
p.intervals.Y | Intervals of detection p-value of Y chromosome probes within each of which the proportion is calculated. Default is the same as that of 'p.intervals.X'
return_with_reference | Logical value to indicate whether or not to return the data derived from references data (Default to FALSE).


### 'plotSexEstimation'
##### Description:
This function generates a scatter plot using 'ggplot' showing the first principal component of two principal component analysis results (pincipal component analysis using the data from X chromosome and the data from Y chromosome)

##### Usage:
```{r}
plotSexEstimation(sex_estimation = NULL, samples = NULL,
  main = "", filename = NULL,
  include_reference = FALSE)
```

##### Arguments:
name | Description
-|-
sex_estimation | result of estimateSex function
samples | the list of samples to be plotted. If NULL (default), all samples in the estimateSex result are plotted.
main | character string to be used as the plot title. Default is an empty string.
filename | When specified, a PDF is generated.
include_reference | Logical value to specify whether to plot the reference data (Default is FALSE).

##### Examples:
```{r}
plotSexEstimation(sex_estimation = sexEst)
```

### 'plotSexDistribution'
##### Description:
This function generates a line chart using 'ggplot' showing the
distribution of beta-values of chrX and chrY and detection
p-values of chrY.

##### Usage:
```{r}
plotSexDistribution(beta.value = NULL, detecP = NULL,
  beta.intervals.X = seq(0, 1, 0.1),
  beta.intervals.Y = seq(0, 1, 0.1),
  p.intervals.X = c(-18, -5, -2, 0),
  p.intervals.Y = c(-18, -5, -2, 0),
  samples = NULL, filename = NULL,
  include_reference = TRUE,
  color = NULL)
```
##### Arguments:
name | description
-|-
beta-value | beta-value matrix
detecP | detection p-value matrix
beta.intervals.X | Intervals of beta-value of X chromosome probes within each of which the proportion is calculated. Default is 'seq(0,1,0.1)', which is 0 to 0.1, 0.1 to 0.2, ..., 0.9 to 1.0.
beta.intervals.Y | Intervals of beta-value of Y chromosome probes within each of which the proportion is calculated. Default is the same as that of 'beta.intervals.X'
p.intervals.X | Intervals of detection p-value of X chromosome probes within each of which the proportion is calculated. Default is -Inf to 1e-18m 1e-18 to 1e-5, 1e-5 to 1e-2, 1e-2 to 0.
p.intervals.Y | Intervals of detection p-value of Y chromosome probes within each of which the proportion is calculated. Default is the same as that of 'p.intervals.X'
samples | the list of samples to be plotted. If NULL (default) all samples in the beta.value matrix are plotted.
filename | When specified, a PDF is generated.
include_reference | Logical value to specify whether to plot reference data (Default is FALSE).
color | Logical value to specify whether to plot the reference data (Default is '#984ea3', one of the shades of purple).
main | character string to be used as the plot title. Default is an empty string.

##### Examples:
```{r}
plotSexDistribution(beta.value=beta, detecP=pval)
```


### 'get.proportion_table'
##### Description:
This function returns the beta-value and detection p-value proportion table within each intervals specified by parameters.

##### Usage:
```{r}
get.proportion_table(beta.value = NULL, detecP = NULL,
   beta.intervals.X = seq(0, 1, 0.1),
   beta.intervals.Y = seq(0, 1, 0.1),
   p.intervals.X = c(-18, -5, -2, 0),
   p.intervals.Y = c(-18, -5, -2, 0),
   samples = NULL)
```

##### Arguments:
name | description
-|-
beta-value | beta-value matrix
detecP | detection p-value matrix
beta.intervals.X |Intervals of beta-value of X chromosome probes within each of which the proportion is calculated. Default is 'seq(0,1,0.1)', which is 0 to 0.1, 0.1 to 0.2, ..., 0.9 to 1.0.
beta.intervals.Y | Intervals of beta-value of Y chromosome probes within each of which the proportion is calculated. Default is the same as that of 'beta.intervals.X'
p.intervals.X | Intervals of detection p-value of X chromosome probes within each of which the proportion is calculated. Default is -Inf to 1e-18m 1e-18 to 1e-5, 1e-5 to 1e-2, 1e-2 to 0.
p.intervals.Y | Intervals of detection p-value of Y chromosome probes within each of which the proportion is calculated. Default is the same as that of 'p.intervals.X'
samples | list of samples to calculate beta/detection p-value proportions. If NULL (default) all samples in the beta-value and detection p-value matrices are used.

##### Examples:
```{r}
proportion_table <- get.proportion_table(beta.value=beta,
  detecP=detecP)
```

# sEst
An R-package for estimating the sex from the Infinium HumanMethylation 450k microarray data

## Installation
### 1) Click on 'sest.tar' file.

### 2) Click on the 'Download' button.

### 3) Install the package from source:

#### a) In command line: 
```
$> R CMD INSTALL '/path/to/sest.tar'
```
#### b) Or, in R environment: 
```{r}
> install.packages(pkgs="/path/to/sest.tar", repos=NULL, type="source")
```


## Example

loading thelibrary and data.
```{r}
library(sest)

load("../data_preparation/beta.bg.XY.rda")
load("../data_preparation/p.XY.rda")
load("../data_preparation/pd.rda")
load("../data_preparation/QC.rda")
```


# sEst example

### loading the library and data (GSE51032)

```{r}
library(sest)
```

### Sample QC, using the per-sample call-rate of chrX probes.

```{r}
samples.GSE51032 <- rownames(pd)
samples.GSE51032.good <- samples.GSE51032[QC[samples.GSE51032, "callrate.X"]>0.95]
```

### Running sex-estimation
```{r}
sex.GSE51032.good <- estimateSex(beta.value=beta.bg.XY[, samples.GSE51032.good], detecP=p.XY[, samples.GSE51032.good], return_with_reference=TRUE)
```

### Contents of the 'estimateSex' result:
'estimateSex' returns a list object containing two data.frames: one is for reference data ('reference') and the other is for the test samples ('test')
```{r}
names(sex.GSE51032.good)
[1] "test"      "reference"
```
```{r}
head(sex.GSE51032.good$reference)
```
| X.PC1 |Y.PC1 | predicted.X | predicted.Y | predicted
-|------|------|-------------|-------------|-----------
male.0001	| -0.2537145 | -0.6624954 |	M	| M	| M
male.0002	| -0.2798323 | -0.6596785 |	M |	M |	M
male.0003	| -0.2161992 | -0.6539044 |	M |	M |	M
male.0004	| -0.2597709 | -0.6551377 |	M |	M |	M
male.0005	| -0.2542591 | -0.6545710 |	M |	M |	M
male.0006	| -0.2634860 | -0.6528368 |	M |	M |	M
... | ... | ... | ... | ... | ...

```{r}
head(sex.GSE51032.good$test)
```
| X.PC1 |Y.PC1 | predicted.X | predicted.Y | predicted
-|------|------|-------------|-------------|-----------
6969568099_R02C02	| 0.1490050	| 0.3504673	| F	| F	| F
6969568052_R02C01	| 0.1546120	| 0.4023443	| F	| F	| F
6969568141_R03C01	| 0.1634140	| 0.2984963	| F	| F	| F
6969568052_R01C02 |	0.1671693	| 0.3757989	| F	| F	| F
7766130010_R05C01 | 0.1484081	| 0.5155788	| F	| F	| F
6969568141_R06C02	| 0.1599462	| 0.4135742	| F	| F	| F
... | ... | ... | ... | ... | ...

### Plot of sex-estimation result by 'plotSexEstimation':
```{r}
plotSexEstimation(sex.GSE51032.good)
```
![plotSexEstimation without reference](figure/unnamed-chunk-6-1.png)

### Plot of sex-estimation result by 'plotSexEstimation', including reference samples:
```{r}
plotSexEstimation(sex.GSE51032.good, include_reference = TRUE)
```

Distribution profile plot of Male-predicted samples by 'plotSexDistribution':
```{r}
samples.GSE51032.good.M <- rownames(sex.GSE51032.good$test)[sex.GSE51032.good$test$predicted == "M"]
plotSexDistribution(beta.value=beta.bg.XY, detecP=p.XY, samples=samples.GSE51032.good.M)
```
![plotSexDistribution, male-predicted samples](figure/unnamed-chunk-7-1.png)

Distribution profile plot of Female-predicted samples by 'plotSexDistribution':
```{r}
samples.GSE51032.good.F <- rownames(sex.GSE51032.good$test)[sex.GSE51032.good$test$predicted == "F"]
plotSexDistribution(beta.value=beta.bg.XY, detecP=p.XY, samples=samples.GSE51032.good.F)
```
![plotSexDistribution, female-predicted samples](figure/unnamed-chunk-8-1.png)

Distribution profile plot of N-predicted samples by 'plotSexDistribution':
```{r}
samples.GSE51032.good.N <- rownames(sex.GSE51032.good$test)[sex.GSE51032.good$test$predicted == "N"]
plotSexDistribution(beta.value=beta.bg.XY, detecP=p.XY, samples=samples.GSE51032.good.N)
```
![plotSexDistribution, N-predicted samples](figure/unnamed-chunk-9-1.png)



Trying the sex-estimation using different beta-value intervals
```{r}
sex.GSE51032.good.5 <- estimateSex(beta.value=beta.bg.XY[, samples.GSE51032.good], detecP=p.XY[, samples.GSE51032.good], return_with_reference=TRUE, beta.intervals.X = seq(0,1,0.2), beta.intervals.Y = seq(0,1,0.2))

plotSexEstimation(sex.GSE51032.good.5)
```
![plotSexEstimation, using wider beta-value intervals](figure/unnamed-chunk-10-1.png)

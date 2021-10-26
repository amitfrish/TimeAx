The TimeAx algorithm allows performing multiple trajectory alignment (MTA) on time-series datasets of individuals each of which is considered as an individual partial trajectory

## TimeAx package installation and code requirements
TimeAx package can be downloaded from github. Please make sure you have the devtools package installed prior to TimeAx installation.

```R
library(devtools)
install_github("amitfrish/TimeAx")
```

## Training a TimeAx model

#GEData: 
A matrix containing profiles (columns) of omics measurments (rows) from multiple individuals and different time points. Profiles for each individual should be ordered by chronological time.
#sampleNames 
A vector containing the individual identity of each sample in the GEData.
#numOfIter: 
Number of consensus trajectories. The default is 100.
#numOfTopGenes: 
Length of the conserved-dynamics-seed of features. The default is 50.
#seed
The conserved-dynamics-seed. If provided, the alignment process will be conducted based on these features. The default is NULL.
#no_cores:
A number for the amount of cores which will be used for the analysis. The defalt (NULL) is total number of cores minus 1.

```R

library(TimeAx)
data(UBCData)

model = modelCreation(DataUBC,UBCSamples)

# Limiting the number of cores used for the training:
model = modelCreation(DataUBC,UBCSamples, no_cores = 2)

```

## The healthy trajectory
We will first characterize the healthy monocytic development in healthy individuals through a developmental trajectory analysis.
The fcs files of the 9 healthy individuals are deposited in the "Data" folder inside this repository ("fcsFilesH6Pops" directory). To obtain cells assigned to only sorted cell populations with known markers' expression patterns, we gated the .fcs files, such that each healthy individual has 6 separate files, each includes cells related to the following cell populations: erythroblasts, mDCs, monoblasts, monocytes, pDCs, and promonocytes. To read these files in a convenient way as an expression set, we wrote the function "readFCSeset" (can be found in the "funcTutorial.R" script) that reads the .fcs data and stores the information about the population name in the phenotypic data.
Specific details about the function and its attributes can be found in the functions R script.


```R
fcsFilesHDir = file.path(dataDir,'fcsFilesH6Pops')
fileNames = c('036','049','050','067','092','132','138','247','277')
esetH = readFCSeset(fileNames, fcsFilesHDir, transform = T)
```

We will first visualize the monocytic developmental axis by projection onto a PCA space. The developmental markers were chosen as those markers with known association to the developmental process of monocytic maturation: CD34, CD117, CD33, CD64, CD14, CD13 and CD11b. Of note, before applying the PCA we first sample the cells with inverse relation to the frequency of their assigned population to obtain a uniform representation of lal the monocytic cell populations.  


```R
#subset the healthy expressionset to include only those cell subsets associated with the monocytic developmental process:
devPops = c('monoblasts','promonocytes','monocytes')
devMarkers = c('CD34','CD117','CD33','CD64','CD13','CD14','CD11b')
esetHDev = esetH[devMarkers,esetH$popName %in% devPops]

#sample cells with inverse relation to the population's frequency:
sampProb = table(esetHDev$popName)/ncol(esetHDev)
esetHDevSamp = esetHDev[,sample(1:ncol(esetHDev), size = 10000, prob = sapply(esetHDev$popName, function(pop){return(1/sampProb[pop])}))]

#apply PCA:
pcaRes = prcomp(t(exprs(esetHDevSamp)), center = T)$x[,c('PC1','PC2')]
pcaResDf = data.frame(pcaRes, popName = esetHDevSamp$popName)
pcaResDf$popName = factor(pcaResDf$popName, levels = c('monoblasts','promonocytes','monocytes'))
ggplot(pcaResDf, aes(x = PC1, y = PC2, color = popName)) + geom_point() + theme_classic() +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  ggtitle('trajectory of sampled data with inverse relation to population frequency')
```
![alternativetext](Figures/Fig4b.PNG)


To obtain the pseudotime values of the healthy cells, we utilize the psupertime algorithm with the cellular assignment into cellular populaitons. The function "BuildDevTrajHealthyPops" located in the *funcTutorial.R* script that runs the psupertime algorithm with the cellular annotations to cell populations. The pseudotime values are stored in the psuperScaled attribute of the healthy expression set. 


```R
esetHDev = BuildDevTrajHealthyPops(esetHDev, markers = devMarkers, orderedPops = c('monoblasts','promonocytes', 'monocytes'))
```

We will now characterize the obtained healthy trajectory in 3 ways (Supp. Fig. 7):

1. Pseudotime values of each cell population along the trajectory. The provided function: "plotDev" provided in the *funcTutorial.R* plots the distribution of each cell population along the pseudotime axis.

2. Markers' expression dynamics along the pseudotime axis. The provided function "cellAlignInter" provided in the *funcTutorial.R* script applies interpolation and scaling of markers' expression levels along the pseudotime axis.

3. Cellular density along the pseudotime trajectory across individuals.


```R
#Plot the developmental populations ordering along the trajectory:
plotDev(esetHDev, method = 'psuperScaled', orderedPops = c('monoblasts','promonocytes', 'monocytes'))

#markers expression dynamics:
interScaledH = cellAlignInter(esetHDev, markers = devMarkers, method = 'psuperScaled')
pheatmap(interScaledH$scaledData, cluster_rows = F, cluster_cols = F,
         main = 'markers expression dynamics along the healthy trajectory',
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100))

#density along the trajectory:
ptXrossInd = data.frame(ind = esetHDev$fileName, pt = esetHDev$psuperScaled)
ggplot() + geom_density(ptXrossInd, mapping = aes(x = pt, y = ..density.., group = ind), color = 'grey') +
  geom_density(ptXrossInd, mapping = aes(x = pt, y = ..density..)) +
  theme_classic() +
  ggtitle('cellular density along the healthy trajectory across patients')

```

![alternativetext](Figures/SuppFig7a.PNG)
![alternativetext](Figures/SuppFig7b.PNG)
![alternativetext](Figures/SuppFig7c.PNG)

## Generation of a developmental classifier for myeloid cell populations
We will now generate the pseudotime trajectories per AML patient in a supervised trajectory approach that first classifies the cells into cell populations by their phenotypic similarity to the healthy populations and then assemble a trajectory that preserves the cell population ordering per sample. 
To classify the AML cells into these populations, we will use the developmental classifier proposed by Good et al which applies nearest neighbour classification using Mahalanobis distance metric. For this, we should first define per healthy population the markers' mean expression along with their covariance matrix across cells using the markers used for classification:


```R
classMarkers = c('CD34','CD123','CD33','CD64','CD11c','CD14','CD13','CD117','CD71','CD11b')

#define popMeans - a list of elements, in each the mean expression value of the markers used for classification:
popMeans = lapply(unique(esetH$popName), function(pop){
  return(apply(exprs(esetH)[classMarkers,esetH$popName == pop],1,mean))
})
names(popMeans) = unique(esetH$popName)

#define popMeans - a list of elements, in each the covariance matrix of the markers used for classification:
popSx = lapply(unique(esetH$popName), function(pop){
  covRes = cov(t(exprs(esetH)[classMarkers,esetH$popName == pop]))
  return(covRes)
})
names(popSx) = unique(esetH$popName)
```

The code for the classifier (function: devMapperAML located in the *funcTutorial.R*) along with details regarding the requied input parameters can be found in the script, and is based on code recieved from the authors of Good et al. 
We will test the performance of our classifier on the healthy data by calculating the confusion matrix using the real and predicted cell types, as obtained by manual gating and the developmental classifier, respectively (Supp. Fig. 9a).


```R
esetH = esetH[classMarkers,]

#classify the cells based on the data:
esetH$predCellTypes = devMapperAML(transEsetBQuery = esetH, refMarkers = classMarkers, popMeans, popSx)
confMatRes = caret::confusionMatrix(data = factor(esetH$predCellTypes, levels = unique(esetH$popName)),
                                    reference = factor(esetH$popName, levels = unique(esetH$popName)))
confMatRes$overall['Accuracy']
confMatRes$table
NormedConfMat = apply(confMatRes$table,2,function(x){return(x/sum(x))})
pheatmap(NormedConfMat, cluster_cols = F, cluster_rows = F,
         main = 'confusion matrix\nreference labels on columns and prediction on rows')
```

![alternativetext](Figures/SuppFig9a.PNG)


## Developmental trajectory assembly for AML samples

The gated AML samples are deposited in this repository. Assembly of the developmental trajectory per AML sample is consisted of 3 steps:
- Read the FCS sample to an expression-set: This functionality is provided by the function "readFCSeset" mentioned above, that appears in the *funcTutorial.R* script.
- Classification into developmental cell populations by using the developmental classifier "devMapperAML" mentioned above, that appears in the *funcTutorial.R* script.
- Using only the developmental cell populations (monoblasts, promonocytes and monocytes) to generate a supervised trajectory using the function: "BuildDevTrajHealthyPops"  mentioned above, that appears in the *funcTutorial.R* script.

Of note, this part is computationally intensive given the size of the AML files that typically include thousands of cells. For convenience, we provide .rds files of the AML expressionsets along with the assembled trajectories generated by this code.


```R
#define markers required for classification and trajectory building, and the developmental cell populations:
classMarkers = c('CD34','CD123','CD33','CD64','CD11c','CD14','CD13','CD117','CD71','CD11b')
devMarkers = c('CD34','CD117','CD33','CD64','CD13','CD14','CD11b')
devPops = c('monoblasts','promonocytes','monocytes')

#Classify each cell to the closest population:
esetAML = readFCSeset(fileNames = sample, fileDir = AMLDir, transform = T)
esetAML$popName = devMapperAML(esetAML, classMarkers, popMeans, popSx)

#Generate a trajectory using the developmental cell populations:
esetAMLDev = esetAML[,esetAML$popName %in% devPops]
esetAMLDev = BuildDevTrajHealthyPops(esetAMLDev, devMarkers, orderedPops = devPops)
```

## The tuMap algorithm

### tuMap algorithm alignment principles
The tuMap algorithm takes as input cancer samples ordered along developmental trajectories assembled for each patient individually and an averaged healthy trajectory that serves as a reference backbone. To overcome the inconsistencies in markers expression that characterize cancer samples, tuMap applies a weighted alignment that attenuates the effect of markers with disrupted expression dynamics relative to the healthy while enhancing the effect of those markers with conserved expression. Specifically, for each cancer sample the tuMap algorithm first quantifies the overall similarity between the expression dynamics of single markers along the cancer relative to the averaged healthy trajectory by using the distance derived from single-marker alignments. Next, tuMap assigns weights for each marker that are inversely associated with this distance by using different transformation functions (Supp. Fig. 4b-c). These weights are used to calculate the pairwise dissimilarity matrix between the cancer sample’s trajectory and the averaged healthy trajectory, on which the optimal alignment of the two trajectories is identified. Finally, tuMap uses the resulting alignment to map single cells from the cancer sample to their approximate location along the healthy development to obtain tuMap pseudotime values which are similarly scaled across cancer samples, allowing for cross-samples comparison. tuMap further provides the typical mapping error along the tuMap trajectory that is inherent to the weighted alignment process (Fig. 2c, Methods). 

### Application of the tuMap algorithm on the AML samples
As mentioned above, for convenience we provide .rds files with the AML expressionsets along with the the obtained pseudotime values.  
As detailed in the manuscript, to calculate the markers' weights, the tuMap first aligns them independently to the averaged healthy trajectory, followed by application of a transformation function which converts the alignment distance to weights. The tuMap algorithm provides 2 transformation functions: signoidal and non sigmoidal and the user should specify the value of a constant parameter that affects the nature of the transformation function. 
For downstream analysis, we will store the AML expression sets along with the tuMap pseudotime values in a new directory: *AMLsamplesDevTrajwTumap*


```R
AMLwTrajDir = file.path(dataDir,'AMLSampleswTraj')
AMLwTrajFiles = list.files(AMLwTrajDir)
tuMapTrajDir = file.path(dataDir,'AMLsamplesDevTrajwTumap')
ifelse(!dir.exists(tuMapTrajDir), dir.create(tuMapTrajDir, FALSE))

markersWeightsXrossPatients = do.call('rbind', lapply(AMLwTrajFiles, function(file){
  print(file)
  esetAMLDev = readRDS(file = file.path(AMLwTrajDir,file))
  sampleID = str_extract(file, '[0-9]+\\_[0-9]+')

  # calculate the weighted dissimilarity matrix between our trajectory and the healthy:
  interScaledAML = cellAlignInter(esetAMLDev, devMarkers, method = 'psuperScaled')

  wDisMat = tumap::WeightedDissimilarityMatrix(x = interScaledAML$scaledData[devMarkers,],
                                               y = interScaledH$scaledData[devMarkers,], const = 1)
                                               
  #global alignment with the weighted distance matrix:
  gAlignWweights = cellAlign::globalAlign(wDisMat$dissimilarity_matrix, step.pattern = symmetric1)
  esetAMLDev$tuMapPt = tumap::pseudotimeScaling(ptCancer = esetAMLDev$psuperScaled, interScaledCancer = interScaledAML,
                                                interScaledHealthy = interScaledH, gAlign = gAlignWweights)

  saveRDS(esetAMLDev, file = file.path(tuMapTrajDir, paste0(sampleID, 'wtuMap.rds')))

  return(data.frame(sampleID = sampleID, t(wDisMat$weights)))
}))
rownames(markersWeightsXrossPatients) = markersWeightsXrossPatients$sampleID

#show the markers' weights as violin plots - Supp. Fig. 9b:
markersWeightsXrossPatients_melt = melt(markersWeightsXrossPatients, id.vars = 'sampleID')
meanWeights = sapply(devMarkers, function(marker){return(mean(markersWeightsXrossPatients_melt$value[markersWeightsXrossPatients_melt$variable == marker]))})
ordMarkers = names(meanWeights)[order(meanWeights)]
markersWeightsXrossPatients_melt$variable = factor(markersWeightsXrossPatients_melt$variable, levels = ordMarkers)
ggplot(markersWeightsXrossPatients_melt, aes(x = variable, y = value, fill = variable)) + geom_jitter(width = 0.1, color = 'black', size = 0.5) + geom_violin(alpha = 0.2) +
  theme_classic() + ggtitle('markers weights across samples')
  
#Organize the healthy data for tuMap alignment error estimation:
hInds = unique(esetHDev$fileName)
intDataHInd = lapply(hInds, function(hInd){
  esetHInd = esetHDev[devMarkers,esetHDev$fileName == hInd]
  esetHInd$psuperOrg = esetHInd$psuperScaled
  esetHInd = BuildDevTrajHealthyPops(esetHInd, devMarkers, orderedPops = c('monoblasts','promonocytes', 'monocytes'))
  return(list(trajInd = esetHInd$psuperScaled,
              trajAvH = esetHInd$psuperOrg,
              interScaled = cellAlignInter(esetHInd, devMarkers, method = 'psuperScaled', scale = T)))
})

#read the expression matrix for one AML sample:
esetAMLDev = readRDS(file = path_to_AML_file)

# calculate the weighted dissimilarity matrix between our trajectory and the healthy:
interScaledAML = cellAlignInter(esetAMLDev, devMarkers, method = 'psuperScaled')

wDisMat = tumap::WeightedDissimilarityMatrix(x = interScaledAML$scaledData[devMarkers,],
                                             y = interScaledH$scaledData[devMarkers,], const = 1)

tuMapError = tumap::estimatetuMapErr(wDisMat, intDataHInd, interScaledH)
```

![alternativetext](Figures/SuppFig9b.PNG)


## Cellular density along the tuMap trajectory analysis
Since the tuMap trajectory is similarly scaled across individuals, we can compare the cellular density along it in ALL patients relative to the healthy trajectory. For this, we will first calculate the cellular density along the healthy and tuMap pseudotime axes in the healthy and AML samples, respectively. Then we will quantify the deviation of the cellular density in AML patients as compared to the healthy using the Earth Movers' Distance metric. Thus, the resulting EMDs reflect the overall cellular maturation degree of the sample with samples with low EMD exhibiting a high abundance of mature cells and vice versa. We thus define the term *stemness index* of a sample as the EMD from the healthy cellular density distribution, allowing for a meaningful characterization of AML samples. We observed that the diagnosis and relapse samples exhibit high *stemness indices*, that were significantly reduced 14 days following induction therapy (paired t-test P = 0.04, 0.04, n = 19, 6 samples with paired day 14 and either diagnosis or relapse samples, respectively), reflecting a response to the induction therapy that was reversed with disease relapse and suggesting potential associations between the density distributions along the tuMap axis and clinical outcome (Fig. 4c). 
For this analysis, we provide the .csv file *samplesData.csv* (located in the clinicalData directory attached to this repository) with the samples' mapping with relevant details for each sample, including the patient's ID and time point.


```R
samplesData = read.csv(file = file.path(dataDir, 'clinicalData', 'samplesData.csv'), stringsAsFactors = F)
tuMapTrajFiles = list.files(tuMapTrajDir)

#calculate the cellular density along the tuMap trajectory across the AML samples:
tumapDensXrossSamples = do.call('rbind', lapply(tuMapTrajFiles, function(tuMapFile){
  esetAML = readRDS(file = file.path(tuMapTrajDir, tuMapFile))
  tumapDens = density(esetAML$tuMapPt, from = 0, to = 1, n = 128, adjust = 2)
  return(tumapDens$y)
}))
rownames(tumapDensXrossSamples) = str_extract(tuMapTrajFiles,'[0-9]+\\_[0-9]+')

#calculate the healthy cellular density along the pseudotime axis:
hDens = density(esetHDev$psuperScaled, from = 0, to = 1, n = 128, adjust = 2)

#calculate per sample the EMD distance of cellular density from the healthy:
EMDtuMapH = sapply(tuMapTrajFiles, function(tuMapFile){
  esetAML = readRDS(file = file.path(tuMapTrajDir, tuMapFile))
  tumapDens = density(esetAML$tuMapPt, from = 0, to = 1, n = 128, adjust = 2)
  return(emdist::emd(cbind(tumapDens$y, tumapDens$x),
                     cbind(hDens$y, hDens$x)))
})
names(EMDtuMapH) = str_extract(names(EMDtuMapH),'[0-9]+\\_[0-9]+')

#show the dynamics of EMDH of patient2 over time:
EMDHDynamicsPt2 = do.call('rbind', lapply(samplesData$SampleID[samplesData$Patient == 2], function(sampleID){
  return(data.frame(patient = 2, sampleID = sampleID,
                    timePoint = samplesData$Timepoint[samplesData$SampleID == sampleID],
                    EMDH = EMDtuMapH[sampleID]))
}))
EMDHDynamicsPt2$timePoint = factor(EMDHDynamicsPt2$timePoint, levels = c('d0','d14','rem','rel'))
ggplot(EMDHDynamicsPt2, aes(x = timePoint, y = EMDH, group = patient)) + geom_line(linetype = 'dashed') + geom_point() +
  theme_classic() + ggtitle('EMDH of longitudinal samples of patient 2')

#order the individuals by the EMD distance from the healthy:
ordInd = names(EMDtuMapH)[order(-EMDtuMapH)]
tumapDensXrossSamplesOrd = tumapDensXrossSamples[ordInd,]

#add annotations per sample describing the time point in which the sample was taken:
annotRows = do.call('rbind', lapply(rownames(tumapDensXrossSamplesOrd), function(sampleID){
  return(data.frame(patient = samplesData$Patient[samplesData$SampleID == sampleID],
                    timePoint = samplesData$Timepoint[samplesData$SampleID == sampleID],
                    EMDH = EMDtuMapH[sampleID]))
}))
rownames(annotRows) = rownames(tumapDensXrossSamplesOrd)

#add the healthy density at the bottom of the table - Supp. Fig. 9c:
tumapDensXrossSamplesOrd = rbind(tumapDensXrossSamplesOrd, hDens$y)
rownames(tumapDensXrossSamplesOrd)[nrow(tumapDensXrossSamplesOrd)] = 'H'
annotRows = rbind(annotRows, data.frame(patient = 'H', timePoint = 'H', EMDH = 0))
rownames(annotRows)[nrow(annotRows)] = 'H'
pheatmap(tumapDensXrossSamplesOrd, cluster_cols = F, cluster_rows = F, border_color = NA,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(100), show_rownames = F, gaps_row = nrow(tumapDensXrossSamplesOrd) - 1,
         main = 'cellular density along the tuMap pseudotime axis', annotation_row = annotRows[,c('timePoint', 'EMDH'), drop = F])

#paired comparison between d0 and d14 healthy EMD - Figure 4c:
EMDHD0D14 = dcast(annotRows, formula = patient ~ timePoint, value.var = 'EMDH')
EMDHD0D14 = EMDHD0D14[!is.na(EMDHD0D14$d14),c('d0','d14')]
t.test(EMDHD0D14$d0, EMDHD0D14$d14, paired = T)

#paired comparison between relapse and d14 healthy EMD:
EMDHD0Rel = dcast(annotRows, formula = patient ~ timePoint, value.var = 'EMDH')
EMDHD0Rel = EMDHD0Rel[!is.na(EMDHD0Rel$d14) & !is.na(EMDHD0Rel$rel), c('d14','rel')]
t.test(EMDHD0Rel$rel, EMDHD0Rel$d14, paired = T)

#boxplot of all the data:
ggplot(subset(annotRows, timePoint %in% c('d0','d14','rel')), aes(x = timePoint, y = EMDH)) + geom_boxplot() + theme_classic() +
  ggtitle('EMD from the healthy distribution along tuMap trajectory')
```
![alternativetext](Figures/SuppFig9c.PNG)

![alternativetext](Figures/Fig4c.PNG)


## Correlation of cellular density along the tuMap trajectory with outcome (overall survival)
To identify associations between the cellular density along the tuMap pseudotime axis and clinical outcome, we calculated a regularized cox proportional hazard model correlating patients’ overall survival with age at AML diagnosis, gender, and cellular density along the tuMap pseudotime axis at diagnosis and 14 days following initiation of induction therapy (Methods). Stratification of patients by the median risk value as predicted by the model yielded two groups whose survival rates significantly differed (Fig. 4d, log-rank test P = 1.78∙10^-2 , n = 19 patients with matched diagnosis and day-14 samples). Apart of age at AML diagnosis, the optimal survival model used cellular density along the tuMap trajectory at day 14 and not day 0 to predict survival, emphasizing the clinical importance in post-treatment patients’ monitoring for outcome prediction. Studying the cellular density distributions along the tuMap pseudotime axis revealed that patients assigned to the low-risk group exhibited an enrichment of mature cells following therapy, in contrast to those patients assigned to the high-risk group, reflecting the desired response to therapy (Supp. Fig. 9d). To compare the predictive performance of the tuMap approach to the one using assignment into discrete cell populations, we calculated a similar survival model using the frequencies of the main monocytic cell populations (monoblasts, promonocytes and mature monocytes). The resulting survival model yielded less accurate stratification into groups (log rank test P = 0.0018, n = 19 patients with matched diagnosis and day-14 samples), highlighting the added value obtained by the continuous characterization of the involved cell lineage in clinical predictions. 
To run this code we provide an additional .csv file named: *patientsData.csv* located in the clinicalData directory that includes the phenotypic data about the patients (age, survival and followup time). 


```R
patientsData = read.csv(file = file.path(dataDir,clinicalData,'patientsData.csv'), stringsAsFactors = F)

#the Surv() function requires the data to be organized as: 0-alive, 1 - dead (https://www.rdocumentation.org/packages/survival/versions/2.11-4/topics/Surv)
patientsData$status.alive = 1 - patientsData$is.alive

#get only those patients with d14 samples:
samplesData = samplesData[samplesData$SampleID %in% rownames(tumapDensXrossSamples),]
d14Patients = samplesData$Patient[samplesData$Timepoint == 'd14']
patientsDataD14 = patientsData[d14Patients,]

#organize predictors data:
densityD0D14 = do.call('rbind', lapply(d14Patients, function(pt){
  #diagnosis age:
  diagAge = patientsData$age.diagnosis[patientsData$PatientID == pt]
  #gender:
  gender = as.numeric(patientsData$Gender[patientsData$PatientID == pt] == 'F')
  #sample ID of the sample taken at day 0 and its cellular density along the tuMap axis:
  d0Sample = samplesData$SampleID[samplesData$Patient == pt & samplesData$Timepoint == 'd0']
  d0Density = tumapDensXrossSamples[d0Sample,]
  #sample ID of the sample taken at day 14 and its cellular density along the tuMap axis:
  d14Sample = samplesData$SampleID[samplesData$Patient == pt & samplesData$Timepoint == 'd14']
  d14Density = tumapDensXrossSamples[d14Sample,]
  #concatenate the data as a vector:
  predPt = c(diagAge, gender, d0Density, d14Density)
  names(predPt) = c('diagAge', 'gender', paste0('d0_',1:length(d0Density)), paste0('d14_',1:length(d14Density)))
  return(predPt)
}))
rownames(densityD0D14) = paste0('pt',d14Patients)

#generate a glmnet regularized COX model:
set.seed(1)
cv.fitTraj <- cv.glmnet(x = densityD0D14,
                        Surv(patientsDataD14[,'followup.alive'], patientsDataD14[,'status.alive']), family="cox",
                        alpha = 1, nfolds = 5)

#choose the best lambda:
best.lambdaTraj = cv.fitTraj$lambda.min

#get the coefficients of the model with the best lambda:
CoefficientsTraj <- coef(cv.fitTraj, s = best.lambdaTraj)
Active.IndexTraj <- which(CoefficientsTraj != 0)
Active.CoefficientsTraj <- CoefficientsTraj[Active.IndexTraj]
coefBestLambdaTraj = data.frame(fet = colnames(densityD0D14)[Active.IndexTraj],
                                coef = CoefficientsTraj[Active.IndexTraj])

#fit a full model with all the data
lasso.modTraj <- glmnet(x = densityD0D14,
                        Surv(patientsDataD14[,'followup.alive'], patientsDataD14[,'status.alive']), family="cox",
                        alpha = 1, lambda = best.lambdaTraj)

#generate a Kaplan Meier survival curve for the best lambda:
lpTraj = predict(lasso.modTraj, newx = densityD0D14, s = best.lambdaTraj)

#put a threshold on the median lpAll:
threshTraj = median(lpTraj[,1])

#calculate log-rank test for 2 groups stratified by the median risk:
patientsDataD14$groupsModelTraj = (lpTraj[,1] <= threshTraj)

# plot the Kaplam Meier plot (Figure 4d)
survDiffObj = survdiff(Surv(patientsDataD14[,'followup.alive'], patientsDataD14[,'status.alive']) ~ patientsDataD14$groupsModelTraj)
pv = 1 - pchisq(survDiffObj$chisq, length(survDiffObj$n) - 1)
fit = survfit(Surv(patientsDataD14[,'followup.alive'], patientsDataD14[,'status.alive']) ~ patientsDataD14$groupsModelTraj)
ggsurvplot(fit, data = patientsDataD14, pval = TRUE)

#show the distribution of d14 cellular density along the trajectory for the low vs. high risk - Supp. Fig. 9d:
densityD14 = densityD0D14[,str_detect(colnames(densityD0D14),'d14')]
densityD14_m = melt(densityD14)
densityD14_m$ind = as.numeric(str_replace(densityD14_m$Var2,'d14_',''))
densityD14_m$traj = (densityD14_m$ind - 1)*(1/127)
densityD14_m$riskGrp = sapply(densityD14_m$Var1, function(ptid){
  return(patientsDataD14$groupsModelTraj[patientsDataD14$PatientID == str_extract(ptid,'[0-9]+')])})
ggplot(densityD14_m, aes(x = traj, y = value, group = Var1, color = riskGrp)) + geom_line() + theme_classic() +
  ggtitle('day 14 density distributions for the high vs. low risk patients')
```

![alternativetext](Figures/Fig4d.PNG)

![alternativetext](Figures/SuppFig9d.PNG)

The TimeAx algorithm allows performing multiple trajectory alignment (MTA) on time-series datasets of individuals each of which is considered as an individual partial trajectory

## TimeAx package installation and code requirements
TimeAx package can be downloaded from github. Please make sure you have the devtools package installed prior to TimeAx installation.

```R
library(devtools)
install_github("amitfrish/TimeAx")
```

## Training a TimeAx model
The user should first train a TimeAx model based on a any kind of logitudinal data of the biological process, with at least 3 samples in each individual trajectory. The model will be later used to infer the pseudotime positions of each sample.

### GEData: 
A matrix containing profiles (columns) of omics measurments (rows) from multiple individuals and different time points. Profiles for each individual should be ordered by chronological time.
### sampleNames 
A vector containing the individual identity of each sample in the GEData.
### numOfIter: 
Number of consensus trajectories. The default is 100.
### numOfTopGenes: 
Length of the conserved-dynamics-seed of features. The default is 50.
### seed
The conserved-dynamics-seed. If provided, the alignment process will be conducted based on these features. The default is NULL.
### no_cores:
A number for the amount of cores which will be used for the analysis. The defalt (NULL) is total number of cores minus 1.

```R

library(TimeAx)
data(UBCData)

model = modelCreation(DataUBC,UBCSamples)

# Limiting the number of cores used for the training:
model = modelCreation(DataUBC,UBCSamples, no_cores = 2)

```

## Inferring pseudotime
Based on the TimeAx model, the user can infer the pseudotime position of each sample, assuming its profile includes the same features as the train data. The output of this step is a list containing the pseudotime positions of each sample (predictions) and it's equivilant certainty score (certainty).

### model:
A TimeAx model.
### GEData: 
A matrix containing profiles (columns) of omics measurments (rows).
### sampleNames 
Used for the robustness analysis. Always keep as NULL.
### no_cores:
A number for the amount of cores which will be used for the analysis. The defalt (NULL) is total number of cores minus 1.
### seed
The conserved-dynamics-seed. If provided, the prediction process will be conducted based on these features. Use the model's seed by keeping the the default value of NULL.
### batchCorrect
Whether to correct the new samples based on the consensus trajectory. The defualt is TRUE.

```R

library(TimeAx)

pseudotimeStats = predictByConsensus(model,DataUBC)
pseudotime = pseudotimeStats$predictions
certainty = pseudotimeStats$certainty

```

## Robustness analysis
Calculates a robustness score for the TimeAx model. High robustness score implies that the model indeed captures a biological process that changes over time. On the other hand, a low robustness score suggests that the model fails to represent a continuous process over time.

### model:
A TimeAx model.
### GEData: 
The matrix containing profiles (columns) of omics measurments (rows), which was used to train the model.
### sampleNames 
A vector containing the individual identity of each sample in the GEData. Same vector as used in the training.
### no_cores:
A number for the amount of cores which will be used for the analysis. The defalt (NULL) is total number of cores minus 1.
### pseudo
The output list of predictByConsensus. If not provided (NULL), pseudotime will be inferred by this function.

```R

library(TimeAx)

pseudotimeStats = predictByConsensus(model,DataUBC)
pseudotime = pseudotimeStats$predictions
certainty = pseudotimeStats$certainty

```

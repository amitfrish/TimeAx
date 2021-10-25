########## Building patient data
#' @keywords internal
dataCreation = function(GEData, sampleNames){
  lapply(unique(sampleNames), function(currSample){
    selectedIndexes = which(sampleNames == currSample)
    selectedTrajectoryBase = 1:length(selectedIndexes)
    trajectory= selectedTrajectoryBase/max(selectedTrajectoryBase)
    currGE = GEData[,selectedIndexes]
    baseData = currGE
    GEDataNorm = t(apply(currGE, 1, function(x){
      minVal = min(x)
      maxVal = max(x)
      (x - minVal)/(maxVal - minVal)
    }))
    list(scaledData = GEDataNorm, traj = trajectory, baseData = baseData, name = currSample, type = "Sample")
  })
}

########## Calculation of consensus trajectory
#' @keywords internal
createTrajectoryFromData = function(sampleData, sampleTraj, trajGenes){
  pcaOfSample = stats::prcomp(t(sampleData[trajGenes,]),scale. = T, center = T)
  distPCA = fields::rdist(pcaOfSample$x[,1:which(summary(pcaOfSample)$importance[3,]>0.9)[1]])
  sampleFixedTraj = cumsum(c(0,unlist(lapply(2:length(sampleTraj),function(i){
    abs(distPCA[i,i-1])
  }))))

  (sampleFixedTraj-min(sampleFixedTraj))/(max(sampleFixedTraj)-min(sampleFixedTraj))
}

########## Data smoothing
#' @keywords internal
computeNewData = function(sampleTraj, trajCond, dataToTransorm,winSz){
  lapply(1:length(sampleTraj), function(i) {
    dist2Others = trajCond - sampleTraj[i]
    weightedData = exp(-(dist2Others^2)/(winSz^2))
    weightedData = weightedData/sum(weightedData)
    dataToTransorm %*% weightedData
  })
}

########## Patient weigths for the alignment
#' @keywords internal
calculateSampleWeights = function(cors){
  relevantCors = diag(cors)
  abs(mean(relevantCors[!is.na(relevantCors)]))
}

########## Pairwise alignment
#' @keywords internal
getAlignment = function(sample1, sample2, trajGenes){
  mutualGenes = intersect(intersect(row.names(sample1$scaledData),row.names(sample2$scaledData)),trajGenes)
    cellAlign::globalAlign(sample1$scaledData[mutualGenes,], sample2$scaledData[mutualGenes,],
                                 scores = list(query = sample1$traj,
                                               ref = sample2$traj),
                                 sigCalc = F, numPerm = 20, verbose = F)
}

########## Model training (inner)
#' @keywords internal
multiAlign = function(listOfSamplesSmall, trajGenes, numOfIter, no_cores){
  #### All pairwise alignments ####
  if(is.null(no_cores)){
    no_cores = max(1, parallel::detectCores() - 1)
  }
  message('Calculating all pairwise alignments:')
  cl<-parallel::makeCluster(no_cores)
  parallel::clusterExport(cl=cl, varlist=c("getAlignment", "listOfSamplesSmall", "trajGenes"), envir=environment())
  doSNOW::registerDoSNOW(cl)
  pb <- utils::txtProgressBar(min = 1, max = length(listOfSamplesSmall), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  `%dopar2%` <- foreach::`%dopar%`
  firstSampleInd = NULL
  allAlignments <- foreach::foreach(firstSampleInd = 1:length(listOfSamplesSmall), .options.snow = opts) %dopar2% {
    firstSample = listOfSamplesSmall[[firstSampleInd]]
    alignmentList = lapply(listOfSamplesSmall, function(secondSample){
      getAlignment(firstSample,secondSample, trajGenes)
    })
    setTxtProgressBar(pb, firstSampleInd)
    alignmentList
  }
  parallel::stopCluster(cl)
  close(pb)

  #### Consensus list ####
  message('Creating consensus list:')
  cl<-parallel::makeCluster(no_cores)
  parallel::clusterExport(cl=cl, varlist=c("calculateSampleWeights", "listOfSamplesSmall","getAlignment", "trajGenes","createTrajectoryFromData","allAlignments"), envir=environment())
  doSNOW::registerDoSNOW(cl)
  pb <- utils::txtProgressBar(min = 1, max = numOfIter, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  `%dopar2%` <- foreach::`%dopar%`
  iterNum = NULL
  consensusList <- foreach::foreach(iterNum = 1:numOfIter, .options.snow = opts) %dopar2% {
    #phylTree = phangorn::upgma(distMatrix)
    phylTree = phangorn::upgma(fields::rdist(sample(1:length(listOfSamplesSmall))))

    #### Couple joining ####
    nodesInTheTree = 1:(phylTree$Nnode+1)
    listOfSamplesForTree = listOfSamplesSmall
    for(i in 1:length(listOfSamplesForTree)){
      listOfSamplesForTree[[i]]$sampleScore = mean(unlist(lapply(1:length(allAlignments[[i]]),function(j){
        currAlign = (allAlignments[[i]])[[j]]
        currAlignmentSteps = currAlign$align[[1]]
        #mutualGenes = intersect(intersect(row.names(listOfSamplesSmall[[i]]$scaledData),row.names(listOfSamplesSmall[[j]]$scaledData)),trajGenes)
        calculateSampleWeights(stats::cor(listOfSamplesSmall[[i]]$scaledData[,currAlignmentSteps$index1],
                                   listOfSamplesSmall[[j]]$scaledData[,currAlignmentSteps$index2]))
        #currAlign$normalizedDistance
      })))
    }

    message('Creating a concensus sample...')
    edgeInformation = cbind(phylTree$edge,phylTree$edge.length)
    allFatherNodes = edgeInformation[,1]

    while(length(unique(nodesInTheTree))>1){
      currEdgeInformation = edgeInformation[edgeInformation[,1] %in% unique(allFatherNodes),]
      fatherOfKnownNodes = unique(allFatherNodes)[which(unlist(lapply(unique(allFatherNodes),function(i){
        length(which(currEdgeInformation[allFatherNodes==i,2] %in% nodesInTheTree))==2
      })))]
      knownNodes = unlist(lapply(fatherOfKnownNodes,function(i){currEdgeInformation[allFatherNodes==i,2]}))
      for(currFather in fatherOfKnownNodes){
        currNodes = currEdgeInformation[allFatherNodes==currFather,2]
        currSample1 = listOfSamplesForTree[[currNodes[1]]]
        currSample2 = listOfSamplesForTree[[currNodes[2]]]

        refSample = currSample1
        secSample = currSample2
        # if(currSample2$sampleScore>currSample1$sampleScore){
        if(length(currSample1$traj)<length(currSample2$traj)){
          refSample = currSample2
          secSample = currSample1
          currNodes = rev(currNodes)
        }

        if(refSample$type=="Comb" | secSample$type=="Comb"){
          currAlignment = getAlignment(secSample,refSample,trajGenes)
        }else{
          currAlignment = (allAlignments[[currNodes[2]]])[[currNodes[1]]]
        }
        currAlignmentSteps = currAlignment$align[[1]]

        ## Weighting expression profiles
        # refSamplePriority = 0.5
        # secSamplePriority = 0.5
        refSamplePriority = refSample$sampleScore/(refSample$sampleScore+secSample$sampleScore)
        secSamplePriority = secSample$sampleScore/(refSample$sampleScore+secSample$sampleScore)

        newExp = refSample$scaledData[,currAlignmentSteps$index2]*refSamplePriority +
          secSample$scaledData[,currAlignmentSteps$index1]*secSamplePriority
        newTraj = seq(0,1,length.out = dim(newExp)[2])

        newBase = refSample$baseData[,currAlignmentSteps$index2]*refSamplePriority +
          secSample$baseData[,currAlignmentSteps$index1]*secSamplePriority

        fixedTraj = newTraj
        alignQuality = calculateSampleWeights(stats::cor(refSample$scaledData[,currAlignmentSteps$index2],
                                                  secSample$scaledData[,currAlignmentSteps$index1]))


        listOfSamplesForTree[[currFather]] = list(scaledData = newExp, traj = fixedTraj, name = currFather, type = "Comb", sampleScore = alignQuality, baseData = newBase)
      }
      nodesInTheTree = c(nodesInTheTree[!(nodesInTheTree %in% knownNodes)],fatherOfKnownNodes)
      allFatherNodes = allFatherNodes[!(allFatherNodes %in% fatherOfKnownNodes)]
    }

    #### Interpreting results ####
    fullAlignedSample = listOfSamplesForTree[[nodesInTheTree]]
    fullAlignedSample$traj = createTrajectoryFromData(fullAlignedSample$scaledData,fullAlignedSample$traj, trajGenes)
    #fullAlignedSample$traj = createTrajectoryFromData(fullAlignedSample$baseData,fullAlignedSample$traj, trajGenes)
    #fullAlignedSample$scaledData = do.call("cbind",computeNewData(fullAlignedSample$traj,fullAlignedSample$traj,fullAlignedSample$scaledData,0.1))

    setTxtProgressBar(pb, iterNum)
    fullAlignedSample
  }
  parallel::stopCluster(cl)
  close(pb)
  message('Model created')
  consensusList
}

#' Training the TimeAx model
#'
#' This function initiate model training using the TimeAx algorithm - performing a multiple trajectory alignment (MTA) on time-series datasets of individuals each of which is considered as an individual partial trajectory.
#'
#' @param GEData A matrix containing profiles (columns) of omics measurments (rows) from multiple individuals and different time points. Profiles for each individual should be ordered by chronological time.
#' @param sampleNames A vector containing the individual identity of each sample in the GEData.
#' @param numOfIter Number of consensus trajectories. The default is 100.
#' @param numOfTopGenes Length of the conserved-dynamics-seed of features. The default is 50.
#' @param seed The conserved-dynamics-seed. If provided, the alignment process will be conducted based on these features. The default is NULL.
#' @param no_cores A number for the amount of cores which will be used for the analysis. The defalt (NULL) is total number of cores minus 1.
#' @return A TimeAx model consists of:
#' \item{consensusList}{List of consensus trajectories.}
#' \item{seed}{The conserved-dynamics-seed that was used in the alignment process.}
#' @references
#' Submitted
#' @examples
#' data(UBCData)
#'
#' # Training the model
#' model = modelCreation(DataUBC,UBCSamples, no_cores = 2)
#'
#' \dontrun{
#'
#' }
#' @export
#' @importFrom "utils" "setTxtProgressBar"
#' @importFrom "stats" "sd" "var"
#' @importFrom "grDevices" "chull"
modelCreation = function(GEData, sampleNames, numOfIter = 100, numOfTopGenes = 50 ,seed = NULL, no_cores = NULL){
  listOfSamples = dataCreation(GEData, sampleNames)
  if(is.null(seed)){
    seed = detectTrajGenes(listOfSamples, sampleNames, numOfTopGenes = numOfTopGenes, no_cores = no_cores)
  }

  listOfSamplesSmall = lapply(listOfSamples,function(currSample){
    currSampleNew = currSample
    currSampleNew$scaledData = currSample$scaledData[seed,]
    currSampleNew$baseData = currSample$baseData[seed,]
    currSampleNew
  })

  consensusList = multiAlign(listOfSamplesSmall, seed, numOfIter = numOfIter, no_cores = no_cores)

  #### Output ####
  list(consensusList = consensusList, seed = seed)
}

#' Selecting conserved-dynamics-seed features
#'
#' @param GEData A matrix containing profiles (columns) of omics measurments (rows) from multiple individuals and different time points. Profiles for each individual should be ordered by chronological time.
#' @param sampleNames A vector containing the individual identity of each sample in the GEData.
#' @param numOfTopGenes Length of the conserved-dynamics-seed of features. The default is 50.
#' @param topGenes Number of initial high variable features to be considered for the seed selection. The default is 4000.
#' @param numOfIterations Number of different random sample selections for the calculation. The default is 20.
#' @param percOfSamples Fraction of samples from each individual, selected in each sample selection. The default is 0.8
#' @param no_cores A number for the amount of cores which will be used for the analysis. The default (NULL) is total number of cores minus 1.
#' @return A list including:
#' The conserved-dynamics-seed. A list of feature, suitable for the model training, ordered from best to worst.
#' @references
#' Submitted
#' @examples
#' data(UBCData)
#'
#' # Selecting conserved-dynamics-seed features
#' seed = detectTrajGenes(DataUBC,UBCSamples, no_cores = 2)
#' @export
detectTrajGenes = function(GEData, sampleNames, numOfTopGenes = 50, topGenes = 4000, numOfIterations = 20, percOfSamples = 0.8, no_cores = NULL){
  if(class(GEData)=='list'){
    listOfSamples = GEData
  }else{
    listOfSamples = dataCreation(GEData, sampleNames)
  }
  baseDataList = lapply(listOfSamples,function(x){x$baseData})
  minSize = round(min(unlist(lapply(baseDataList,function(x){dim(x)[2]})))*percOfSamples)

  message('Initial feature selection')
  mutualGeneNames = row.names(baseDataList[[1]])
  for(currSample in baseDataList){
    mutualGeneNames = intersect(mutualGeneNames,row.names(currSample))
  }
  if(length(mutualGeneNames)<(2*topGenes)){topGenes = floor(length(mutualGeneNames)/2)}
  dataForHighVar = as.data.frame(t(do.call(cbind,lapply(baseDataList,function(currSampleData){
    baseOriginal = currSampleData[mutualGeneNames,]
    t(apply(baseOriginal, 1, function(x){
      minVal = min(x)
      maxVal = max(x)
      (x - minVal)/(maxVal - minVal)
    }))
  }))))

  maxUniques = dim(dataForHighVar)[1]+2-2*length(baseDataList)
  minAllowedUniques = 0.8*maxUniques
  genesUniques = lengths(lapply(dataForHighVar, unique))
  selectedGeneIndexes = which(genesUniques>minAllowedUniques)
  if(length(selectedGeneIndexes)<2*topGenes){
    selectedGeneIndexes = order(genesUniques,decreasing = T)[1:(2*topGenes)]
  }
  genesSD = apply(dataForHighVar[,selectedGeneIndexes],2,sd)
  selectedGeneIndexes = selectedGeneIndexes[order(genesSD,decreasing = T)[1:topGenes]]
  mutualGeneNames = mutualGeneNames[selectedGeneIndexes]

  baseDataList = lapply(baseDataList,function(x){x[mutualGeneNames,]})

  message('Choosing conserved-dynamics-seed:')
  if(is.null(no_cores)){
    no_cores = min(numOfIterations,max(1, parallel::detectCores() - 1))
  }
  cl<-parallel::makeCluster(no_cores)
  parallel::clusterExport(cl=cl, varlist=c("mutualGeneNames", "minSize","baseDataList"), envir=environment())
  doSNOW::registerDoSNOW(cl)
  pb <- utils::txtProgressBar(min = 1, max = numOfIterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  `%dopar2%` <- foreach::`%dopar%`
  iteration = NULL
  geneScoreList = foreach::foreach(iteration = 1:numOfIterations, .options.snow = opts) %dopar2% {
    sampledData = t(do.call(rbind,lapply(baseDataList, function(currSample){
      selectedCells = sort(sample(1:dim(currSample)[2],minSize))
      newData = currSample[,selectedCells]
      colnames(newData) = 1:dim(newData)[2]
      newData
    })))

    setTxtProgressBar(pb, iteration)

    sapply(1:length(mutualGeneNames), function(i){
      currMatrix = sampledData[,seq(i,dim(sampledData)[2],length(mutualGeneNames))]
      mean(stats::cor(currMatrix,method = "spearman"))
    })
  }
  parallel::stopCluster(cl)
  close(pb)

  geneScore = colMeans(do.call(rbind,geneScoreList))
  mutualGeneNames[order(geneScore,decreasing = T)[1:numOfTopGenes]]
}

#' Infer pseudotime for new samples bassed on the TimeAx model
#'
#' @param model A TimeAx model.
#' @param GEData A matrix containing profiles (columns) of omics measurments (rows).
#' @param sampleNames Used for the robustness analysis. Always keep as NULL.
#' @param no_cores A number for the amount of cores which will be used for the analysis. The default (NULL) is total number of cores minus 1.
#' @param seed The conserved-dynamics-seed. If provided, the prediction process will be conducted based on these features. Use the model's seed by keeping the the default value of NULL.
#' @param batchCorrect Whether to correct the new samples based on the consensus trajectory. The defualt is TRUE.
#' @return A prediction list consists of:
#' \item{predictions}{The final pseudotime position for each sample.}
#' \item{certainty}{A certainty score for each position. Lower scores means higher certainty.}
#' @references
#' Submitted
#' @examples
#' data(UBCData)
#'
#' # Training the model
#' model = modelCreation(DataUBC,UBCSamples, no_cores = 2)
#'
#' # Inferring pseudotime positions
#' pseudotimeStats = predictByConsensus(model,DataUBC, no_cores = 2)
#' pseudotime = pseudotimeStats$predictions
#' certainty = pseudotimeStats$certainty
#' @export
predictByConsensus = function(model, GEData, sampleNames = NULL, no_cores = NULL, seed = NULL, batchCorrect = T){
  if(is.null(seed)){
    seed = intersect(model$seed, row.names(GEData))
  }
  # if(is.null(sampleNames)){
  #   sampleNames = 1:dim(GEData)[2]
  # }
  if(is.null(no_cores)){
    no_cores = max(1, parallel::detectCores() - 1)
  }

  cleanModel = lapply(model$consensusList,function(x){
    list(data = x$baseData[seed,], traj = x$traj)
  })
  GEData = GEData[seed,,drop = F]

  if(is.null(sampleNames)){
    message('Predicting samples pseudotime positions:')
  }else{
    message('Predicting robustness pseudotime positions:')
  }
  cl<-parallel::makeCluster(no_cores)
  parallel::clusterExport(cl=cl, varlist=c("seed", "GEData","sampleNames","computeNewData","cleanModel"), envir=environment())
  doSNOW::registerDoSNOW(cl)
  pb <- utils::txtProgressBar(min = 1, max = length(cleanModel), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  `%dopar2%` <- foreach::`%dopar%`
  consensusInd = NULL
  predictionStats <- foreach::foreach(consensusInd = 1:length(cleanModel), .options.snow = opts) %dopar2% {
    sampleConsensus = cleanModel[[consensusInd]]
    ref = sampleConsensus$data
    combinedData = cbind(ref, GEData)
    typeVector = c(rep(0, dim(ref)[2]),rep(1,dim(GEData)[2]))
    dataAfterCombat = combinedData
    if(dim(GEData)[2]>2 & batchCorrect){
      dataAfterCombat = sva::ComBat(combinedData, batch = typeVector)
    }
    refNorm = dataAfterCombat[,typeVector==0]
    testDataNorm = dataAfterCombat[,typeVector==1]

    if(is.null(sampleNames)){
      corMatrix = stats::cor(refNorm,testDataNorm,method = "spearman")
      corMatrix = apply(as.matrix(corMatrix),2,function(x){
        unlist(computeNewData(sampleConsensus$traj,sampleConsensus$traj,t(as.matrix(x)),0.1))
      })
      currMaxIndexes = apply(corMatrix,2,which.max)
      prediction = sampleConsensus$traj[currMaxIndexes]
    }else{
      prediction = unlist(lapply(1:length(unique(sampleNames)), function(firstSampleInd){
        currSample = unique(sampleNames)[firstSampleInd]
        testSample = testDataNorm[,currSample==sampleNames,drop=F]
        corMatrix = t(stats::cor(as.matrix(testSample)[seed,],refNorm[seed,],method = "spearman"))
        corMatrix = 1 - apply(as.matrix(corMatrix),2,function(x){
          unlist(computeNewData(sampleConsensus$traj,sampleConsensus$traj,t(as.matrix(x)),0.1))
        })

        startMatrix = cbind(0,t(as.matrix(corMatrix[,1])))
        matrixList = NULL
        if(dim(corMatrix)[2]>1){
          matrixList = lapply(2:dim(corMatrix)[2], function(i){
            currConsMatrix = t(matrix(corMatrix[,i],nrow = dim(corMatrix)[1],ncol = dim(corMatrix)[1]))
            currConsMatrix[lower.tri(currConsMatrix)] <- 0
            currConsMatrix
          })
        }
        endMatrix = rbind(matrix(1, ncol = 1, nrow = dim(corMatrix)[1]),0)

        namesForBigMatrix = c("Start", paste("Cons",1:dim(corMatrix)[1],1,sep="_"))
        if(dim(corMatrix)[2]>1){
          namesForBigMatrix = c(namesForBigMatrix, unlist(lapply(2:dim(corMatrix)[2], function(i){
            paste("Cons",1:dim(corMatrix)[1],i,sep="_")
          })))
        }
        namesForBigMatrix = c(namesForBigMatrix,"End")

        matrixList = c(list(startMatrix),matrixList,list(endMatrix))
        bigMatrix = Matrix::.bdiag(matrixList)
        colnames(bigMatrix) = row.names(bigMatrix) = namesForBigMatrix

        dfForGraph = data.frame(from = namesForBigMatrix[bigMatrix@i+1], to = namesForBigMatrix[bigMatrix@j+1], cost = bigMatrix@x)
        modelGraph = cppRouting::makegraph(dfForGraph,directed = T)
        shortestPath<-cppRouting::get_path_pair(modelGraph,from="Start",to="End")
        pathNodes = rev(shortestPath$Start_End)
        finalIndexesForPrediction = as.numeric(sapply(pathNodes[2:(length(pathNodes)-1)],function(x){unlist(strsplit(x,"_"))[2]}))
        currPrediction = sampleConsensus$traj[finalIndexesForPrediction]
        #certentiy = mean(1-diag(as.matrix(corMatrix[finalIndexesForPrediction,])))

        # modelGraph = igraph::graph_from_adjacency_matrix(bigMatrix,weighted = T)
        # shortestPath = shortest_paths(modelGraph, from = "Start", to = "End")
        # pathNodes = shortestPath$vpath[[1]]
        # sampleConsensus$traj[as.numeric(sapply(pathNodes$name[2:(length(pathNodes)-1)],function(x){unlist(strsplit(x,"_"))[2]}))]

        setTxtProgressBar(pb, consensusInd)
        #list(prediction = prediction, certentiy = certentiy)
        currPrediction
      }))
    }
    prediction
  }
  parallel::stopCluster(cl)
  close(pb)

  predictionMatrix = do.call(rbind,predictionStats)
  finalPredictions = colMeans(predictionMatrix)
  sampleCertainty = apply(predictionMatrix,2,sd)
  #finalPredictions = unlist(lapply(predictionStats, function(currPredictions){currPredictions$prediction}))
  #sampleCertainty = unlist(lapply(predictionStats, function(currPredictions){currPredictions$certentiy}))
  list(predictions = finalPredictions, certainty = sampleCertainty)
}

#' Calculate a robustness score for the TimeAx model
#'
#' @param model A TimeAx model.
#' @param GEData The matrix containing profiles (columns) of omics measurments (rows), which was used to train the model.
#' @param sampleNames A vector containing the individual identity of each sample in the GEData. Same vector as used in the training.
#' @param no_cores A number for the amount of cores which will be used for the analysis. The default (NULL) is total number of cores minus 1.
#' @param pseudo The output list of predictByConsensus. If not provided (NULL), pseudotime will be inferred by this function.
#' @return A robustness list consists of:
#' \item{robustnessPseudo}{Robustness pseudotime positions for all samples.}
#' \item{score}{TimeAx robustness score for the model.}
#' @references
#' Submitted
#' @examples
#' data(UBCData)
#'
#' # Training the model
#' model = modelCreation(DataUBC,UBCSamples,no_cores = 2)
#'
#' # Inferring pseudotime positions
#' robustnessStats = robustness(model,DataUBC,UBCSamples,no_cores = 2)
#' robustnessPseudo = robustnessStats$robustnessPseudo
#' robustnessScore = robustnessStats$score
#' @export
robustness = function(model, GEData, sampleNames, pseudo = NULL, no_cores = NULL){
  if(is.null(pseudo)){
    pseudo = predictByConsensus(model, GEData, no_cores = no_cores)$predictions
  }else{
    pseudo = pseudo$predictions
  }
  pseudoRobust = predictByConsensus(model, GEData, sampleNames, no_cores = no_cores)$predictions
  list(robustnessPseudo = pseudoRobust, score = stats::cor(pseudoRobust,pseudo))
}


#' UBC RNA-seq data.
#'
#' @format A matrix with 27 profiles (columns) of 5000 genes.
"DataUBC"

#' UBC sample labels.
#'
#' @format A vector with 27 sample labels for 5 individuals.
"UBCSamples"

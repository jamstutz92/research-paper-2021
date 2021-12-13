# VARIABLES
Rh1DataWithCluster <- data.frame
DGRPDataWithCluster <- data.frame
clustNumList <- vector()
combinedClusterData <- list()
initialListData <- list()

###

# METHOD - PREPARE K-MEANS DATA USING CLUSTERING AND PLOT RESULTS FOR STRAINS
kMeansAnalysis <- function(numClusters, eyeOrExpress, fileName) {
  # Create filtered data list
  if(eyeOrExpress == "Eye") {
    filteredInitData <- data.frame(initStrainList, initAverageEyeSize)
    names(filteredInitData) <- c("Strain", "Mean_Eye_Size")
    rownames(filteredInitData) <- filteredInitData$Strain
  } else {
    filteredInitData <- data.frame(filteredGeneList)
    rownames(filteredInitData) <- initStrainList
  }

  # Retrieve optimal/recommended number of clusters via factoextra graph
  # Run with silhouette analysis
  print(fviz_nbclust(filteredInitData, kmeans, method = "silhouette"))
  if(eyeOrExpress == "Eye") { 
    ggsave("RH1G69D Data/Other Plots/kMeansEyeRecClusters.jpeg")   
  } else {
    ggsave("RH1G69D Data/Other Plots/kMeansExpressRecClusters.jpeg")   
  }
  
  # Producing K-means result
  kMeansResult <- kmeans(filteredInitData, numClusters, nstart = 25)
  
  # Organizing data from K-means result
  aggregate(filteredInitData, by=list(cluster=kMeansResult$cluster), mean)
  if(eyeOrExpress == "Eye") { 
    Rh1DataWithCluster <<- cbind(filteredInitData, cluster = kMeansResult$cluster)  
  } else {
    DGRPDataWithCluster <<- cbind(filteredInitData, cluster = kMeansResult$cluster)   
  }

  # Visualizing K-means result and creating cluster groups
  if(eyeOrExpress == "Eye") { 
    print(fviz_cluster(kMeansResult, filteredInitData, choose.vars = c("Mean_Eye_Size", "Strain"), show.clust.cent = TRUE, 
      ellipse.type = "convex", labelsize = 7, ylab = "", main = "Cluster Plot for Rh1G69D Strains"))
    ggsave("RH1G69D Data/Other Plots/kMeansEyeClusters.jpeg")   
  } else {
   print(fviz_cluster(kMeansResult, filteredInitData, choose.vars = names(filteredInitData), show.clust.cent = TRUE,
      ellipse.type = "convex", labelsize = 7, main = "Cluster Plot for Averaged Expressions"))
    ggsave(paste("RH1G69D Data/Other Plots/", fileName, ".jpeg", sep = ""))
  }
  
  # (Only for expressions) Get list of expression clusters with eye sizes attached
  if(eyeOrExpress == "Express") {
    vectorList <- list()
    vectorNames <- c()
    
    for(i in 1:max(DGRPDataWithCluster$cluster)) {
      currVector <- vector()
      vectorNames <- append(vectorNames, paste("clust", i, sep = ""))
      
      for(j in 1:length(DGRPDataWithCluster$cluster)) {
        if(i == DGRPDataWithCluster$cluster[[j]]) {
          currVector <- append(currVector, Rh1DataWithCluster$Mean_Eye_Size[[j]])
        }
      }
      vectorList[[i]] <- currVector
    }
    names(vectorList) <- vectorNames
    combinedClusterData <<- vectorList
  }

  # (Only for eye size) Show cluster with eye size array
  if(eyeOrExpress == "Eye") { 
    print(ggplot(Rh1DataWithCluster, aes(initAverageEyeSize,initStrainList,color=cluster)) + geom_point() + 
      labs(title="Sorted K-Means Data for Rh1G69D Strains", x="Mean Eye Size (Pixels)", y="Strain #") + 
      theme(axis.text.x = element_text(angle = 90), axis.text.y = element_blank(), axis.title.y = element_blank()))
    ggsave("RH1G69D Data/Other Plots/kMeansRh1ExpressionStrains.jpeg")
  }
}

# METHOD - PLOT K-MEANS CATEGORIZATION FOR Rh1 DATA
kMeansRh1Categorization <- function(fileName, titleName, sentData) {
  # Get ordered version of Rh1 cluster data by cluster
  clustNumList <<- vector()
  sentRh1Data <- sentData
  orderedRh1ClusterData <- sentRh1Data[order(sentRh1Data$Mean_Eye_Size),]
  
  # Get arrangement of clusters by number
  dupes <- duplicated(orderedRh1ClusterData$cluster)
  for(i in 1:length(orderedRh1ClusterData$cluster)) {
    if(dupes[i] == FALSE) {
      clustNumList <<- c(clustNumList, orderedRh1ClusterData$cluster[[i]])
    }
  }
  
  # Produce separated K-Means lists by cluster and get max length of them
  sepKMeansList <- c()
  inc <- 1
  maxLen <- 0
  for(i in clustNumList) {
    currClust <- sort(orderedRh1ClusterData$Mean_Eye_Size[orderedRh1ClusterData$cluster == i], decreasing = FALSE)
    sepKMeansList[[inc]] <- currClust
    inc <- inc + 1
    maxLen <- maxLen + length(currClust)
  }
  emptyList <- c(1:maxLen)
  emptyList[!is.na(emptyList)] <- NA
  
  # Produce combined K-Means lists
  combKMeansList <- c()
  inc <- 1
  secondInc <- 1
  for(cluster in sepKMeansList) {
    currCombList <- emptyList
    for(element in cluster) {
      currCombList[secondInc] <- element
      secondInc <- secondInc + 1
    }
    combKMeansList[[inc]] <- currCombList
    inc <- inc + 1
  }

  # Create names for resulting cluster list
  inc <- 1
  for(i in clustNumList) {
    clustNumList[inc] <<- paste("Cluster", toString(i))
    inc <- inc + 1
  }
  names(combKMeansList) <- clustNumList
  
  combinedListData <- data.frame()
  if(titleName == "(Lowest Expression Cluster)") {
    combinedListData <- combKMeansList
    combinedListData <- c(combinedListData, initialListData)
  } else if(titleName == "(Highest Expression Cluster)") {
    initialListData <<- combKMeansList
  }
  
  # Plot combined K-Means list
  if(titleName == "(Lowest Expression Cluster)") {
    jpeg(paste("RH1G69D Data/Other Plots/", "Rh1SortedCombinedExpressStrains.jpeg", sep = ""))
    plot(unlist(combinedListData),type="p",xlim=c(1,max(sapply(combinedListData,length))), 
         main = "Sorted Rh1G69D Strains", xlab = "Index Position", ylab = "Mean Eye Size (Pixels)")
    mapply(points,combinedListData,col=seq_along(combinedListData),lty=2)
    legend("bottomright",names(combinedListData),lty=2,bty="n",cex=0.75,col=seq_along(combinedListData))
    dev.off()
  }
  combinedListData <- combKMeansList
  jpeg(paste("RH1G69D Data/Other Plots/", fileName, sep = ""))
  plot(unlist(combinedListData),type="p",xlim=c(1,max(sapply(combinedListData,length))), 
       main = paste("Sorted Rh1G69D Strains", titleName, sep = " "), xlab = "Index Position", ylab = "Mean Eye Size (Pixels)")
  mapply(points,combinedListData,col=seq_along(combinedListData),lty=2)
  legend("bottomright",names(combinedListData),lty=2,bty="n",cex=0.75,col=seq_along(combinedListData))
  dev.off()
}

###

# METHOD - Determine suspected candidate modifiers by cross-referencing gene set of clusters
# to find gene names associated with the lowest and highest average clusters
getSuspectedKMeansModifiers <- function(fileName, eyeOrExpress, sentData) {
  extremeClusterData <- sentData
  
  # Make dividers
  geneDividerValue <- (highestExpression - lowestExpression) / length(clustNumList)
  geneLowestMarker <- lowestExpression + geneDividerValue
  geneHighestMarker <- highestExpression - geneDividerValue
  
  # Create extreme quadrants for expression values
  finalAverageGeneExpress <- vector()
  combSuspectedList <<- list()
  
  if(eyeOrExpress == "Eye") {
    for(i in filteredGeneList) {
      finalAverageGeneExpress <- append(finalAverageGeneExpress, (sum(i) / length(i)))
    }  
  } else {
    for(fetchedGene in filteredGeneList) {
      currVector <- c()
      for(i in 1:length(fetchedGene)) {
        for(j in 1:length(extremeClusterData$Strain)) {
          if(i == j) { 
            currVector <- append(currVector, fetchedGene[i])
          }
        }
      }
      finalAverageGeneExpress <- append(finalAverageGeneExpress, sum(currVector) / length(currVector))
    }
  }
  
  lowestClust <- vector(); highestClust <- vector()
  for(i in finalAverageGeneExpress) {
    if(i <= geneLowestMarker) {
      lowestClust <- append(lowestClust, i)
    } else if(i >= geneHighestMarker) {
      highestClust <- append(highestClust, i)
    } else { }
  }
  lowestClust <- sort(lowestClust, decreasing = FALSE); 
  highestClust <- sort(highestClust, decreasing = FALSE); 

  # Filter out more values from lowest and highest clusters so only extremes remain
  finalLowList <- vector(); finalHighList <- vector()
  lowestClust <- head(lowestClust, 10)
  highestClust <- tail(highestClust, 10)
  lastPlotData <- vector()

  for(i in 1:length(names(filteredGeneList))) { 
    if(eyeOrExpress == "Eye") { 
      currAverage <- (sum(filteredGeneList[[i]]) / length(filteredGeneList[[i]]))
    } else {
      currAverage <- finalAverageGeneExpress[[i]]
    }
    for(j in lowestClust) {
      if(currAverage == j) {
        finalLowList <- append(finalLowList, names(filteredGeneList)[[i]])
        lastPlotData <- append(lastPlotData, j)
      }
    }
  }
  for(i in 1:length(names(filteredGeneList))) { 
    if(eyeOrExpress == "Eye") { 
      currAverage <- (sum(filteredGeneList[[i]]) / length(filteredGeneList[[i]]))
    } else {
      currAverage <- finalAverageGeneExpress[[i]]
    }
    for(j in highestClust) {
      if(currAverage == j) {
        finalHighList <- append(finalHighList, names(filteredGeneList)[[i]])
        lastPlotData <- append(lastPlotData, j)
      }
    }
  }

  # Make final plot between suspected candidate modifier genes (x-axis) and average gene expression (y-axis)
  combSuspectedList <<- c(finalLowList, finalHighList)
  finalDataframe <- data.frame(lastPlotData, row.names = combSuspectedList)
  
  finalDataframe$Quadrant <- c(clustNumList[1],clustNumList[1],clustNumList[1],clustNumList[1],clustNumList[1],clustNumList[1],clustNumList[1],clustNumList[1],clustNumList[1],clustNumList[1],
                               clustNumList[length(clustNumList)],clustNumList[length(clustNumList)],clustNumList[length(clustNumList)],clustNumList[length(clustNumList)],clustNumList[length(clustNumList)],
                               clustNumList[length(clustNumList)],clustNumList[length(clustNumList)],clustNumList[length(clustNumList)],clustNumList[length(clustNumList)],clustNumList[length(clustNumList)])
  
  ggplot(finalDataframe, aes(rownames(finalDataframe), lastPlotData, color=Quadrant)) + geom_point() + 
    labs(title="Potential candidate modifiers - K-Means", x="Genes", y="Relative Expression") +
    scale_y_continuous(limits=c(2, 15)) + theme(axis.text.x = element_text(angle = 90))
  ggsave(paste("RH1G69D Data/Other Plots/", fileName, sep = ""))
}

###

# METHOD - FILTER OUT STRAINS INTO DATA FRAME THAT ONLY INCLUDES 2 EXTREME CLUSTERS
getExtremeClusteredData <- function(numPartitions, eyeOrExpress) {
  # Run two separate methods depending on whether using eye size or average expression data
  if(eyeOrExpress == "Eye") { 
    # EYE SIZE METHOD
    # Get highest and lowest cluster numbers
    highestClust <- c()
    lowestClust <- c()
    highestEyeSize <- max(Rh1DataWithCluster$Mean_Eye_Size)
    lowestEyeSize <- min(Rh1DataWithCluster$Mean_Eye_Size)
    
    for(i in 1:length(Rh1DataWithCluster$Mean_Eye_Size)) {
      if(Rh1DataWithCluster$Mean_Eye_Size[i] == highestEyeSize) {
        highestClust <- Rh1DataWithCluster$cluster[i]
      }
      if(Rh1DataWithCluster$Mean_Eye_Size[i] == lowestEyeSize) {
        lowestClust <- Rh1DataWithCluster$cluster[i]
      }
    }
    
    # Create filtered Rh1 cluster data frame
    lowEndData <- Rh1DataWithCluster[Rh1DataWithCluster$cluster==lowestClust,]
    highEndData <- Rh1DataWithCluster[Rh1DataWithCluster$cluster==highestClust,]
    
    extremeUpperEyeClusterData <<- highEndData 
    extremeLowerEyeClusterData <<- lowEndData
    extremeCombinedTest <<- rbind(highEndData,lowEndData)
  } else {
    # AVERAGE EXPRESSION METHOD
    # Initial variables
    sortedRh1ClusterData <- Rh1DataWithCluster[order(Rh1DataWithCluster$Mean_Eye_Size),]
    meanDividerValue <- (max(sortedRh1ClusterData$Mean_Eye_Size) - min(sortedRh1ClusterData$Mean_Eye_Size)) / numPartitions
    highThresholdValue <- max(sortedRh1ClusterData$Mean_Eye_Size) - meanDividerValue
    lowThresholdValue <- min(sortedRh1ClusterData$Mean_Eye_Size) + meanDividerValue

    highestClustEyeSize <- c(); lowestClustEyeSize <- c();
    highestClustStrains <- c(); lowestClustStrains <- c(); 
    highestClustNum <- 0; lowestClustNum <- 0;
    
    # Retrieve highest and lowest value clusters and their mean eye sizes
    getHighLowClustValues <- function(highOrLow) {
      prevClustEyeSize <- c()
      decision <- highOrLow
      
      for(i in 1:length(combinedClusterData)) {
        currClustEyeSize <- c()
        for(j in combinedClusterData[[i]]) {
          if(decision == "Low") {
            if(j < lowThresholdValue) {
              currClustEyeSize <- append(currClustEyeSize, j)
            }          
          } else {
            if(j > highThresholdValue) {
              currClustEyeSize <- append(currClustEyeSize, j)
            }          
          }
        }
        if(length(currClustEyeSize) > length(prevClustEyeSize)) {
          prevClustEyeSize <- currClustEyeSize
          if(decision == "Low") {
            lowestClustEyeSize <<- combinedClusterData[[i]]
            lowestClustNum <<- i        
          } else {
            highestClustEyeSize <<- combinedClusterData[[i]]
            highestClustNum <<- i           
          }
        }
      } 
    }
    getHighLowClustValues("Low")
    getHighLowClustValues("High")
    
    # Retrieve highest and lowest strain numbers of extreme clusters
    for(i in 1:length(Rh1DataWithCluster$Mean_Eye_Size)) {
      for(j in 1:length(highestClustEyeSize)) {
        if(Rh1DataWithCluster$Mean_Eye_Size[[i]] == highestClustEyeSize[[j]]) { 
          highestClustStrains <- append(highestClustStrains, Rh1DataWithCluster$Strain[[i]])
        }
      }
    }
    for(i in 1:length(Rh1DataWithCluster$Mean_Eye_Size)) {
      for(j in 1:length(lowestClustEyeSize)) {
        if(Rh1DataWithCluster$Mean_Eye_Size[[i]] == lowestClustEyeSize[[j]]) { 
          lowestClustStrains <- append(lowestClustStrains, Rh1DataWithCluster$Strain[[i]])
        }
      }
    }

    # Create filtered Rh1 cluster data frame
    lowEndData <- data.frame(lowestClustStrains, lowestClustEyeSize, lowestClustNum)
    colnames(lowEndData) <- c("Strain", "Mean_Eye_Size", "cluster")
    highEndData <- data.frame(highestClustStrains, highestClustEyeSize, highestClustNum)
    colnames(highEndData) <- c("Strain", "Mean_Eye_Size", "cluster")

    extremeUpperClusterData <<- highEndData 
    extremeLowerClusterData <<- lowEndData
    extremeCombinedTest <<- rbind(highEndData,lowEndData)
  }
}

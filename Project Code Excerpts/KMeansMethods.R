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
         main = "Sorted Rh1G69D Strains", xlab = "Index Position", ylab = "Mean Eye Size (Pixels)", xaxt='n')
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

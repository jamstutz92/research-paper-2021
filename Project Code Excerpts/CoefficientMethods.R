# VARIABLES
combDataFrame <- data.frame

# METHOD - Retrieve correlation coefficients and p-values for all genes and get highest/lowest extreme values for them
getAllCoeffValues <- function(pathName, eyeOrExpress, titleName, sentData, shouldComb) {
  extremeClusterData <- sentData
  
  # Needed variables
  correlationCoefficients <- c()
  corrCoeffGeneNames <- c()
  pValueList <- c()
  pValueSignificanceList <- c()
  
  # Produce all correlation coefficients
  inc <- 1
  for(i in names(filteredGeneList)) {
    currVector <- filteredGeneList[[i]]
    
    if(eyeOrExpress == "Express") { 
      fetchedGene <- c()
      for(k in 1:length(currVector)) {
        for(l in 1:length(extremeClusterData$Strain)) {
          if(k == l) { 
            fetchedGene <- append(fetchedGene, currVector[k])
          }
        }
      }
    } else {
      fetchedGene <- currVector
    }
  
    correlationCoefficients[inc] <- cor(fetchedGene, sortedAverageEyeSizes, method = correlationMethod)
    corrCoeffGeneNames[inc] <- names(filteredGeneList[i])
    
    inc <- inc + 1
  }
  
  # Plot out combined plot of highest and lowest coefficient values
  coeffDataFrame <- data.frame(coeff = correlationCoefficients, gene = corrCoeffGeneNames)
  coeffDataFrame <- coeffDataFrame[order(coeffDataFrame$coeff),]
  lowestCoeff <- head(coeffDataFrame$coeff, 10)
  highestCoeff <- tail(coeffDataFrame$coeff, 10)
  lowestNames <- head(coeffDataFrame$gene, 10)
  highestNames <- tail(coeffDataFrame$gene, 10)
  
  combCoeff <- c(lowestCoeff, highestCoeff)
  combCoeffNames <- c(lowestNames, highestNames)
  
  correlationCoefficients <- c()
  corrCoeffGeneNames <- c()
  
  # Plot the correlated gene expressions and average eye sizes of extreme coefficient genes
  inc <- 1
  for(i in names(filteredGeneList)) {
    for(j in combCoeffNames) {
      if(j == i) {
        Clusters <- sortedRh1ClusterData$cluster
        currVector <- filteredGeneList[[i]]
        
        if(eyeOrExpress == "Express") { 
          fetchedGene <- c()
          for(k in 1:length(currVector)) {
            for(l in 1:length(extremeClusterData$Strain)) {
              if(k == l) { 
                fetchedGene <- append(fetchedGene, currVector[k])
              }
            }
          }
        } else {
          fetchedGene <- currVector
        }
        
        tempDataFrame <- data.frame(sortedAverageEyeSizes, fetchedGene)
        
        # Calculate p-value and attach to two vectors for later
        pValue <- cor.test(fetchedGene, sortedAverageEyeSizes, method = correlationMethod)
        pValue <- round(pValue$p.value, digits = 5)
        pValueList <- c(pValueList, pValue)
        if(pValue > 0.05) {
          pValueSignificanceList <- c(pValueSignificanceList, "Not Significant")
        } else {
          pValueSignificanceList <- c(pValueSignificanceList, "Significant") 
        }
        
        # Calculate correlation coefficients from current vectors
        correlationCoefficients[inc] <- cor(fetchedGene, sortedAverageEyeSizes, method = correlationMethod)
        corrCoeffGeneNames[inc] <- names(filteredGeneList[i])
        
        mid<-mean(Clusters)
        print(ggplot(tempDataFrame, aes(x=sortedAverageEyeSizes, y=fetchedGene, label=sortedRh1ClusterData$Strain, color=Clusters)) + 
                labs(x="Average Eye Size (Pixels)", y="Average Genetic Expression", title=paste("Gene", names(filteredGeneList[i]), ", P-value:", pValue)) +
                geom_point() + geom_text(aes(label=sortedRh1ClusterData$Strain), hjust=0, vjust=0,size=2) + scale_color_gradient2(midpoint=mid, low="blue", mid="red", high="orange"))
        ggsave(paste(pathName, "Gene", inc, ".jpeg", sep = ""))
        
        inc <- inc + 1
        Sys.sleep(0)
      }
    }
  }
  tempDataFrame <- data.frame(coeffNames = corrCoeffGeneNames, coeff = correlationCoefficients, pList = pValueList, pValues = pValueSignificanceList)
  cat(corrCoeffGeneNames, file = "geneNames.txt", sep = "\n", append = TRUE)
  
  # If upper and lower bounds of clusters treated as separate frames, combine them for later to plot
  Cluster <- c()
  for(i in 1:length(tempDataFrame$coeffNames)) {
    if((substr(titleName, 0, 2) == "(L") || substr(titleName, 0, 2) == "(A") {
      Cluster <- append(Cluster, "Lowest Eye Size Cluster")
    } else {
      Cluster <- append(Cluster, "Highest Eye Size Cluster")
    }
  }
  if(shouldComb == TRUE) {
    tempDataFrame <- cbind(tempDataFrame, Cluster)
    combDataFrame <<- rbind(combDataFrame, tempDataFrame)
  } else {
    tempDataFrame <- cbind(tempDataFrame, Cluster)
    combDataFrame <<- tempDataFrame
  }

  # Plot correlation coefficients per suspected gene
  print(ggplot(tempDataFrame, aes(x = coeffNames, y = coeff, label = pList, color=pValues)) + geom_point() + 
    labs(title=paste("Suspected Genes", titleName, sep = " "), x="Genes", y="Coefficient Value") +
    theme(axis.text.x = element_text(angle = 90)) + 
    geom_text_repel(aes(label=pList), hjust=0, vjust=2, size=2)) +
    scale_x_discrete(expand=c(0.1, 0)) + 
    scale_colour_manual(values = c("Significant" = "dark red", "Not Significant" = "dark blue"))
  ggsave(paste(pathName, "CoefficientPlot.jpeg", sep = ""))
}

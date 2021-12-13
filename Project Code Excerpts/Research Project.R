# James Amstutz

# Install package - ggplot2 helps with visualization in graphs
install.packages("ggplot2")
# factoextra produces cluster graphs
install.packages("factoextra")
# For testing concordant/discordant pairs when using Kendall's correlation analysis
install.packages("DescTools")

# Enable packages/libraries
# require("seqinr") - Keep if writing output files used in future
require("ggplot2")
require("factoextra")
require("ggrepel")
require("stringr")
require("DescTools")

# Set working directory
setwd("/Users/jamesamstutz/repos/Bioinformatics/Thesis")

# VARIABLES
finalAverageEyeSize <- vector()  # Average eye sizes per strain
initAverageEyeSize <- vector() # Initial, unsorted average eye sizes
initStrainList <- vector() # Initial, unsorted strain list

averageDGRPData <- list() # The average genetic expression of each filtered DGRP gene
filteredGeneList <- list() # Combination of filtered DGRP gene Data and average eye sizes for each gene
combSuspectedList <- list() # Combined list of suspected candidate modifier genes

dgrpData <- data.frame # Base DGRP data retrieved from text file
filteredDGRPData <- data.frame # DGRP data with strains not featured in RH1 data filtered out
extremeClusterData <- data.frame # Data that only contains extreme clusters
extremeCombinedTest <- data.frame # Combines extreme cluster data frames to plot combination correlation graph
oldRh1DataWithCluster <- data.frame # Maintain Rh1DataWithCluster before extreme methods run

highestExpression <- 0 # Highest possible gene expression found
lowestExpression <- 9999 # Lowest possible gene expression found

correlationMethod <- "pearson" # Works as 'pearson', 'kendall' or 'spearman'

# METHOD - FETCH RH1 AND DGRP DATA FROM TEXT FILES
fetchData <- function() {
  # Read Rh1G69D.txt file and create strain list numbers
  Rh1Data <- read.delim("Rh1G69D.txt", header = TRUE, sep = "", dec = ".")
  
  # Removing this from data selection due to it not having two annotated lines from DGRP data
  for(i in 1:length(Rh1Data$Strain)) {
    if(Rh1Data$Strain[i] == "RAL513") {
      Rh1Data <- Rh1Data[-c(i),]
      break
    }
  }

  # Remove RAL prefix for Rh1 strains
  for(i in 1:length(Rh1Data$Strain)) {
    substr <- substring(Rh1Data$Strain[i], 4, 6)
    if(substring(substr, 0, 1) == "0") {
      substr <- substring(substr, 2, 3)
    }
    Rh1Data$Strain[i] <- substr
  }
  
  # Plot scattered Rh1 data as-is, before categorization
  ggplot(Rh1Data, aes(Mean_Eye_Size, Strain)) + geom_point() +
    labs(x = "Mean Eye Size (pixels)") +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_text(size = 16), plot.title = element_text(size = 16)) +
    ggtitle("Rh1G69D Expression Strains") + 
    scale_x_continuous(limits = c(14000,27500), expand = c(0.1, 0))
  ggsave("RH1G69D Data/Other Plots/Rh1ExpressionStrains.jpeg", width=12, height=8)

  # Read dgrp.array.exp.female.txt file
  dgrpData <<- read.delim("dgrp.array.exp.female.txt", header = TRUE, sep = "", dec = ".")
  
  # Return fetched Rh1 data for use in next function
  return(Rh1Data)
}

###

# METHOD - GET AVERAGE EYE SIZES AND FILTER GENETIC DATA
filterData <- function(Rh1Data) {
  Rh1StrainList <- sort(Rh1Data$Strain, decreasing = FALSE)

  # 1: CREATE DGRP STRAIN LIST
  # Create list of line/strain numbers for later use from DGRP data, excluding gene column
  DGRPStrainList <- vector()
  dgrpColNames <- colnames(dgrpData)
  
  for(i in 1:length(dgrpColNames)) {
    if(dgrpColNames[i] != "gene") {
      subStr <- substring(dgrpColNames[i], 6, 8)
      
      if(grepl(".", subStr, fixed = TRUE) == TRUE) {
        subStr <- substring(dgrpColNames[i], 6, 7)
      }
      DGRPStrainList[i] <- subStr 
    }
  }
  
  DGRPStrainList <- strtoi(DGRPStrainList[!is.na(DGRPStrainList)])
  DGRPStrainList <- sort(DGRPStrainList, decreasing = FALSE)
  
  # 3: CREATE FILTERED STRAIN/LINE LISTS AND AVERAGE EYE SIZES
  # Init strain data is used later for k-means plotting
  initStrainList <<- intersect(Rh1Data$Strain, DGRPStrainList)
  
  # Cross-reference two strain/line lists to get matching values
  finalStrainList <- sort(intersect(Rh1StrainList, DGRPStrainList), decreasing = FALSE)
  
  # Retrieve the relevant average eye sizes from the cross-referenced strain list
  for(i in 1:length(Rh1Data$Strain)) {
    for(j in 1:length(finalStrainList)) {
      if(Rh1Data$Strain[i] == finalStrainList[j]) {
        finalAverageEyeSize <<- append(finalAverageEyeSize, Rh1Data$Mean_Eye_Size[i])
      }  
    }
  }
  
  # Init and sorted average eye size data is used later for k-means plotting
  initAverageEyeSize <<- finalAverageEyeSize
  finalAverageEyeSize <<- sort(finalAverageEyeSize, decreasing = FALSE)
  
  # Filter out DGRP strains/lines that are not found in Rh1 data
  drops <- vector()
  filteredData <<- dgrpData
  comparisonList <- is.element(DGRPStrainList, finalStrainList)

  for(i in 1:length(comparisonList)) {
    if(comparisonList[i] == "FALSE") {
      drops <- append(drops, names(filteredData)[i+1])
    }
  }
  
  #drops <- drops[!drops %in% drops[5:6]] - Uncomment if using p53 data
  filteredDGRPData <<- filteredData[ , !(names(filteredData) %in% drops)]

  ###
  
  # 4: GET AVERAGE OF TWO LINE EXPRESSIONS FOR DGRP DATA AND PLOT EXAMPLE
  # Get average RNA expression from lines 1-2 for each gene
  count <- 0
  inc <- 1
  
  # For loop
  for(i in filteredDGRPData) {
    # Exclude gene name column from loop
    if(i[1] != "FBgn0000014") {
      if(count == 0) {
        # Fetch first line expression (e.g. line_21:1)
        line1 <- i
        count <- count+1
        
        # If last line expression, grab line_913.2 since it is otherwise skipped
        if(names(filteredDGRPData[inc]) == "line_913.1") {
          line2 <- i
          count <- 0
          combineList <- c(1:length(line1))
          
          # Get average of first two lines fetched
          for(j in 1:length(line1)) {
            combineList[j] <- (line1[j] + line2[j]) / 2
          }
          
          # Append list with average of first two line expressions
          listForData <- list(combineList)
          averageDGRPData <<- c(averageDGRPData, listForData) 
          inc <- inc + 1
        }
        
        inc <- inc + 1
      } else {
        # Fetch second line expression (e.g. line_21:2) and reset count
        line2 <- i
        count <- 0
        combineList <- c(1:length(line1))

        # Get average of first two lines fetched
        for(j in 1:length(line1)) {
          combineList[j] <- (line1[j] + line2[j]) / 2
        }
        
        # Append list with average of first two line expressions
        listForData <- list(combineList)
        averageDGRPData <<- c(averageDGRPData, listForData)
        inc <- inc + 1
      }
    }
  }

  ###
  
  # 5: Retrieve filtered DGRP expressions associated with Rh1 strains and combine
  initFilteredGeneList <- list()
  
  for(i in 1:length(filteredDGRPData$gene)) {
    currGene <- vector()
    
    for(j in 1:length(averageDGRPData)) {
      currGene <- append(currGene, averageDGRPData[[j]][i])
    } 
    initFilteredGeneList[[i]] <- currGene
  }
  names(initFilteredGeneList) <- filteredDGRPData$gene
  filteredGeneList <<- initFilteredGeneList
}

###

# METHOD - Reset all global variables before running extremeTesting
resetData <- function() {
  # Reset most data values
  finalAverageEyeSize <<- vector()
  initAverageEyeSize <<- vector()
  initStrainList <<- vector()
  averageDGRPData <<- list()
  filteredGeneList <<- list()
  combSuspectedList <<- list()
  filteredDGRPData <<- data.frame
  highestExpression <<- 0
  lowestExpression <<- 9999
  clustNumList <<- vector()  
}

###

# METHOD - Print combined coefficient/p-list plots
printCombData <- function(pathName, titleName) {
  # Save and display combined data - Title deprecated
  print(ggplot(combDataFrame, aes(x = coeffNames, y = coeff, label = pList, color=Cluster)) + geom_point() + 
    labs(x="Genes", y="Coefficient Value") +
    theme(axis.text.x = element_text(angle = 90), legend.title = element_blank(), legend.position = "bottom") + 
    geom_text_repel(aes(label=pList), hjust=0, vjust=2, size=2)) +
    scale_x_discrete(expand=c(0.1, 0)) + 
    scale_colour_manual(values = c("Lowest Eye Size Cluster" = "dark red", "Highest Eye Size Cluster" = "dark blue"))
  ggsave(paste("RH1G69D Data/", correlationMethod, " clusters/CoefficientPlot", pathName, ".jpeg", sep = ""))
}

###

# Get user input (typically use 6 for numClust and 4 for numPartitions)
numClusters <- strtoi(readline(prompt="Enter # clusters to be produced: "))
numPartitions <- strtoi(readline(prompt="Enter # partitions for express cluster evaluation: "))
correlationMethod <<- readline(prompt="Enter correlation method (pearson, kendall or spearman): ")

# Create directories for results
dir.create("RH1G69D Data")
dir.create("RH1G69D Data/Other Plots")
dir.create(paste("RH1G69D Data/", correlationMethod, " clusters", sep=""))
dir.create(paste("RH1G69D Data/", correlationMethod, " clusters/All Eye Clusters", sep=""))
dir.create(paste("RH1G69D Data/", correlationMethod, " clusters/All Upper Express Clusters", sep=""))
dir.create(paste("RH1G69D Data/", correlationMethod, " clusters/All Lower Express Clusters", sep=""))
dir.create(paste("RH1G69D Data/", correlationMethod, " clusters/Extreme Upper Eye Clusters", sep=""))
dir.create(paste("RH1G69D Data/", correlationMethod, " clusters/Extreme Lower Eye Clusters", sep=""))
dir.create(paste("RH1G69D Data/", correlationMethod, " clusters/Extreme Upper Express Clusters", sep=""))
dir.create(paste("RH1G69D Data/", correlationMethod, " clusters/Extreme Lower Express Clusters", sep=""))

# Run initial program
initRh1Data <- fetchData()
filterData(initRh1Data)

# Create sample data frame to display filtered DGRP results for strains 21 and 26
firstSampleSet <- data.frame(name = "Strain 21", sampleGenes = filteredDGRPData$gene[1:9], sampleAverages = averageDGRPData[[1]][1:9])
secondSampleSet <- data.frame(name = "Strain 26", sampleGenes = filteredDGRPData$gene[1:9], sampleAverages = averageDGRPData[[2]][1:9])
combSampleSet <- rbind(firstSampleSet, secondSampleSet)

print(ggplot(combSampleSet, aes(x = sampleGenes, y = sampleAverages, color=name)) + geom_point() + 
  labs(title="Sample DGRP Array for Strains 21 & 26", x="Gene #", y="Average RNA Expression") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90))) +
  scale_colour_manual(values = c("Strain 21" = "dark red", "Strain 26" = "dark blue"))
ggsave("RH1G69D Data/Other Plots/DGRPStrainExample.jpeg", width=8, height=4)

# Call KMeansMethods.R file and plot Rh1 strain data and candidate modifiers
source("KMeansMethods.R")
kMeansAnalysis(numClusters, "Eye", "")
kMeansAnalysis(numClusters, "Express", "kMeansExpressClusters")
getExtremeClusteredData(numPartitions, "Express")

# Retrieve suspected candidate modifiers based on extremes of values
kMeansRh1Categorization("Rh1SortedEyeStrains.jpeg", "(All Eye Clusters)", Rh1DataWithCluster)
kMeansRh1Categorization("Rh1SortedUpperExpressStrains.jpeg", "(Highest Expression Cluster)", extremeUpperClusterData)
kMeansRh1Categorization("Rh1SortedLowerExpressStrains.jpeg", "(Lowest Expression Cluster)", extremeLowerClusterData)

# Get correlation coefficients for suspected candidate modifiers using all K-Means clusters
source("CoefficientMethods.R")
# For eye size
sortedRh1ClusterData <- Rh1DataWithCluster[order(Rh1DataWithCluster$Strain),]
sortedAverageEyeSizes <- sortedRh1ClusterData$Mean_Eye_Size
getAllCoeffValues(paste("RH1G69D Data/", correlationMethod, " clusters/All Eye Clusters/", sep=""), "Eye", "(All Eye Clusters)", NULL, FALSE)

# For expression values
sortedRh1ClusterData <- Rh1DataWithCluster[order(extremeUpperClusterData$Strain),]
sortedAverageEyeSizes <- sortedRh1ClusterData$Mean_Eye_Size
combDataFrame <<- data.frame
getAllCoeffValues(paste("RH1G69D Data/", correlationMethod, " clusters/All Upper Express Clusters/", sep=""), "Express", "(Highest Expression Cluster)", extremeUpperClusterData, FALSE)
sortedRh1ClusterData <- Rh1DataWithCluster[order(extremeLowerClusterData$Strain),]
sortedAverageEyeSizes <- sortedRh1ClusterData$Mean_Eye_Size
getAllCoeffValues(paste("RH1G69D Data/", correlationMethod, " clusters/All Lower Express Clusters/", sep=""), "Express", "(Lowest Expression Cluster)", extremeLowerClusterData, TRUE)
printCombData("CombExpress", "(Combined Expression Clusters)")

###

# Run second half of algorithm for genes using extreme clusters, for eye size and expressions
# (Upper bound eye size cluster)
oldRh1DataWithCluster <<- Rh1DataWithCluster
resetData()

# Filter out Rh1 strains/average eye sizes into data frame that only includes extreme 2 clusters
Rh1DataWithCluster <<- oldRh1DataWithCluster
getExtremeClusteredData(numPartitions, "Eye")
filterData(extremeUpperEyeClusterData)
Rh1DataWithCluster <<- extremeUpperEyeClusterData  
  
# Repeat methods/variable sets used earlier, this time with extreme cluster data
kMeansRh1Categorization("extremeSortedUpperEyeStrains.jpeg", "(Highest Eye Clusters w/ Filtering)", Rh1DataWithCluster)
sortedRh1ClusterData <- extremeUpperEyeClusterData[order(extremeUpperEyeClusterData$Strain),]
sortedAverageEyeSizes <- sortedRh1ClusterData$Mean_Eye_Size

# Repeat gathering coefficient values, this time with extreme cluster data
combDataFrame <<- data.frame
getAllCoeffValues(paste("RH1G69D Data/", correlationMethod, " clusters/Extreme Upper Eye Clusters/", sep=""), "Eye", "(Highest Eye Cluster w/ Filtering)", NULL, FALSE) 

# (Lower bound eye size cluster)
resetData()

# Filter out Rh1 strains/average eye sizes into data frame that only includes extreme 2 clusters
Rh1DataWithCluster <<- oldRh1DataWithCluster
getExtremeClusteredData(numPartitions, "Eye")
filterData(extremeLowerEyeClusterData)
Rh1DataWithCluster <<- extremeLowerEyeClusterData   
  
# Repeat methods/variable sets used earlier, this time with extreme cluster data
kMeansRh1Categorization("extremeSortedLowerEyeStrains.jpeg", "(Lowest Eye Cluster w/ Filtering)", Rh1DataWithCluster)
sortedRh1ClusterData <- extremeLowerEyeClusterData[order(extremeLowerEyeClusterData$Strain),]
sortedAverageEyeSizes <- sortedRh1ClusterData$Mean_Eye_Size

# Repeat gathering coefficient values, this time with extreme cluster data
getAllCoeffValues(paste("RH1G69D Data/", correlationMethod, " clusters/Extreme Lower Eye Clusters/", sep=""), "Eye", "(Lowest Eye Cluster w/ Filtering)", NULL, TRUE) 
printCombData("CombFiltEyeSize", "(Combined Filtered Eye Clusters)")

###

# (Upper bound expression cluster)
resetData()

# Filter out Rh1 strains/average eye sizes into data frame that only includes extreme 2 clusters
Rh1DataWithCluster <<- oldRh1DataWithCluster
getExtremeClusteredData(numPartitions, "Express")
filterData(extremeUpperClusterData)
Rh1DataWithCluster <<- extremeUpperClusterData

# Repeat methods/variable sets used earlier, this time with extreme cluster data
kMeansRh1Categorization("extremeSortedUpperExpressStrains.jpeg", "(Highest Expression Cluster w/ Filtering)", Rh1DataWithCluster)
sortedRh1ClusterData <- extremeUpperClusterData[order(extremeUpperClusterData$Strain),]
sortedAverageEyeSizes <- sortedRh1ClusterData$Mean_Eye_Size

# Repeat gathering coefficient values, this time with extreme cluster data
combDataFrame <<- data.frame
getAllCoeffValues(paste("RH1G69D Data/", correlationMethod, " clusters/Extreme Upper Express Clusters/", sep=""), "Express", "(Highest Expression Cluster w/ Filtering)", extremeUpperClusterData, FALSE) 

# (Lower bound expression cluster)
resetData()

# Filter out Rh1 strains/average eye sizes into data frame that only includes extreme 2 clusters
Rh1DataWithCluster <<- oldRh1DataWithCluster
getExtremeClusteredData(numPartitions, "Express")
filterData(extremeLowerClusterData)
Rh1DataWithCluster <<- extremeLowerClusterData

# Repeat methods/variable sets used earlier, this time with extreme cluster data
kMeansRh1Categorization("extremeSortedLowerExpressStrains.jpeg", "(Lowest Expression Cluster w/ Filtering)", Rh1DataWithCluster)
sortedRh1ClusterData <- extremeLowerClusterData[order(extremeLowerClusterData$Strain),]
sortedAverageEyeSizes <- sortedRh1ClusterData$Mean_Eye_Size

# Repeat gathering coefficient values, this time with extreme cluster data
getAllCoeffValues(paste("RH1G69D Data/", correlationMethod, " clusters/Extreme Lower Express Clusters/", sep=""), "Express", "(Lowest Expression Cluster w/ Filtering)", extremeLowerClusterData, TRUE) 
printCombData("CombFiltExpress", "(Combined Filtered Expression Clusters)")

###

# Old lines of code go here - can re-implement if needed

# (Use before K-Means clustering)
# Call QuadrantMethods.R file and plot Rh1 strain data and candidate modifiers
#source("QuadrantMethods.R")
#quadrantRh1Categorization() 

# (Use between getting suspected candidates for clustering)
#getSuspectedQuadrantModifiers()
#getSuspectedKMeansModifiers("kMeansEyeCandidates.jpeg", "Eye", extremeUpperClusterData)
#getSuspectedKMeansModifiers("kMeansExpressCandidates.jpeg", "Express", extremeUpperClusterData)
#getSuspectedKMeansModifiers("kMeansExpressCandidates.jpeg", "Express", extremeLowerClusterData)

# (Use before getting all correlation coefficient values)
#getSuspectedCoeffValues("RH1G69D Data/All Eye Clusters/", "Eye")
#getSuspectedCoeffValues("RH1G69D Data/All Express Clusters/", "Express")

# (Use during second half of project for extreme values)
#getSuspectedKMeansModifiers("extremeKMeansEyeCandidates.jpeg", "Eye", extremeUpperClusterData)
#getSuspectedCoeffValues("RH1G69D Data/Extreme Eye Clusters/", "Eye")
#getSuspectedKMeansModifiers("extremeKMeansExpressCandidates.jpeg", "Express", extremeUpperClusterData)
#getSuspectedCoeffValues(paste("RH1G69D Data/Extreme", "Express", "Clusters/", sep = " "), "Express")

# (Use to determine max and min gene expression values for deprecated code)
# Max
#for(i in filteredGeneList) {
  #if(max(i) > highestExpression) {
    #highestExpression <<- max(i)
  #}
#}
# Min
#for(i in filteredGeneList) {
  #if(min(i) < lowestExpression) {
    #lowestExpression <<- min(i)
  #}
#}


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
require("scales")

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

# Only run if producing sample cluster plot
testFrame <- data.frame(replicate(10,sample(0:1000,100,rep=TRUE)))
kMeansResult <- kmeans(testFrame, 3, nstart = 25)
print(fviz_cluster(kMeansResult, testFrame, choose.vars = names(testFrame), show.clust.cent = TRUE,
                   ellipse.type = "convex", labelsize = 7, main = NULL))
ggsave(paste("RH1G69D Data/Other Plots/sampleCluster.jpeg", sep = ""))

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
    scale_x_continuous(label = comma, limits = c(14000,27500), expand = c(0.1, 0))
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

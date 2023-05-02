#Breaking_up_merged_peaks_max_2kb_bins_parallel_2021October20.R
# Author: Belle Moyers

################################################################################
#If you want to run things at the command line, run the following command:
#Rscript /path/to/this/script.R $merged_Bed_file $unmerged_Bed_file 
#Rscript /cluster/home/bmoyers/Scripts/BrainTF_AlleleSpecificBinding/for_jacob/Breaking_up_merged_peaks_max_2kb_bins_2021October19.R /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/Jacob_Merged_Peaks/union_DLPFC.rds /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/Jacob_Merged_Peaks/unlist_DLPFC.rds
#If you want to make the regions even smaller than 2kb, then you can also add the
#$maxLength argument to the end.
#
#the "merged_Bed_file" should correspond to the location of your "union_TISSUE.rds" and 
#the "unmerged_Bed_file" should correspond to the location of your "unlist_TISSUE.rds".
#the names of these files might seem a little confusing at first.  This is mostly because 
#I made this script by modifying a prior script, and didn't want to change variable names.
#
#If you are running this in R studio, you'll be copy-pasting things.  In that case, 
#don't copy-paste the lines below which read 
#args <- commandArgs(trailingOnly=T)
#merged_Bed_file <- args[1]
#unmerged_Bed_file <- args[2]
#
#Instead, enter the paths to your files manually.
#
#Depending on how many regions are >2kb long, this will take 2 hours or more to run.
#When it's done, you should have a file in the same directory as your "merged_Bed_file"
#with a slightly different name.
################################################################################

################################################################################
################################################################################
#If running this on the cluster, make sure you've loaded the proper module first.
################################################################################
################################################################################


#module load cluster/R/4.0.0



################################################################################
################################################################################
#Load Libraries and set options.
################################################################################
################################################################################

library(Rsamtools)
library(GenomicRanges)
library(parallel)
options("scipen"=100)


################################################################################
################################################################################
#Define any functions necessary for the Script.
################################################################################
################################################################################





################################################################################
#This is the main function that will break up each region of >maxLength length 
#into smaller regions, and return a data.frame of those smaller regions.
#it calls the above function, breakBins, for each line of bad_Bed.
#
#bad_Bed is a GRanges object containing all regions of length > maxLength
#unmerged_Regions is the original "unlist_TISSUE" genomic ranges object.
#maxLength is the maximum length you want regions to be, default 2kb.
#
#I've rewritten this function in a few ways that are necessary for parallelization.
################################################################################

mainFunction <- function(theChromosome, bad_Bed, unmerged_Bed, maxLength) {
  
  library(Rsamtools)
  library(GenomicRanges)
  options("scipen"=100)

  ################################################################################
  #The following function takes a region whose length is > maxLength  and breaks 
  #it into smaller regions with a maximum length of maxLength and with the start 
  #and ends overlapping perfectly with the start and end of the original region. 
  #it also determines the number and identity of unique TFs bound in each subregion.
  #
  #thisLine is a GRanges object with a single entry (length>maxLength)
  #unmerged_Regions is the original "unlist_TISSUE" genomic ranges object.
  #maxLength is the maximum length you want regions to be, default 2kb.
  ################################################################################
  
  breakBins <- function(thisLine, unmerged_Bed, maxLength) {
    
    
    ##############################################################################
    #Intersect the umerged_Bed with the region of interest to pull out ONLY the peaks 
    #we're interested in.
    ##############################################################################
    theIntersect <- as.data.frame(findOverlaps(thisLine, unmerged_Bed))
    thesePeaks <- unmerged_Bed[unique(theIntersect[,2])]
    
    ##############################################################################
    #Identify the maximum bp we're interested in.
    ##############################################################################
    the_last <- end(thisLine)+(maxLength-1)
    
    ##############################################################################
    #Identify all possible starting locations for breaks, defined as 
    #every base pair from (regionStart)-(maxLength-1) to (regionStart).  
    #For every one of those [maxLength] number of possible start locations, 
    #identify exactly where each breakpoint will be.
    #Turn this into a dataframe.
    #the last column of this dataframe stores which start location corresponds to 
    # a given line.
    ##############################################################################
    all_starts <- (start(thisLine)-(maxLength-1)):(start(thisLine))
    starts_df <- c()
    for (j in 1:length(all_starts)) {
      theseBreaks <- seq(all_starts[j], the_last, by=maxLength)
      thisTable <- cbind(rep(as.character(seqnames(thisLine)), length(theseBreaks)), theseBreaks, theseBreaks, rep(all_starts[j], length(theseBreaks)))
      starts_df <- rbind(starts_df, thisTable)
    }
    starts_df <- as.data.frame(starts_df)
    starts_df[,2] <- as.numeric(as.character(starts_df[,2]))
    starts_df[,3] <- as.numeric(as.character(starts_df[,3]))
    starts_df[,4] <- as.numeric(as.character(starts_df[,4]))
    colnames(starts_df) <- c("chr", "start", "end", "all_starts")
    
    ##############################################################################
    #make a GRanges object out of your starts data frame.  Intersect it with the 
    #peaks.  This gives you an intersection of every possible break point 
    #with all of the peaks.  Essentially this tells you, "If I choose this starting
    #point, how many peaks am I going to split over multiple regions?"
    ##############################################################################
    starts_df_gr <- makeGRangesFromDataFrame(starts_df, keep.extra.columns=TRUE, ignore.strand=TRUE)
    theIntersect <- as.data.frame(findOverlaps(starts_df_gr, thesePeaks))
    theIntersect <- theIntersect[order(theIntersect[,1]),]
    
    ##############################################################################
    #For every single break point, identify the number of peaks being split.
    ##############################################################################
    splitPeaksCount <- rep(0, nrow(starts_df))
    for (j in 1:nrow(starts_df)) {
      splitPeaksCount[j] <- nrow(theIntersect[theIntersect[,1]==j,])
    }
    starts_df$splitCount <- splitPeaksCount
    
    ##############################################################################
    #For every possible start location, count the total number of peak splits 
    #it creates.  Beause each start point corresponds to multiple break points, we 
    #have to sum multiple break points for each start to determine how "bad" a start
    #that is.
    ##############################################################################
    final_splitCounts <- rep(0, length(all_starts))
    for (j in 1:length(all_starts)) {
      totalSplit <- sum(starts_df[starts_df$all_starts==all_starts[j],"splitCount"])
      final_splitCounts[j] <- totalSplit
    }
    
    ##############################################################################
    #Identify the start points which create the minimum number of peak splits.
    #if there are multiple, choose the first one that occurs.
    ##############################################################################
    bestStarts <- all_starts[final_splitCounts==min(final_splitCounts)]
    minLoc <- bestStarts[1]
    
    ##############################################################################
    #Now that we've identified the best start, use that start location and all 
    #of the break points it creates to form the poorly-named variable "trueLines"
    #this is a data frame that is close to our final regions.
    ##############################################################################
    trueStarts <- seq(minLoc, the_last, by=maxLength)
    trueStarts <- c(max(trueStarts[1]-maxLength,0), trueStarts)
    trueEnds <- trueStarts + (maxLength-1)
    trueLines <- cbind(rep(as.character(seqnames(thisLine)), length(trueStarts)), trueStarts, trueEnds)
    trueLines <- as.data.frame(trueLines)
    trueLines[,1] <- as.character(trueLines[,1])
    trueLines[,2] <- as.numeric(as.character(trueLines[,2]))
    trueLines[,3] <- as.numeric(as.character(trueLines[,3]))
    colnames(trueLines) <- c("chr", "start", "end")
    
    ##############################################################################
    #Bound categories, as defined by Jacob.
    ##############################################################################
    bound_categories <- c("1-5", "6-20", "21-50", "51+")
    
    ##############################################################################
    #broken_Lines will actually be our final, returned list of regions.
    #why I defined a whole new variable for this is anybody's guess.  
    #I might have been coding tipsy.  I did this years ago.
    #
    #For each line of "trueLines", intersect that region with the peaks of interest.  
    #determine total number of unique TFs bound in that region and their 
    #identity.  Store this in broken_Lines.
    #convert broken_Lines to a data frame and give it proper column names.
    ##############################################################################
    broken_Lines <- c()
    for (j in 1:nrow(trueLines)) {
      true_gr <- makeGRangesFromDataFrame(trueLines[j,], keep.extra.columns=TRUE, ignore.strand=TRUE)
      theIntersect <- as.data.frame(findOverlaps(true_gr, thesePeaks))
      final_peaks <- thesePeaks[theIntersect[,2]]
      thisNum <- length(final_peaks$name)
      thisTFSet <- paste(unique(final_peaks$name), collapse=",")
      if(thisNum <=5 ) { thisCat <- bound_categories[1] }
      if(thisNum > 5 && thisNum <=20 ) { thisCat <- bound_categories[2] }
      if(thisNum > 20 && thisNum <= 50 ) { thisCat <- bound_categories[3] }
      if(thisNum > 50  ) { thisCat <- bound_categories[4] }
      broken_Lines <- rbind(broken_Lines, c(trueLines[j,], thisNum, thisCat, thisTFSet))
    }
    broken_Lines <- as.data.frame(broken_Lines)
    colnames(broken_Lines) <- c("seqnames", "start", "end", "count", "bound", "TF")
    
    ##############################################################################
    #Because of the way I defined the starts and ends, it's possible to have a weird 
    #situation where I've accidentally created regions to which nothing is bound, entirely 
    #outside of our peaks of interest.  This was necessary to make sure the 
    #start and end points of our broken regions lined up perfectly with the start 
    #and end points of the original region.  
    #If there are any regions in broken_Lines with nothing bound, remove them.
    ##############################################################################
    broken_Lines <- broken_Lines[broken_Lines$count!=0,]
    
    ##############################################################################
    #If the start of our first broken_Lines region is less than the start of our 
    #original region, discard the bps before our region starts.
    #If the end of our last broken_Lines region is less than the end of our original 
    #region, discard the bps after our region starts.
    #redetermine the lengths after this change.
    ##############################################################################
    broken_Lines$start <- as.numeric(as.character(broken_Lines$start))
    broken_Lines$end <- as.numeric(as.character(broken_Lines$end))  
    if(broken_Lines[1,2]<as.numeric(start(thisLine))) { broken_Lines[1,2] <- as.numeric(start(thisLine)) }
    if(broken_Lines[nrow(broken_Lines),3]>as.numeric(end(thisLine))) { broken_Lines[nrow(broken_Lines),3] <- as.numeric(end(thisLine)) }
    
    broken_Lines$length <- broken_Lines$end - broken_Lines$start + 1
    
    ##############################################################################
    #These were just some quality control tests to make sure none of the above 
    #screws anything up.  I want to make sure I didn't accidentally make a region with
    #negative length, and I wanted to make sure that the total length of our broken 
    #regions exactly matched the length of our original regions.  This worked for all 
    #regions in the DLPFC.  If you ever see the script print "PROBLEM", something went 
    #wrong.  You may or may not be on your own.
    ##############################################################################
    if(broken_Lines[1,2] > broken_Lines[1,3] ) {print("PROBLEM")}
    if(broken_Lines[1,2] > broken_Lines[2,2] ) {print("PROBLEM")}
    if(broken_Lines[nrow(broken_Lines),3] < broken_Lines[nrow(broken_Lines),2])  {print("PROBLEM")}
    if(broken_Lines[nrow(broken_Lines),3] < broken_Lines[nrow(broken_Lines)-1,3])  {print("PROBLEM")}
    if(sum(broken_Lines$length) != as.numeric(width(thisLine))) {print("PROBLEM")}
    
    ##############################################################################
    #return the broken-up regions.
    ##############################################################################
    return(broken_Lines)
  }
  
  
  
  ##############################################################################
  #restrict to the chromosome of interest
  ##############################################################################
  bad_Bed <- bad_Bed[seqnames(bad_Bed)==theChromosome]
  
  
  ##############################################################################
  #Set up the object that you'll eventually return.
  ##############################################################################
  theRetLines <- c()
  
  ##############################################################################
  #For each region in bad_Bed, break it up.
  #add the broken up regions to your return object.
  ##############################################################################
  theBreaks <- c(floor(length(bad_Bed)*seq(0.01,1, by=0.01)))
  #print(paste("Breaking regions of length > ", maxLength, " and determining TFs bound to subregions", sep=""))
  for (i in 1:length(bad_Bed)) {
    theseLines <- breakBins(bad_Bed[i], unmerged_Bed, maxLength)
    theRetLines <- rbind(theRetLines, theseLines)
    #if(i%in%theBreaks) {
    #  print(paste("Processed ", which(theBreaks==i), "% of ", length(bad_Bed), " regions", sep=""))
    #}
  }
  
  ##############################################################################
  #return the final object.
  ##############################################################################
  return(theRetLines)
  
}




################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################


################################################################################
#Load in the arguments from the command line.
#If you are not running from the command line, skip the next 5 lines
################################################################################
args <- commandArgs(trailingOnly=T)
merged_Bed_file <- args[1]
unmerged_Bed_file <- args[2]
if(length(args)>2) { maxLength <- as.numeric(as.character(args[3])) }
if(length(args)==2) { maxLength <- 2000 }



################################################################################
#If you're running this in Rstudio, use the following 3 lines (uncomment them), 
#modifying them to point to the correct paths or have the correct lengths.
################################################################################
#merged_Bed_file <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/Jacob_Merged_Peaks/union_DLPFC.rds"
#unmerged_Bed_file <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/Jacob_Merged_Peaks/unlist_DLPFC.rds"
#maxLength <- 2000


################################################################################
#This defines the final output file.
################################################################################
saveFile <- gsub(".rds", paste("_", maxLength, "_max_length.rds", sep=""), merged_Bed_file, fixed=T)


################################################################################
#Read in your RDS files.
################################################################################
merged_Bed <- readRDS(merged_Bed_file)
unmerged_Bed <- readRDS(unmerged_Bed_file)


################################################################################
#Identify the widths of your regions.
################################################################################
theWidths <- width(merged_Bed)


################################################################################
#Identify the regions that are smaller than or equal to your maxLength (2kb default).
################################################################################
good_Bed <- merged_Bed[theWidths<=maxLength]

#Ignore this line-- was used to speed up tests.
#good_Bed <- good_Bed[1:4000,]

################################################################################
#Determine the number and identify of the TFs bound in those regions. 
#There are more elegant and faster ways to do this. I use a for loop 
#because it lets me be sure of what's happening at each step and not 
#second-guess myself about whether I've done this right.
#takes a few minutes to run.
################################################################################
overlaps<-as.data.frame(findOverlaps(good_Bed, unmerged_Bed))
TF_list <- rep(NA, length(good_Bed))
nBound <- rep(0, length(good_Bed))
totalLength <- length(good_Bed)
theBreaks <- c(floor(totalLength*seq(0.1,1, by=0.1)))
print(paste("Determining the unique TFs bound to regions of length <= ", maxLength, sep=""))
for (i in 1:length(good_Bed)) {
  if(i%in%theBreaks) {
    print(paste("Processed ", which(theBreaks==i), "0 percent of ", totalLength, " regions", sep=""))
  }
  thisSet <- unmerged_Bed$name[overlaps[overlaps[,1]==i,2]]
  thisSet <- unique(thisSet)
  nBound[i] <- length(thisSet)
  thisSet <- paste(thisSet, collapse=",")
  TF_list[i] <- thisSet
}

bound <- rep("1-5", length(nBound))

bound[nBound>5] <- "6-20"
bound[nBound>20] <- "21-50"
bound[nBound>50] <- "51+"

good_Bed$count <- nBound 
good_Bed$bound <- bound
good_Bed$TF <- TF_list




################################################################################
#Identify the regions that are larger than your maxLength (2kb default).
################################################################################
bad_Bed <- merged_Bed[theWidths>maxLength]

#Ignore this line.  Was used to speed up tests.
#bad_Bed <- bad_Bed[1:100,]

################################################################################
#run the main function on these regions.  This will split them up into a 
#set of bins whose maximum size is 2kb, and which perfectly overlap with the 
#original start and end points of each region.  It'll identify the number 
#and identity of TFs in each split region.  
#This takes a couple hours to run on DLPFC, possibly longer for other tissues.
#If you make maxLength shorter, it will likely take significantly longer.
#It returns a data frame.
################################################################################


theChroms <- unique(as.character(as.data.frame(bad_Bed)[,1]))
#
#new_bad_Bed <- c()
#for (i in 1:length(theChroms)) {
#  thisSet <- as.data.frame(bad_Bed)[as.data.frame(bad_Bed)[,1]==theChroms[i],]
#  thisSet <- thisSet[1:20,]
#  new_bad_Bed <- rbind(new_bad_Bed, thisSet)
#}
#
#
numCores <- detectCores()
theClust <- makeCluster(numCores)

split_Bed_List <- clusterApply(theClust, x=theChroms, fun <- mainFunction, bad_Bed <- bad_Bed, unmerged_Bed <- unmerged_Bed, maxLength <- maxLength)

#split_Bed <- mainFunction(bad_Bed, unmerged_Bed, maxLength)

split_Bed <-   do.call("rbind", split_Bed_List)


################################################################################
#Merge the split_Bed and good_Bed back together.
#Since split_Bed is already a dataframe, I convery the good_Bed to a dataframe 
#to make it easy to combine and sort them.  The time difference here is minimal,
#and this way I knew what I was doing.
################################################################################
good_Bed_df <- as.data.frame(good_Bed)
good_Bed_df <- good_Bed_df[,c(1:3,6:8)]

final_Bed <- rbind(good_Bed_df, split_Bed[,1:6])
final_Bed <- final_Bed[order(as.character(final_Bed[,1]), as.numeric(as.character(final_Bed[,2]))),]

################################################################################
#convert that final_Bed back to a GenomicRanges object.
################################################################################
final_Bed_gr <- makeGRangesFromDataFrame(final_Bed, ignore.strand=T, keep.extra.columns=T)

################################################################################
#save the final Granges object.
################################################################################
saveRDS(final_Bed_gr, file = saveFile )








#
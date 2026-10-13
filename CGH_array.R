library(ggplot2)
library(dplyr)
require(DNAcopy)
library(data.table)

#Required Functions
#readCghAgilent : read CGH Agilent .txt files
#path : files path
readCghAgilent = function(path){
  
  #import library marray(need limma)
  require(marray)
  
  #ls	
  list_files = list.files(path,full.names=T)
  
  #read data	
  rawData = read.maimages(list_files, source = "agilent", columns = list(R = "rProcessedSignal", G = "gProcessedSignal"), annotation = c("Row","Col","ControlType", "ProbeName","SystematicName"), names = basename(list_files))
  
  #dye assignment
  
  #return rawData
  return(rawData)
}

readmRNAAgilentTxt = function(path){
  
  require(marray)
  
  list_files = list.files(path,full.names=T,pattern=".txt")
  
  #read data	
  rawData = read.maimages(list_files, source = "agilent", columns = list(R = "rProcessedSignal", G = "gProcessedSignal", LogRatio = "LogRatio",pVal="PValueLogRatio"), annotation = c("Row","Col","ControlType", "ProbeName","SystematicName","GeneName","Description"), names = basename(list_files))
  
  return(rawData)
  
}

readCghAgilentTxt = function(path){
  
  #import library marray(need limma)
  require(marray)
  
  #read data	
  rawData = read.maimages(path, source = "agilent", columns = list(R = "rProcessedSignal", G = "gProcessedSignal"), annotation = c("Row","Col","ControlType", "ProbeName","SystematicName","rBGSubSignal","gBGSubSignal","gMedianSignal","rMedianSignal","gBGMedianSignal","rBGMedianSignal"), names = basename(path))
  
  #dye assignment
  
  #return rawData
  return(rawData)
}

readCghAgilent_with_tkprogressBar = function(path){
  
  #import library marray(need limma)
  require(marray)
  
  rawData = list()
  
  #ls	
  list_files = list.files(path,full.names=T)[1:2]
  
  
  #progressBar
  pb = tkProgressBar(title="Loading data",label="Loading Data...",min=0,max=length(list_files))
  
  #read data
  for(i in 1:length(list_files)){	
    
    if(i==1){
      rawData = read.maimages(list_files[1], source = "agilent", columns = list(R = "rProcessedSignal", G = "gProcessedSignal"), annotation = c("Row","Col","ControlType", "ProbeName","GeneName","SystematicName"), names = basename(list_files[1]))
    }
    else if(i != 1){
      tmp = read.maimages(list_files[i], source = "agilent", columns = list(R = "rProcessedSignal", G = "gProcessedSignal"), annotation = c("Row","Col","ControlType", "ProbeName","GeneName","SystematicName"), names = basename(list_files[i]))
      rawData = cbind(rawData,tmp)
    }
    setTkProgressBar(pb,i)
  }
  
  
  close(pb)
  
  #return rawData
  return(rawData)
}




#normalize_cgh : background correction and normalization within arrays
#rawData : RGList object 
#ma:MAList object
normalize_cgh = function(rawData,method="none",backgroundCorrect=F){
  
  if(backgroundCorrect==T){
    ma <- normalizeWithinArrays(backgroundCorrect(rawData, method = "minimum"),method = method)
  }
  else {
    ma <- normalizeWithinArrays(rawData,method = method)	
  }
  return(ma)
  
}

main = function(segData,seuil=5){
  result = vector()
  id_sample = unique(segData$output$ID)
  for(i in 1:length(id_sample)){
    for(j in 1:24){
      segment.sample.chr = segData$output[which(segData$output$ID==id_sample[i] & segData$output$chrom==j),]	
      res = DNAcopy_to_matrix4.0(segment.sample.chr,seuil)
      result = rbind(result,res)
    }
  }
  result_matrix = seg_to_matrix(result)
  return(result_matrix)
}

# Convert SegFile to matrix 
#DNAcopy : an object of type DNAcopy(cf DNAcopy Package)
#CNA : an object of type CNA (cf DNAcopy Package)
#seuil : segments < will be merged
seg_to_matrix = function(segment){
  
  #get sample id from DNAcopy
  id_sample = unique(segment$ID)
  
  #matrix initialization
  segment_matrix = matrix(0,nrow=sum(segment[which(segment$ID==id_sample[1]),"num.mark"]),ncol=length(id_sample))
  
  #foreach samples
  for(i in 1:length(id_sample)){
    segment.sample = segment[which(segment$ID==id_sample[i]),]
    start = 1 
    end = 0
    for(j in 1:nrow(segment.sample)){
      end = end + segment.sample[j,"num.mark"]
      seg_mean = segment.sample[j,"seg.mean"]
      #cat(start,"\t",end,"\t",seg_mean,"\n")
      segment_matrix[start:end,i]=seg_mean
      start = end + 1
    }
  }
  
  #return matrix
  colnames(segment_matrix)=id_sample
  return(segment_matrix)
}

#segments with a size less than threshold are merge to the closest (in terms of distance), and a new segment value is calculated as the mean of the two merge segment.
DNAcopy_to_matrix4.0 = function(segment.sample.chr,seuil){
  segment.merge = segment.sample.chr
  res = vector()
  merge = TRUE
  i = 1
  while(i <= nrow(segment.sample.chr)){
    
    
    numMark = as.numeric(as.character(segment.sample.chr$num.mark))[i]
    if(numMark < seuil & i == 1){
      segment.sample.chr = merge_segment(segment.sample.chr,i,"next")
      i = 1
    }
    else if(numMark < seuil & i == nrow(segment.sample.chr)){
      segment.sample.chr = merge_segment(segment.sample.chr,i,"last")
      i = 1
    }
    else if(numMark < seuil){
      d_left = abs(as.numeric(as.character(segment.sample.chr$loc.end))[i-1] - as.numeric(as.character(segment.sample.chr$loc.start))[i])
      d_right = abs(as.numeric(as.character(segment.sample.chr$loc.start))[i+1] - as.numeric(as.character(segment.sample.chr$loc.end))[i])
      if(d_left <= d_right){
        segment.sample.chr = merge_segment(segment.sample.chr,i,"last")
      }
      else{
        segment.sample.chr = merge_segment(segment.sample.chr,i,"next")
      }
      i = 1
    }
    else{
      i = i + 1
    }
  }
  return(segment.sample.chr)
}

merge_segment = function(segment,i,merge="last"){
  if(merge=="last"){
    segment_merge = segment
    segment_merge[i-1,"loc.end"]=segment_merge[i,"loc.end"]
    segment_merge[i-1,"num.mark"]=segment_merge[i,"num.mark"]+segment_merge[i-1,"num.mark"]
    segment_merge[i-1,"seg.mean"] = (segment_merge[i-1,"seg.mean"]*segment_merge[i-1,"num.mark"] + segment_merge[i,"seg.mean"]*segment_merge[i,"num.mark"])/(segment_merge[i,"num.mark"]+segment_merge[i-1,"num.mark"])
    segment_merge = segment_merge[-i,]
  }
  else if(merge=="next"){
    segment_merge = segment
    segment_merge[i+1,"loc.end"]=segment_merge[i,"loc.end"]
    segment_merge[i+1,"num.mark"]=segment_merge[i,"num.mark"]+segment_merge[i+1,"num.mark"]
    segment_merge[i+1,"seg.mean"] = (segment_merge[i+1,"seg.mean"]*segment_merge[i+1,"num.mark"] + segment_merge[i,"seg.mean"]*segment_merge[i,"num.mark"])/(segment_merge[i,"num.mark"]+segment_merge[i+1,"num.mark"])
    segment_merge = segment_merge[-i,]
  }
  return(segment_merge)
}

#segments with a size less than threshold are merge to the closest (in terms of distance), and a new segment value is calculated as the mean of the two merge segment.
DNAcopy_to_matrix5.0 = function(segment.sample.chr,seuil){
  segment.merge = segment.sample.chr
  res = vector()
  merge = TRUE
  i = 1
  while(i <= nrow(segment.sample.chr)){
    
    
    numMark = as.numeric(as.character(segment.sample.chr$num.mark))[i]
    if(numMark < seuil & i == 1){
      #print(segment.sample.chr[i,])
      segment.sample.chr = merge_segment(segment.sample.chr,i,"next")
      i = 1
    }
    else if(numMark < seuil & i == nrow(segment.sample.chr)){
      #print(segment.sample.chr[i,])	
      segment.sample.chr = merge_segment(segment.sample.chr,i,"last")
      i = 1
    }
    else if(numMark < seuil){
      d_left = abs(as.numeric(as.character(segment.sample.chr$seg.mean))[i-1] - as.numeric(as.character(segment.sample.chr$seg.mean))[i])
      d_right = abs(as.numeric(as.character(segment.sample.chr$seg.mean))[i+1] - as.numeric(as.character(segment.sample.chr$seg.mean))[i])
      if(d_left <= d_right){
        segment.sample.chr = merge_segment2.0(segment.sample.chr,i,"last")
      }
      else{
        #print(segment.sample.chr[i,])
        segment.sample.chr = merge_segment2.0(segment.sample.chr,i,"next")
      }
      i = 1
    }
    else{
      i = i + 1
    }
  }
  return(segment.sample.chr)
}

merge_segment2.0 = function(segment,i,merge="last"){
  if(merge=="last"){
    segment_merge = segment
    segment_merge[i-1,"loc.end"]=segment_merge[i,"loc.end"]
    segment_merge[i-1,"num.mark"]=segment_merge[i,"num.mark"]+segment_merge[i-1,"num.mark"]
    #segment_merge[i-1,"seg.mean"] = segment_merge[i-1,"seg.mean"]
    segment_merge = segment_merge[-i,]
  }
  else if(merge=="next"){
    segment_merge = segment
    segment_merge[i+1,"loc.end"]=segment_merge[i,"loc.end"]
    segment_merge[i+1,"num.mark"]=segment_merge[i,"num.mark"]+segment_merge[i+1,"num.mark"]
    #segment_merge[i+1,"seg.mean"] = segment_merge[i+1,"seg.mean"]
    segment_merge = segment_merge[-i,]
  }
  return(segment_merge)
}

main2.0 = function(segData,seuil=5){
  result = vector()
  id_sample = unique(segData$output$ID)
  for(i in 1:length(id_sample)){
    cat(i,"\n")
    for(j in 1:23){
      segment.sample.chr = segData$output[which(segData$output$ID==id_sample[i] & segData$output$chrom==j),]	
      print(segment.sample.chr)
      res = DNAcopy_to_matrix5.0(segment.sample.chr,seuil)
      result = rbind(result,res)
    }
  }
  result_matrix = seg_to_matrix(result)
  return(result_matrix)
}

DNAcopy_mergeLevels <- function(cna,segments,pv.thres=0.0001,ansari.sign=0.05,thresMin=0.05,thresMax=0.5,verbose=1,scale=TRUE)
{
  library(aCGH)
  # ratios et segments sont 2 matrices probes*tumeurs
  
  
  mergedSegments <- c();
  
  for (i in  3:ncol(cna))
  {	cat(i,"\n")
    mergedValues <- mergeLevels(cna[,i],segments[,i-2],pv.thres=pv.thres,ansari.sign=ansari.sign,thresMin=thresMin,thresMax=thresMax,verbose=verbose,scale=scale)
    mergedSegments <- cbind(mergedSegments,mergedValues$vecMerged)
  }
  
  return(mergedSegments)
  
}


calcul_cytoband = function(Chip,file){
  #read file with cytoband info
  hg_table = read.table(file,sep="\t",header=T,comment.char="*")
  
  #chr_name
  chr_name = c(1:24)
  real_chr_name = c(1:22,"X","Y")
  
  #init var
  Chip$Cytoband = NA
  
  for( i in 1:length(chr_name)){
    #get cytoband for chr i
    cytoband_chr_i =hg_table[which(hg_table$chrom==paste("chr",real_chr_name[i],sep="")),]
    cat(i,"\n")
    #foreach cytoband of chr i
    for (j in 1:nrow(cytoband_chr_i)){
      #get start and end position
      start_cytoband = cytoband_chr_i$chromStart[j]
      end_cytoband = cytoband_chr_i$chromEnd[j]
      cat(j,"\n")
      #get index which are in cytoband j of chr i
      ind=which(Chip$ChrName==chr_name[i] & (as.numeric(as.matrix(Chip$Start)) >= start_cytoband & as.numeric(as.matrix(Chip$Start)) <= end_cytoband)	)
      #print unfound index
      #if(length(ind)==0){			
      #	print(cytoband_chr_i[j,])
      #}
      #store result in matrix
      Chip$Cytoband[ind]=as.character(as.matrix(cytoband_chr_i$name[j]))
    }
  }
  return(Chip)
}

segment_discretisation = function(ma_filtered,segment_matrix){
  source("function_R/analyse_genomique/segment_discretisation.R")
  source("function_R/analyse_genomique/sigma.R")
  matrice_discretise = matrix(nrow=nrow(segment_matrix),ncol=ncol(segment_matrix))
  chrom_name = gsub("chr([0-9XY]+):.*", "\\1", ma_filtered$genes[, "SystematicName"])
  ind_X_Y = which(chrom_name=="X" | chrom_name=="Y")
  for( i in 1:ncol(segment_matrix)){
    sigma = sigma(ma_filtered$M[,i])
    seuil_gain = sigma*1.75
    seuil_amp = sigma*7
    seuil_perte = sigma*-1.75
    
    ind_gain = which(segment_matrix[,i]>=seuil_gain &  segment_matrix[,i]<seuil_amp)
    ind_amp = which(segment_matrix[,i]>=seuil_amp)	
    ind_perte = which(segment_matrix[,i]<=seuil_perte)	
    
    matrice_discretise[,i]=0
    matrice_discretise[ind_perte,i]=-1
    matrice_discretise[ind_amp,i]=10
    matrice_discretise[ind_gain,i]=1
  }
  
  return(matrice_discretise)
  
}


#Load data

rawData = readCghAgilent(path)

ma = normalize_cgh(rawData)

chrom_name = gsub("chr([0-9XY]+):.*", "\\1", ma$genes[, "SystematicName"])
chrom_start = as.numeric(gsub(".*:([0-9]+)-.*", "\\1",ma$genes[, "SystematicName"]))
chrom_stop = as.numeric(gsub(".*-([0-9]+$)","\\1",ma$genes[,"SystematicName"]))
drop_index = c(which(!chrom_name %in% c(1:22,"X","Y")),which(ma$genes[,"ControlType"] != 0))
rawData_filtered = rawData[-drop_index,]
ma_filtered = ma[-drop_index,]
chrom_name=chrom_name[-drop_index]
chrom_start=chrom_start[-drop_index]
chrom_stop=chrom_stop[-drop_index]

duplicate = which(duplicated(rawData_filtered$genes$ProbeName))
rawData_filtered = rawData_filtered[-duplicate,]
ma_filtered = ma_filtered[-duplicate,]
chrom_start=chrom_start[-duplicate]
chrom_name=chrom_name[-duplicate]
chrom_stop=chrom_stop[-duplicate]

chrom_name[which(chrom_name=="X")]=23
chrom_name[which(chrom_name=="Y")]=24
chrom_name=as.numeric(chrom_name)

cna <- CNA(ma_filtered$M, chrom_name, chrom_start, data.type = "logratio",sampleid=colnames(ma_filtered$M))

segData <- segment(smooth.CNA(cna))

output <- segData$output

output$ID <- recode(output$ID,
                    "H.H.txt"="Unconverted",
                    "H.M.txt"="Converted to 2-M2high",
                    "H.T.txt"="Converted to 1-M2high",
                    "MCAM.txt"="2-M2high",
                    "MMEhigh.txt"="M2low",
                    "THY1.txt"="1-M2high")
seg_17 <- output[output$chrom == 17,]
seg_17$ID <- factor(seg_17$ID, levels = c("1-M2high", "2-M2high", "M2low", 
                                          "Unconverted", "Converted to 2-M2high", "Converted to 1-M2high"))


cna_plot <- cna
colnames(cna_plot) <- c("chrom", "maploc", "Unconverted", "Converted to 2-M2high", "Converted to 1-M2high", "2-M2high", "M2low", "1-M2high")
data_17 <- cna_plot[cna_plot$chrom==17,]
data_17 <- melt.data.table(as.data.table(data_17), id.vars = c("chrom", "maploc"), measure.vars = c("Unconverted", "Converted to 2-M2high", "Converted to 1-M2high",
                                                                                                    "2-M2high", "M2low", "1-M2high"))
colnames(data_17)[3] <- "ID"
data_17$ID <- factor(data_17$ID, levels = c("1-M2high", "2-M2high", "M2low", 
                                            "Unconverted", "Converted to 2-M2high", "Converted to 1-M2high"))

ggplot(seg_17, mapping = aes(x=maploc, y=value)) +
  geom_point(data=data_17, mapping = aes(x=maploc, y=value), size=0.1, alpha=0.1) +
  geom_segment(aes(x=loc.start, xend=loc.end, y=seg.mean, yend=seg.mean, color = ID), linewidth = 1) +
  facet_grid(ID~chrom, space = "free", scales = "free_x") + 
  theme_classic() +
  theme(panel.spacing = unit(0, 'mm'), 
        panel.border = element_rect(size = 0.5, color = 'grey10', fill = NA),
        legend.position = "none",
        axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        strip.background = element_blank(), strip.text = element_text(size=10),
        axis.line = element_line(size = .8), axis.ticks = element_line(size = .8),
        plot.title = element_text(size = 12)) + 
  ylim(c(-2, 2)) + 
  xlab('Genome position (bp)') +
  ylab("Log2 Ratio") +
  ggtitle("Chromosome 17") +
  geom_hline(yintercept = c(-1, 0, 1), linetype = "dashed")






seg_17_ctrl <- seg_17[seg_17$ID %in% c("M2low", "1-M2high", "2-M2high"),]
data_17_ctrl <- data_17[data_17$ID %in% c("M2low", "1-M2high", "2-M2high"),]

ggplot(seg_17_ctrl, mapping = aes(x=maploc, y=value)) +
  geom_point(data=data_17_ctrl, mapping = aes(x=maploc, y=value), size=0.2, alpha=0.1) +
  geom_segment(aes(x=loc.start, xend=loc.end, y=seg.mean, yend=seg.mean, color = ID), linewidth = 1) +
  facet_grid(ID~., space = "free", scales = "free_x") + 
  theme_classic() +
  theme(panel.spacing = unit(0, 'mm'), 
        panel.border = element_rect(size = 0.5, color = 'grey10', fill = NA),
        legend.position = "none",
        axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        strip.background = element_blank(), strip.text = element_text(size=10),
        axis.line = element_line(size = .8), axis.ticks = element_line(size = .8),
        plot.title = element_text(size = 12)) + 
  ylim(c(-2, 2)) + 
  xlab('Genome position (bp)') +
  ylab("Log2 Ratio") +
  ggtitle("Chromosome 17") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", size=0.5) +
  scale_color_manual(values = c("#fab6B6","#fc4c4c","#576db1"))






seg_17_5aza <- seg_17[seg_17$ID %in% c("Unconverted", "Converted to 2-M2high", "Converted to 1-M2high"),]
data_17_5aza <- data_17[data_17$ID %in% c("Unconverted", "Converted to 2-M2high", "Converted to 1-M2high"),]

ggplot(seg_17_5aza, mapping = aes(x=maploc, y=value)) +
  geom_point(data=data_17_5aza, mapping = aes(x=maploc, y=value), size=0.2, alpha=0.1) +
  geom_segment(aes(x=loc.start, xend=loc.end, y=seg.mean, yend=seg.mean, color = ID), linewidth = 1) +
  facet_grid(ID~chrom, space = "free", scales = "free_x") + 
  theme_classic() +
  theme(panel.spacing = unit(0, 'mm'), 
        panel.border = element_rect(size = 0.5, color = 'grey10', fill = NA),
        legend.position = "none",
        axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        strip.background = element_blank(), strip.text = element_text(size=10),
        axis.line = element_line(size = .8), axis.ticks = element_line(size = .8),
        plot.title = element_text(size = 12)) + 
  ylim(c(-2, 2)) + 
  xlab('Genome position (bp)') +
  ylab("Log2 Ratio") +
  ggtitle("Chromosome 17") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", size=0.5) +
  scale_color_manual(values = c("#576db1","#fab6B6","#fc4c4c"))

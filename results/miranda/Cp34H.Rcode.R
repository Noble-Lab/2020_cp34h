# R code documentation for Cp34H Long-term Extreme Subzero, High Salinity, Low Nutrients Incubations
# Author: Miranda Mudge, Erin Firth
# Date: 11/4/2020

####################################################################################################
####################################################################################################
####################################################################################################
# Figures
####################################################################################################

# Figure 1
# Number of Culturable Cells vs Time

# Data Files
# December Data.xlsx

library(tidyverse)
library(readxl)
library(ggthemes)
library(cowplot)
library(palmerpenguins)

#Pull data from Excel
alldata <- as_tibble(read_excel("December Data.xlsx", 
                                sheet = "Sheet1", col_types = c("text", "numeric", "numeric", 
                                                                "text", "text", "numeric", "text", 
                                                                "text", "text", "text", "numeric", 
                                                                "numeric", "numeric", "numeric", 
                                                                "numeric", "numeric")))
avgdata <- as_tibble(read_excel("December Data.xlsx", 
                                sheet = "avglines",  col_types = c("text", "numeric", "numeric", 
                                                                   "text", "text", "numeric", "text", 
                                                                   "text", "text", "text", "numeric", 
                                                                   "numeric", "numeric", "numeric", 
                                                                   "numeric", "numeric")))

#Subset data according to how MES had done (only "full_leu" and "no", unfrozen)
adata<-subset(alldata, 
              nutr=="full_leu" & tempnum!=-1 & froze=="no" | 
                nutr=="no" & tempnum!=-1 & froze=="no"     #| 
              #nutr=="full_thy" & tempnum!=-1 & froze=="no"
)
avdata<-subset(avgdata, 
               nutr=="full_leu" & tempnum!=-1 | 
                 nutr=="no" & tempnum!=-1       #| 
               #nutr=="full_thy" & tempnum!=-1
)

#Subset out data used for first black dots
inidata<-subset(alldata,tempnum==-1)
inidata2<-subset(alldata,tempnum==-1)
inidata$tempnum<-"-5"
inidata2$tempnum<-"-10"

#Factor and order data for use in facet grid
adata$tempnum<-factor(adata$tempnum,levels=c("-5","-10"))
avdata$tempnum<-factor(avdata$tempnum,levels=c("-5","-10"))
inidata$tempnum<-factor(inidata$tempnum,levels=c("-5","-10"))
inidata2$tempnum<-factor(inidata2$tempnum,levels=c("-5","-10"))

# New facet label names 
salinity.labs <- c("ASW", "Brine")
names(salinity.labs) <- c("ASW", "Equilibrium Brine")
temp.labs <- c("-5°C", "-10°C")
names(temp.labs) <- c("-5", "-10")

#Plot, labeled as "p" for the ability to use ggdraw
ggplot(data=adata)+
  #Log scale (MES data fixes all Zeroes to 0.1) with a little extra height for plot labels
  scale_y_log10()+
  coord_cartesian(ylim=c(1E-1,1E11),expand= TRUE)+
  #Jitter the data points slightly as requested, keeping same color and fill as MES)
  geom_jitter(data=filter(adata, nutr=="no" & salinity=="ASW"),mapping=aes(x=time_d,y=mpn),size=3.2, shape=1,color="dodgerblue",width=4,height=.5)+
  geom_jitter(data=filter(adata, nutr=="no" & salinity=="Equilibrium Brine"),mapping=aes(x=time_d,y=mpn),size=3.2,shape=1,color="orange",width=4,height=.5)+
  geom_jitter(data=filter(adata, nutr=="full_leu" & salinity=="ASW"),mapping=aes(x=time_d,y=mpn),size=3.2,shape=16,color="dodgerblue",width=4,height=.5)+
  geom_jitter(data=filter(adata, nutr=="full_leu" & salinity=="Equilibrium Brine"),mapping=aes(x=time_d,y=mpn),size=3.2,shape=16,color="orange",width=4,height=.5)+
  #Plot the black dots
  geom_point(data=inidata,mapping=aes(x=time_d,y=mpn),color="black",size=3.2)+
  geom_point(data=inidata2,mapping=aes(x=time_d,y=mpn),color="black",size=3.2)+
  #Plot the average lines
  geom_path(data=filter(avdata, nutr=="no"),mapping=aes(x=time_d,y=mpn),color="grey")+
  geom_path(data=filter(avdata, nutr=="full_leu"),mapping=aes(x=time_d,y=mpn),color="black")+
  #Grid (plus adding degree symbol to temperature)
  facet_grid(salinity~tempnum,labeller=labeller(salinity =salinity.labs,tempnum=temp.labs))+
  #Base theme plus adjustment to include lines in facet grid
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black",size=1),             
        strip.background = element_blank(),        
        strip.text = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15)) +
  labs(y=expression("Number of Culturable Cells (cells mL"^"-1"*")"), x="Time (days)")

####################################################################################################
####################################################################################################
####################################################################################################

# Figure 2
# NMDS of Sub-populations, Temperature, Salinity

# Data Files
# cp34h.eigan.envirfit.csv
# check.2019_April_30_NASA_Cp34H_STN_ABACUS_output.csv

library("BioStatR")
source('biostats.R')
library(dplyr)
library(tidyverse)
library(vegan)
library(methods)
library("grDevices")

abacus<-read.csv('check.2019_April_30_NASA_Cp34H_STN_ABACUS_output.csv', header=T, row.names=1) #load file
#subset file to ADJNSAF columns
nsaf<-select(abacus, contains('ADJNSAF')) #pull adjnsaf columns from abacus file
#keep only proteins with at least 2 unique peptides
numspec.uniq<-cbind(nsaf, abacus$ALL_NUMPEPSUNIQ) #grab numpepsuniq and bind to numspecadj
twopeps<-subset(numspec.uniq, select=X2019_APRIL_30_NASA_CP34H_ST18N_30_ADJNSAF:X2019_AUG_26_CP34H_STN_KCONT_10_ADJNSAF, numspec.uniq[,47]>1)
#keep only proteins that are not contaminants
twopeps$protein<-row.names(twopeps)
prot<-subset(twopeps, grepl(paste('Q', collapse="|"), twopeps$protein))
#shorten column names
colnames(prot)<-sub('X2019_APRIL_30_NASA_CP34H_ST', "", colnames(prot)) #remove first part from column name
colnames(prot)<-sub("_ADJNSAF", "", colnames(prot)) #remove end part from column name
#change typo column name and more modifications
names(prot)[names(prot)=='18N_30']<-'N18_30'
names(prot)[names(prot)=='N11_12']<-'N19_12'
colnames(prot)<-sub('X2019_AUG_26_CP34H_STN_', "", colnames(prot))
colnames(prot)<-sub('N', "", colnames(prot))
prot<-subset(prot, select=-protein) #remove column
prot2 <- prot[, c(15, 26:27, 31:39, 2:14, 1, 16:25, 28:30, 40:47)] #reorder to relevancy
prot2 <- prot2[, c(1:7, 9:15, 17:26, 16, 27:47)]
#transform data for plotting - uses vegan and source:biostats
prot.t<-t(prot2[,1:46]) 
prot.tra<-(prot.t+1)
prot.tra<-data.trans(prot.tra, method='log', plot=F)
nmds1<-metaMDS(prot.tra, distance='bray', k=2, trymax=100, autotransform=F)
#plot points with names
NMDS.adjNSAF<-ordiplot(nmds1, choices=c(1,2), type='text', display='sites', cex=0.5)
meta <- read.csv("cp34h.eigan.envirfit.csv", header = T)
meta.1<-meta[,-1]
rownames(meta.1)<-meta[,1]

# Plot Figure 2A
nmds_color<-ordiplot(nmds1,choices = c(1,2), type = 'none', display = 'sites', cex=0.5) 
points(nmds_color,'sites', 
       col=c(rep('dodgerblue',6), rep('orange',3), 
                                 rep('dodgerblue',9), rep('orange',4), 
                                 rep('dodgerblue',4), rep('orange',4), 
                                 rep('dodgerblue',4), rep('orange',4),rep('black',8)),
       bg=c(rep('gray100',6), 
            rep('gray100',3), 
            rep('dodgerblue',9), 
            rep('orange',4), 
            rep('gray100',4), 
            rep('gray100',4), 
            rep('dodgerblue',4), 
            rep('orange',4),
            rep('black',8)), 
       lwd=c(rep(2, 49)), 
       pch=c(rep(21,22),
             rep(25,16), 
             rep(21,8))) 

legend("bottomleft", legend=c("ASW -5°C","Brine -5°C","ASW -5°C Nutr","Brine -5°C Nutr","ASW -10°C","Brine -10°C","ASW -10°C Nutr","Brine -10°C Nutr","Controls"), 
       pch = c(21,21,21,21,25,25,25,25,21),
       col=c("dodgerblue","orange","dodgerblue","orange","dodgerblue","orange","dodgerblue","orange","black") ,
       pt.bg=c("gray100","gray100","dodgerblue","orange","gray100","gray100","dodgerblue","orange","black") ,
       cex = 1.0, inset = 0.015, pt.lwd = 2)

ordihull(
  nmds1,
  meta.1$Specific,
  kind = "ehull",
  display = "sites",
  draw = c("polygon"),
  col = c("gray99", "dodgerblue", "gray99", "orange", "black", "gray85"),
  alpha = 50,
  border = c("dodgerblue", "dodgerblue", "orange", "orange", "black", "gray40"),
  lty = c(2,2,2,2,2,2),
  lwd = 2.5
)

text(-0.109,0.074,labels="A", col = "Black", font = 2,  cex = 1.7)
text(-0.095,0.015,labels="Controls", col = "Black", font = 2,  cex = 1.0)
text(0.04,-0.013,labels="ASW -5°C Nutr", col = "Black", font = 2,  cex = 1.0)
text(0.045,-0.06,labels="ASW -5°C", col = "Black", font = 2,  cex = 1.0)
text(-0.022,0.012,labels="Brine -5°C", col = "Black", font = 2,  cex = 1.0)
text(0.098,0.012,labels="Brine -5°C Nutr", col = "Black", font = 2,  cex = 1.0)
text(0.035,0.058,labels="All -10°C", col = "Black", font = 2,  cex = 1.0)
text(0.12,-0.086,labels="Sub-populations", col = "Black", font = 2,  cex = 1.8)

# Plot Figure 2B
nmds_color<-ordiplot(nmds1,choices = c(1,2), type = 'none', display = 'sites', cex=1) 
points(nmds_color,'sites', 
       col=c(rep('dodgerblue',6), rep('orange',3), 
             rep('dodgerblue',9), rep('orange',4), 
             rep('dodgerblue',4), rep('orange',4), 
             rep('dodgerblue',4), rep('orange',4),rep('black',8)),
       bg=c(rep('gray100',6), 
            rep('gray100',3), 
            rep('dodgerblue',9), 
            rep('orange',4), 
            rep('gray100',4), 
            rep('gray100',4), 
            rep('dodgerblue',4), 
            rep('orange',4),
            rep('black',8)), 
       lwd=c(rep(2, 49)), 
       pch=c(rep(21,22),
             rep(25,16), 
             rep(21,8)),
       cex=c(rep(1.9,49)) )

ordihull(
  nmds1,
  meta.1$Temperature,
  display = "sites",
  draw = c("polygon"),
  col = c("#08306B", "#9ECAE1", "#DEEBF7"),
  #alpha = c(90, 70, 80), #1 #to achieve color, run ordihull 3x with alpha in sequence
  #alpha = c(50, 5, 5), #2
  #alpha = c(30, 5, 5), #3
  lty = c(2,2,2),
  lwd = 2.5
)

text(-0.103,0.097,labels="B", col = "Black", font = 2,  cex = 2.9)
text(-0.095,0.02,labels="-1°C", col = "Black", font = 2,  cex = 1.7)
text(0.025,0.07,labels="-10°C", col = "Black", font = 2,  cex = 1.7)
text(0.1,0.03,labels="-5°C", col = "Black", font = 2,  cex = 1.7)
text(0.095,-0.11,labels="Temperature", col = "Black", font = 2,  cex = 3.2)

# Plot Figure 2C
nmds_color<-ordiplot(nmds1,choices = c(1,2), type = 'none', display = 'sites', cex=0.5) 
points(nmds_color,'sites', 
       col=c(rep('dodgerblue',6), rep('orange',3), 
                                 rep('dodgerblue',9), rep('orange',4), 
                                 rep('dodgerblue',4), rep('orange',4), 
                                 rep('dodgerblue',4), rep('orange',4),rep('black',8)),
       bg=c(rep('gray100',6), 
            rep('gray100',3), 
            rep('dodgerblue',9), 
            rep('orange',4), 
            rep('gray100',4), 
            rep('gray100',4), 
            rep('dodgerblue',4), 
            rep('orange',4),
            rep('black',8)), 
       lwd=c(rep(2, 49)), 
       pch=c(rep(21,22),
             rep(25,16), 
             rep(21,8)),
       cex=c(rep(1.9,49)) )

ordihull(
  nmds1,
  meta.1$Salinity,
  display = "sites",
  draw = c("polygon"),
  #col = NULL,
  col = c("dodgerblue","orange"),
  alpha = 30,
  border = c("dodgerblue", "orange"),
  lty = c(2,2),
  lwd = 2.5
)

text(-0.103,0.097,labels="C", col = "Black", font = 2,  cex = 2.9)
text(-0.065,0.04,labels="ASW", col = "Black", font = 2,  cex = 1.8)
text(0.095,0.03,labels="Brine", col = "Black", font = 2,  cex = 1.8)
text(0.12,-0.11,labels="Salinity", col = "Black", font = 2,  cex = 3.2)

####################################################################################################
####################################################################################################
####################################################################################################

# Figure 3
# Upset and Venn diagram plots - Unique and shared proteins

# Data Files
# cp34h.upset.data_matrix.uniq_prot.csv #also "mat.prot"

# Perform code from line 112 - 131

library("BioStatR")
source('biostats.R')
library(dplyr)
library(tidyverse)
library(vegan)
library(methods)
library("grDevices")
library(UpSetR)

#make individual file for each condition
prot.A5 <-prot2[, c(1:6)]
prot.B5 <- prot2[, c(7:9)]
prot.A5N <-prot2[, c(10:18)]
prot.B5N <-prot2[, c(19:22)]
prot.A10 <-prot2[, c(23:26)]
prot.B10 <-prot2[, c(27:30)]
prot.A10N <-prot2[, c(31:34)]
prot.B10N <-prot2[, c(35:38)]
prot.KC <-prot2[, c(39:46)]

#filter each condition for at least 1 peptide in every sample
#filter out rows with all zeros
prot.A5 <- prot.A5[!rowSums(prot.A5[] == 0) >= 1,]  #1139
prot.A5N <- prot.A5N[!rowSums(prot.A5N[] == 0) >= 1,]    # 1165
prot.B5 <- prot.B5[!rowSums(prot.B5[] == 0) >= 1,]  #1238
prot.B5N <- prot.B5N[!rowSums(prot.B5N[] == 0) >= 1,]  #980
prot.A10 <- prot.A10[!rowSums(prot.A10[] == 0) >= 1,]  #1140
prot.A10N <- prot.A10N[!rowSums(prot.A10N[] == 0) >= 1,]    # 1178
prot.B10 <- prot.B10[!rowSums(prot.B10[] == 0) >= 1,]  #1203
prot.B10N <- prot.B10N[!rowSums(prot.B10N[] == 0) >= 1,]  #1172
prot.KC <- prot.KC[!rowSums(prot.KC[] == 0) >= 1,]  #1297

# Figure 3A

#make column of rownames = Protein column
prot.A5 <- tibble::rownames_to_column(prot.A5, "Protein")
prot.B5 <- tibble::rownames_to_column(prot.B5, "Protein")
prot.A5N <- tibble::rownames_to_column(prot.A5N, "Protein")
prot.B5N <- tibble::rownames_to_column(prot.B5N, "Protein")
prot.KC <- tibble::rownames_to_column(prot.KC, "Protein")

#make vectors of individual protein lists using protein column
prot.A5.prot <- dplyr::pull(prot.A5, Protein)
prot.B5.prot <- dplyr::pull(prot.B5, Protein)
prot.A5N.prot <- dplyr::pull(prot.A5N, Protein)
prot.B5N.prot <- dplyr::pull(prot.B5N, Protein)
prot.KC.prot <- dplyr::pull(prot.KC, Protein)

# Reference for code section below
# https://stackoverflow.com/questions/53575257/make-binary-presence-absence-data-matrix-from-multiple-lists-in-r
#########
# create a master list
master_list <- unique(c(prot.A5.prot, prot.B5.prot, prot.A5N.prot, prot.B5N.prot, prot.KC.prot))
prot.A5.prot <- c(prot.A5.prot, rep('ll', length(master_list) - length(prot.A5.prot)))
prot.B5.prot <- c(prot.B5.prot, rep('ll', length(master_list) - length(prot.B5.prot)))
prot.A5N.prot <- c(prot.A5N.prot, rep('ll', length(master_list) - length(prot.A5N.prot)))
prot.B5N.prot <- c(prot.B5N.prot, rep('ll', length(master_list) - length(prot.B5N.prot)))
prot.KC.prot <- c(prot.KC.prot, rep('ll', length(master_list) - length(prot.KC.prot)))
mat.prot <- matrix(c(as.integer(master_list %in% prot.A5.prot),
                     as.integer(master_list %in% prot.B5.prot),
                     as.integer(master_list %in% prot.A5N.prot),
                     as.integer(master_list %in% prot.B5N.prot),
                     as.integer(master_list %in% prot.KC.prot)),
                   nrow = length(master_list), 
                   dimnames = list(master_list))
#########

#created column names manually
colnames(mat.prot, do.NULL = FALSE)
colnames(mat.prot) <- c('ASW -5°C', "Brine -5°C", "ASW -5°C Nutr", "Brine -5°C Nutr", "Control" ) 
#turn matrix into data frame
mat.prot <- as.data.frame(mat.prot)
mat.prot <- tibble::rownames_to_column(mat.prot, "Protein")
#can remove protein column, doesn't matter
mat.prot.1<-mat.prot[,-1]
rownames(mat.prot.1)<-mat.prot[,1]

# Create plot
upset(mat.prot, 
      main.bar.color = "gray65", #top bar color
      matrix.color = "gray75", #points in grid color
      order.by = "freq", #order top bars by decreasing size - freq
      point.size=5,
      text.scale = 2,
      line.size = 1.2,
      shade.alpha = 0.1, #decrease background color in grid
      keep.order = TRUE,
      mainbar.y.label = "Number of Unique Shared Proteins Detected",
      sets = c("Brine -5°C","Brine -5°C Nutr","ASW -5°C","ASW -5°C Nutr", "Control"), #read bottom to top = backwards
      sets.bar.color=c("darkgoldenrod3","goldenrod1","deepskyblue2","lightskyblue","black"),
      
      queries = list(list(query = intersects, params = list( "ASW -5°C", "Brine -5°C"), color = "darkgreen", active = T), 
                     list(query = intersects, params = list( "ASW -5°C Nutr", "Brine -5°C Nutr"), color = "springgreen4", active = T), 
                     
                     list(query = intersects, params = list("Control"), color = "black", active = T), 
                     
                     list(query = intersects, params = list("ASW -5°C", "Brine -5°C","ASW -5°C Nutr", "Brine -5°C Nutr" ), color = "firebrick4", active = T), 
                     
                     list(query = intersects, params = list( "ASW -5°C","ASW -5°C Nutr"), color = "royalblue1", active = T), 
                     list(query = intersects, params = list("ASW -5°C Nutr"), color = "lightskyblue", active = T),
                     list(query = intersects, params = list("ASW -5°C"), color = "deepskyblue2", active = T),
                     
                     list(query = intersects, params = list( "Brine -5°C","Brine -5°C Nutr"), color = "tan1", active = T), 
                     list(query = intersects, params = list("Brine -5°C Nutr"), color = "goldenrod1", active = T),
                     list(query = intersects, params = list("Brine -5°C"), color = "darkgoldenrod3", active = T)
      ) 
) 

# Figure 3B
library(VennDiagram)

# Reference code
#https://stackoverflow.com/questions/51251689/impossible-error-message-in-r-for-venndiagram
venn.plot <- draw.quintuple.venn(
  area.vector = c(68, 7, 9, 4, 12, 8, 45, 4, 16, 3, 2, 13, 1, 3, 5, 3, 30, 8, 3, 39, 15, 9, 3, 2, 3, 16, 12, 31, 134, 8, 874),
  category = c("Control", "ASW -5°C Nutr", "ASW -5°C", "Brine -5°C Nutr", "Brine -5°C"),
  lty = "blank",
  fill = c("gray30", "lightskyblue1", "deepskyblue3", "goldenrod1", "darkgoldenrod4"),
  cat.col = c("gray30", "lightskyblue1", "deepskyblue3", "goldenrod1", "darkgoldenrod4"),
  direct.area = T,
  lwd = 0,
  cat.dist = 0.34,
  margin = 0.3,
  cat.cex = 2, 
  cat.just = rep(list(c(0.6, -0.9)), 5), 
  cex = c(2, 2, 2, 2, 2, 1.5, 1.3, 1.5, 1.3, 1.5, 1.3, 1.5, 1.3, 1.5, 1.3,
          1.5, 1.3, 1.5, 1.3, 1.5, 1.3, 1.5, 1.3, 1.5, 1.3, 1.5, 1.5, 1.5, 1.5, 1.5, 2)
  
)

####################################################################################################
####################################################################################################
####################################################################################################

# Figure 4
# Heatmap

library(tidyverse)
library(dplyr)
library(pheatmap)
library(RColorBrewer)


#read in files
annotations <- read.csv("34h_annotations.csv", header = TRUE)
nsaf <- read.csv("nsaf.csv", header = T)
prot.eggnog.fullset <- read.csv('volcano.eggnog.csv', header = T)

#grab relevant columns from dataset
eggnog<-subset(prot.eggnog.fullset, select=c('Protein', 'LogFoldChange.ASWneg5', 'LogFoldChange.ASWneg5Nutr', 
                                             'LogFoldChange.Brineneg5', 'LogFoldChange.Brineneg5Nutr'))
# select LFC +/- 2 using eggnog dataset of only LFC values at -5C, converts others to zero
lfc.2 <- eggnog %>% 
  mutate(LogFoldChange.ASWneg5 = replace(LogFoldChange.ASWneg5,between(LogFoldChange.ASWneg5, -2, 2), 0))
lfc.2 <- lfc.2 %>% 
  mutate(LogFoldChange.ASWneg5Nutr = replace(LogFoldChange.ASWneg5Nutr,between(LogFoldChange.ASWneg5Nutr, -2, 2), 0)) 
lfc.2 <- lfc.2 %>% 
  mutate(LogFoldChange.Brineneg5 = replace(LogFoldChange.Brineneg5, between(LogFoldChange.Brineneg5, -2, 2), 0)) 
lfc.2 <- lfc.2 %>% 
  mutate(LogFoldChange.Brineneg5Nutr = replace(LogFoldChange.Brineneg5Nutr, between(LogFoldChange.Brineneg5Nutr, -2, 2), 0)) 

##eliminate first column of proteins and set rownames to be protein names
lfc.2.0<-lfc.2[,-1]              
rownames(lfc.2.0)<-lfc.2[,1]     
#remove rows with only zeros
lfc.2.0.cut <- lfc.2.0[apply(lfc.2.0[,-1], 1, function(x) !all(x==0)),]  ##make sure no columns with words in dataframe
# make Protein column from rownames
lfc.2.0.prot <- cbind(rownames(lfc.2.0.cut), data.frame(lfc.2.0.cut, row.names=NULL))
names(lfc.2.0.prot)[1] <- "Protein"
lfc.2.0.annotations <- merge(lfc.2.0.prot, nsaf, by = "Protein", all.x = TRUE)
lfc.2.0.merge <- merge(lfc.2.0.prot, eggnog, by= "Protein", all.x = TRUE)
lfc.2.0.mer <- subset(lfc.2.0.merge, select = c("Protein", 'LogFoldChange.ASWneg5.y', 'LogFoldChange.ASWneg5Nutr.y', 
                                                'LogFoldChange.Brineneg5.y', 'LogFoldChange.Brineneg5Nutr.y'))
lfc.2.0.print <- merge(lfc.2.0.mer, annotations, by= "Protein", all.x = TRUE)

#use NMDS code to simplify column names
nsaf.lfc2.0<-subset(lfc.2.0.annotations, select=X2019_APRIL_30_NASA_CP34H_ST18N_30_ADJNSAF:X2019_AUG_26_CP34H_STN_KCONT_10_ADJNSAF)
colnames(nsaf.lfc2.0)<-sub('X2019_APRIL_30_NASA_CP34H_ST', "", colnames(nsaf.lfc2.0)) #remove first part from column name
colnames(nsaf.lfc2.0)<-sub("_ADJNSAF", "", colnames(nsaf.lfc2.0)) #remove end part from column name
#change typo column name and more modifications
names(nsaf.lfc2.0)[names(nsaf.lfc2.0)=='18N_30']<-'N18_30'
names(nsaf.lfc2.0)[names(nsaf.lfc2.0)=='N11_12']<-'N19_12'
colnames(nsaf.lfc2.0)<-sub('X2019_AUG_26_CP34H_STN_', "", colnames(nsaf.lfc2.0))
colnames(nsaf.lfc2.0)<-sub('N', "", colnames(nsaf.lfc2.0))
nsaf.lfc2.0.order <- nsaf.lfc2.0[, c(15, 26:27, 31:39, 2:14, 1, 16:25, 28:30, 40:47)] #reorder to relevancy
nsaf.lfc2.0.order <- nsaf.lfc2.0.order[, c(1:7, 9:15, 17:26, 16, 27:47)]
rownames(nsaf.lfc2.0.order)<-lfc.2.0.annotations[,1]     #set rownames to be protein names
#rename
nsaf.means.lfc.2.0 <- nsaf.lfc2.0.order

##### Make NSAF row means for each group
nsaf.means.lfc.2.0$AvgA5 = rowMeans(nsaf.means.lfc.2.0[,grepl('1_103|2_105|2_16	3_101|3_25|4_26', names(nsaf.means.lfc.2.0),)],na.rm=FALSE)
#
nsaf.means.lfc.2.0$AvgA5N = rowMeans(nsaf.means.lfc.2.0[,grepl('9_104|9_42|10_09|10_102|11_106|11_36|12_10_107|12_110|12_18', names(nsaf.means.lfc.2.0),)],na.rm=FALSE)
#
nsaf.means.lfc.2.0$AvgB5 = rowMeans(nsaf.means.lfc.2.0[,grepl('7_109|8_108|8_17', names(nsaf.means.lfc.2.0),)],na.rm=FALSE)
#
nsaf.means.lfc.2.0$AvgB5N = rowMeans(nsaf.means.lfc.2.0[,grepl('13_45|14_41|15_10|16_19', names(nsaf.means.lfc.2.0),)],na.rm=FALSE)
#
nsaf.means.lfc.2.0$AvgC = rowMeans(nsaf.means.lfc.2.0[,grepl('KCOT_03|KCOT_04|KCOT_05|KCOT_06|KCOT_07|KCOT_08|KCOT_09|KCOT_10', names(nsaf.means.lfc.2.0),)],na.rm=FALSE)

nsaf.avg.lfc2.0<-nsaf.means.lfc.2.0[,grepl('AvgA5|AvgA5N|AvgB5|AvgB5N|AvgC', names(nsaf.means.lfc.2.0),)]

#Normalize each value by the row mean
nsaf.means.lfc.2.0$normAvgA5 = (nsaf.means.lfc.2.0$AvgA5)/rowMeans(nsaf.means.lfc.2.0[,grepl('AvgA5|AvgA5N|AvgB5|AvgB5N|AvgC', 
                                                                                             names(nsaf.means.lfc.2.0),)], 
                                                                   na.rm =FALSE)
nsaf.means.lfc.2.0$normAvgA5N = (nsaf.means.lfc.2.0$AvgA5N)/rowMeans(nsaf.means.lfc.2.0[,grepl('AvgA5|AvgA5N|AvgB5|AvgB5N|AvgC', 
                                                                                               names(nsaf.means.lfc.2.0),)], 
                                                                     na.rm =FALSE)
nsaf.means.lfc.2.0$normAvgB5 = (nsaf.means.lfc.2.0$AvgB5)/rowMeans(nsaf.means.lfc.2.0[,grepl('AvgA5|AvgA5N|AvgB5|AvgB5N|AvgC', 
                                                                                             names(nsaf.means.lfc.2.0),)], 
                                                                   na.rm =FALSE)
nsaf.means.lfc.2.0$normAvgB5N = (nsaf.means.lfc.2.0$AvgB5N)/rowMeans(nsaf.means.lfc.2.0[,grepl('AvgA5|AvgA5N|AvgB5|AvgB5N|AvgC', 
                                                                                               names(nsaf.means.lfc.2.0),)], 
                                                                     na.rm =FALSE)
nsaf.means.lfc.2.0$normAvgC = (nsaf.means.lfc.2.0$AvgC)/rowMeans(nsaf.means.lfc.2.0[,grepl('AvgA5|AvgA5N|AvgB5|AvgB5N|AvgC', 
                                                                                           names(nsaf.means.lfc.2.0),)], 
                                                                 na.rm =FALSE)

norm.avg.lfc2.0<-nsaf.means.lfc.2.0[,grepl('normAvgA5|normAvgA5N|normAvgB5|normAvgB5N|normAvgC', names(nsaf.means.lfc.2.0),)]

#merge with annotations to get new rownames: gene id
header.2.0.end <- cbind(rownames(norm.avg.lfc2.0), data.frame(norm.avg.lfc2.0, row.names=NULL))
names(header.2.0.end)[1] <- "Protein"
name.cntrl.heat2.0 <- merge(header.2.0.end, annotations, by= "Protein", all.x = TRUE)

# Corrected file names
# Can start here
nsaf.means.lfc.2.0 <- read.csv("new.nsaf.heatmap.lfc2.0.eucclust12.cntrlsep.118prot.csv", header=T)
heat.set.2.0<-nsaf.means.lfc.2.0[,-1]              
rownames(heat.set.2.0)<-nsaf.means.lfc.2.0[,1] 
labels_row <- dplyr::pull(heat.set.2.0, heat.annotate)
#change column names
names(heat.set.2.0)[names(heat.set.2.0)=='normAvgC']<-'Control'
names(heat.set.2.0)[names(heat.set.2.0)=='normAvgA5N']<-'ASW -5°C Nutrients'
names(heat.set.2.0)[names(heat.set.2.0)=='normAvgA5']<-'ASW -5°C'
names(heat.set.2.0)[names(heat.set.2.0)=='normAvgB5N']<-'Brine -5°C Nutrients'
names(heat.set.2.0)[names(heat.set.2.0)=='normAvgB5']<-'Brine -5°C'
#rotate column names
draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

out<-pheatmap(heat.set.2.0[1:5], 
              color = colorRampPalette((brewer.pal(n = 5, name ="Blues")))(100),
              border_color = NA,
              cellwidth = 55, 
              cellheight = 10.5,
              angle_col="0",
              labels_row = labels_row,
              show_rownames=F, 
              fontsize_row = 10,
              fontsize_col = 10,
              show_colnames=T,
              clustering_distance_rows="euclidean",
              cluster_rows=T, 
              cluster_cols= T, 
              fontsize=5, 
              treeheight_row = 30,
              cutree_rows = 12, 
              cutree_cols = 2, 
              legend = T
)

# Generate Supplemental Data File
Clust_List.2.0.named<-heat.set.2.0[c(out$tree_row[['order']]), out$tree_col[['order']]]
clust.2.0.end <- cbind(rownames(Clust_List.2.0.named), data.frame(Clust_List.2.0.named, row.names=NULL))
names(clust.2.0.end)[1] <- "Protein"
clust.end.lfc2.0 <- left_join(clust.2.0.end, annotations, by = "Protein", copy = FALSE)
#write.csv(clust.end.lfc2.0, "nsaf.heatmap.lfc2.0.eucclust12.cntrlsep.118prot.csv")

####################################################################################################
####################################################################################################
####################################################################################################

# Figure 5
# Volcano plots

library(ggrepel)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(cowplot)
library(ggplot2)

merge.ann <- read.csv("volcano_plot.annotations.figure.csv", header = T)

#making eggNog colors
eggnog.names <- c(
  'grey90',
  "#332288",
  "forestgreen",
  #"#44AA99" ,
  "purple",
  "#88CCEE",
  "darkgoldenrod3",
  "#CC6677",
  #"#AA4499", 
  "plum2",
  "black"
  
)

c10 <- c(
  'Other' = "grey90",
  'Amino acid transport and metabolism'="#332288" , 
  'Energy production and conversion'="forestgreen" ,
  'Coenzyme transport and metabolism'="purple", 
  'Carbohydrate transport and metabolism'="#88CCEE" ,
  'Lipid transport and metabolism'="darkgoldenrod3",
  'Intracellular trafficking, secretion, and vesicular transport'="#CC6677",  
  'Cell motility'="plum2"   ,
  "Non" = "black"
  
)

breaks<-c(
  'Other',
  'Amino acid transport and metabolism',
  'Energy production and conversion',
  'Coenzyme transport and metabolism' ,
  'Carbohydrate transport and metabolism',
  'Lipid transport and metabolism',
  'Intracellular trafficking, secretion, and vesicular transport' , #
  'Cell motility', 
  "Non"
  
)

# Figure 5A
#ASWneg5Nutr 
merge.ann.A5N <- read.csv("volcano_plot.annotations.figure.full.A5N.csv", header = T)
merge.ann.A5N <- merge.ann.A5N[order(-as.numeric(factor(merge.ann.A5N$Gene.Ontology))),]
A5N.only <- read.csv("volcano_plot.annotations.figure.A5N.csv", header = TRUE)

high = merge.ann.A5N$LogFoldChange.ASWneg5Nutr>=0.5 & abs(merge.ann.A5N$Zstatistic.ASWneg5Nutr) >=2
low =  merge.ann.A5N$LogFoldChange.ASWneg5Nutr<=-0.5 & abs(merge.ann.A5N$Zstatistic.ASWneg5Nutr) >=2
volcano_ASWneg5Nutr <- ggplot(merge.ann.A5N, aes(x = LogFoldChange.ASWneg5Nutr, y = abs(Zstatistic.ASWneg5Nutr),  col=Gene.Ontology)) #fill keeps legend to right, but not colored legend
volcano_ASWneg5Nutr  + geom_point(aes(color = ifelse(high |low , eggnog.names, "black"), shape = 16))+  #keeps legend correct
  geom_point(size = 3)+ #point size, leave aes out, need line to add points
  scale_shape_identity()+
  scale_y_continuous(breaks=seq(0,40,10), limits = c(0,40))+
  scale_x_continuous(breaks=seq(-4,4,2), limits = c(-4,4.5))+
  scale_color_manual(values = ifelse(high|low , c10, "black"), limits = breaks)+ 
  xlab("Log2 Fold Change") +
  ylab("Z-score")+
  ggtitle("Control vs ASW -5°C Nutrients")+
  theme_classic(base_size = 10) +
  theme(legend.position = "none", text = element_text(size=18)) +
  geom_hline(yintercept = 2, linetype= "dashed")+
  geom_vline (xintercept = 0.5, linetype='dashed')+
  geom_vline (xintercept = -0.5, linetype='dashed') +  
  geom_label_repel(data=A5N.only, #separate dataset for specific labels
                   aes(label=A5N, size = 10), #call protein name, font size?
                   force        = 5, #how far point is from bubble
                   nudge_x      = 0.2, #push along x axis
                   direction    = "both", #bubbles move x and y directions
                   hjust        = 1, # orientation to the right
                   segment.size = 0.5, #line to point thickness
                   colour = A5N.only$color, #color of bubble
                   label.size = 0.5, #bubble edge thickness
                   label.r = 0.5, #shape of bubble roundness
                   box.padding = 0.8) #repel boxes from each other

# Figure 5B
# ASWneg5 
merge.ann.A5 <- read.csv("volcano_plot.annotations.figure.full.A5.csv", header = T)
merge.ann.A5 <- merge.ann.A5[order(-as.numeric(factor(merge.ann.A5$Functional.Annotation))),]
A5.only <- read.csv("volcano_plot.annotations.figure.A5.edits.csv", header = TRUE)

high = merge.ann.A5$LogFoldChange.ASWneg5>=0.5 & abs(merge.ann.A5$Zstatistic.ASWneg5) >=2
low =  merge.ann.A5$LogFoldChange.ASWneg5<=-0.5 & abs(merge.ann.A5$Zstatistic.ASWneg5) >=2
volcano_ASWneg5 <- ggplot(merge.ann.A5, aes(x = LogFoldChange.ASWneg5, y = abs(Zstatistic.ASWneg5),  col=Functional.Annotation)) #fill keeps legend to right, but not colored legend
volcano_ASWneg5  + geom_point(aes(color = ifelse(high |low , eggnog.names, "black"), shape = 16))+  #keeps legend correct
  geom_point(size = 3)+ #point size, leave aes out, need line to add points
  scale_shape_identity()+
  scale_y_continuous(breaks=seq(0,40,10), limits = c(0,40))+
  scale_x_continuous(breaks=seq(-4,4,2), limits = c(-4.6,4.5))+
  scale_color_manual(values = ifelse(high|low , c10, "black"), limits = breaks)+ 
  xlab("Log2 Fold Change") +
  ylab("Z-score")+
  ggtitle("Control vs ASW -5°C")+
  theme_classic(base_size = 10) +
  theme(legend.position = "right", text = element_text(size=18)) + 
  geom_hline(yintercept = 2, linetype= "dashed")+
  geom_vline (xintercept = 0.5, linetype='dashed')+
  geom_vline (xintercept = -0.5, linetype='dashed') +  
  geom_label_repel(data=A5.only, #separate dataset for specific labels
                   aes(label=A5, size = 10), #call protein name, font size?
                   force        = 5, #how far point is from bubble
                   nudge_x      = 0.2, #push along x axis
                   direction    = "both", #bubbles move x and y directions
                   hjust        = 1, # orientation to the right
                   segment.size = 0.5, #line to point thickness
                   colour = A5.only$color, #color of bubble
                   label.size = 0.5, #bubble edge thickness
                   label.r = 0.5, #shape of bubble roundness
                   box.padding = 0.8) #repel boxes from each other

# Figure 5C
# Brineneg5Nutr 
merge.ann.B5N <- read.csv("volcano_plot.annotations.figure.full.B5N.csv", header = T)
merge.ann.B5N <- merge.ann.B5N[order(-as.numeric(factor(merge.ann.B5N$Gene.Ontology))),]
B5N.only <- read.csv("volcano_plot.annotations.figure.B5N.csv", header = TRUE)
B5N.pos <- read.csv("volcano_plot.annotations.figure.B5N.pos.csv", header = TRUE)
B5N.neg <- read.csv("volcano_plot.annotations.figure.B5N.neg.csv", header = TRUE)

high = merge.ann.B5N$LogFoldChange.Brineneg5Nutr>=0.5 & abs(merge.ann.B5N$Zstatistic.Brineneg5Nutr) >=2
low =  merge.ann.B5N$LogFoldChange.Brineneg5Nutr<=-0.5 & abs(merge.ann.B5N$Zstatistic.Brineneg5Nutr) >=2
volcano_Brineneg5Nutr <- ggplot(merge.ann.B5N, aes(x = LogFoldChange.Brineneg5Nutr, y = abs(Zstatistic.Brineneg5Nutr),  col=Gene.Ontology)) #fill keeps legend to right, but not colored legend
volcano_Brineneg5Nutr  + geom_point(aes(color = ifelse(high |low , eggnog.names, "black"), shape = 16))+  #keeps legend correct
  geom_point(size = 3)+ #point size, leave aes out, need line to add points
  scale_shape_identity()+
  scale_y_continuous(breaks=seq(0,40,10), limits = c(0,40))+
  scale_x_continuous(breaks=seq(-4,4,2), limits = c(-4,4))+
  scale_color_manual(values = ifelse(high|low , c10, "black"), limits = breaks)+ 
  xlab("Log2 Fold Change") +
  ylab("Z-score")+
  ggtitle("Control vs Brine -5°C Nutrients")+
  theme_classic(base_size = 10) +
  theme(legend.position = "none", text = element_text(size=18)) + 
  geom_hline(yintercept = 2, linetype= "dashed")+
  geom_vline (xintercept = 0.5, linetype='dashed')+
  geom_vline (xintercept = -0.5, linetype='dashed') +  
  #separate pos and neg labels for easier visualization
  geom_label_repel(data=B5N.pos, #separate dataset for specific labels
                   aes(label=B5N, size = 10), #call protein name, font size?
                   force        = 3, #how far point is from bubble
                   nudge_x      = 0.5, #push along x axis
                   nudge_y = 0.5,
                   direction    = "both", #bubbles move x and y directions
                   hjust        = 1, # orientation to the right
                   vjust = -1,
                   segment.size = 0.5, #line to point thickness
                   colour = B5N.pos$color, #color of bubble
                   label.size = 0.5, #bubble edge thickness
                   label.r = 0.5, #shape of bubble roundness
                   box.padding = 0.7)  + #repel boxes from each other 
  geom_label_repel(data=B5N.neg, #separate dataset for specific labels
                   aes(label=B5N, size = 10), #call protein name, font size?
                   force        = 5, #how far point is from bubble
                   nudge_x      = 0.5, #push along x axis
                   nudge_y = 0.5,
                   direction    = "y", #bubbles move x and y directions
                   vjust = -1,
                   segment.size = 0.5, #line to point thickness
                   colour = B5N.neg$color, #color of bubble
                   label.size = 0.5, #bubble edge thickness
                   label.r = 0.5, #shape of bubble roundness
                   box.padding = 0.2) #repel boxes from each other

# Figure 5D
# Brineneg5 
merge.ann.B5 <- read.csv("volcano_plot.annotations.figure.full.B5.csv", header = T)
merge.ann.B5 <- merge.ann.B5[order(-as.numeric(factor(merge.ann.B5$Gene.Ontology))),]
B5.only <- read.csv("volcano_plot.annotations.figure.B5.csv", header = TRUE)

high = merge.ann.B5$LogFoldChange.Brineneg5>=0.5 & abs(merge.ann.B5$Zstatistic.Brineneg5) >=2
low =  merge.ann.B5$LogFoldChange.Brineneg5<=-0.5 & abs(merge.ann.B5$Zstatistic.Brineneg5) >=2
volcano_Brineneg5 <- ggplot(merge.ann.B5, aes(x = LogFoldChange.Brineneg5, y = abs(Zstatistic.Brineneg5),  col=Gene.Ontology)) #fill keeps legend to right, but not colored legend
volcano_Brineneg5  + geom_point(aes(color = ifelse(high |low , eggnog.names, "black"), shape = 16))+  #keeps legend correct
  geom_point(size = 3)+ #point size, leave aes out, need line to add points
  scale_shape_identity()+
  scale_y_continuous(breaks=seq(0,40,10), limits = c(0,40))+
  scale_x_continuous(breaks=seq(-4,4,2), limits = c(-4,4))+
  scale_color_manual(values = ifelse(high|low , c10, "black"), limits = breaks)+ 
  xlab("Log2 Fold Change") +
  ylab("Z-score")+
  ggtitle("Control vs Brine -5°C")+
  theme_classic(base_size = 10) +
  theme(legend.position = "none", text = element_text(size=18)) + 
  geom_hline(yintercept = 2, linetype= "dashed")+
  geom_vline (xintercept = 0.5, linetype='dashed')+
  geom_vline (xintercept = -0.5, linetype='dashed') +  
  geom_label_repel(data=B5.only, #separate dataset for specific labels
                   aes(label=B5, size = 10), #call protein name, font size?
                   force        = 10, #how far point is from bubble
                   nudge_x      = 0.8, #push along x axis
                   direction    = "both", #bubbles move x and y directions
                   hjust        = 1, # orientation to the right
                   segment.size = 0.5, #line to point thickness
                   colour = B5.only$color, #color of bubble
                   label.size = 0.5, #bubble edge thickness
                   label.r = 0.5, #shape of bubble roundness
                   box.padding = 0.6) #repel boxes from each other

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
# Supplemental Figures
####################################################################################################
# Supplemental Figure S2
# Cell abundance over time per condition

#Pull in necessary packages
library(tidyverse)
library(readxl)
library(ggthemes)
library(cowplot)

#Pull data from Excel
alldata <- as_tibble(read_excel("December Data.xlsx", 
                                sheet = "Sheet1", col_types = c("text", "numeric", "numeric", 
                                                                "text", "text", "numeric", "text", 
                                                                "text", "text", "text", "numeric", 
                                                                "numeric", "numeric", "numeric", 
                                                                "numeric", "numeric")))
avgdata <- as_tibble(read_excel("December Data.xlsx", 
                                sheet = "avglines",  col_types = c("text", "numeric", "numeric", 
                                                                   "text", "text", "numeric", "text", 
                                                                   "text", "text", "text", "numeric", 
                                                                   "numeric", "numeric", "numeric", 
                                                                   "numeric", "numeric")))

#Subset data according to how MES had done (only "full_leu" and "no", unfrozen)
adata<-subset(alldata, nutr=="full_leu" & tempnum!=-1 & froze=="no" | nutr=="no" & tempnum!=-1 & froze=="no")
avdata<-subset(avgdata, nutr=="full_leu" & tempnum!=-1 | nutr=="no" & tempnum!=-1)

#Subset out data used for first black dots
inidata<-subset(alldata,tempnum==-1)
inidata2<-subset(alldata,tempnum==-1)
inidata$tempnum<-"-5"
inidata2$tempnum<-"-10"

#Factor and order data for use in facet grid
adata$tempnum<-factor(adata$tempnum,levels=c("-5","-10"))
avdata$tempnum<-factor(avdata$tempnum,levels=c("-5","-10"))
inidata$tempnum<-factor(inidata$tempnum,levels=c("-5","-10"))
inidata2$tempnum<-factor(inidata2$tempnum,levels=c("-5","-10"))

# New facet label names 
salinity.labs <- c("ASW", "Brine")
names(salinity.labs) <- c("ASW", "Equilibrium Brine")
temp.labs <- c("-5°C", "-10°C")
names(temp.labs) <- c("-5", "-10")

#Plot, labeled as "p" for the ability to use ggdraw
ggplot(data=adata)+
  #Log scale (MES data fixes all Zeroes to 0.1) with a little extra height for plot labels
  scale_y_log10()+
  coord_cartesian(ylim=c(1E5,1E9),expand= TRUE)+
  #Jitter the data points slightly as requested, keeping same color and fill as MES)
  geom_jitter(data=filter(adata, nutr=="no" & salinity=="ASW"),mapping=aes(x=time_d,y=cellnum),size=3.2, shape=1,color="dodgerblue",width=4,height=.5)+
  geom_jitter(data=filter(adata, nutr=="no" & salinity=="Equilibrium Brine"),mapping=aes(x=time_d,y=cellnum),size=3.2,shape=1,color="orange",width=4,height=.5)+
  geom_jitter(data=filter(adata, nutr=="full_leu" & salinity=="ASW"),mapping=aes(x=time_d,y=cellnum),size=3.2,shape=16,color="dodgerblue",width=4,height=.5)+
  geom_jitter(data=filter(adata, nutr=="full_leu" & salinity=="Equilibrium Brine"),mapping=aes(x=time_d,y=cellnum),size=3.2,shape=16,color="orange",width=4,height=.5)+
  #geom_jitter(data=filter(adata, nutr=="full_thy" & salinity=="ASW"),mapping=aes(x=time_d,y=cellnum),size=2,shape=16,color="green",width=4,height=.5)+
  #geom_jitter(data=filter(adata, nutr=="full_thy" & salinity=="Equilibrium Brine"),mapping=aes(x=time_d,y=cellnum),size=2,shape=16,color="green",width=4,height=.5)+
  #Plot the black dots
  geom_point(data=inidata,mapping=aes(x=time_d,y=cellnum),color="black",size=3.2)+
  geom_point(data=inidata2,mapping=aes(x=time_d,y=cellnum),color="black",size=3.2)+
  #Plot the average lines
  geom_path(data=filter(avdata, nutr=="no"),mapping=aes(x=time_d,y=cellnum),color="grey")+
  geom_path(data=filter(avdata, nutr=="full_leu"),mapping=aes(x=time_d,y=cellnum),color="black")+
  #geom_path(data=filter(avdata, nutr=="full_thy"),mapping=aes(x=time_d,y=cellnum),color="green")+
  #Grid (plus adding degree symbol to temperature)
  facet_grid(salinity~tempnum,labeller=labeller(salinity =salinity.labs,tempnum=temp.labs))+
  #Base theme plus adjustment to include lines in facet grid
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black",size=1),             
        strip.background = element_blank(),        
        strip.text = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15)) +
  labs(y=expression("Cell Abundance (cells mL"^"-1"*")"), x="Time (days)")

####################################################################################################
####################################################################################################
####################################################################################################

# Supplemental Figure S3
# Number culturable cells over time

library(tidyverse)
library(readxl)
library(cowplot)

#Change the path here to match the location of the file in your machine
tb<-as_tibble(read_excel("Master Growth Curve Analysis.xlsx",
                         sheet = "34H Dec RData"))
#dec is -1C in DecASW + Glu+YE
dec<-filter(tb,Nutr=="Glu+YE")
#asw is 6C in ASW + GLVT
asw<-filter(tb,Nutr=="GLVT")

#hline and rect are hardcoded T0 count values +/- stdev
ggplot(dec,aes(Time,MPN))+
  coord_cartesian(xlim=c(0,20),ylim=c(1e4,1e9))+
  scale_y_log10()+
  geom_hline(aes(yintercept=1.96e6),color="blue")+
  geom_hline(aes(yintercept=1.36e6),color="blue")+
  geom_rect(aes(xmin=-20,xmax=40,ymin=1.36e6-6.96e5,ymax=1.96e6+6.91e5),fill="blue",alpha=0.03)+
  geom_point()+
  geom_errorbar(mapping=aes(ymin=minstdev,ymax=maxstdev))+
  theme_bw() +
  labs(y=expression("Number of Culturable Cells (cells mL"^"-1"*")"), x="Time (days)")

####################################################################################################
####################################################################################################
####################################################################################################
# Supplemental Figure S4
# Counts of proteins in each condition sorted by eggNog terms

library(RColorBrewer)
library(pheatmap)

# Rotate x axis labels
# https://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

# Code for number of proteins increased
inc <- read.csv("counts.eggnog.inc.csv", header = T)
inc <- subset(inc, select = c(1:9))

names(inc)[names(inc)=='A5.inc']<-'ASW -5°C'
names(inc)[names(inc)=='A5N.inc']<-'ASW -5°C Nutrients'
names(inc)[names(inc)=='B5.inc']<-'Brine -5°C'
names(inc)[names(inc)=='B5N.inc']<-'Brine -5°C Nutrients'
names(inc)[names(inc)=='A10.inc']<-'ASW -10°C'
names(inc)[names(inc)=='A10N.inc']<-'ASW -10°C Nutrients'
names(inc)[names(inc)=='B10.inc']<-'Brine -10°C'
names(inc)[names(inc)=='B10N.inc']<-'Brine -10°C Nutrients'

inc.plot<-inc[,-1]
rownames(inc.plot)<-inc[,1]
format(inc,scientific=0)

hm.col.inc<-brewer.pal(9,'Reds')
heat1 <- pheatmap(inc.plot, display_numbers = T, 
                  cluster_rows = F, cluster_cols = F, 
                  color = hm.col.inc, show_rownames = F, 
                  treeheight_row = 0, treeheight_col = 0, legend = T, 
                  cellwidth=30, cellheight=15, number_format = "%.0f",
                  number_color = "grey65")

# Code for number of proteins decreased
dec <- read.csv("counts.eggnog.dec.csv", header = T)
dec <- subset(dec, select = c(1:9))
names(dec)[names(dec)=='A5.dec']<-'ASW -5°C'
names(dec)[names(dec)=='A5N.dec']<-'ASW -5°C Nutrients'
names(dec)[names(dec)=='B5.dec']<-'Brine -5°C'
names(dec)[names(dec)=='B5N.dec']<-'Brine -5°C Nutrients'
names(dec)[names(dec)=='A10.dec']<-'ASW -10°C'
names(dec)[names(dec)=='A10N.dec']<-'ASW -10°C Nutrients'
names(dec)[names(dec)=='B10.dec']<-'Brine -10°C'
names(dec)[names(dec)=='B10N.dec']<-'Brine -10°C Nutrients'

dec.plot<-dec[,-1]
rownames(dec.plot)<-dec[,1]

hm.col.dec<-brewer.pal(9,'Blues')
heat2 <- pheatmap(dec.plot, display_numbers = T, color = hm.col.dec, 
                  cluster_rows = F, cluster_cols = F,
                  treeheight_row = 0, treeheight_col = 0, legend = T, 
                  cellwidth=30, cellheight=15, number_format = "%.0f")



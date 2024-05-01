#code for generating elemental depth records from micro-XRF generated maps
#Dynamic Time Warping (DTW) has been applied on these depth records
#This code is for the Steinbruch Schmidt Usseln Limestone (ULS)
#ULS consists of two sections, ULS 1 and ULS 2

#steps within this code:
#1) removing areas without data using k-means clustering
#2) plot all individual segments together as a composite map
#3) determine overlap between segments and combine them (ULS1 and ULS2)
#4) average all rows in order to generate a composite depth record

library(tidyverse)
library(astrochron)
library(devtools)
library(ptw)
library(RColorBrewer)

setwd("insert wd")

#--------------Initial data checking and formatting--------------------------------

#Load in data:

ULS1_Ca<-read.csv("ULS1_Ca.csv", header=T)
ULS1_Fe<-read.csv("ULS1_Fe.csv", header=T)
ULS1_Ti<-read.csv("ULS1_Ti.csv", header=T)
ULS1_K<-read.csv("ULS1_K.csv", header=T)
ULS1_Mn<-read.csv("ULS1_Mn.csv", header=T)
ULS1_Si<-read.csv("ULS1_Si.csv", header=T)

ULS2_Ca<-read.csv("ULS2_Ca.csv", header=T)
ULS2_Fe<-read.csv("ULS2_Fe.csv", header=T)
ULS2_Ti<-read.csv("ULS2_Ti.csv", header=T)
ULS2_K<-read.csv("ULS2_K.csv", header=T)
ULS2_Mn<-read.csv("ULS2_Mn.csv", header=T)
ULS2_Si<-read.csv("ULS2_Si.csv", header=T)


#plot as heatmaps for inspection (optional):
heatmap(as.matrix(ULS1_Ca), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULS1_Fe), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULS1_K), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULS1_Ti), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULS1_Si), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULS1_Mn), Rowv=NA, Colv=NA, asp=1)

heatmap(as.matrix(ULS2_Ca), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULS2_Fe), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULS2_K), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULS2_Ti), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULS2_Si), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULS2_Mn), Rowv=NA, Colv=NA, asp=1)



#-------------------------K-means clustering--------------------------------------- 
#sum and sqrt all different elemental maps to increase the contrast between rock and non-rock area
ULS1_sum<-sqrt(ULS1_Ca+ULS1_Fe+ULS1_K+ULS1_Ti+ULS1_Si+ULS1_Mn)
ULS2_sum<-sqrt(ULS2_Ca+ULS2_Fe+ULS2_K+ULS2_Ti+ULS2_Si+ULS2_Mn)

#visually inspect the contrast level (optional):
heatmap(as.matrix(ULS1_sum), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULS2_sum), Rowv=NA, Colv=NA, asp=1)




#-------------------------------ULS1----------------------------------------------
#Gather dataframe into single column
ULS1_g<-gather(ULS1_sum, col, val, 1:1170, factor_key=T)
ULS1_g[is.na(ULS1_g)]<-0

#carry out kmeans clustering:
ULS1_g_km <- kmeans(ULS1_g[,2],2,iter.max = 100, nstart = 20)

#spread dataset again with just clusters to check whether clustering worked and
#which cluster nr corresponds to the rock surface
S_1<-cbind(ULS1_g, ULS1_g_km$cluster)
colnames(S_1)[3]<-"cluster"
St_1<-S_1[,c(1,3)] #check whether the clustering works. If it does, turn selected value into NA
Sp_1<-pivot_wider(St_1, names_from=col, values_from=cluster)
Sp_1<-unnest(Sp_1) 
heatmap(as.matrix(Sp_1), Rowv=NA, Colv=NA, asp=1)
view(Sp_1) #cluster 1 = rock surface, cluster 2 = non-rock
#NB pay attention to this when generating your own data, might switch around due to the iterative nature of kmeans!



#create new dataframe with original data for rock surface and NA for non-rock surface
#for each elemental map:
#Ca:
ULS1_Ca_g<-gather(ULS1_Ca, col, val, 1:1170, factor_key=T)
ULS1_Ca_g[is.na(ULS1_Ca_g)]<-0
ULS1_Ca_g<-cbind(ULS1_Ca_g, S_1$cluster)
colnames(ULS1_Ca_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULS1_Ca_g$col))

for (i in 1:length(ULS1_Ca_g$col))
{
  if (ULS1_Ca_g$cluster[i]==1){
    Sv_1[i]<-ULS1_Ca_g$val[i]}
  else if (ULS1_Ca_g$cluster[i]==2){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULS1_Ca_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULS1_Ca_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULS1_Ca_f<-unnest(ULS1_Ca_f) 
heatmap(as.matrix(ULS1_Ca_f), Rowv=NA, Colv=NA, asp=1)


#Fe:
ULS1_Fe_g<-gather(ULS1_Fe, col, val, 1:1170, factor_key=T)
ULS1_Fe_g[is.na(ULS1_Fe_g)]<-0
ULS1_Fe_g<-cbind(ULS1_Fe_g, S_1$cluster)
colnames(ULS1_Fe_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULS1_Fe_g$col))

for (i in 1:length(ULS1_Fe_g$col))
{
  if (ULS1_Fe_g$cluster[i]==1){
    Sv_1[i]<-ULS1_Fe_g$val[i]}
  else if (ULS1_Fe_g$cluster[i]==2){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULS1_Fe_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULS1_Fe_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULS1_Fe_f<-unnest(ULS1_Fe_f) 
heatmap(as.matrix(ULS1_Fe_f), Rowv=NA, Colv=NA, asp=1)


#K:
ULS1_K_g<-gather(ULS1_K, col, val, 1:1170, factor_key=T)
ULS1_K_g[is.na(ULS1_K_g)]<-0
ULS1_K_g<-cbind(ULS1_K_g, S_1$cluster)
colnames(ULS1_K_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULS1_K_g$col))

for (i in 1:length(ULS1_K_g$col))
{
  if (ULS1_K_g$cluster[i]==1){
    Sv_1[i]<-ULS1_K_g$val[i]}
  else if (ULS1_K_g$cluster[i]==2){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULS1_K_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULS1_K_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULS1_K_f<-unnest(ULS1_K_f) 
heatmap(as.matrix(ULS1_K_f), Rowv=NA, Colv=NA, asp=1)


#Ti:
ULS1_Ti_g<-gather(ULS1_Ti, col, val, 1:1170, factor_key=T)
ULS1_Ti_g[is.na(ULS1_Ti_g)]<-0
ULS1_Ti_g<-cbind(ULS1_Ti_g, S_1$cluster)
colnames(ULS1_Ti_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULS1_Ti_g$col))

for (i in 1:length(ULS1_Ti_g$col))
{
  if (ULS1_Ti_g$cluster[i]==1){
    Sv_1[i]<-ULS1_Ti_g$val[i]}
  else if (ULS1_Ti_g$cluster[i]==2){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULS1_Ti_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULS1_Ti_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique Tiey-value combos
ULS1_Ti_f<-unnest(ULS1_Ti_f) 
heatmap(as.matrix(ULS1_Ti_f), Rowv=NA, Colv=NA, asp=1)


#Si:
ULS1_Si_g<-gather(ULS1_Si, col, val, 1:1170, factor_key=T)
ULS1_Si_g[is.na(ULS1_Si_g)]<-0
ULS1_Si_g<-cbind(ULS1_Si_g, S_1$cluster)
colnames(ULS1_Si_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULS1_Si_g$col))

for (i in 1:length(ULS1_Si_g$col))
{
  if (ULS1_Si_g$cluster[i]==1){
    Sv_1[i]<-ULS1_Si_g$val[i]}
  else if (ULS1_Si_g$cluster[i]==2){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULS1_Si_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULS1_Si_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique Siey-value combos
ULS1_Si_f<-unnest(ULS1_Si_f) 
heatmap(as.matrix(ULS1_Si_f), Rowv=NA, Colv=NA, asp=1)


#Mn:
ULS1_Mn_g<-gather(ULS1_Mn, col, val, 1:1170, factor_key=T)
ULS1_Mn_g[is.na(ULS1_Mn_g)]<-0
ULS1_Mn_g<-cbind(ULS1_Mn_g, S_1$cluster)
colnames(ULS1_Mn_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULS1_Mn_g$col))

for (i in 1:length(ULS1_Mn_g$col))
{
  if (ULS1_Mn_g$cluster[i]==1){
    Sv_1[i]<-ULS1_Mn_g$val[i]}
  else if (ULS1_Mn_g$cluster[i]==2){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULS1_Mn_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULS1_Mn_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique Mney-value combos
ULS1_Mn_f<-unnest(ULS1_Mn_f) 
heatmap(as.matrix(ULS1_Mn_f), Rowv=NA, Colv=NA, asp=1)




#Write everything to csv:
write.csv(ULS1_Ca_f, "ULS1_Ca_NA.csv")
write.csv(ULS1_Fe_f, "ULS1_Fe_NA.csv")
write.csv(ULS1_K_f, "ULS1_K_NA.csv")
write.csv(ULS1_Ti_f, "ULS1_Ti_NA.csv")
write.csv(ULS1_Si_f, "ULS1_Si_NA.csv")
write.csv(ULS1_Mn_f, "ULS1_Mn_NA.csv")



#-------------------------------ULS2----------------------------------------------
#Gather dataframe into single column
ULS2_g<-gather(ULS2_sum, col, val, 1:1160, factor_key=T)
ULS2_g[is.na(ULS2_g)]<-0

#carry out kmeans clustering:
ULS2_g_km <- kmeans(ULS2_g[,2],2,iter.max = 100, nstart = 20)

#spread dataset again with just clusters to check whether clustering worked and
#which cluster nr corresponds to the rock surface
S_1<-cbind(ULS2_g, ULS2_g_km$cluster)
colnames(S_1)[3]<-"cluster"
St_1<-S_1[,c(1,3)] #check whether the clustering works. If it does, turn selected value into NA
Sp_1<-pivot_wider(St_1, names_from=col, values_from=cluster) 
Sp_1<-unnest(Sp_1) 
heatmap(as.matrix(Sp_1), Rowv=NA, Colv=NA, asp=1)
view(Sp_1) #cluster 2 = rock surface, cluster 1 = non-rock
#NB pay attention to this when generating your own data, might switch around due to the iterative nature of kmeans!



#create new dataframe with original data for rock surface and NA for non-rock surface
#for each elemental map:
#Ca:
ULS2_Ca_g<-gather(ULS2_Ca, col, val, 1:1160, factor_key=T)
ULS2_Ca_g[is.na(ULS2_Ca_g)]<-0
ULS2_Ca_g<-cbind(ULS2_Ca_g, S_1$cluster)
colnames(ULS2_Ca_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULS2_Ca_g$col))

for (i in 1:length(ULS2_Ca_g$col))
{
  if (ULS2_Ca_g$cluster[i]==2){
    Sv_1[i]<-ULS2_Ca_g$val[i]}
  else if (ULS2_Ca_g$cluster[i]==1){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULS2_Ca_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULS2_Ca_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULS2_Ca_f<-unnest(ULS2_Ca_f) 
heatmap(as.matrix(ULS2_Ca_f), Rowv=NA, Colv=NA, asp=1)


#Fe:
ULS2_Fe_g<-gather(ULS2_Fe, col, val, 1:1160, factor_key=T)
ULS2_Fe_g[is.na(ULS2_Fe_g)]<-0
ULS2_Fe_g<-cbind(ULS2_Fe_g, S_1$cluster)
colnames(ULS2_Fe_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULS2_Fe_g$col))

for (i in 1:length(ULS2_Fe_g$col))
{
  if (ULS2_Fe_g$cluster[i]==2){
    Sv_1[i]<-ULS2_Fe_g$val[i]}
  else if (ULS2_Fe_g$cluster[i]==1){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULS2_Fe_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULS2_Fe_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULS2_Fe_f<-unnest(ULS2_Fe_f) 
heatmap(as.matrix(ULS2_Fe_f), Rowv=NA, Colv=NA, asp=1)


#K:
ULS2_K_g<-gather(ULS2_K, col, val, 1:1160, factor_key=T)
ULS2_K_g[is.na(ULS2_K_g)]<-0
ULS2_K_g<-cbind(ULS2_K_g, S_1$cluster)
colnames(ULS2_K_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULS2_K_g$col))

for (i in 1:length(ULS2_K_g$col))
{
  if (ULS2_K_g$cluster[i]==2){
    Sv_1[i]<-ULS2_K_g$val[i]}
  else if (ULS2_K_g$cluster[i]==1){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULS2_K_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULS2_K_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULS2_K_f<-unnest(ULS2_K_f) 
heatmap(as.matrix(ULS2_K_f), Rowv=NA, Colv=NA, asp=1)


#Ti:
ULS2_Ti_g<-gather(ULS2_Ti, col, val, 1:1160, factor_key=T)
ULS2_Ti_g[is.na(ULS2_Ti_g)]<-0
ULS2_Ti_g<-cbind(ULS2_Ti_g, S_1$cluster)
colnames(ULS2_Ti_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULS2_Ti_g$col))

for (i in 1:length(ULS2_Ti_g$col))
{
  if (ULS2_Ti_g$cluster[i]==2){
    Sv_1[i]<-ULS2_Ti_g$val[i]}
  else if (ULS2_Ti_g$cluster[i]==1){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULS2_Ti_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULS2_Ti_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique Tiey-value combos
ULS2_Ti_f<-unnest(ULS2_Ti_f) 
heatmap(as.matrix(ULS2_Ti_f), Rowv=NA, Colv=NA, asp=1)


#Si:
ULS2_Si_g<-gather(ULS2_Si, col, val, 1:1160, factor_key=T)
ULS2_Si_g[is.na(ULS2_Si_g)]<-0
ULS2_Si_g<-cbind(ULS2_Si_g, S_1$cluster)
colnames(ULS2_Si_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULS2_Si_g$col))

for (i in 1:length(ULS2_Si_g$col))
{
  if (ULS2_Si_g$cluster[i]==2){
    Sv_1[i]<-ULS2_Si_g$val[i]}
  else if (ULS2_Si_g$cluster[i]==1){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULS2_Si_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULS2_Si_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique Siey-value combos
ULS2_Si_f<-unnest(ULS2_Si_f) 
heatmap(as.matrix(ULS2_Si_f), Rowv=NA, Colv=NA, asp=1)


#Mn:
ULS2_Mn_g<-gather(ULS2_Mn, col, val, 1:1160, factor_key=T)
ULS2_Mn_g[is.na(ULS2_Mn_g)]<-0
ULS2_Mn_g<-cbind(ULS2_Mn_g, S_1$cluster)
colnames(ULS2_Mn_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULS2_Mn_g$col))

for (i in 1:length(ULS2_Mn_g$col))
{
  if (ULS2_Mn_g$cluster[i]==2){
    Sv_1[i]<-ULS2_Mn_g$val[i]}
  else if (ULS2_Mn_g$cluster[i]==1){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULS2_Mn_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULS2_Mn_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique Mney-value combos
ULS2_Mn_f<-unnest(ULS2_Mn_f) 
heatmap(as.matrix(ULS2_Mn_f), Rowv=NA, Colv=NA, asp=1)




#Write everything to csv:
write.csv(ULS2_Ca_f, "ULS2_Ca_NA.csv")
write.csv(ULS2_Fe_f, "ULS2_Fe_NA.csv")
write.csv(ULS2_K_f, "ULS2_K_NA.csv")
write.csv(ULS2_Ti_f, "ULS2_Ti_NA.csv")
write.csv(ULS2_Si_f, "ULS2_Si_NA.csv")
write.csv(ULS2_Mn_f, "ULS2_Mn_NA.csv")


#-------------------------Composite maps------------------------------------------
#Read in data; can't automatically use the above due to the iterative nature of kmeans.
#Therefore, use pre-saved data. 
#this data consists of the (when necessary) rotated, flipped, and cropped maps
#(see code above)
ULS1_Ca<-read.csv("ULS1_Ca_NA.csv")
ULS1_Fe<-read.csv("ULS1_Fe_NA.csv")
ULS1_K<-read.csv("ULS1_K_NA.csv")
ULS1_Ti<-read.csv("ULS1_Ti_NA.csv")
ULS1_Si<-read.csv("ULS1_Si_NA.csv")
ULS1_Mn<-read.csv("ULS1_Mn_NA.csv")

ULS2_Ca<-read.csv("ULS2_Ca_NA.csv")
ULS2_Fe<-read.csv("ULS2_Fe_NA.csv")
ULS2_K<-read.csv("ULS2_K_NA.csv")
ULS2_Ti<-read.csv("ULS2_Ti_NA.csv")
ULS2_Si<-read.csv("ULS2_Si_NA.csv")
ULS2_Mn<-read.csv("ULS2_Mn_NA.csv")



#-----------------------------formatting------------------------------------------
#remove row index that R adds when saving csv
#(column index is added too, but this is recognized)
ULS1_Ca<-ULS1_Ca[,2:ncol(ULS1_Ca)]
ULS1_Fe<-ULS1_Fe[,2:ncol(ULS1_Fe)]
ULS1_K<-ULS1_K[,2:ncol(ULS1_K)]
ULS1_Ti<-ULS1_Ti[,2:ncol(ULS1_Ti)]
ULS1_Si<-ULS1_Si[,2:ncol(ULS1_Si)]
ULS1_Mn<-ULS1_Mn[,2:ncol(ULS1_Mn)]

ULS2_Ca<-ULS2_Ca[,2:ncol(ULS2_Ca)]
ULS2_Fe<-ULS2_Fe[,2:ncol(ULS2_Fe)]
ULS2_K<-ULS2_K[,2:ncol(ULS2_K)]
ULS2_Ti<-ULS2_Ti[,2:ncol(ULS2_Ti)]
ULS2_Si<-ULS2_Si[,2:ncol(ULS2_Si)]
ULS2_Mn<-ULS2_Mn[,2:ncol(ULS2_Mn)]

#----------------------------determine overlap------------------------------------
#require the diagonal overlap to obtain complete record
#determine displacement of ULS2 (done with the help of thin-sections and rock photos)
#pad both matrices with zeroes so they have the same length and 
#ULS2 is aligned correctly, then combine both
#Then average this record (with na.rm=T) and add depth scale

#ULS1: overlaps with ULS2 from ca. row 870 (83 mm) to end (row 1170, 113 mm)
#Full length of ULS1 is 1170 rows - ULS2 needs to be padded with 1170-300 = 870 rows to left

#Then, new length of ULS2 is 1160+870=2030
#So ULS1 needs to be padded with 2030-1170 = 860 rows to right

#--------------------------Ca----------------------------------------------------
#Padding & matching colnames
ULS1_Ca_padded<-padzeros(ULS1_Ca, nzeros=860, side="right")
ULS2_Ca_padded<-padzeros(ULS2_Ca, nzeros=870, side="left")

colnames(ULS1_Ca_padded)<-c(1:2030)
colnames(ULS2_Ca_padded)<-c(1:2030)

#Joing matrices by rows
ULS_Ca_padded<-rbind(ULS2_Ca_padded, ULS1_Ca_padded)

#Convert zeroes to NAs and check
ULS_Ca_padded[ULS_Ca_padded==0]<-NA
heatmap(as.matrix(ULS_Ca_padded[,30:2010]), Rowv=NA, Colv=NA, asp=1)


#--------------------------Fe----------------------------------------------------
#Padding & matching colnames
ULS1_Fe_padded<-padzeros(ULS1_Fe, nzeros=860, side="right")
ULS2_Fe_padded<-padzeros(ULS2_Fe, nzeros=870, side="left")

colnames(ULS1_Fe_padded)<-c(1:2030)
colnames(ULS2_Fe_padded)<-c(1:2030)

#Joing matrices by rows
ULS_Fe_padded<-rbind(ULS2_Fe_padded, ULS1_Fe_padded)

#Convert zeroes to NAs and check
ULS_Fe_padded[ULS_Fe_padded==0]<-NA
heatmap(as.matrix(ULS_Fe_padded[,30:2010]), Rowv=NA, Colv=NA, asp=1)


#--------------------------K----------------------------------------------------
#Padding & matching colnames
ULS1_K_padded<-padzeros(ULS1_K, nzeros=860, side="right")
ULS2_K_padded<-padzeros(ULS2_K, nzeros=870, side="left")

colnames(ULS1_K_padded)<-c(1:2030)
colnames(ULS2_K_padded)<-c(1:2030)

#Joing matrices by rows
ULS_K_padded<-rbind(ULS2_K_padded, ULS1_K_padded)

#Convert zeroes to NAs and check
ULS_K_padded[ULS_K_padded==0]<-NA
heatmap(as.matrix(ULS_K_padded[,30:2010]), Rowv=NA, Colv=NA, asp=1)


#--------------------------Ti----------------------------------------------------
#Padding & matching colnames
ULS1_Ti_padded<-padzeros(ULS1_Ti, nzeros=860, side="right")
ULS2_Ti_padded<-padzeros(ULS2_Ti, nzeros=870, side="left")

colnames(ULS1_Ti_padded)<-c(1:2030)
colnames(ULS2_Ti_padded)<-c(1:2030)

#Joing matrices by rows
ULS_Ti_padded<-rbind(ULS2_Ti_padded, ULS1_Ti_padded)

#Convert zeroes to NAs and check
ULS_Ti_padded[ULS_Ti_padded==0]<-NA
heatmap(as.matrix(ULS_Ti_padded[,30:2010]), Rowv=NA, Colv=NA, asp=1)


#--------------------------Si----------------------------------------------------
#Padding & matching colnames
ULS1_Si_padded<-padzeros(ULS1_Si, nzeros=860, side="right")
ULS2_Si_padded<-padzeros(ULS2_Si, nzeros=870, side="left")

colnames(ULS1_Si_padded)<-c(1:2030)
colnames(ULS2_Si_padded)<-c(1:2030)

#Joing matrices by rows
ULS_Si_padded<-rbind(ULS2_Si_padded, ULS1_Si_padded)

#Convert zeroes to NAs and check
ULS_Si_padded[ULS_Si_padded==0]<-NA
heatmap(as.matrix(ULS_Si_padded[,30:2010]), Rowv=NA, Colv=NA, asp=1)


#--------------------------Mn----------------------------------------------------
#Padding & matching colnames
ULS1_Mn_padded<-padzeros(ULS1_Mn, nzeros=860, side="right")
ULS2_Mn_padded<-padzeros(ULS2_Mn, nzeros=870, side="left")

colnames(ULS1_Mn_padded)<-c(1:2030)
colnames(ULS2_Mn_padded)<-c(1:2030)

#Joing matrices by rows
ULS_Mn_padded<-rbind(ULS2_Mn_padded, ULS1_Mn_padded)

#Convert zeroes to NAs and check
ULS_Mn_padded[ULS_Mn_padded==0]<-NA
heatmap(as.matrix(ULS_Mn_padded[,30:2010]), Rowv=NA, Colv=NA, asp=1)



#write all to csv:
write.csv(ULS_Ca_padded, "ULS_comp_map_Ca.csv")
write.csv(ULS_Fe_padded, "ULS_comp_map_Fe.csv")
write.csv(ULS_K_padded, "ULS_comp_map_K.csv")
write.csv(ULS_Ti_padded, "ULS_comp_map_Ti.csv")
write.csv(ULS_Si_padded, "ULS_comp_map_Si.csv")
write.csv(ULS_Mn_padded, "ULS_comp_map_Mn.csv")


#------------------------Composite records-----------------------------------------
#read in data so the whole script doesn't have to be run (time-consuming)
ULS_Ca<-read.csv("ULS_comp_map_Ca.csv")
ULS_Fe<-read.csv("ULS_comp_map_Fe.csv")
ULS_K<-read.csv("ULS_comp_map_K.csv")
ULS_Ti<-read.csv("ULS_comp_map_Ti.csv")
ULS_Si<-read.csv("ULS_comp_map_Si.csv")
ULS_Mn<-read.csv("ULS_comp_map_Mn.csv")

#Once again remove row index
ULS_Ca<-ULS_Ca[,2:ncol(ULS_Ca)]
ULS_Fe<-ULS_Fe[,2:ncol(ULS_Fe)]
ULS_K<-ULS_K[,2:ncol(ULS_K)]
ULS_Ti<-ULS_Ti[,2:ncol(ULS_Ti)]
ULS_Si<-ULS_Si[,2:ncol(ULS_Si)]
ULS_Mn<-ULS_Mn[,2:ncol(ULS_Mn)]


#------------------------Ca-------------------------------------------------------
#Average data:
ULS_Ca_mean<-colMeans(ULS_Ca, na.rm=T)

#add depth column, step size was 100 um = 0.1 mm
ULS_Ca_d<-seq(0, length(ULS_Ca_mean)-1, 1)*0.1 
ULS_Ca_depth<-cbind(ULS_Ca_d, ULS_Ca_mean)

#colnames
colnames(ULS_Ca_depth)<-c("depth_mm", "Ca")

#to dataframe
ULS_Ca_depth<-as.data.frame(ULS_Ca_depth)

#cut off start and end which are empty (see heatmap)
ULS_Ca_depth<-ULS_Ca_depth[30:2010,]

#plot
ggplot()+
  geom_line(data=ULS_Ca_depth, aes(x=depth_mm, y=Ca), col="navy")+
  theme_classic()


#------------------------Fe-------------------------------------------------------
#Average data:
ULS_Fe_mean<-colMeans(ULS_Fe, na.rm=T)

#add depth column, step size was 100 um = 0.1 mm
ULS_Fe_d<-seq(0, length(ULS_Fe_mean)-1, 1)*0.1 
ULS_Fe_depth<-cbind(ULS_Fe_d, ULS_Fe_mean)

#colnames
colnames(ULS_Fe_depth)<-c("depth_mm", "Fe")

#to dataframe
ULS_Fe_depth<-as.data.frame(ULS_Fe_depth)

#cut off start and end which are empty (see heatmap)
ULS_Fe_depth<-ULS_Fe_depth[30:2010,]

#plot
ggplot()+
  geom_line(data=ULS_Fe_depth, aes(x=depth_mm, y=Fe), col="navy")+
  theme_classic()


#------------------------K-------------------------------------------------------
#Average data:
ULS_K_mean<-colMeans(ULS_K, na.rm=T)

#add depth column, step size was 100 um = 0.1 mm
ULS_K_d<-seq(0, length(ULS_K_mean)-1, 1)*0.1 
ULS_K_depth<-cbind(ULS_K_d, ULS_K_mean)

#colnames
colnames(ULS_K_depth)<-c("depth_mm", "K")

#to dataframe
ULS_K_depth<-as.data.frame(ULS_K_depth)

#cut off start and end which are empty (see heatmap)
ULS_K_depth<-ULS_K_depth[30:2010,]

#plot
ggplot()+
  geom_line(data=ULS_K_depth, aes(x=depth_mm, y=K), col="navy")+
  theme_classic()


#------------------------Ti-------------------------------------------------------
#Average data:
ULS_Ti_mean<-colMeans(ULS_Ti, na.rm=T)

#add depth column, step size was 100 um = 0.1 mm
ULS_Ti_d<-seq(0, length(ULS_Ti_mean)-1, 1)*0.1 
ULS_Ti_depth<-cbind(ULS_Ti_d, ULS_Ti_mean)

#colnames
colnames(ULS_Ti_depth)<-c("depth_mm", "Ti")

#to dataframe
ULS_Ti_depth<-as.data.frame(ULS_Ti_depth)

#cut off start and end which are empty (see heatmap)
ULS_Ti_depth<-ULS_Ti_depth[30:2010,]

#plot
ggplot()+
  geom_line(data=ULS_Ti_depth, aes(x=depth_mm, y=Ti), col="navy")+
  theme_classic()


#------------------------Si-------------------------------------------------------
#Average data:
ULS_Si_mean<-colMeans(ULS_Si, na.rm=T)

#add depth column, step size was 100 um = 0.1 mm
ULS_Si_d<-seq(0, length(ULS_Si_mean)-1, 1)*0.1 
ULS_Si_depth<-cbind(ULS_Si_d, ULS_Si_mean)

#colnames
colnames(ULS_Si_depth)<-c("depth_mm", "Si")

#to dataframe
ULS_Si_depth<-as.data.frame(ULS_Si_depth)

#cut off start and end which are empty (see heatmap)
ULS_Si_depth<-ULS_Si_depth[30:2010,]

#plot
ggplot()+
  geom_line(data=ULS_Si_depth, aes(x=depth_mm, y=Si), col="navy")+
  theme_classic()


#------------------------Mn-------------------------------------------------------
#Average data:
ULS_Mn_mean<-colMeans(ULS_Mn, na.rm=T)

#add depth column, step size was 100 um = 0.1 mm
ULS_Mn_d<-seq(0, length(ULS_Mn_mean)-1, 1)*0.1 
ULS_Mn_depth<-cbind(ULS_Mn_d, ULS_Mn_mean)

#colnames
colnames(ULS_Mn_depth)<-c("depth_mm", "Mn")

#to dataframe
ULS_Mn_depth<-as.data.frame(ULS_Mn_depth)

#cut off start and end which are empty (see heatmap)
ULS_Mn_depth<-ULS_Mn_depth[30:2010,]

#plot
ggplot()+
  geom_line(data=ULS_Mn_depth, aes(x=depth_mm, y=Mn), col="navy")+
  theme_classic()



#write everything to csv
write.csv(ULS_Ca_depth, "ULS_Ca_depth.csv")
write.csv(ULS_Fe_depth, "ULS_Fe_depth.csv")
write.csv(ULS_K_depth, "ULS_K_depth.csv")
write.csv(ULS_Ti_depth, "ULS_Ti_depth.csv")
write.csv(ULS_Si_depth, "ULS_Si_depth.csv")
write.csv(ULS_Mn_depth, "ULS_Mn_depth.csv")



#------------------------Quick plotting of records-------------------------------
#read in data:
ULS_Ca_depth<-read.csv("ULS_Ca_depth.csv")
ULS_Fe_depth<-read.csv("ULS_Fe_depth.csv")
ULS_K_depth<-read.csv("ULS_K_depth.csv")
ULS_Ti_depth<-read.csv("ULS_Ti_depth.csv")
ULS_Mn_depth<-read.csv("ULS_Mn_depth.csv")
ULS_Si_depth<-read.csv("ULS_Si_depth.csv")

#formatting: remove first column
ULS_Ca_depth<-ULS_Ca_depth[,2:3]
ULS_Fe_depth<-ULS_Fe_depth[,2:3]
ULS_K_depth<-ULS_K_depth[,2:3]
ULS_Ti_depth<-ULS_Ti_depth[,2:3]
ULS_Mn_depth<-ULS_Mn_depth[,2:3]
ULS_Si_depth<-ULS_Si_depth[,2:3]

#Plot individual records:
ggplot()+
  geom_line(data=ULS_Ca_depth, aes(x=depth_mm, y=Ca), col="navy")+
  theme_classic()

ggplot()+
  geom_line(data=ULS_Fe_depth, aes(x=depth_mm, y=Fe), col="navy")+
  theme_classic()

ggplot()+
  geom_line(data=ULS_K_depth, aes(x=depth_mm, y=K), col="navy")+
  theme_classic()

ggplot()+
  geom_line(data=ULS_Ti_depth, aes(x=depth_mm, y=Ti), col="navy")+
  theme_classic()

ggplot()+
  geom_line(data=ULS_Si_depth, aes(x=depth_mm, y=Si), col="navy")+
  theme_classic()

ggplot()+
  geom_line(data=ULS_Mn_depth, aes(x=depth_mm, y=Mn), col="navy")+
  theme_classic()



#Plot together:
Ca_p<-ggplot()+
  geom_line(data=ULS_Ca_depth, aes(x=depth_mm, y=Ca), col="navy")+
  coord_flip()+
  theme_classic()

Fe_p<-ggplot()+
  geom_line(data=ULS_Fe_depth, aes(x=depth_mm, y=Fe), col="navy")+
  coord_flip()+
  theme_classic()

K_p<-ggplot()+
  geom_line(data=ULS_K_depth, aes(x=depth_mm, y=K), col="navy")+
  coord_flip()+
  theme_classic()

Ti_p<-ggplot()+
  geom_line(data=ULS_Ti_depth, aes(x=depth_mm, y=Ti), col="navy")+
  coord_flip()+
  theme_classic()

Si_p<-ggplot()+
  geom_line(data=ULS_Si_depth, aes(x=depth_mm, y=Si), col="navy")+
  coord_flip()+
  theme_classic()

Mn_p<-ggplot()+
  geom_line(data=ULS_Mn_depth, aes(x=depth_mm, y=Mn), col="navy")+
  coord_flip()+
  theme_classic()


cowplot::plot_grid(Ca_p, Fe_p, K_p, Ti_p, Si_p, Mn_p, ncol=6)

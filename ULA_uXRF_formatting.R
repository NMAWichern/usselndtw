#code for generating elemental depth records from micro-XRF generated maps
#Dynamic Time Warping (DTW) has been applied on these depth records
#This code is for the Arfeld Usseln Limestone (ULA)
#ULA consists of five sections, ULA 1 to 5

#steps within this code:
#1) removing areas without data using k-means clustering
#2) plot all individual segments together as a composite map
#3) determine overlap between segments and combine them (ULA 1 to 5)
#4) average all rows in order to generate a composite depth record


library(tidyverse)
library(astrochron)
library(devtools)
library(ptw)
library(RColorBrewer)

setwd("insert wd")

#--------------Initial data checking and formatting-------------------------------

#Load in data:
ULA1_Ca<-read.csv("ULA1_Ca.csv", header=T)
ULA1_Fe<-read.csv("ULA1_Fe.csv", header=T)
ULA1_Ti<-read.csv("ULA1_Ti.csv", header=T)
ULA1_K<-read.csv("ULA1_K.csv", header=T)
ULA1_Mn<-read.csv("ULA1_Mn.csv", header=T)
ULA1_Si<-read.csv("ULA1_Si.csv", header=T)

ULA2_Ca<-read.csv("ULA2_Ca.csv", header=T)
ULA2_Fe<-read.csv("ULA2_Fe.csv", header=T)
ULA2_Ti<-read.csv("ULA2_Ti.csv", header=T)
ULA2_K<-read.csv("ULA2_K.csv", header=T)
ULA2_Mn<-read.csv("ULA2_Mn.csv", header=T)
ULA2_Si<-read.csv("ULA2_Si.csv", header=T)

ULA3_Ca<-read.csv("ULA3_Ca.csv", header=T)
ULA3_Fe<-read.csv("ULA3_Fe.csv", header=T)
ULA3_Ti<-read.csv("ULA3_Ti.csv", header=T)
ULA3_K<-read.csv("ULA3_K.csv", header=T)
ULA3_Mn<-read.csv("ULA3_Mn.csv", header=T)
ULA3_Si<-read.csv("ULA3_Si.csv", header=T)

ULA4_Ca<-read.csv("ULA4_Ca.csv", header=T)
ULA4_Fe<-read.csv("ULA4_Fe.csv", header=T)
ULA4_Ti<-read.csv("ULA4_Ti.csv", header=T)
ULA4_K<-read.csv("ULA4_K.csv", header=T)
ULA4_Mn<-read.csv("ULA4_Mn.csv", header=T)
ULA4_Si<-read.csv("ULA4_Si.csv", header=T)

ULA5_Ca<-read.csv("ULA5_Ca_rot.csv", header=T) #rotated version so bedding is mostly vertical
ULA5_Fe<-read.csv("ULA5_Fe_rot.csv", header=T)
ULA5_Ti<-read.csv("ULA5_Ti_rot.csv", header=T)
ULA5_K<-read.csv("ULA5_K_rot.csv", header=T)
ULA5_Mn<-read.csv("ULA5_Mn_rot.csv", header=T)
ULA5_Si<-read.csv("ULA5_Si_rot.csv", header=T)


#flip horizontally: 1, 5
ULA1_Ca<-ULA1_Ca[350:1,1001:1]
ULA1_Fe<-ULA1_Fe[350:1,1001:1]
ULA1_K<-ULA1_K[350:1,1001:1]
ULA1_Ti<-ULA1_Ti[350:1,1001:1]
ULA1_Si<-ULA1_Si[350:1,1001:1]
ULA1_Mn<-ULA1_Mn[350:1,1001:1]

colnames(ULA1_Ca)<-c(1:1001) #because it's confusing when they're in descending order
rownames(ULA1_Ca)<-c(1:350)
colnames(ULA1_Fe)<-c(1:1001) 
rownames(ULA1_Fe)<-c(1:350)
colnames(ULA1_K)<-c(1:1001) 
rownames(ULA1_K)<-c(1:350)
colnames(ULA1_Ti)<-c(1:1001) 
rownames(ULA1_Ti)<-c(1:350)
colnames(ULA1_Si)<-c(1:1001) 
rownames(ULA1_Si)<-c(1:350)
colnames(ULA1_Mn)<-c(1:1001) 
rownames(ULA1_Mn)<-c(1:350)

ULA5_Ca<-ULA5_Ca[674:1,968:1]
ULA5_Fe<-ULA5_Fe[674:1,968:1]
ULA5_K<-ULA5_K[674:1,968:1]
ULA5_Ti<-ULA5_Ti[674:1,968:1]
ULA5_Si<-ULA5_Si[674:1,968:1]
ULA5_Mn<-ULA5_Mn[674:1,968:1]

colnames(ULA5_Ca)<-c(1:968)
rownames(ULA5_Ca)<-c(1:674)
colnames(ULA5_Fe)<-c(1:968)
rownames(ULA5_Fe)<-c(1:674)
colnames(ULA5_K)<-c(1:968)
rownames(ULA5_K)<-c(1:674)
colnames(ULA5_Ti)<-c(1:968)
rownames(ULA5_Ti)<-c(1:674)
colnames(ULA5_Si)<-c(1:968)
rownames(ULA5_Si)<-c(1:674)
colnames(ULA5_Mn)<-c(1:968)
rownames(ULA5_Mn)<-c(1:674)


#REMOVE EXCESS AREAS
#ULA1: flip (done), make square
ULA1_Ca<-ULA1_Ca[107:287,37:982]
ULA1_Fe<-ULA1_Fe[107:287,37:982]
ULA1_K<-ULA1_K[107:287,37:982]
ULA1_Ti<-ULA1_Ti[107:287,37:982]
ULA1_Si<-ULA1_Si[107:287,37:982]
ULA1_Mn<-ULA1_Mn[107:287,37:982]

#ULA2: shorten a bit
ULA2_Ca<-ULA2_Ca[,11:1102]
ULA2_Fe<-ULA2_Fe[,11:1102]
ULA2_K<-ULA2_K[,11:1102]
ULA2_Ti<-ULA2_Ti[,11:1102]
ULA2_Si<-ULA2_Si[,11:1102]
ULA2_Mn<-ULA2_Mn[,11:1102]


#ULA3: shorten a bit
ULA3_Ca<-ULA3_Ca[,28:1171]
ULA3_Fe<-ULA3_Fe[,28:1171]
ULA3_K<-ULA3_K[,28:1171]
ULA3_Ti<-ULA3_Ti[,28:1171]
ULA3_Si<-ULA3_Si[,28:1171]
ULA3_Mn<-ULA3_Mn[,28:1171]

#ULA4: shorten a bit
ULA4_Ca<-ULA4_Ca[1:223,4:918]
ULA4_Fe<-ULA4_Fe[1:223,4:918]
ULA4_K<-ULA4_K[1:223,4:918]
ULA4_Ti<-ULA4_Ti[1:223,4:918]
ULA4_Si<-ULA4_Si[1:223,4:918]
ULA4_Mn<-ULA4_Mn[1:223,4:918]

#ULA5: make square, overlap with ULA4 is good. 
ULA5_Ca<-ULA5_Ca[80:621,100:900]
ULA5_Fe<-ULA5_Fe[80:621,100:900]
ULA5_K<-ULA5_K[80:621,100:900]
ULA5_Ti<-ULA5_Ti[80:621,100:900]
ULA5_Si<-ULA5_Si[80:621,100:900]
ULA5_Mn<-ULA5_Mn[80:621,100:900]


#plot as heatmaps:
heatmap(as.matrix(ULA1_Ca), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA1_Fe), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA1_K), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA1_Ti), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA1_Si), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA1_Mn), Rowv=NA, Colv=NA, asp=1)

heatmap(as.matrix(ULA2_Ca), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA2_Fe), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA2_K), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA2_Ti), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA2_Si), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA2_Mn), Rowv=NA, Colv=NA, asp=1)

heatmap(as.matrix(ULA3_Ca), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA3_Fe), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA3_K), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA3_Ti), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA3_Si), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA3_Mn), Rowv=NA, Colv=NA, asp=1)

heatmap(as.matrix(ULA4_Ca), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA4_Fe), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA4_K), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA4_Ti), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA4_Si), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA4_Mn), Rowv=NA, Colv=NA, asp=1)

heatmap(as.matrix(ULA5_Ca), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA5_Fe), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA5_K), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA5_Ti), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA5_Si), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA5_Mn), Rowv=NA, Colv=NA, asp=1)




#-------------------------K-means clustering--------------------------------------- 
#sum and sqrt all different elemental maps to increase the contrast between rock and non-rock area
ULA1_sum<-sqrt(ULA1_Ca+ULA1_Fe+ULA1_K+ULA1_Ti+ULA1_Si+ULA1_Mn)
ULA2_sum<-sqrt(ULA2_Ca+ULA2_Fe+ULA2_K+ULA2_Ti+ULA2_Si+ULA2_Mn)
ULA3_sum<-sqrt(ULA3_Ca+ULA3_Fe+ULA3_K+ULA3_Ti+ULA3_Si+ULA3_Mn)
ULA4_sum<-sqrt(ULA4_Ca+ULA4_Fe+ULA4_K+ULA4_Ti+ULA4_Si+ULA4_Mn)
ULA5_sum<-sqrt(ULA5_Ca+ULA5_Fe+ULA5_K+ULA5_Ti+ULA5_Si+ULA5_Mn)

heatmap(as.matrix(ULA1_sum), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA2_sum), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA3_sum), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA4_sum), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULA5_sum), Rowv=NA, Colv=NA, asp=1)
#visually the contrast seems good in these





#-------------------------------ULA1----------------------------------------------
#Gather dataframe into single column
ULA1_g<-gather(ULA1_sum, col, val, 1:946, factor_key=T)
ULA1_g[is.na(ULA1_g)]<-0

#carry out kmeans clustering:
ULA1_g_km <- kmeans(ULA1_g[,2],2,iter.max = 100, nstart = 20)

#spread dataset again with just clusters to check whether clustering worked and
#which cluster nr corresponds to the rock surface
S_1<-cbind(ULA1_g, ULA1_g_km$cluster)
colnames(S_1)[3]<-"cluster"
St_1<-S_1[,c(1,3)] #check whether the clustering works. If it does, turn selected value into NA
Sp_1<-pivot_wider(St_1, names_from=col, values_from=cluster)
Sp_1<-unnest(Sp_1) 
heatmap(as.matrix(Sp_1), Rowv=NA, Colv=NA, asp=1)
view(Sp_1) #cluster 1 = rock surface, cluster 2 = non-rock
#NB pay attention to this when generating your own data, might switch around due to the iterative nature of kmeans!



#create new dataframe with original data for rock surface and NA for non-rock surface
#for each elemental map
#Ca:
ULA1_Ca_g<-gather(ULA1_Ca, col, val, 1:946, factor_key=T)
ULA1_Ca_g[is.na(ULA1_Ca_g)]<-0
ULA1_Ca_g<-cbind(ULA1_Ca_g, S_1$cluster)
colnames(ULA1_Ca_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULA1_Ca_g$col))

for (i in 1:length(ULA1_Ca_g$col))
{
  if (ULA1_Ca_g$cluster[i]==1){
    Sv_1[i]<-ULA1_Ca_g$val[i]}
  else if (ULA1_Ca_g$cluster[i]==2){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULA1_Ca_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULA1_Ca_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA1_Ca_f<-unnest(ULA1_Ca_f) 
heatmap(as.matrix(ULA1_Ca_f), Rowv=NA, Colv=NA, asp=1)


#Fe:
ULA1_Fe_g<-gather(ULA1_Fe, col, val, 1:946, factor_key=T)
ULA1_Fe_g[is.na(ULA1_Fe_g)]<-0
ULA1_Fe_g<-cbind(ULA1_Fe_g, S_1$cluster)
colnames(ULA1_Fe_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULA1_Fe_g$col))

for (i in 1:length(ULA1_Fe_g$col))
{
  if (ULA1_Fe_g$cluster[i]==1){
    Sv_1[i]<-ULA1_Fe_g$val[i]}
  else if (ULA1_Fe_g$cluster[i]==2){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULA1_Fe_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULA1_Fe_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA1_Fe_f<-unnest(ULA1_Fe_f) 
heatmap(as.matrix(ULA1_Fe_f), Rowv=NA, Colv=NA, asp=1)



#K:
ULA1_K_g<-gather(ULA1_K, col, val, 1:946, factor_key=T)
ULA1_K_g[is.na(ULA1_K_g)]<-0
ULA1_K_g<-cbind(ULA1_K_g, S_1$cluster)
colnames(ULA1_K_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULA1_K_g$col))

for (i in 1:length(ULA1_K_g$col))
{
  if (ULA1_K_g$cluster[i]==1){
    Sv_1[i]<-ULA1_K_g$val[i]}
  else if (ULA1_K_g$cluster[i]==2){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULA1_K_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULA1_K_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA1_K_f<-unnest(ULA1_K_f) 
heatmap(as.matrix(ULA1_K_f), Rowv=NA, Colv=NA, asp=1)



#Ti:
ULA1_Ti_g<-gather(ULA1_Ti, col, val, 1:946, factor_key=T)
ULA1_Ti_g[is.na(ULA1_Ti_g)]<-0
ULA1_Ti_g<-cbind(ULA1_Ti_g, S_1$cluster)
colnames(ULA1_Ti_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULA1_Ti_g$col))

for (i in 1:length(ULA1_Ti_g$col))
{
  if (ULA1_Ti_g$cluster[i]==1){
    Sv_1[i]<-ULA1_Ti_g$val[i]}
  else if (ULA1_Ti_g$cluster[i]==2){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULA1_Ti_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULA1_Ti_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA1_Ti_f<-unnest(ULA1_Ti_f) 
heatmap(as.matrix(ULA1_Ti_f), Rowv=NA, Colv=NA, asp=1)



#Si:
ULA1_Si_g<-gather(ULA1_Si, col, val, 1:946, factor_key=T)
ULA1_Si_g[is.na(ULA1_Si_g)]<-0
ULA1_Si_g<-cbind(ULA1_Si_g, S_1$cluster)
colnames(ULA1_Si_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULA1_Si_g$col))

for (i in 1:length(ULA1_Si_g$col))
{
  if (ULA1_Si_g$cluster[i]==1){
    Sv_1[i]<-ULA1_Si_g$val[i]}
  else if (ULA1_Si_g$cluster[i]==2){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULA1_Si_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULA1_Si_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA1_Si_f<-unnest(ULA1_Si_f) 
heatmap(as.matrix(ULA1_Si_f), Rowv=NA, Colv=NA, asp=1)



#Mn:
ULA1_Mn_g<-gather(ULA1_Mn, col, val, 1:946, factor_key=T)
ULA1_Mn_g[is.na(ULA1_Mn_g)]<-0
ULA1_Mn_g<-cbind(ULA1_Mn_g, S_1$cluster)
colnames(ULA1_Mn_g)[3]<-"cluster"
Sv_1<-array(data=NA, dim=length(ULA1_Mn_g$col))

for (i in 1:length(ULA1_Mn_g$col))
{
  if (ULA1_Mn_g$cluster[i]==1){
    Sv_1[i]<-ULA1_Mn_g$val[i]}
  else if (ULA1_Mn_g$cluster[i]==2){
    Sv_1[i]<-NA}
}

Sv_1<-cbind(ULA1_Mn_g$col, Sv_1)
colnames(Sv_1)<-c("col", "val")

ULA1_Mn_f<-pivot_wider(as.data.frame(Sv_1), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA1_Mn_f<-unnest(ULA1_Mn_f) 
heatmap(as.matrix(ULA1_Mn_f), Rowv=NA, Colv=NA, asp=1)

write.csv(ULA1_Ca_f, "ULA1_Ca_NA.csv")
write.csv(ULA1_Fe_f, "ULA1_Fe_NA.csv")
write.csv(ULA1_K_f, "ULA1_K_NA.csv")
write.csv(ULA1_Ti_f, "ULA1_Ti_NA.csv")
write.csv(ULA1_Si_f, "ULA1_Si_NA.csv")
write.csv(ULA1_Mn_f, "ULA1_Mn_NA.csv")








#-------------------------------ULA2----------------------------------------------
#Gather dataframe into single column
ULA2_g<-gather(ULA2_sum, col, val, 1:1092, factor_key=T)
ULA2_g[is.na(ULA2_g)]<-0

#carry out kmeans clustering:
ULA2_g_km <- kmeans(ULA2_g[,2],2,iter.max = 100, nstart = 20)

#spread dataset again with just clusters to check whether clustering worked and
#which cluster nr corresponds to the rock surface
S_2<-cbind(ULA2_g, ULA2_g_km$cluster)
colnames(S_2)[3]<-"cluster"
St_2<-S_2[,c(1,3)] #check whether the clustering works. If it does, turn selected value into NA
Sp_2<-pivot_wider(St_2, names_from=col, values_from=cluster)
Sp_2<-unnest(Sp_2) 
heatmap(as.matrix(Sp_2), Rowv=NA, Colv=NA, asp=1)
view(Sp_2) #cluster 1 = rock surface, cluster 2 = non-rock
#NB pay attention to this when generating your own data, might switch around due to the iterative nature of kmeans!



#create new dataframe with original data for rock surface and NA for non-rock surface
#for each elemental map
#Ca:
ULA2_Ca_g<-gather(ULA2_Ca, col, val, 1:1092, factor_key=T)
ULA2_Ca_g[is.na(ULA2_Ca_g)]<-0
ULA2_Ca_g<-cbind(ULA2_Ca_g, S_2$cluster)
colnames(ULA2_Ca_g)[3]<-"cluster"
Sv_2<-array(data=NA, dim=length(ULA2_Ca_g$col))

for (i in 2:length(ULA2_Ca_g$col))
{
  if (ULA2_Ca_g$cluster[i]==2){
    Sv_2[i]<-ULA2_Ca_g$val[i]}
  else if (ULA2_Ca_g$cluster[i]==1){
    Sv_2[i]<-NA}
}

Sv_2<-cbind(ULA2_Ca_g$col, Sv_2)
colnames(Sv_2)<-c("col", "val")

ULA2_Ca_f<-pivot_wider(as.data.frame(Sv_2), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA2_Ca_f<-unnest(ULA2_Ca_f) 
heatmap(as.matrix(ULA2_Ca_f), Rowv=NA, Colv=NA, asp=1)


#Fe:
ULA2_Fe_g<-gather(ULA2_Fe, col, val, 1:1092, factor_key=T)
ULA2_Fe_g[is.na(ULA2_Fe_g)]<-0
ULA2_Fe_g<-cbind(ULA2_Fe_g, S_2$cluster)
colnames(ULA2_Fe_g)[3]<-"cluster"
Sv_2<-array(data=NA, dim=length(ULA2_Fe_g$col))

for (i in 1:length(ULA2_Fe_g$col))
{
  if (ULA2_Fe_g$cluster[i]==2){
    Sv_2[i]<-ULA2_Fe_g$val[i]}
  else if (ULA2_Fe_g$cluster[i]==1){
    Sv_2[i]<-NA}
}

Sv_2<-cbind(ULA2_Fe_g$col, Sv_2)
colnames(Sv_2)<-c("col", "val")

ULA2_Fe_f<-pivot_wider(as.data.frame(Sv_2), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA2_Fe_f<-unnest(ULA2_Fe_f) 
heatmap(as.matrix(ULA2_Fe_f), Rowv=NA, Colv=NA, asp=1)



#K:
ULA2_K_g<-gather(ULA2_K, col, val, 1:1092, factor_key=T)
ULA2_K_g[is.na(ULA2_K_g)]<-0
ULA2_K_g<-cbind(ULA2_K_g, S_2$cluster)
colnames(ULA2_K_g)[3]<-"cluster"
Sv_2<-array(data=NA, dim=length(ULA2_K_g$col))

for (i in 1:length(ULA2_K_g$col))
{
  if (ULA2_K_g$cluster[i]==2){
    Sv_2[i]<-ULA2_K_g$val[i]}
  else if (ULA2_K_g$cluster[i]==1){
    Sv_2[i]<-NA}
}

Sv_2<-cbind(ULA2_K_g$col, Sv_2)
colnames(Sv_2)<-c("col", "val")

ULA2_K_f<-pivot_wider(as.data.frame(Sv_2), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA2_K_f<-unnest(ULA2_K_f) 
heatmap(as.matrix(ULA2_K_f), Rowv=NA, Colv=NA, asp=1)



#Ti:
ULA2_Ti_g<-gather(ULA2_Ti, col, val, 1:1092, factor_key=T)
ULA2_Ti_g[is.na(ULA2_Ti_g)]<-0
ULA2_Ti_g<-cbind(ULA2_Ti_g, S_2$cluster)
colnames(ULA2_Ti_g)[3]<-"cluster"
Sv_2<-array(data=NA, dim=length(ULA2_Ti_g$col))

for (i in 1:length(ULA2_Ti_g$col))
{
  if (ULA2_Ti_g$cluster[i]==2){
    Sv_2[i]<-ULA2_Ti_g$val[i]}
  else if (ULA2_Ti_g$cluster[i]==1){
    Sv_2[i]<-NA}
}

Sv_2<-cbind(ULA2_Ti_g$col, Sv_2)
colnames(Sv_2)<-c("col", "val")

ULA2_Ti_f<-pivot_wider(as.data.frame(Sv_2), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA2_Ti_f<-unnest(ULA2_Ti_f) 
heatmap(as.matrix(ULA2_Ti_f), Rowv=NA, Colv=NA, asp=1)



#Si:
ULA2_Si_g<-gather(ULA2_Si, col, val, 1:1092, factor_key=T)
ULA2_Si_g[is.na(ULA2_Si_g)]<-0
ULA2_Si_g<-cbind(ULA2_Si_g, S_2$cluster)
colnames(ULA2_Si_g)[3]<-"cluster"
Sv_2<-array(data=NA, dim=length(ULA2_Si_g$col))

for (i in 1:length(ULA2_Si_g$col))
{
  if (ULA2_Si_g$cluster[i]==2){
    Sv_2[i]<-ULA2_Si_g$val[i]}
  else if (ULA2_Si_g$cluster[i]==1){
    Sv_2[i]<-NA}
}

Sv_2<-cbind(ULA2_Si_g$col, Sv_2)
colnames(Sv_2)<-c("col", "val")

ULA2_Si_f<-pivot_wider(as.data.frame(Sv_2), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA2_Si_f<-unnest(ULA2_Si_f) 
heatmap(as.matrix(ULA2_Si_f), Rowv=NA, Colv=NA, asp=1)



#Mn:
ULA2_Mn_g<-gather(ULA2_Mn, col, val, 1:1092, factor_key=T)
ULA2_Mn_g[is.na(ULA2_Mn_g)]<-0
ULA2_Mn_g<-cbind(ULA2_Mn_g, S_2$cluster)
colnames(ULA2_Mn_g)[3]<-"cluster"
Sv_2<-array(data=NA, dim=length(ULA2_Mn_g$col))

for (i in 1:length(ULA2_Mn_g$col))
{
  if (ULA2_Mn_g$cluster[i]==2){
    Sv_2[i]<-ULA2_Mn_g$val[i]}
  else if (ULA2_Mn_g$cluster[i]==1){
    Sv_2[i]<-NA}
}

Sv_2<-cbind(ULA2_Mn_g$col, Sv_2)
colnames(Sv_2)<-c("col", "val")

ULA2_Mn_f<-pivot_wider(as.data.frame(Sv_2), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA2_Mn_f<-unnest(ULA2_Mn_f) 
heatmap(as.matrix(ULA2_Mn_f), Rowv=NA, Colv=NA, asp=1)

write.csv(ULA2_Ca_f, "ULA2_Ca_NA.csv")
write.csv(ULA2_Fe_f, "ULA2_Fe_NA.csv")
write.csv(ULA2_K_f, "ULA2_K_NA.csv")
write.csv(ULA2_Ti_f, "ULA2_Ti_NA.csv")
write.csv(ULA2_Si_f, "ULA2_Si_NA.csv")
write.csv(ULA2_Mn_f, "ULA2_Mn_NA.csv")





#-------------------------------ULA3----------------------------------------------
#Gather dataframe into single column
ULA3_g<-gather(ULA3_sum, col, val, 1:1144, factor_key=T)
ULA3_g[is.na(ULA3_g)]<-0

#carry out kmeans clustering:
ULA3_g_km <- kmeans(ULA3_g[,2],2,iter.max = 100, nstart = 20)

#spread dataset again with just clusters to check whether clustering worked and
#which cluster nr corresponds to the rock surface
S_3<-cbind(ULA3_g, ULA3_g_km$cluster)
colnames(S_3)[3]<-"cluster"
St_3<-S_3[,c(1,3)] #check whether the clustering works. If it does, turn selected value into NA
Sp_3<-pivot_wider(St_3, names_from=col, values_from=cluster) 
Sp_3<-unnest(Sp_3) 
heatmap(as.matrix(Sp_3), Rowv=NA, Colv=NA, asp=1)
view(Sp_3) #cluster 2 = rock surface, cluster 1 = non-rock
#NB pay attention to this when generating your own data, might switch around due to the iterative nature of kmeans!



#create new dataframe with original data for rock surface and NA for non-rock surface
#for each elemental map
#Ca:
ULA3_Ca_g<-gather(ULA3_Ca, col, val, 1:1144, factor_key=T)
ULA3_Ca_g[is.na(ULA3_Ca_g)]<-0
ULA3_Ca_g<-cbind(ULA3_Ca_g, S_3$cluster)
colnames(ULA3_Ca_g)[3]<-"cluster"
Sv_3<-array(data=NA, dim=length(ULA3_Ca_g$col))

for (i in 3:length(ULA3_Ca_g$col))
{
  if (ULA3_Ca_g$cluster[i]==2){
    Sv_3[i]<-ULA3_Ca_g$val[i]}
  else if (ULA3_Ca_g$cluster[i]==1){
    Sv_3[i]<-NA}
}

Sv_3<-cbind(ULA3_Ca_g$col, Sv_3)
colnames(Sv_3)<-c("col", "val")

ULA3_Ca_f<-pivot_wider(as.data.frame(Sv_3), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA3_Ca_f<-unnest(ULA3_Ca_f) 
heatmap(as.matrix(ULA3_Ca_f), Rowv=NA, Colv=NA, asp=1)


#Fe:
ULA3_Fe_g<-gather(ULA3_Fe, col, val, 1:1144, factor_key=T)
ULA3_Fe_g[is.na(ULA3_Fe_g)]<-0
ULA3_Fe_g<-cbind(ULA3_Fe_g, S_3$cluster)
colnames(ULA3_Fe_g)[3]<-"cluster"
Sv_3<-array(data=NA, dim=length(ULA3_Fe_g$col))

for (i in 1:length(ULA3_Fe_g$col))
{
  if (ULA3_Fe_g$cluster[i]==2){
    Sv_3[i]<-ULA3_Fe_g$val[i]}
  else if (ULA3_Fe_g$cluster[i]==1){
    Sv_3[i]<-NA}
}

Sv_3<-cbind(ULA3_Fe_g$col, Sv_3)
colnames(Sv_3)<-c("col", "val")

ULA3_Fe_f<-pivot_wider(as.data.frame(Sv_3), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA3_Fe_f<-unnest(ULA3_Fe_f) 
heatmap(as.matrix(ULA3_Fe_f), Rowv=NA, Colv=NA, asp=1)



#K:
ULA3_K_g<-gather(ULA3_K, col, val, 1:1144, factor_key=T)
ULA3_K_g[is.na(ULA3_K_g)]<-0
ULA3_K_g<-cbind(ULA3_K_g, S_3$cluster)
colnames(ULA3_K_g)[3]<-"cluster"
Sv_3<-array(data=NA, dim=length(ULA3_K_g$col))

for (i in 1:length(ULA3_K_g$col))
{
  if (ULA3_K_g$cluster[i]==2){
    Sv_3[i]<-ULA3_K_g$val[i]}
  else if (ULA3_K_g$cluster[i]==1){
    Sv_3[i]<-NA}
}

Sv_3<-cbind(ULA3_K_g$col, Sv_3)
colnames(Sv_3)<-c("col", "val")

ULA3_K_f<-pivot_wider(as.data.frame(Sv_3), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA3_K_f<-unnest(ULA3_K_f) 
heatmap(as.matrix(ULA3_K_f), Rowv=NA, Colv=NA, asp=1)



#Ti:
ULA3_Ti_g<-gather(ULA3_Ti, col, val, 1:1144, factor_key=T)
ULA3_Ti_g[is.na(ULA3_Ti_g)]<-0
ULA3_Ti_g<-cbind(ULA3_Ti_g, S_3$cluster)
colnames(ULA3_Ti_g)[3]<-"cluster"
Sv_3<-array(data=NA, dim=length(ULA3_Ti_g$col))

for (i in 1:length(ULA3_Ti_g$col))
{
  if (ULA3_Ti_g$cluster[i]==2){
    Sv_3[i]<-ULA3_Ti_g$val[i]}
  else if (ULA3_Ti_g$cluster[i]==1){
    Sv_3[i]<-NA}
}

Sv_3<-cbind(ULA3_Ti_g$col, Sv_3)
colnames(Sv_3)<-c("col", "val")

ULA3_Ti_f<-pivot_wider(as.data.frame(Sv_3), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA3_Ti_f<-unnest(ULA3_Ti_f) 
heatmap(as.matrix(ULA3_Ti_f), Rowv=NA, Colv=NA, asp=1)



#Si:
ULA3_Si_g<-gather(ULA3_Si, col, val, 1:1144, factor_key=T)
ULA3_Si_g[is.na(ULA3_Si_g)]<-0
ULA3_Si_g<-cbind(ULA3_Si_g, S_3$cluster)
colnames(ULA3_Si_g)[3]<-"cluster"
Sv_3<-array(data=NA, dim=length(ULA3_Si_g$col))

for (i in 1:length(ULA3_Si_g$col))
{
  if (ULA3_Si_g$cluster[i]==2){
    Sv_3[i]<-ULA3_Si_g$val[i]}
  else if (ULA3_Si_g$cluster[i]==1){
    Sv_3[i]<-NA}
}

Sv_3<-cbind(ULA3_Si_g$col, Sv_3)
colnames(Sv_3)<-c("col", "val")

ULA3_Si_f<-pivot_wider(as.data.frame(Sv_3), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA3_Si_f<-unnest(ULA3_Si_f) 
heatmap(as.matrix(ULA3_Si_f), Rowv=NA, Colv=NA, asp=1)



#Mn:
ULA3_Mn_g<-gather(ULA3_Mn, col, val, 1:1144, factor_key=T)
ULA3_Mn_g[is.na(ULA3_Mn_g)]<-0
ULA3_Mn_g<-cbind(ULA3_Mn_g, S_3$cluster)
colnames(ULA3_Mn_g)[3]<-"cluster"
Sv_3<-array(data=NA, dim=length(ULA3_Mn_g$col))

for (i in 1:length(ULA3_Mn_g$col))
{
  if (ULA3_Mn_g$cluster[i]==2){
    Sv_3[i]<-ULA3_Mn_g$val[i]}
  else if (ULA3_Mn_g$cluster[i]==1){
    Sv_3[i]<-NA}
}

Sv_3<-cbind(ULA3_Mn_g$col, Sv_3)
colnames(Sv_3)<-c("col", "val")

ULA3_Mn_f<-pivot_wider(as.data.frame(Sv_3), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA3_Mn_f<-unnest(ULA3_Mn_f) 
heatmap(as.matrix(ULA3_Mn_f), Rowv=NA, Colv=NA, asp=1)

write.csv(ULA3_Ca_f, "ULA3_Ca_NA.csv")
write.csv(ULA3_Fe_f, "ULA3_Fe_NA.csv")
write.csv(ULA3_K_f, "ULA3_K_NA.csv")
write.csv(ULA3_Ti_f, "ULA3_Ti_NA.csv")
write.csv(ULA3_Si_f, "ULA3_Si_NA.csv")
write.csv(ULA3_Mn_f, "ULA3_Mn_NA.csv")









#-----------------------------ULA4--------------------------------------------------

#Gather dataframe into single column
ULA4_g<-gather(ULA4_sum, col, val, 1:915, factor_key=T)
ULA4_g[is.na(ULA4_g)]<-0

#carry out kmeans clustering:
ULA4_g_km <- kmeans(ULA4_g[,2],2,iter.max = 100, nstart = 20)

#spread dataset again with just clusters to check whether clustering worked and
#which cluster nr corresponds to the rock surface
S_4<-cbind(ULA4_g, ULA4_g_km$cluster)
colnames(S_4)[3]<-"cluster"
St_4<-S_4[,c(1,3)] #check whether the clustering works. If it does, turn selected value into NA
Sp_4<-pivot_wider(St_4, names_from=col, values_from=cluster)
Sp_4<-unnest(Sp_4) 
heatmap(as.matrix(Sp_4), Rowv=NA, Colv=NA, asp=1)
view(Sp_4) #cluster 1 = rock surface, cluster 2 = non-rock
#NB pay attention to this when generating your own data, might switch around due to the iterative nature of kmeans!



#create new dataframe with original data for rock surface and NA for non-rock surface
#for each elemental map
#Ca:
ULA4_Ca_g<-gather(ULA4_Ca, col, val, 1:915, factor_key=T)
ULA4_Ca_g[is.na(ULA4_Ca_g)]<-0
ULA4_Ca_g<-cbind(ULA4_Ca_g, S_4$cluster)
colnames(ULA4_Ca_g)[3]<-"cluster"
Sv_4<-array(data=NA, dim=length(ULA4_Ca_g$col))

for (i in 1:length(ULA4_Ca_g$col))
{
  if (ULA4_Ca_g$cluster[i]==1){
    Sv_4[i]<-ULA4_Ca_g$val[i]}
  else if (ULA4_Ca_g$cluster[i]==2){
    Sv_4[i]<-NA}
}

Sv_4<-cbind(ULA4_Ca_g$col, Sv_4)
colnames(Sv_4)<-c("col", "val")

ULA4_Ca_f<-pivot_wider(as.data.frame(Sv_4), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA4_Ca_f<-unnest(ULA4_Ca_f) 
heatmap(as.matrix(ULA4_Ca_f), Rowv=NA, Colv=NA, asp=1)



#Fe:
ULA4_Fe_g<-gather(ULA4_Fe, col, val, 1:915, factor_key=T)
ULA4_Fe_g[is.na(ULA4_Fe_g)]<-0
ULA4_Fe_g<-cbind(ULA4_Fe_g, S_4$cluster)
colnames(ULA4_Fe_g)[3]<-"cluster"
Sv_4<-array(data=NA, dim=length(ULA4_Fe_g$col))

for (i in 1:length(ULA4_Fe_g$col))
{
  if (ULA4_Fe_g$cluster[i]==1){
    Sv_4[i]<-ULA4_Fe_g$val[i]}
  else if (ULA4_Fe_g$cluster[i]==2){
    Sv_4[i]<-NA}
}

Sv_4<-cbind(ULA4_Fe_g$col, Sv_4)
colnames(Sv_4)<-c("col", "val")

ULA4_Fe_f<-pivot_wider(as.data.frame(Sv_4), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA4_Fe_f<-unnest(ULA4_Fe_f) 
heatmap(as.matrix(ULA4_Fe_f), Rowv=NA, Colv=NA, asp=1)



#K:
ULA4_K_g<-gather(ULA4_K, col, val, 1:915, factor_key=T)
ULA4_K_g[is.na(ULA4_K_g)]<-0
ULA4_K_g<-cbind(ULA4_K_g, S_4$cluster)
colnames(ULA4_K_g)[3]<-"cluster"
Sv_4<-array(data=NA, dim=length(ULA4_K_g$col))

for (i in 1:length(ULA4_K_g$col))
{
  if (ULA4_K_g$cluster[i]==1){
    Sv_4[i]<-ULA4_K_g$val[i]}
  else if (ULA4_K_g$cluster[i]==2){
    Sv_4[i]<-NA}
}

Sv_4<-cbind(ULA4_K_g$col, Sv_4)
colnames(Sv_4)<-c("col", "val")

ULA4_K_f<-pivot_wider(as.data.frame(Sv_4), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA4_K_f<-unnest(ULA4_K_f) 
heatmap(as.matrix(ULA4_K_f), Rowv=NA, Colv=NA, asp=1)



#Ti:
ULA4_Ti_g<-gather(ULA4_Ti, col, val, 1:915, factor_key=T)
ULA4_Ti_g[is.na(ULA4_Ti_g)]<-0
ULA4_Ti_g<-cbind(ULA4_Ti_g, S_4$cluster)
colnames(ULA4_Ti_g)[3]<-"cluster"
Sv_4<-array(data=NA, dim=length(ULA4_Ti_g$col))

for (i in 1:length(ULA4_Ti_g$col))
{
  if (ULA4_Ti_g$cluster[i]==1){
    Sv_4[i]<-ULA4_Ti_g$val[i]}
  else if (ULA4_Ti_g$cluster[i]==2){
    Sv_4[i]<-NA}
}

Sv_4<-cbind(ULA4_Ti_g$col, Sv_4)
colnames(Sv_4)<-c("col", "val")

ULA4_Ti_f<-pivot_wider(as.data.frame(Sv_4), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA4_Ti_f<-unnest(ULA4_Ti_f) 
heatmap(as.matrix(ULA4_Ti_f), Rowv=NA, Colv=NA, asp=1)



#Si:
ULA4_Si_g<-gather(ULA4_Si, col, val, 1:915, factor_key=T)
ULA4_Si_g[is.na(ULA4_Si_g)]<-0
ULA4_Si_g<-cbind(ULA4_Si_g, S_4$cluster)
colnames(ULA4_Si_g)[3]<-"cluster"
Sv_4<-array(data=NA, dim=length(ULA4_Si_g$col))

for (i in 1:length(ULA4_Si_g$col))
{
  if (ULA4_Si_g$cluster[i]==1){
    Sv_4[i]<-ULA4_Si_g$val[i]}
  else if (ULA4_Si_g$cluster[i]==2){
    Sv_4[i]<-NA}
}

Sv_4<-cbind(ULA4_Si_g$col, Sv_4)
colnames(Sv_4)<-c("col", "val")

ULA4_Si_f<-pivot_wider(as.data.frame(Sv_4), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA4_Si_f<-unnest(ULA4_Si_f) 
heatmap(as.matrix(ULA4_Si_f), Rowv=NA, Colv=NA, asp=1)



#Mn:
ULA4_Mn_g<-gather(ULA4_Mn, col, val, 1:915, factor_key=T)
ULA4_Mn_g[is.na(ULA4_Mn_g)]<-0
ULA4_Mn_g<-cbind(ULA4_Mn_g, S_4$cluster)
colnames(ULA4_Mn_g)[3]<-"cluster"
Sv_4<-array(data=NA, dim=length(ULA4_Mn_g$col))

for (i in 1:length(ULA4_Mn_g$col))
{
  if (ULA4_Mn_g$cluster[i]==1){
    Sv_4[i]<-ULA4_Mn_g$val[i]}
  else if (ULA4_Mn_g$cluster[i]==2){
    Sv_4[i]<-NA}
}

Sv_4<-cbind(ULA4_Mn_g$col, Sv_4)
colnames(Sv_4)<-c("col", "val")

ULA4_Mn_f<-pivot_wider(as.data.frame(Sv_4), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA4_Mn_f<-unnest(ULA4_Mn_f) 
heatmap(as.matrix(ULA4_Mn_f), Rowv=NA, Colv=NA, asp=1)

write.csv(ULA4_Ca_f, "ULA4_Ca_NA.csv")
write.csv(ULA4_Fe_f, "ULA4_Fe_NA.csv")
write.csv(ULA4_K_f, "ULA4_K_NA.csv")
write.csv(ULA4_Ti_f, "ULA4_Ti_NA.csv")
write.csv(ULA4_Si_f, "ULA4_Si_NA.csv")
write.csv(ULA4_Mn_f, "ULA4_Mn_NA.csv")








#-----------------------------ULA5--------------------------------------------------
#Gather dataframe into single column
ULA5_g<-gather(ULA5_sum, col, val, 1:801, factor_key=T)
ULA5_g[is.na(ULA5_g)]<-0

#carry out kmeans clustering:
ULA5_g_km <- kmeans(ULA5_g[,2],2,iter.max = 100, nstart = 20)

#spread dataset again with just clusters to check whether clustering worked and
#which cluster nr corresponds to the rock surface
S_5<-cbind(ULA5_g, ULA5_g_km$cluster)
colnames(S_5)[3]<-"cluster"
St_5<-S_5[,c(1,3)] #check whether the clustering works. If it does, turn selected value into NA
Sp_5<-pivot_wider(St_5, names_from=col, values_from=cluster)
Sp_5<-unnest(Sp_5) 
heatmap(as.matrix(Sp_5), Rowv=NA, Colv=NA, asp=1)
view(Sp_5) #cluster 2 = rock surface, cluster 1 = non-rock
#NB pay attention to this when generating your own data, might switch around due to the iterative nature of kmeans!


#create new dataframe with original data for rock surface and NA for non-rock surface
#for each elemental map
#Ca:
ULA5_Ca_g<-gather(ULA5_Ca, col, val, 1:801, factor_key=T)
ULA5_Ca_g[is.na(ULA5_Ca_g)]<-0
ULA5_Ca_g<-cbind(ULA5_Ca_g, S_5$cluster)
colnames(ULA5_Ca_g)[3]<-"cluster"
Sv_5<-array(data=NA, dim=length(ULA5_Ca_g$col))

for (i in 1:length(ULA5_Ca_g$col))
{
  if (ULA5_Ca_g$cluster[i]==1){
    Sv_5[i]<-ULA5_Ca_g$val[i]}
  else if (ULA5_Ca_g$cluster[i]==2){
    Sv_5[i]<-NA}
}

Sv_5<-cbind(ULA5_Ca_g$col, Sv_5)
colnames(Sv_5)<-c("col", "val")

ULA5_Ca_f<-pivot_wider(as.data.frame(Sv_5), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA5_Ca_f<-unnest(ULA5_Ca_f) 
heatmap(as.matrix(ULA5_Ca_f), Rowv=NA, Colv=NA, asp=1)


#Fe:
ULA5_Fe_g<-gather(ULA5_Fe, col, val, 1:801, factor_key=T)
ULA5_Fe_g[is.na(ULA5_Fe_g)]<-0
ULA5_Fe_g<-cbind(ULA5_Fe_g, S_5$cluster)
colnames(ULA5_Fe_g)[3]<-"cluster"
Sv_5<-array(data=NA, dim=length(ULA5_Fe_g$col))

for (i in 1:length(ULA5_Fe_g$col))
{
  if (ULA5_Fe_g$cluster[i]==1){
    Sv_5[i]<-ULA5_Fe_g$val[i]}
  else if (ULA5_Fe_g$cluster[i]==2){
    Sv_5[i]<-NA}
}

Sv_5<-cbind(ULA5_Fe_g$col, Sv_5)
colnames(Sv_5)<-c("col", "val")

ULA5_Fe_f<-pivot_wider(as.data.frame(Sv_5), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA5_Fe_f<-unnest(ULA5_Fe_f) 
heatmap(as.matrix(ULA5_Fe_f), Rowv=NA, Colv=NA, asp=1)



#K:
ULA5_K_g<-gather(ULA5_K, col, val, 1:801, factor_key=T)
ULA5_K_g[is.na(ULA5_K_g)]<-0
ULA5_K_g<-cbind(ULA5_K_g, S_5$cluster)
colnames(ULA5_K_g)[3]<-"cluster"
Sv_5<-array(data=NA, dim=length(ULA5_K_g$col))

for (i in 1:length(ULA5_K_g$col))
{
  if (ULA5_K_g$cluster[i]==1){
    Sv_5[i]<-ULA5_K_g$val[i]}
  else if (ULA5_K_g$cluster[i]==2){
    Sv_5[i]<-NA}
}

Sv_5<-cbind(ULA5_K_g$col, Sv_5)
colnames(Sv_5)<-c("col", "val")

ULA5_K_f<-pivot_wider(as.data.frame(Sv_5), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA5_K_f<-unnest(ULA5_K_f) 
heatmap(as.matrix(ULA5_K_f), Rowv=NA, Colv=NA, asp=1)



#Ti:
ULA5_Ti_g<-gather(ULA5_Ti, col, val, 1:801, factor_key=T)
ULA5_Ti_g[is.na(ULA5_Ti_g)]<-0
ULA5_Ti_g<-cbind(ULA5_Ti_g, S_5$cluster)
colnames(ULA5_Ti_g)[3]<-"cluster"
Sv_5<-array(data=NA, dim=length(ULA5_Ti_g$col))

for (i in 1:length(ULA5_Ti_g$col))
{
  if (ULA5_Ti_g$cluster[i]==1){
    Sv_5[i]<-ULA5_Ti_g$val[i]}
  else if (ULA5_Ti_g$cluster[i]==2){
    Sv_5[i]<-NA}
}

Sv_5<-cbind(ULA5_Ti_g$col, Sv_5)
colnames(Sv_5)<-c("col", "val")

ULA5_Ti_f<-pivot_wider(as.data.frame(Sv_5), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA5_Ti_f<-unnest(ULA5_Ti_f) 
heatmap(as.matrix(ULA5_Ti_f), Rowv=NA, Colv=NA, asp=1)



#Si:
ULA5_Si_g<-gather(ULA5_Si, col, val, 1:801, factor_key=T)
ULA5_Si_g[is.na(ULA5_Si_g)]<-0
ULA5_Si_g<-cbind(ULA5_Si_g, S_5$cluster)
colnames(ULA5_Si_g)[3]<-"cluster"
Sv_5<-array(data=NA, dim=length(ULA5_Si_g$col))

for (i in 1:length(ULA5_Si_g$col))
{
  if (ULA5_Si_g$cluster[i]==1){
    Sv_5[i]<-ULA5_Si_g$val[i]}
  else if (ULA5_Si_g$cluster[i]==2){
    Sv_5[i]<-NA}
}

Sv_5<-cbind(ULA5_Si_g$col, Sv_5)
colnames(Sv_5)<-c("col", "val")

ULA5_Si_f<-pivot_wider(as.data.frame(Sv_5), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA5_Si_f<-unnest(ULA5_Si_f) 
heatmap(as.matrix(ULA5_Si_f), Rowv=NA, Colv=NA, asp=1)



#Mn:
ULA5_Mn_g<-gather(ULA5_Mn, col, val, 1:801, factor_key=T)
ULA5_Mn_g[is.na(ULA5_Mn_g)]<-0
ULA5_Mn_g<-cbind(ULA5_Mn_g, S_5$cluster)
colnames(ULA5_Mn_g)[3]<-"cluster"
Sv_5<-array(data=NA, dim=length(ULA5_Mn_g$col))

for (i in 1:length(ULA5_Mn_g$col))
{
  if (ULA5_Mn_g$cluster[i]==1){
    Sv_5[i]<-ULA5_Mn_g$val[i]}
  else if (ULA5_Mn_g$cluster[i]==2){
    Sv_5[i]<-NA}
}

Sv_5<-cbind(ULA5_Mn_g$col, Sv_5)
colnames(Sv_5)<-c("col", "val")

ULA5_Mn_f<-pivot_wider(as.data.frame(Sv_5), names_from=col, values_from=val) #can't use spread() with non-unique key-value combos
ULA5_Mn_f<-unnest(ULA5_Mn_f) 
heatmap(as.matrix(ULA5_Mn_f), Rowv=NA, Colv=NA, asp=1)

write.csv(ULA5_Ca_f, "ULA5_Ca_NA.csv")
write.csv(ULA5_Fe_f, "ULA5_Fe_NA.csv")
write.csv(ULA5_K_f, "ULA5_K_NA.csv")
write.csv(ULA5_Ti_f, "ULA5_Ti_NA.csv")
write.csv(ULA5_Si_f, "ULA5_Si_NA.csv")
write.csv(ULA5_Mn_f, "ULA5_Mn_NA.csv")





#-------------------------Composite maps------------------------------------------
#Read in data; can't automatically use the above due to the iterative nature of kmeans.
#Therefore, use pre-saved data. 
#this data consists of the (when necessary) rotated, flipped, and cropped maps
#(see code above)

ULA1_Ca<-read.csv("ULA1_Ca_NA.csv")
ULA1_Fe<-read.csv("ULA1_Fe_NA.csv")
ULA1_K<-read.csv("ULA1_K_NA.csv")
ULA1_Ti<-read.csv("ULA1_Ti_NA.csv")
ULA1_Si<-read.csv("ULA1_Si_NA.csv")
ULA1_Mn<-read.csv("ULA1_Mn_NA.csv")

ULA2_Ca<-read.csv("ULA2_Ca_NA.csv")
ULA2_Fe<-read.csv("ULA2_Fe_NA.csv")
ULA2_K<-read.csv("ULA2_K_NA.csv")
ULA2_Ti<-read.csv("ULA2_Ti_NA.csv")
ULA2_Si<-read.csv("ULA2_Si_NA.csv")
ULA2_Mn<-read.csv("ULA2_Mn_NA.csv")

ULA3_Ca<-read.csv("ULA3_Ca_NA.csv")
ULA3_Fe<-read.csv("ULA3_Fe_NA.csv")
ULA3_K<-read.csv("ULA3_K_NA.csv")
ULA3_Ti<-read.csv("ULA3_Ti_NA.csv")
ULA3_Si<-read.csv("ULA3_Si_NA.csv")
ULA3_Mn<-read.csv("ULA3_Mn_NA.csv")

ULA4_Ca<-read.csv("ULA4_Ca_NA.csv")
ULA4_Fe<-read.csv("ULA4_Fe_NA.csv")
ULA4_K<-read.csv("ULA4_K_NA.csv")
ULA4_Ti<-read.csv("ULA4_Ti_NA.csv")
ULA4_Si<-read.csv("ULA4_Si_NA.csv")
ULA4_Mn<-read.csv("ULA4_Mn_NA.csv")

ULA5_Ca<-read.csv("ULA5_Ca_NA.csv")
ULA5_Fe<-read.csv("ULA5_Fe_NA.csv")
ULA5_K<-read.csv("ULA5_K_NA.csv")
ULA5_Ti<-read.csv("ULA5_Ti_NA.csv")
ULA5_Si<-read.csv("ULA5_Si_NA.csv")
ULA5_Mn<-read.csv("ULA5_Mn_NA.csv")



#-----------------------------formatting------------------------------------------
#remove row index that R adds when saving csv
#(column index is added too, but this is recognized)
ULA1_Ca<-ULA1_Ca[,2:ncol(ULA1_Ca)]
ULA1_Fe<-ULA1_Fe[,2:ncol(ULA1_Fe)]
ULA1_K<-ULA1_K[,2:ncol(ULA1_K)]
ULA1_Ti<-ULA1_Ti[,2:ncol(ULA1_Ti)]
ULA1_Si<-ULA1_Si[,2:ncol(ULA1_Si)]
ULA1_Mn<-ULA1_Mn[,2:ncol(ULA1_Mn)]

ULA2_Ca<-ULA2_Ca[,2:ncol(ULA2_Ca)]
ULA2_Fe<-ULA2_Fe[,2:ncol(ULA2_Fe)]
ULA2_K<-ULA2_K[,2:ncol(ULA2_K)]
ULA2_Ti<-ULA2_Ti[,2:ncol(ULA2_Ti)]
ULA2_Si<-ULA2_Si[,2:ncol(ULA2_Si)]
ULA2_Mn<-ULA2_Mn[,2:ncol(ULA2_Mn)]

ULA3_Ca<-ULA3_Ca[,2:ncol(ULA3_Ca)]
ULA3_Fe<-ULA3_Fe[,2:ncol(ULA3_Fe)]
ULA3_K<-ULA3_K[,2:ncol(ULA3_K)]
ULA3_Ti<-ULA3_Ti[,2:ncol(ULA3_Ti)]
ULA3_Si<-ULA3_Si[,2:ncol(ULA3_Si)]
ULA3_Mn<-ULA3_Mn[,2:ncol(ULA3_Mn)]

ULA4_Ca<-ULA4_Ca[,2:ncol(ULA4_Ca)]
ULA4_Fe<-ULA4_Fe[,2:ncol(ULA4_Fe)]
ULA4_K<-ULA4_K[,2:ncol(ULA4_K)]
ULA4_Ti<-ULA4_Ti[,2:ncol(ULA4_Ti)]
ULA4_Si<-ULA4_Si[,2:ncol(ULA4_Si)]
ULA4_Mn<-ULA4_Mn[,2:ncol(ULA4_Mn)]

ULA5_Ca<-ULA5_Ca[,2:ncol(ULA5_Ca)]
ULA5_Fe<-ULA5_Fe[,2:ncol(ULA5_Fe)]
ULA5_K<-ULA5_K[,2:ncol(ULA5_K)]
ULA5_Ti<-ULA5_Ti[,2:ncol(ULA5_Ti)]
ULA5_Si<-ULA5_Si[,2:ncol(ULA5_Si)]
ULA5_Mn<-ULA5_Mn[,2:ncol(ULA5_Mn)]



#----------------------------determine overlap-------------------------------------
#--------------------------Ca----------------------------------------------------

#ULA1-ULA2:
heatmap(as.matrix(ULA1_Ca[,1:915]), Rowv=NA, Colv=NA, asp=1) #at 915, 31 left
heatmap(as.matrix(ULA2_Ca[,45:1092]), Rowv=NA, Colv=NA, asp=1) #at 45, 45 left
#31+45=76 overlap, (1092+946)-76=1965 length for both
ULA1_Ca_padded<-padzeros(data=ULA1_Ca, nzeros=(1965-946), side="right")
ULA2_Ca_padded<-padzeros(data=ULA2_Ca, nzeros=(1965-1092), side="left")
colnames(ULA1_Ca_padded)<-c(1:1965)
colnames(ULA2_Ca_padded)<-c(1:1965)
ULA12_Ca_padded<-rbind(ULA2_Ca_padded, ULA1_Ca_padded)
heatmap(as.matrix(ULA12_Ca_padded), Rowv=NA, Colv=NA, asp=1) 
#Ok. Don't worry about the zeroes, convert everything to NA in the end.

#ULA2-3:
heatmap(as.matrix(ULA12_Ca_padded[,1:1910]), Rowv=NA, Colv=NA, asp=1)  #at 1910, 55 left
heatmap(as.matrix(ULA3_Ca[,60:1144]), Rowv=NA, Colv=NA, asp=1) #at 60, 60 left
#60+55=115 overlap, (1965+1144)-115=2994 length for both
ULA12_Ca_padded<-padzeros(data=ULA12_Ca_padded, nzeros=(2994-1965), side="right")
ULA3_Ca_padded<-padzeros(data=ULA3_Ca, nzeros=(2994-1144), side="left")
colnames(ULA12_Ca_padded)<-c(1:2994)
colnames(ULA3_Ca_padded)<-c(1:2994)
ULA123_Ca_padded<-rbind(ULA3_Ca_padded, ULA12_Ca_padded)
heatmap(as.matrix(ULA123_Ca_padded), Rowv=NA, Colv=NA, asp=1) 

#ULA3-4
heatmap(as.matrix(ULA123_Ca_padded[,1:2820]), Rowv=NA, Colv=NA, asp=1)  #at 2820, 174 left
heatmap(as.matrix(ULA4_Ca), Rowv=NA, Colv=NA, asp=1) #at 0, 0 left
#0+174=174 overlap, (2994+915)-174=3735 length for both
ULA123_Ca_padded<-padzeros(data=ULA123_Ca_padded, nzeros=(3735-2994), side="right")
ULA4_Ca_padded<-padzeros(data=ULA4_Ca, nzeros=(3735-915), side="left")
colnames(ULA123_Ca_padded)<-c(1:3735)
colnames(ULA4_Ca_padded)<-c(1:3735)
ULA1234_Ca_padded<-rbind(ULA123_Ca_padded, ULA4_Ca_padded)
heatmap(as.matrix(ULA1234_Ca_padded), Rowv=NA, Colv=NA, asp=1) 

#ULA4-5
heatmap(as.matrix(ULA1234_Ca_padded[,1:3665]), Rowv=NA, Colv=NA, asp=1)  #at 3665, 70 left
heatmap(as.matrix(ULA5_Ca), Rowv=NA, Colv=NA, asp=1) #at 0, 0 left
#0+70=70 overlap, (3735+802)-70=4467 length for both
ULA1234_Ca_padded<-padzeros(data=ULA1234_Ca_padded, nzeros=(4467-3735), side="right")
ULA5_Ca_padded<-padzeros(data=ULA5_Ca, nzeros=(4467-801), side="left")
colnames(ULA1234_Ca_padded)<-c(1:4467)
colnames(ULA5_Ca_padded)<-c(1:4467)
ULA12345_Ca_padded<-rbind(ULA5_Ca_padded, ULA1234_Ca_padded)
heatmap(as.matrix(ULA12345_Ca_padded), Rowv=NA, Colv=NA, asp=1) 


#Rename full version
ULA_Ca_padded<-ULA12345_Ca_padded

#Convert zeroes to NA as they would disrupt averaging
ULA_Ca_padded[ULA_Ca_padded==0]<-NA
heatmap(as.matrix(ULA_Ca_padded), Rowv=NA, Colv=NA, asp=1) 









#--------------------------Fe----------------------------------------------------

#ULA1-ULA2:
ULA1_Fe_padded<-padzeros(data=ULA1_Fe, nzeros=(1965-946), side="right")
ULA2_Fe_padded<-padzeros(data=ULA2_Fe, nzeros=(1965-1092), side="left")
colnames(ULA1_Fe_padded)<-c(1:1965)
colnames(ULA2_Fe_padded)<-c(1:1965)
ULA12_Fe_padded<-rbind(ULA2_Fe_padded, ULA1_Fe_padded)

#ULA2-3:
ULA12_Fe_padded<-padzeros(data=ULA12_Fe_padded, nzeros=(2994-1965), side="right")
ULA3_Fe_padded<-padzeros(data=ULA3_Fe, nzeros=(2994-1144), side="left")
colnames(ULA12_Fe_padded)<-c(1:2994)
colnames(ULA3_Fe_padded)<-c(1:2994)
ULA123_Fe_padded<-rbind(ULA3_Fe_padded, ULA12_Fe_padded)

#ULA3-4
ULA123_Fe_padded<-padzeros(data=ULA123_Fe_padded, nzeros=(3735-2994), side="right")
ULA4_Fe_padded<-padzeros(data=ULA4_Fe, nzeros=(3735-915), side="left")
colnames(ULA123_Fe_padded)<-c(1:3735)
colnames(ULA4_Fe_padded)<-c(1:3735)
ULA1234_Fe_padded<-rbind(ULA123_Fe_padded, ULA4_Fe_padded)

#ULA4-5
ULA1234_Fe_padded<-padzeros(data=ULA1234_Fe_padded, nzeros=(4467-3735), side="right")
ULA5_Fe_padded<-padzeros(data=ULA5_Fe, nzeros=(4467-801), side="left")
colnames(ULA1234_Fe_padded)<-c(1:4467)
colnames(ULA5_Fe_padded)<-c(1:4467)
ULA12345_Fe_padded<-rbind(ULA5_Fe_padded, ULA1234_Fe_padded)

#Rename full version
ULA_Fe_padded<-ULA12345_Fe_padded

#Convert zeroes to NA as they would disrupt averaging
ULA_Fe_padded[ULA_Fe_padded==0]<-NA
heatmap(as.matrix(ULA_Fe_padded), Rowv=NA, Colv=NA, asp=1) 








#--------------------------K----------------------------------------------------

#ULA1-ULA2:
ULA1_K_padded<-padzeros(data=ULA1_K, nzeros=(1965-946), side="right")
ULA2_K_padded<-padzeros(data=ULA2_K, nzeros=(1965-1092), side="left")
colnames(ULA1_K_padded)<-c(1:1965)
colnames(ULA2_K_padded)<-c(1:1965)
ULA12_K_padded<-rbind(ULA2_K_padded, ULA1_K_padded)

#ULA2-3:
ULA12_K_padded<-padzeros(data=ULA12_K_padded, nzeros=(2994-1965), side="right")
ULA3_K_padded<-padzeros(data=ULA3_K, nzeros=(2994-1144), side="left")
colnames(ULA12_K_padded)<-c(1:2994)
colnames(ULA3_K_padded)<-c(1:2994)
ULA123_K_padded<-rbind(ULA3_K_padded, ULA12_K_padded)

#ULA3-4
ULA123_K_padded<-padzeros(data=ULA123_K_padded, nzeros=(3735-2994), side="right")
ULA4_K_padded<-padzeros(data=ULA4_K, nzeros=(3735-915), side="left")
colnames(ULA123_K_padded)<-c(1:3735)
colnames(ULA4_K_padded)<-c(1:3735)
ULA1234_K_padded<-rbind(ULA123_K_padded, ULA4_K_padded)

#ULA4-5
ULA1234_K_padded<-padzeros(data=ULA1234_K_padded, nzeros=(4467-3735), side="right")
ULA5_K_padded<-padzeros(data=ULA5_K, nzeros=(4467-801), side="left")
colnames(ULA1234_K_padded)<-c(1:4467)
colnames(ULA5_K_padded)<-c(1:4467)
ULA12345_K_padded<-rbind(ULA5_K_padded, ULA1234_K_padded)


#Rename full version
ULA_K_padded<-ULA12345_K_padded

#Convert zeroes to NA as they would disrupt averaging
ULA_K_padded[ULA_K_padded==0]<-NA
heatmap(as.matrix(ULA_K_padded), Rowv=NA, Colv=NA, asp=1) 






#-------------------------Ti----------------------------------------------------

#ULA1-ULA2:
ULA1_Ti_padded<-padzeros(data=ULA1_Ti, nzeros=(1965-946), side="right")
ULA2_Ti_padded<-padzeros(data=ULA2_Ti, nzeros=(1965-1092), side="left")
colnames(ULA1_Ti_padded)<-c(1:1965)
colnames(ULA2_Ti_padded)<-c(1:1965)
ULA12_Ti_padded<-rbind(ULA2_Ti_padded, ULA1_Ti_padded)

#ULA2-3:
ULA12_Ti_padded<-padzeros(data=ULA12_Ti_padded, nzeros=(2994-1965), side="right")
ULA3_Ti_padded<-padzeros(data=ULA3_Ti, nzeros=(2994-1144), side="left")
colnames(ULA12_Ti_padded)<-c(1:2994)
colnames(ULA3_Ti_padded)<-c(1:2994)
ULA123_Ti_padded<-rbind(ULA3_Ti_padded, ULA12_Ti_padded)

#ULA3-4
ULA123_Ti_padded<-padzeros(data=ULA123_Ti_padded, nzeros=(3735-2994), side="right")
ULA4_Ti_padded<-padzeros(data=ULA4_Ti, nzeros=(3735-915), side="left")
colnames(ULA123_Ti_padded)<-c(1:3735)
colnames(ULA4_Ti_padded)<-c(1:3735)
ULA1234_Ti_padded<-rbind(ULA123_Ti_padded, ULA4_Ti_padded)

#ULA4-5
ULA1234_Ti_padded<-padzeros(data=ULA1234_Ti_padded, nzeros=(4467-3735), side="right")
ULA5_Ti_padded<-padzeros(data=ULA5_Ti, nzeros=(4467-801), side="left")
colnames(ULA1234_Ti_padded)<-c(1:4467)
colnames(ULA5_Ti_padded)<-c(1:4467)
ULA12345_Ti_padded<-rbind(ULA5_Ti_padded, ULA1234_Ti_padded)


#Rename full version
ULA_Ti_padded<-ULA12345_Ti_padded

#Convert zeroes to NA as they would disrupt averaging
ULA_Ti_padded[ULA_Ti_padded==0]<-NA
heatmap(as.matrix(ULA_Ti_padded), Rowv=NA, Colv=NA, asp=1) 





#-------------------------Si----------------------------------------------------

#ULA1-ULA2:
ULA1_Si_padded<-padzeros(data=ULA1_Si, nzeros=(1965-946), side="right")
ULA2_Si_padded<-padzeros(data=ULA2_Si, nzeros=(1965-1092), side="left")
colnames(ULA1_Si_padded)<-c(1:1965)
colnames(ULA2_Si_padded)<-c(1:1965)
ULA12_Si_padded<-rbind(ULA2_Si_padded, ULA1_Si_padded)

#ULA2-3:
ULA12_Si_padded<-padzeros(data=ULA12_Si_padded, nzeros=(2994-1965), side="right")
ULA3_Si_padded<-padzeros(data=ULA3_Si, nzeros=(2994-1144), side="left")
colnames(ULA12_Si_padded)<-c(1:2994)
colnames(ULA3_Si_padded)<-c(1:2994)
ULA123_Si_padded<-rbind(ULA3_Si_padded, ULA12_Si_padded)

#ULA3-4
ULA123_Si_padded<-padzeros(data=ULA123_Si_padded, nzeros=(3735-2994), side="right")
ULA4_Si_padded<-padzeros(data=ULA4_Si, nzeros=(3735-915), side="left")
colnames(ULA123_Si_padded)<-c(1:3735)
colnames(ULA4_Si_padded)<-c(1:3735)
ULA1234_Si_padded<-rbind(ULA123_Si_padded, ULA4_Si_padded)

#ULA4-5
ULA1234_Si_padded<-padzeros(data=ULA1234_Si_padded, nzeros=(4467-3735), side="right")
ULA5_Si_padded<-padzeros(data=ULA5_Si, nzeros=(4467-801), side="left")
colnames(ULA1234_Si_padded)<-c(1:4467)
colnames(ULA5_Si_padded)<-c(1:4467)
ULA12345_Si_padded<-rbind(ULA5_Si_padded, ULA1234_Si_padded)


#Rename full version
ULA_Si_padded<-ULA12345_Si_padded

#Convert zeroes to NA as they would disrupt averaging
ULA_Si_padded[ULA_Si_padded==0]<-NA
heatmap(as.matrix(ULA_Si_padded), Rowv=NA, Colv=NA, asp=1) 






#-------------------------Mn----------------------------------------------------

#ULA1-ULA2:
ULA1_Mn_padded<-padzeros(data=ULA1_Mn, nzeros=(1965-946), side="right")
ULA2_Mn_padded<-padzeros(data=ULA2_Mn, nzeros=(1965-1092), side="left")
colnames(ULA1_Mn_padded)<-c(1:1965)
colnames(ULA2_Mn_padded)<-c(1:1965)
ULA12_Mn_padded<-rbind(ULA2_Mn_padded, ULA1_Mn_padded)

#ULA2-3:
ULA12_Mn_padded<-padzeros(data=ULA12_Mn_padded, nzeros=(2994-1965), side="right")
ULA3_Mn_padded<-padzeros(data=ULA3_Mn, nzeros=(2994-1144), side="left")
colnames(ULA12_Mn_padded)<-c(1:2994)
colnames(ULA3_Mn_padded)<-c(1:2994)
ULA123_Mn_padded<-rbind(ULA3_Mn_padded, ULA12_Mn_padded)

#ULA3-4
ULA123_Mn_padded<-padzeros(data=ULA123_Mn_padded, nzeros=(3735-2994), side="right")
ULA4_Mn_padded<-padzeros(data=ULA4_Mn, nzeros=(3735-915), side="left")
colnames(ULA123_Mn_padded)<-c(1:3735)
colnames(ULA4_Mn_padded)<-c(1:3735)
ULA1234_Mn_padded<-rbind(ULA123_Mn_padded, ULA4_Mn_padded)

#ULA4-5
ULA1234_Mn_padded<-padzeros(data=ULA1234_Mn_padded, nzeros=(4467-3735), side="right")
ULA5_Mn_padded<-padzeros(data=ULA5_Mn, nzeros=(4467-801), side="left")
colnames(ULA1234_Mn_padded)<-c(1:4467)
colnames(ULA5_Mn_padded)<-c(1:4467)
ULA12345_Mn_padded<-rbind(ULA5_Mn_padded, ULA1234_Mn_padded)


#Rename full verMnon
ULA_Mn_padded<-ULA12345_Mn_padded

#Convert zeroes to NA as they would disrupt averaging
ULA_Mn_padded[ULA_Mn_padded==0]<-NA
heatmap(as.matrix(ULA_Mn_padded), Rowv=NA, Colv=NA, asp=1) 







#write everything to csv:
write.csv(ULA_Ca_padded, "ULA_comp_map_Ca.csv") 
write.csv(ULA_Fe_padded, "ULA_comp_map_Fe.csv") 
write.csv(ULA_K_padded, "ULA_comp_map_K.csv") 
write.csv(ULA_Ti_padded, "ULA_comp_map_Ti.csv") 
write.csv(ULA_Si_padded, "ULA_comp_map_Si.csv")
write.csv(ULA_Mn_padded, "ULA_comp_map_Mn.csv")




#------------------------Composite records-----------------------------------------
#read in data so the whole script doesn't have to be run
ULA_Ca<-read.csv("ULA_comp_map_Ca.csv")
ULA_Fe<-read.csv("ULA_comp_map_Fe.csv")
ULA_K<-read.csv("ULA_comp_map_K.csv")
ULA_Ti<-read.csv("ULA_comp_map_Ti.csv")
ULA_Si<-read.csv("ULA_comp_map_Si.csv")
ULA_Mn<-read.csv("ULA_comp_map_Mn.csv")

#Once again remove row index
ULA_Ca<-ULA_Ca[,2:ncol(ULA_Ca)]
ULA_Fe<-ULA_Fe[,2:ncol(ULA_Fe)]
ULA_K<-ULA_K[,2:ncol(ULA_K)]
ULA_Ti<-ULA_Ti[,2:ncol(ULA_Ti)]
ULA_Si<-ULA_Si[,2:ncol(ULA_Si)]
ULA_Mn<-ULA_Mn[,2:ncol(ULA_Mn)]


#------------------------Ca-------------------------------------------------------
#Average data:
ULA_Ca_mean<-colMeans(ULA_Ca, na.rm=T)

#add depth column; slightly complicated because of varying step sizes
#step size was 100 um = 0.1 mm for ULA2, ULA3, ULA4, 946-3665 2719
#95 um=0.095 mm for ULA1, 1-946 946
#133 um= 0.133 mm for ULA5, 3665-4518 853
#ignore overlaps - we are looking at cm, not micro- or mm scale changes. 

#add depth columns
ULA_Ca_d1<-seq(from=0, to=89.775, by=0.095)
ULA_Ca_d234<-seq(from=max(ULA_Ca_d1)+0.1, to=271.8+max(ULA_Ca_d1)+0.1, by=0.1)
ULA_Ca_d5<-seq(from=max(ULA_Ca_d234)+0.133, to=113.316+max(ULA_Ca_d234)+0.133, by=0.133)
ULA_Ca_d<-c(ULA_Ca_d1, ULA_Ca_d234, ULA_Ca_d5)

ULA_Ca_depth<-cbind(ULA_Ca_d, ULA_Ca_mean)
ULA_Ca_depth<-as.data.frame(ULA_Ca_depth[1:4516,])
colnames(ULA_Ca_depth)<-c("depth_mm", "Ca")

#plot
ggplot()+
  geom_line(data=ULA_Ca_depth, aes(x=depth_mm, y=Ca), col="navy")+
  theme_classic()


#------------------------Fe-------------------------------------------------------
#Average data:
ULA_Fe_mean<-colMeans(ULA_Fe, na.rm=T)

#add depth columns
ULA_Fe_d1<-seq(from=0, to=89.775, by=0.095)
ULA_Fe_d234<-seq(from=max(ULA_Fe_d1)+0.1, to=271.8+max(ULA_Fe_d1)+0.1, by=0.1)
ULA_Fe_d5<-seq(from=max(ULA_Fe_d234)+0.133, to=113.316+max(ULA_Fe_d234)+0.133, by=0.133)
ULA_Fe_d<-c(ULA_Fe_d1, ULA_Fe_d234, ULA_Fe_d5)

ULA_Fe_depth<-cbind(ULA_Fe_d, ULA_Fe_mean)
ULA_Fe_depth<-as.data.frame(ULA_Fe_depth[1:4516,])
colnames(ULA_Fe_depth)<-c("depth_mm", "Fe")

#plot
ggplot()+
  geom_line(data=ULA_Fe_depth, aes(x=depth_mm, y=Fe), col="navy")+
  theme_classic()



#------------------------K-------------------------------------------------------
#Average/Sum data:
ULA_K_mean<-colMeans(ULA_K, na.rm=T)

#add depth columns
ULA_K_d1<-seq(from=0, to=89.775, by=0.095)
ULA_K_d234<-seq(from=max(ULA_K_d1)+0.1, to=271.8+max(ULA_K_d1)+0.1, by=0.1)
ULA_K_d5<-seq(from=max(ULA_K_d234)+0.133, to=113.316+max(ULA_K_d234)+0.133, by=0.133)
ULA_K_d<-c(ULA_K_d1, ULA_K_d234, ULA_K_d5)

ULA_K_depth<-cbind(ULA_K_d, ULA_K_mean)
ULA_K_depth<-as.data.frame(ULA_K_depth[1:4516,])
colnames(ULA_K_depth)<-c("depth_mm", "K")

#plot
ggplot()+
  geom_line(data=ULA_K_depth, aes(x=depth_mm, y=K), col="navy")+
  theme_classic()



#------------------------Ti-------------------------------------------------------
#Average/Sum data:
ULA_Ti_mean<-colMeans(ULA_Ti, na.rm=T)

#add depth columns
ULA_Ti_d1<-seq(from=0, to=89.775, by=0.095)
ULA_Ti_d234<-seq(from=max(ULA_Ti_d1)+0.1, to=271.8+max(ULA_Ti_d1)+0.1, by=0.1)
ULA_Ti_d5<-seq(from=max(ULA_Ti_d234)+0.133, to=113.316+max(ULA_Ti_d234)+0.133, by=0.133)
ULA_Ti_d<-c(ULA_Ti_d1, ULA_Ti_d234, ULA_Ti_d5)

ULA_Ti_depth<-cbind(ULA_Ti_d, ULA_Ti_mean)
ULA_Ti_depth<-as.data.frame(ULA_Ti_depth[1:4516,])
colnames(ULA_Ti_depth)<-c("depth_mm", "Ti")

#plot
ggplot()+
  geom_line(data=ULA_Ti_depth, aes(x=depth_mm, y=Ti), col="navy")+
  theme_classic()



#------------------------Si-------------------------------------------------------
#Average/Sum data:
ULA_Si_mean<-colMeans(ULA_Si, na.rm=T)

#add depth columns
ULA_Si_d1<-seq(from=0, to=89.775, by=0.095)
ULA_Si_d234<-seq(from=max(ULA_Si_d1)+0.1, to=271.8+max(ULA_Si_d1)+0.1, by=0.1)
ULA_Si_d5<-seq(from=max(ULA_Si_d234)+0.133, to=113.316+max(ULA_Si_d234)+0.133, by=0.133)
ULA_Si_d<-c(ULA_Si_d1, ULA_Si_d234, ULA_Si_d5)

ULA_Si_depth<-cbind(ULA_Si_d, ULA_Si_mean)
ULA_Si_depth<-as.data.frame(ULA_Si_depth[1:4516,])
colnames(ULA_Si_depth)<-c("depth_mm", "Si")

#plot
ggplot()+
  geom_line(data=ULA_Si_depth, aes(x=depth_mm, y=Si), col="navy")+
  theme_classic()



#------------------------Mn-------------------------------------------------------
#Average/Sum data:
ULA_Mn_mean<-colMeans(ULA_Mn, na.rm=T)

#add depth columns
ULA_Mn_d1<-seq(from=0, to=89.775, by=0.095)
ULA_Mn_d234<-seq(from=max(ULA_Mn_d1)+0.1, to=271.8+max(ULA_Mn_d1)+0.1, by=0.1)
ULA_Mn_d5<-seq(from=max(ULA_Mn_d234)+0.133, to=113.316+max(ULA_Mn_d234)+0.133, by=0.133)
ULA_Mn_d<-c(ULA_Mn_d1, ULA_Mn_d234, ULA_Mn_d5)

ULA_Mn_depth<-cbind(ULA_Mn_d, ULA_Mn_mean)
ULA_Mn_depth<-as.data.frame(ULA_Mn_depth[1:4516,])
colnames(ULA_Mn_depth)<-c("depth_mm", "Mn")

#plot
ggplot()+
  geom_line(data=ULA_Mn_depth, aes(x=depth_mm, y=Mn), col="navy")+
  theme_classic()




#write everything to csv:
write.csv(ULA_Ca_depth, "ULA_Ca_depth.csv")
write.csv(ULA_Fe_depth, "ULA_Fe_depth.csv")
write.csv(ULA_K_depth, "ULA_K_depth.csv")
write.csv(ULA_Ti_depth, "ULA_Ti_depth.csv")
write.csv(ULA_Si_depth, "ULA_Si_depth.csv")
write.csv(ULA_Mn_depth, "ULA_Mn_depth.csv")



#------------------------Quick plotting of records-------------------------------
#read in data:
ULA_Ca_depth<-read.csv("ULA_Ca_depth.csv")
ULA_Fe_depth<-read.csv("ULA_Fe_depth.csv")
ULA_K_depth<-read.csv("ULA_K_depth.csv")
ULA_Ti_depth<-read.csv("ULA_Ti_depth.csv")
ULA_Mn_depth<-read.csv("ULA_Mn_depth.csv")
ULA_Si_depth<-read.csv("ULA_Si_depth.csv")

#formatting: remove first column
ULA_Ca_depth<-ULA_Ca_depth[,2:3]
ULA_Fe_depth<-ULA_Fe_depth[,2:3]
ULA_K_depth<-ULA_K_depth[,2:3]
ULA_Ti_depth<-ULA_Ti_depth[,2:3]
ULA_Mn_depth<-ULA_Mn_depth[,2:3]
ULA_Si_depth<-ULA_Si_depth[,2:3]

#Plot individual records:
ggplot()+
  geom_line(data=ULA_Ca_depth, aes(x=depth_mm, y=Ca), col="navy")+
  theme_classic()

ggplot()+
  geom_line(data=ULA_Fe_depth, aes(x=depth_mm, y=Fe), col="navy")+
  theme_classic()

ggplot()+
  geom_line(data=ULA_K_depth, aes(x=depth_mm, y=K), col="navy")+
  theme_classic()

ggplot()+
  geom_line(data=ULA_Ti_depth, aes(x=depth_mm, y=Ti), col="navy")+
  theme_classic()

ggplot()+
  geom_line(data=ULA_Si_depth, aes(x=depth_mm, y=Si), col="navy")+
  theme_classic()

ggplot()+
  geom_line(data=ULA_Mn_depth, aes(x=depth_mm, y=Mn), col="navy")+
  theme_classic()



#Plot together:
Ca_p<-ggplot()+
  geom_line(data=ULA_Ca_depth, aes(x=depth_mm, y=Ca), col="navy")+
  coord_flip()+
  theme_classic()

Fe_p<-ggplot()+
  geom_line(data=ULA_Fe_depth, aes(x=depth_mm, y=Fe), col="navy")+
  coord_flip()+
  theme_classic()

K_p<-ggplot()+
  geom_line(data=ULA_K_depth, aes(x=depth_mm, y=K), col="navy")+
  coord_flip()+
  theme_classic()

Ti_p<-ggplot()+
  geom_line(data=ULA_Ti_depth, aes(x=depth_mm, y=Ti), col="navy")+
  coord_flip()+
  theme_classic()

Si_p<-ggplot()+
  geom_line(data=ULA_Si_depth, aes(x=depth_mm, y=Si), col="navy")+
  coord_flip()+
  theme_classic()

Mn_p<-ggplot()+
  geom_line(data=ULA_Mn_depth, aes(x=depth_mm, y=Mn), col="navy")+
  coord_flip()+
  theme_classic()


cowplot::plot_grid(Ca_p, Fe_p, K_p, Ti_p, Si_p, Mn_p, ncol=6)

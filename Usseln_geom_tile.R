#This script contains code to generate the elemental heatmaps as shown in the paper
#using ggplot's geom_tile function

library(tidyverse)
library(viridis)

setwd("insert wd")
#----------------------------WINSENBERG (ULW)-------------------------------------------

#------------------------------load data-----------------------------------------
ulw_a1_Ca<-read.csv("ULWa1_Ca.csv")
ulw_a1_Fe<-read.csv("ULWa1_Fe.csv")
ulw_a1_Ti<-read.csv("ULWa1_Ti.csv")
ulw_a1_K<-read.csv("ULWa1_K.csv")

#check (optional)
heatmap(as.matrix(ulw_a1_Ca), Rowv=NA, Colv=NA, asp=1)

#------------------------------formatting----------------------------------------
#change colnames; keep first column as it is needed later
colnames(ulw_a1_Ca)<-c(1:dim(ulw_a1_Ca)[2])
colnames(ulw_a1_Fe)<-c(1:dim(ulw_a1_Fe)[2])
colnames(ulw_a1_Ti)<-c(1:dim(ulw_a1_Ti)[2])
colnames(ulw_a1_K)<-c(1:dim(ulw_a1_K)[2])

#add prefix, easier to refer to
colnames(ulw_a1_Ca)[2:527] <- paste('X', colnames(ulw_a1_Ca[2:527]), sep="") 
colnames(ulw_a1_Ca)[1]<-"rows"

colnames(ulw_a1_Fe)[2:527] <- paste('X', colnames(ulw_a1_Fe[2:527]), sep="") 
colnames(ulw_a1_Fe)[1]<-"rows"

colnames(ulw_a1_K)[2:527] <- paste('X', colnames(ulw_a1_K[2:527]), sep="") 
colnames(ulw_a1_K)[1]<-"rows"

colnames(ulw_a1_Ti)[2:527] <- paste('X', colnames(ulw_a1_Ti[2:527]), sep="") 
colnames(ulw_a1_Ti)[1]<-"rows"

#to long format
ulw_a1_Ca<-gather(ulw_a1_Ca, cols, val, X2:X527, factor_key = T)
ulw_a1_Fe<-gather(ulw_a1_Fe, cols, val, X2:X527, factor_key = T)
ulw_a1_Ti<-gather(ulw_a1_Ti, cols, val, X2:X527, factor_key = T)
ulw_a1_K<-gather(ulw_a1_K, cols, val, X2:X527, factor_key = T)


#create (K+Ti)/Ca record
kt<-(ulw_a1_K$val+ulw_a1_Ti$val)/ulw_a1_Ca$val
ulw_a1_r<-cbind(ulw_a1_Ca$rows, ulw_a1_Ca$cols, kt)
colnames(ulw_a1_r)<-c("rows", "cols", "val")
ulw_a1_r<-as.data.frame(ulw_a1_r)


#-----------------------------plot heatmaps----------------------------------------
#Ca
ggplot()+
  geom_tile(data=ulw_a1_Ca, aes(x=cols, y=rows, fill=val)) + 
  scale_fill_viridis(na.value=NA, option="D")+
  coord_fixed()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#Fe
ggplot()+
  geom_tile(data=ulw_a1_Fe, aes(x=cols, y=rows, fill=val)) + 
  scale_fill_viridis(na.value=NA, option="D")+
  coord_fixed()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#(K+Ti)/Ca
ggplot()+
  geom_tile(data=ulw_a1_r, aes(x=cols, y=rows, fill=val)) + 
  scale_fill_viridis(na.value=NA, option="F", limits=c(0, 0.5))+
  coord_fixed()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())






#---------------------------ARFELD (ULA)-----------------------------------------------

#------------------------------load data-----------------------------------------
ula_2_Ca<-read.csv("ULA2_Ca_NA.csv")
ula_2_Fe<-read.csv("ULA2_Fe_NA.csv")
ula_2_Ti<-read.csv("ULA2_Ti_NA.csv")
ula_2_K<-read.csv("ULA2_K_NA.csv")


#check (optional)
heatmap(as.matrix(ula_2_Ca), Rowv=NA, Colv=NA, asp=1)

#------------------------------formatting----------------------------------------
#change colnames; keep first column as it is needed later
colnames(ula_2_Ca)<-c(1:dim(ula_2_Ca)[2])
colnames(ula_2_Fe)<-c(1:dim(ula_2_Fe)[2])
colnames(ula_2_Ti)<-c(1:dim(ula_2_Ti)[2])
colnames(ula_2_K)<-c(1:dim(ula_2_K)[2])

#add prefix, easier to refer to
colnames(ula_2_Ca)[2:1093] <- paste('X', colnames(ula_2_Ca[2:1093]), sep="") 
colnames(ula_2_Ca)[1]<-"rows"

colnames(ula_2_Fe)[2:1093] <- paste('X', colnames(ula_2_Fe[2:1093]), sep="") 
colnames(ula_2_Fe)[1]<-"rows"

colnames(ula_2_K)[2:1093] <- paste('X', colnames(ula_2_K[2:1093]), sep="") 
colnames(ula_2_K)[1]<-"rows"

colnames(ula_2_Ti)[2:1093] <- paste('X', colnames(ula_2_Ti[2:1093]), sep="") 
colnames(ula_2_Ti)[1]<-"rows"

#to long format
ula_2_Ca<-gather(ula_2_Ca, cols, val, X2:X1093, factor_key = T)
ula_2_Fe<-gather(ula_2_Fe, cols, val, X2:X1093, factor_key = T)
ula_2_Ti<-gather(ula_2_Ti, cols, val, X2:X1093, factor_key = T)
ula_2_K<-gather(ula_2_K, cols, val, X2:X1093, factor_key = T)


#create (K+Ti)/Ca record
kt<-(ula_2_K$val+ula_2_Ti$val)/ula_2_Ca$val
ula_2_r<-cbind(ula_2_Ca$rows, ula_2_Ca$cols, kt)
colnames(ula_2_r)<-c("rows", "cols", "val")
ula_2_r<-as.data.frame(ula_2_r)

#-----------------------------plot heatmaps----------------------------------------
#Ca
ggplot()+
  geom_tile(data=ula_2_Ca, aes(x=cols, y=rows, fill=val)) + 
  scale_fill_viridis(na.value=NA, option="D")+
  coord_fixed()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#Fe
ggplot()+
  geom_tile(data=ula_2_Fe, aes(x=cols, y=rows, fill=val)) + 
  scale_fill_viridis(na.value=NA, option="D")+
  coord_fixed()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#(K+Ti)/Ca
ggplot()+
  geom_tile(data=ula_2_r, aes(x=cols, y=rows, fill=val)) + 
  scale_fill_viridis(na.value=NA, option="F", limits=c(0, 0.1))+
  coord_fixed()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())






#-----------------------STEINBRUCH SCHMIDT (ULS)-----------------------------------

#------------------------------load data-----------------------------------------
uls_2_Ca<-read.csv("ULS2_Ca_NA.csv")
uls_2_Fe<-read.csv("ULS2_Fe_NA.csv")
uls_2_Ti<-read.csv("ULS2_Ti_NA.csv")
uls_2_K<-read.csv("ULS2_K_NA.csv")

#check (optional)
heatmap(as.matrix(uls_2_Ca), Rowv=NA, Colv=NA, asp=1)

#------------------------------formatting----------------------------------------
#change colnames; keep first column as it is needed later
colnames(uls_2_Ca)<-c(1:dim(uls_2_Ca)[2])
colnames(uls_2_Fe)<-c(1:dim(uls_2_Fe)[2])
colnames(uls_2_Ti)<-c(1:dim(uls_2_Ti)[2])
colnames(uls_2_K)<-c(1:dim(uls_2_K)[2])

#add prefix, easier to refer to
colnames(uls_2_Ca)[2:1161] <- paste('X', colnames(uls_2_Ca[2:1161]), sep="") 
colnames(uls_2_Ca)[1]<-"rows"

colnames(uls_2_Fe)[2:1161] <- paste('X', colnames(uls_2_Fe[2:1161]), sep="") 
colnames(uls_2_Fe)[1]<-"rows"

colnames(uls_2_K)[2:1161] <- paste('X', colnames(uls_2_K[2:1161]), sep="") 
colnames(uls_2_K)[1]<-"rows"

colnames(uls_2_Ti)[2:1161] <- paste('X', colnames(uls_2_Ti[2:1161]), sep="") 
colnames(uls_2_Ti)[1]<-"rows"

#to long format
uls_2_Ca<-gather(uls_2_Ca, cols, val, X2:X1161, factor_key = T)
uls_2_Fe<-gather(uls_2_Fe, cols, val, X2:X1161, factor_key = T)
uls_2_Ti<-gather(uls_2_Ti, cols, val, X2:X1161, factor_key = T)
uls_2_K<-gather(uls_2_K, cols, val, X2:X1161, factor_key = T)


#create (K+Ti)/Ca record
kt<-(uls_2_K$val+uls_2_Ti$val)/uls_2_Ca$val
uls_2_r<-cbind(uls_2_Ca$rows, uls_2_Ca$cols, kt)
colnames(uls_2_r)<-c("rows", "cols", "val")
uls_2_r<-as.data.frame(uls_2_r)

#-----------------------------plot heatmaps----------------------------------------
#Ca
ggplot()+
  geom_tile(data=uls_2_Ca, aes(x=cols, y=rows, fill=val)) + 
  scale_fill_viridis(na.value=NA, option="D")+
  coord_fixed()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#Fe
ggplot()+
  geom_tile(data=uls_2_Fe, aes(x=cols, y=rows, fill=val)) + 
  scale_fill_viridis(na.value=NA, option="D")+
  coord_fixed()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#(K+Ti)/Ca
ggplot()+
  geom_tile(data=uls_2_r, aes(x=cols, y=rows, fill=val)) + 
  scale_fill_viridis(na.value=NA, option="F", limits=c(0, 0.5))+
  coord_fixed()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

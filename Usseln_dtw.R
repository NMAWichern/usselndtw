#This code applies dynamic time warping (dtw) to the Usseln Limestone depth records
#consisting of uXRF-generated elemental data.
#Usseln Limestone specimens: Steinbruch Schmidt (ULS), Winsenberg (ULW), Arfeld (ULA)

setwd("insert wd")

#necessary packages:
library(tidyverse)
library(astrochron)
library(cowplot)
library(dtw)

set.seed(10) #dtw is iterative, therefore set a seed for consistent results.

#------------------------------Loading data:-----------------------------------------
ULA_Ca<-read.csv("ULA_Ca_depth.csv", header=T)
ULA_Fe<-read.csv("ULA_Fe_depth.csv", header=T)
ULA_K<-read.csv("ULA_K_depth.csv", header=T)
ULA_Ti<-read.csv("ULA_Ti_depth.csv", header=T)
ULA_Mn<-read.csv("ULA_Mn_depth.csv", header=T)

ULS_Ca<-read.csv("ULS_Ca_depth.csv", header=T)
ULS_Fe<-read.csv("ULS_Fe_depth.csv", header=T)
ULS_K<-read.csv("ULS_K_depth.csv", header=T)
ULS_Ti<-read.csv("ULS_Ti_depth.csv", header=T)
ULS_Mn<-read.csv("ULS_Mn_depth.csv", header=T)

ULW_Ca<-read.csv("ULw_Ca_depth.csv", header=T)
ULW_Fe<-read.csv("ULw_Fe_depth.csv", header=T)
ULW_K<-read.csv("ULw_K_depth.csv", header=T)
ULW_Ti<-read.csv("ULw_Ti_depth.csv", header=T)
ULW_Mn<-read.csv("ULW_Mn_depth.csv", header=T)


#----------------------------Formatting data------------------------------------------
#remove first column
ULA_Ca<-ULA_Ca[,2:3]
ULA_Fe<-ULA_Fe[,2:3]
ULA_K<-ULA_K[,2:3]
ULA_Ti<-ULA_Ti[,2:3]
ULA_Mn<-ULA_Mn[,2:3]

ULS_Ca<-ULS_Ca[,2:3]
ULS_Fe<-ULS_Fe[,2:3]
ULS_K<-ULS_K[,2:3]
ULS_Ti<-ULS_Ti[,2:3]
ULS_Mn<-ULS_Mn[,2:3]

ULW_Ca<-ULW_Ca[,2:3]
ULW_Fe<-ULW_Fe[,2:3]
ULW_K<-ULW_K[,2:3]
ULW_Ti<-ULW_Ti[,2:3]
ULW_Mn<-ULW_Mn[,2:3]


#Scale to cm (all in mm)
ULA_Ca$depth_mm<-ULA_Ca$depth_mm*0.1
ULA_Fe$depth_mm<-ULA_Fe$depth_mm*0.1
ULA_K$depth_mm<-ULA_K$depth_mm*0.1
ULA_Ti$depth_mm<-ULA_Ti$depth_mm*0.1
ULA_Mn$depth_mm<-ULA_Mn$depth_mm*0.1

colnames(ULA_Ca)<-c("depth_cm", "Ca")
colnames(ULA_Fe)<-c("depth_cm", "Fe")
colnames(ULA_K)<-c("depth_cm", "K")
colnames(ULA_Ti)<-c("depth_cm", "Ti")
colnames(ULA_Mn)<-c("depth_cm", "Mn")

ULS_Ca$depth_mm<-ULS_Ca$depth_mm*0.1
ULS_Fe$depth_mm<-ULS_Fe$depth_mm*0.1
ULS_K$depth_mm<-ULS_K$depth_mm*0.1
ULS_Ti$depth_mm<-ULS_Ti$depth_mm*0.1
ULS_Mn$depth_mm<-ULS_Mn$depth_mm*0.1

colnames(ULS_Ca)<-c("depth_cm", "Ca")
colnames(ULS_Fe)<-c("depth_cm", "Fe")
colnames(ULS_K)<-c("depth_cm", "K")
colnames(ULS_Ti)<-c("depth_cm", "Ti")
colnames(ULS_Mn)<-c("depth_cm", "Mn")

ULW_Ca$depth_mm<-ULW_Ca$depth_mm*0.1
ULW_Fe$depth_mm<-ULW_Fe$depth_mm*0.1
ULW_K$depth_mm<-ULW_K$depth_mm*0.1
ULW_Ti$depth_mm<-ULW_Ti$depth_mm*0.1
ULW_Mn$depth_mm<-ULW_Mn$depth_mm*0.1

colnames(ULW_Ca)<-c("depth_cm", "Ca")
colnames(ULW_Fe)<-c("depth_cm", "Fe")
colnames(ULW_K)<-c("depth_cm", "K")
colnames(ULW_Ti)<-c("depth_cm", "Ti")
colnames(ULW_Mn)<-c("depth_cm", "Mn")


#Make combined records, easier for plotting
ULA<-cbind(ULA_Ca$depth_cm, ULA_Ca$Ca, ULA_Fe$Fe, ULA_K$K, ULA_Ti$Ti, ULA_Mn$Mn)
ULS<-cbind(ULS_Ca$depth_cm, ULS_Ca$Ca, ULS_Fe$Fe, ULS_K$K, ULS_Ti$Ti, ULS_Mn$Mn)
ULW<-cbind(ULW_Ca$depth_cm, ULW_Ca$Ca, ULW_Fe$Fe, ULW_K$K, ULW_Ti$Ti, ULW_Mn$Mn)

colnames(ULA)<-c("depth_cm", "Ca", "Fe", "K", "Ti", "Mn")
colnames(ULS)<-c("depth_cm", "Ca", "Fe", "K", "Ti", "Mn")
colnames(ULW)<-c("depth_cm", "Ca", "Fe", "K", "Ti", "Mn")

ULA<-as.data.frame(ULA)
ULS<-as.data.frame(ULS)
ULW<-as.data.frame(ULW)



#-------------------------Plotting - single elements (optional)-------------------
#Single elements -Arfeld
A_Ca_p<-ggplot()+
  geom_line(data=ULA, aes(x=depth_cm, y=Ca), col="navy")+
  labs(x="depth (cm)", y="Ca (counts)", title="Arfeld Ca")+
  coord_flip()+
  theme_classic()

A_Fe_p<-ggplot()+
  geom_line(data=ULA, aes(x=depth_cm, y=Fe), col="darkred")+
  labs(x="depth (cm)", y="Fe (counts)", title="Arfeld Fe")+
  coord_flip()+
  theme_classic()

A_K_p<-ggplot()+
  geom_line(data=ULA, aes(x=depth_cm, y=K), col="purple")+
  labs(x="depth (cm)", y="K (counts)", title="Arfeld K")+
  coord_flip()+
  theme_classic()

A_Ti_p<-ggplot()+
  geom_line(data=ULA, aes(x=depth_cm, y=Ti), col="darkgreen")+
  labs(x="depth (cm)", y="Ti (counts)", title="Arfeld Ti")+
  coord_flip()+
  theme_classic()

cowplot::plot_grid(A_Ca_p, A_Fe_p, A_K_p, A_Ti_p, ncol=4)



#Single elements - SBS
S_Ca_p<-ggplot()+
  geom_line(data=ULS, aes(x=depth_cm, y=Ca), col="navy")+
  labs(x="depth (cm)", y="Ca (counts)", title="SBS Ca")+
  coord_flip()+
  theme_classic()

S_Fe_p<-ggplot()+
  geom_line(data=ULS, aes(x=depth_cm, y=Fe), col="darkred")+
  labs(x="depth (cm)", y="Fe (counts)", title="SBS Fe")+
  coord_flip()+
  theme_classic()

S_K_p<-ggplot()+
  geom_line(data=ULS, aes(x=depth_cm, y=K), col="purple")+
  labs(x="depth (cm)", y="K (counts)", title="SBS K")+
  coord_flip()+
  theme_classic()

S_Ti_p<-ggplot()+
  geom_line(data=ULS, aes(x=depth_cm, y=Ti), col="darkgreen")+
  labs(x="depth (cm)", y="Ti (counts)", title="SBS Ti")+
  coord_flip()+
  theme_classic()

cowplot::plot_grid(S_Ca_p, S_Fe_p, S_K_p, S_Ti_p, ncol=4)






#Single elements - Winsenberg
W_Ca_p<-ggplot()+
  geom_line(data=ULW, aes(x=depth_cm, y=Ca), col="navy")+
  labs(x="depth (cm)", y="Ca (counts)", title="Winsenberg Ca")+
  coord_flip()+
  theme_classic()

W_Fe_p<-ggplot()+
  geom_line(data=ULW, aes(x=depth_cm, y=Fe), col="darkred")+
  labs(x="depth (cm)", y="Fe (counts)", title="Winsenberg Fe")+
  coord_flip()+
  theme_classic()

W_K_p<-ggplot()+
  geom_line(data=ULW, aes(x=depth_cm, y=K), col="purple")+
  labs(x="depth (cm)", y="K (counts)", title="Winsenberg K")+
  coord_flip()+
  theme_classic()

W_Ti_p<-ggplot()+
  geom_line(data=ULW, aes(x=depth_cm, y=Ti), col="darkgreen")+
  labs(x="depth (cm)", y="Ti (counts)", title="Winsenberg Ti")+
  coord_flip()+
  theme_classic()

cowplot::plot_grid(W_Ca_p, W_Fe_p, W_K_p, W_Ti_p, ncol=4)




#All sites - Mn:
W_Mn_p<-ggplot()+
  geom_line(data=ULW, aes(x=depth_cm, y=Mn), col="darkred")+
  labs(x="depth (cm)", y="Mn (counts)", title="Winsenberg Mn")+
  xlim(0,75)+
  coord_flip()+
  theme_classic()

A_Mn_p<-ggplot()+
  geom_line(data=ULA, aes(x=depth_cm, y=Mn), col="darkred")+
  labs(x="depth (cm)", y="Mn (counts)", title="Arfeld Mn")+
  xlim(0,75)+
  coord_flip()+
  theme_classic()

S_Mn_p<-ggplot()+
  geom_line(data=ULS, aes(x=depth_cm, y=Mn), col="darkred")+
  labs(x="depth (cm)", y="Mn (counts)", title="SBS Mn")+
  xlim(0,75)+
  coord_flip()+
  theme_classic()

#remove spike at end of SBS:
library(astrochron)
ULS_Mn_t<-trim(ULS_Mn, c=10)

S_Mnt_p<-ggplot()+
  geom_line(data=ULS_Mn_t, aes(x=X1, y=X2), col="darkred")+
  labs(x="depth (cm)", y="Mn (counts)", title="SBS Mn - trimmed")+
  xlim(0,75)+
  coord_flip()+
  theme_classic()

cowplot::plot_grid(W_Mn_p, A_Mn_p, S_Mn_p, S_Mnt_p, ncol=4)



#-----------------------Plotting - ratios (optional)--------------------------------
#Ratios -Arfeld
A_TiKCa_p<-ggplot()+
  geom_line(data=ULA, aes(x=depth_cm, y=(Ti+K)/Ca), col="navy")+
  labs(x="depth (cm)", y="(Ti+K/Ca) [-]", title="Arfeld (Ti+K)/Ca")+
  coord_flip()+
  theme_classic()

A_FeCa_p<-ggplot()+
  geom_line(data=ULA, aes(x=depth_cm, y=Fe/Ca), col="darkred")+
  labs(x="depth (cm)", y="Fe/Ca [-]", title="Arfeld Fe/Ca")+
  coord_flip()+
  theme_classic()

A_KCa_p<-ggplot()+
  geom_line(data=ULA, aes(x=depth_cm, y=K/Ca), col="purple")+
  labs(x="depth (cm)", y="K/Ca [-]", title="Arfeld K/Ca")+
  coord_flip()+
  theme_classic()

A_KTi_p<-ggplot()+
  geom_line(data=ULA, aes(x=depth_cm, y=K/Ti), col="darkgreen")+
  labs(x="depth (cm)", y="K/Ti [-]", title="Arfeld K/Ti")+
  coord_flip()+
  theme_classic()

cowplot::plot_grid(A_TiKCa_p, A_FeCa_p, A_KCa_p, A_KTi_p, ncol=4)





#Ratios - SBS
#without log conversion: limited variability, slightly improved by log10 scaling
ggplot()+
  geom_line(data=ULS, aes(x=depth_cm, y=(Ti+K)/log10(Ca)), col="navy")+
  labs(x="depth (cm)", y="(Ti+K/Ca) [-]", title="SBS (Ti+K)/log(Ca)")+
  coord_flip()+
  theme_classic()

S_TiKCa_p<-ggplot()+
  geom_line(data=ULS, aes(x=depth_cm, y=(Ti+K)/log10(Ca)), col="navy")+
  labs(x="depth (cm)", y="(Ti+K/Ca) [-]", title="SBS (Ti+K)/log(Ca)")+
  coord_flip()+
  theme_classic()

S_FeCa_p<-ggplot()+
  geom_line(data=ULS, aes(x=depth_cm, y=Fe/log10(Ca)), col="darkred")+
  labs(x="depth (cm)", y="Fe/Ca [-]", title="SBS Fe/log(Ca)")+
  coord_flip()+
  theme_classic()

S_KCa_p<-ggplot()+
  geom_line(data=ULS, aes(x=depth_cm, y=K/log10(Ca)), col="purple")+
  labs(x="depth (cm)", y="K/Ca [-]", title="SBS K/log(Ca)")+
  coord_flip()+
  theme_classic()

S_KTi_p<-ggplot()+
  geom_line(data=ULS, aes(x=depth_cm, y=K/Ti), col="darkgreen")+
  labs(x="depth (cm)", y="K/Ti [-]", title="SBS K/Ti")+
  coord_flip()+
  theme_classic()

cowplot::plot_grid(S_TiKCa_p, S_FeCa_p, S_KCa_p, S_KTi_p, ncol=4)






#Ratios - Winsenberg
W_TiKCa_p<-ggplot()+
  geom_line(data=ULW, aes(x=depth_cm, y=((Ti+K)/Ca)), col="navy")+
  labs(x="depth (cm)", y="(Ti+K/Ca) [-]", title="Winsenberg (Ti+K)/Ca")+
  coord_flip()+
  theme_classic()

W_FeCa_p<-ggplot()+
  geom_line(data=ULW, aes(x=depth_cm, y=Fe/Ca), col="darkred")+
  labs(x="depth (cm)", y="Fe/Ca [-]", title="Winsenberg Fe/Ca")+
  coord_flip()+
  theme_classic()

W_KCa_p<-ggplot()+
  geom_line(data=ULW, aes(x=depth_cm, y=K/Ca), col="purple")+
  labs(x="depth (cm)", y="K/Ca [-]", title="Winsenberg K/Ca")+
  coord_flip()+
  theme_classic()

W_KTi_p<-ggplot()+
  geom_line(data=ULW, aes(x=depth_cm, y=K/Ti), col="darkgreen")+
  labs(x="depth (cm)", y="K/Ti [-]", title="Winsenberg K/Ti")+
  coord_flip()+
  theme_classic()

cowplot::plot_grid(W_TiKCa_p, W_FeCa_p, W_KCa_p, W_KTi_p, ncol=4)





#For paper:

cowplot::plot_grid(W_Ca_p, W_TiKCa_p, ncol=2)
cowplot::plot_grid(A_Ca_p, A_TiKCa_p, ncol=2)
cowplot::plot_grid(S_Ca_p, S_TiKCa_p, ncol=2)


#--------------Create ratio-records and normalize datasets-----------------------
#Remove NAs from ULW_K record:
ULW_K[is.na(ULW_K)]<-0
#elemental records:
ULA_Ca_n<-cbind(ULA_Ca$depth_cm, ((ULA_Ca$Ca-mean(ULA_Ca$Ca))/sd(ULA_Ca$Ca)))
ULA_Fe_n<-cbind(ULA_Fe$depth_cm, ((ULA_Fe$Fe-mean(ULA_Fe$Fe))/sd(ULA_Fe$Fe)))
ULA_K_n<-cbind(ULA_K$depth_cm, ((ULA_K$K-mean(ULA_K$K))/sd(ULA_K$K)))
ULA_Ti_n<-cbind(ULA_Ti$depth_cm, ((ULA_Ti$Ti-mean(ULA_Ti$Ti))/sd(ULA_Ti$Ti)))

ULW_Ca_n<-cbind(ULW_Ca$depth_cm, ((ULW_Ca$Ca-mean(ULW_Ca$Ca))/sd(ULW_Ca$Ca)))
ULW_Fe_n<-cbind(ULW_Fe$depth_cm, ((ULW_Fe$Fe-mean(ULW_Fe$Fe))/sd(ULW_Fe$Fe)))
ULW_K_n<-cbind(ULW_K$depth_cm, ((ULW_K$K-mean(ULW_K$K))/sd(ULW_K$K)))
ULW_Ti_n<-cbind(ULW_Ti$depth_cm, ((ULW_Ti$Ti-mean(ULW_Ti$Ti))/sd(ULW_Ti$Ti)))

ULS_Ca_n<-cbind(ULS_Ca$depth_cm, ((ULS_Ca$Ca-mean(ULS_Ca$Ca))/sd(ULS_Ca$Ca)))
ULS_Fe_n<-cbind(ULS_Fe$depth_cm, ((ULS_Fe$Fe-mean(ULS_Fe$Fe))/sd(ULS_Fe$Fe)))
ULS_K_n<-cbind(ULS_K$depth_cm, ((ULS_K$K-mean(ULS_K$K))/sd(ULS_K$K)))
ULS_Ti_n<-cbind(ULS_Ti$depth_cm, ((ULS_Ti$Ti-mean(ULS_Ti$Ti))/sd(ULS_Ti$Ti)))


colnames(ULA_Ca_n)<-c("depth_cm", "Ca")
colnames(ULA_Fe_n)<-c("depth_cm", "Fe")
colnames(ULA_K_n)<-c("depth_cm", "K")
colnames(ULA_Ti_n)<-c("depth_cm", "Ti")

colnames(ULW_Ca_n)<-c("depth_cm", "Ca")
colnames(ULW_Fe_n)<-c("depth_cm", "Fe")
colnames(ULW_K_n)<-c("depth_cm", "K")
colnames(ULW_Ti_n)<-c("depth_cm", "Ti")

colnames(ULS_Ca_n)<-c("depth_cm", "Ca")
colnames(ULS_Fe_n)<-c("depth_cm", "Fe")
colnames(ULS_K_n)<-c("depth_cm", "K")
colnames(ULS_Ti_n)<-c("depth_cm", "Ti")



ULA_Ca_n<-as.data.frame(ULA_Ca_n)
ULA_Fe_n<-as.data.frame(ULA_Fe_n)
ULA_K_n<-as.data.frame(ULA_K_n)
ULA_Ti_n<-as.data.frame(ULA_Ti_n)

ULW_Ca_n<-as.data.frame(ULW_Ca_n)
ULW_Fe_n<-as.data.frame(ULW_Fe_n)
ULW_K_n<-as.data.frame(ULW_K_n)
ULW_Ti_n<-as.data.frame(ULW_Ti_n)

ULS_Ca_n<-as.data.frame(ULS_Ca_n)
ULS_Fe_n<-as.data.frame(ULS_Fe_n)
ULS_K_n<-as.data.frame(ULS_K_n)
ULS_Ti_n<-as.data.frame(ULS_Ti_n)



#ratios:
ULA_TKC<-cbind(ULA_Ti$depth_cm, (ULA_Ti$Ti+ULA_K$K)/ULA_Ca$Ca)
ULA_TKC<-as.data.frame(ULA_TKC)
ULA_TKCn<-cbind(ULA_TKC$V1, ((ULA_TKC$V2-mean(ULA_TKC$V2))/sd(ULA_TKC$V2)))
colnames(ULA_TKCn)<-c("depth_cm", "ratio")
ULA_TKCn<-as.data.frame(ULA_TKCn)

ULA_KT<-cbind(ULA_K$depth_cm, ULA_K$K/ULA_Ti$Ti)
ULA_KT<-as.data.frame(ULA_KT)
ULA_KTn<-cbind(ULA_KT$V1, ((ULA_KT$V2-mean(ULA_KT$V2))/sd(ULA_KT$V2)))
colnames(ULA_KTn)<-c("depth_cm", "ratio")
ULA_KTn<-as.data.frame(ULA_KTn)


ULW_TKC<-cbind(ULW_Ti$depth_cm, (ULW_Ti$Ti+ULW_K$K)/ULW_Ca$Ca)
ULW_TKC<-as.data.frame(ULW_TKC)
ULW_TKCn<-cbind(ULW_TKC$V1, ((ULW_TKC$V2-mean(ULW_TKC$V2))/sd(ULW_TKC$V2)))
colnames(ULW_TKCn)<-c("depth_cm", "ratio")
ULW_TKCn<-as.data.frame(ULW_TKCn)

ULW_KT<-cbind(ULW_K$depth_cm, ULW_K$K/ULW_Ti$Ti)
ULW_KT<-as.data.frame(ULW_KT)
ULW_KTn<-cbind(ULW_KT$V1, ((ULW_KT$V2-mean(ULW_KT$V2))/sd(ULW_KT$V2)))
colnames(ULW_KTn)<-c("depth_cm", "ratio")
ULW_KTn<-as.data.frame(ULW_KTn)

ULS_TKC<-cbind(ULS_Ti$depth_cm, (ULS_Ti$Ti+ULS_K$K)/ULS_Ca$Ca)
ULS_TKC<-as.data.frame(ULS_TKC)
ULS_TKCn<-cbind(ULS_TKC$V1, ((ULS_TKC$V2-mean(ULS_TKC$V2))/sd(ULS_TKC$V2)))
colnames(ULS_TKCn)<-c("depth_cm", "ratio")
ULS_TKCn<-as.data.frame(ULS_TKCn)

ULS_KT<-cbind(ULS_K$depth_cm, ULS_K$K/ULS_Ti$Ti)
ULS_KT<-as.data.frame(ULS_KT)
ULS_KTn<-cbind(ULS_KT$V1, ((ULS_KT$V2-mean(ULS_KT$V2))/sd(ULS_KT$V2)))
colnames(ULS_KTn)<-c("depth_cm", "ratio")
ULS_KTn<-as.data.frame(ULS_KTn)


#----------------interpolate to same amount of datapoints------------------------
#Winsenberg (ULW) is the longest, base ULA and ULS off of that
#use astrochron's linterp(), dt=length/datapoints
ULA_Ca_ni<-linterp(ULA_Ca_n, dt=max(ULA_Ca_n$depth_cm)/length(ULW_Ca_n$depth_cm))
ULA_Fe_ni<-linterp(ULA_Fe_n, dt=max(ULA_Fe_n$depth_cm)/length(ULW_Ca_n$depth_cm))
ULA_Ti_ni<-linterp(ULA_Ti_n, dt=max(ULA_Ti_n$depth_cm)/length(ULW_Ca_n$depth_cm))
ULA_K_ni<-linterp(ULA_K_n, dt=max(ULA_K_n$depth_cm)/length(ULW_Ca_n$depth_cm))
ULA_KT_ni<-linterp(ULA_KTn, dt=max(ULA_KTn$depth_cm)/length(ULW_Ca_n$depth_cm))
ULA_TKC_ni<-linterp(ULA_TKCn, dt=max(ULA_TKCn$depth_cm)/length(ULW_Ca_n$depth_cm))

ULS_Ca_ni<-linterp(ULS_Ca_n, dt=max(ULS_Ca_n$depth_cm)/length(ULW_Ca_n$depth_cm))
ULS_Fe_ni<-linterp(ULS_Fe_n, dt=max(ULS_Fe_n$depth_cm)/length(ULW_Ca_n$depth_cm))
ULS_Ti_ni<-linterp(ULS_Ti_n, dt=max(ULS_Ti_n$depth_cm)/length(ULW_Ca_n$depth_cm))
ULS_K_ni<-linterp(ULS_K_n, dt=max(ULS_K_n$depth_cm)/length(ULW_Ca_n$depth_cm))
ULS_KT_ni<-linterp(ULS_KTn, dt=max(ULS_KTn$depth_cm)/length(ULW_Ca_n$depth_cm))
ULS_TKC_ni<-linterp(ULS_TKCn, dt=max(ULS_TKCn$depth_cm)/length(ULW_Ca_n$depth_cm))



#--------------------------------DTW---------------------------------------------

#-------------------------WINSENBERG-ARFELD--------------------------------------
#Winsenberg = W, Arfeld = A

#Ca:
AW_Ca<-dtw(ULA_Ca_ni$Ca, ULW_Ca_n$Ca,open.begin=T, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T, window.type="none")

dtwPlotThreeWay(AW_Ca)
plot(AW_Ca)
dtwPlotTwoWay(AW_Ca, offset=8)

a<-ULA_Ca_ni[AW_Ca$index1,1]
w<-ULW_Ca_n[AW_Ca$index2,1]
aw<-cbind(a,w)
aw<-aw[!duplicated(aw[,2]), ]
aw<-aw[!duplicated(aw[,1]), ]
ULA_Ca_ni_t<-tune(ULA_Ca_ni, aw)

ggplot()+
  geom_line(data=ULW_Ca_n, aes(x=depth_cm, y=Ca), col="black")+
  geom_line(data=ULA_Ca_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

AW_Ca$normalizedDistance
AW_Ca$distance



#Fe:
AW_Fe<-dtw(ULA_Fe_ni$Fe, ULW_Fe_n$Fe, open.begin=T, open.end=T, step.pattern=asymmetricP05,
           window.type="none", keep.internals=T)

dtwPlotThreeWay(AW_Fe)
dtwPlotTwoWay(AW_Fe, offset=8)

a<-ULA_Fe_ni[AW_Fe$index1,1]
w<-ULW_Fe_n[AW_Fe$index2,1]
aw<-cbind(a,w)
aw<-aw[!duplicated(aw[,2]), ]
aw<-aw[!duplicated(aw[,1]), ]
ULA_Fe_ni_t<-tune(ULA_Fe_ni, aw)

ggplot()+
  geom_line(data=ULW_Fe_n, aes(x=depth_cm, y=Fe), col="black")+
  geom_line(data=ULA_Fe_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

AW_Fe$distance
AW_Fe$normalizedDistance




#K:
AW_K<-dtw(ULA_K_ni$K, ULW_K_n$K, open.begin=T, open.end=T, step.pattern=asymmetricP05,
           window.type="none", keep.internals=T)

dtwPlotThreeWay(AW_K)
#dtwPlotTwoWay(AW_K, offset=8)

a<-ULA_K_ni[AW_K$index1,1]
w<-ULW_K_n[AW_K$index2,1]
aw<-cbind(a,w)
aw<-aw[!duplicated(aw[,2]), ]
aw<-aw[!duplicated(aw[,1]), ]
ULA_K_ni_t<-tune(ULA_K_ni, aw)

ggplot()+
  geom_line(data=ULW_K_n, aes(x=depth_cm, y=K), col="black")+
  geom_line(data=ULA_K_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

AW_K$distance
AW_K$normalizedDistance





#Ti:
AW_Ti<-dtw(ULA_Ti_ni$Ti, ULW_Ti_n$Ti, open.begin=T, open.end=T, step.pattern=asymmetricP05,
           window.type="none", keep.internals=T)

dtwPlotThreeWay(AW_Ti)
#dtwPlotTwoWay(AW_Ti, offset=8)

a<-ULA_Ti_ni[AW_Ti$index1,1]
w<-ULW_Ti_n[AW_Ti$index2,1]
aw<-cbind(a,w)
aw<-aw[!duplicated(aw[,2]), ]
aw<-aw[!duplicated(aw[,1]), ]
ULA_Ti_ni_t<-tune(ULA_Ti_ni, aw)

ggplot()+
  geom_line(data=ULW_Ti_n, aes(x=depth_cm, y=Ti), col="black")+
  geom_line(data=ULA_Ti_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

AW_Ti$distance
AW_Ti$normalizedDistance






#(Ti+K)/Ca:
AW_tkc<-dtw(ULA_TKC_ni$ratio, ULW_TKCn$ratio, open.begin=T, open.end=F, step.pattern=asymmetricP05,
            window.type="none", keep.internals=T)

dtwPlotThreeWay(AW_tkc)
#dtwPlotTwoWay(AW_tkc, offset=8)

a<-ULA_TKC_ni[AW_tkc$index1,1]
w<-ULW_TKCn[AW_tkc$index2,1]
aw<-cbind(a,w)
aw<-aw[!duplicated(aw[,2]), ]
aw<-aw[!duplicated(aw[,1]), ]
ULA_tkc_ni_t<-tune(ULA_TKC_ni, aw)

ggplot()+
  geom_line(data=ULW_TKCn, aes(x=depth_cm, y=ratio), col="black")+
  geom_line(data=ULA_tkc_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

AW_tkc$distance
AW_tkc$normalizedDistance




#K/Ti:
AW_KT<-dtw(ULA_KT_ni$ratio, ULW_KTn$ratio, open.begin=T, open.end=T, step.pattern=asymmetricP05,
           keep.internals=T)

dtwPlotThreeWay(AW_KT)
#dtwPlotTwoWay(AW_KT, offset=8)

a<-ULA_KT_ni[AW_KT$index1,1]
w<-ULW_KTn[AW_KT$index2,1]
aw<-cbind(a,w)
aw<-aw[!duplicated(aw[,2]), ]
aw<-aw[!duplicated(aw[,1]), ]
ULA_KT_ni_t<-tune(ULA_KT_ni, aw)

ggplot()+
  geom_line(data=ULW_KTn, aes(x=depth_cm, y=ratio), col="black")+
  geom_line(data=ULA_KT_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

AW_KT$distance
AW_KT$normalizedDistance


#---------------------------alignment plots A-W----------------------------------------

#-------------------------------Open-open-----------------------------------------
AW_Ca<-dtw(ULA_Ca_ni$Ca, ULW_Ca_n$Ca,open.begin=T, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
AW_Fe<-dtw(ULA_Fe_ni$Fe, ULW_Fe_n$Fe,open.begin=T, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
AW_Ti<-dtw(ULA_Ti_ni$Ti, ULW_Ti_n$Ti,open.begin=T, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
AW_K<-dtw(ULA_K_ni$K, ULW_K_n$K,open.begin=T, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
AW_KT<-dtw(ULA_KT_ni$ratio, ULW_KTn$ratio, open.begin=T, open.end=T, 
           step.pattern=asymmetricP05, keep.internals=T)
AW_tkc<-dtw(ULA_TKC_ni$ratio, ULW_TKCn$ratio, open.begin=T, open.end=T,
            step.pattern=asymmetricP05,keep.internals=T)


plot(AW_Ca, col='#fde725',ylab="Winsenberg", xlab="Arfeld", main="Open - Open")
lines(x=AW_Ca$index1, y=AW_Ca$index2, col="#fde725", lw=3)
lines(x=AW_Fe$index1, y=AW_Fe$index2, col="#7ad151", lw=3)
lines(x=AW_Ti$index1, y=AW_Ti$index2, col="#22a884", lw=3)
lines(x=AW_K$index1, y=AW_K$index2, col="#2a788e", lw=3)
lines(x=AW_KT$index1, y=AW_KT$index2, col="#414487", lw=3)
lines(x=AW_tkc$index1, y=AW_tkc$index2, col="#440154", lw=3)
legend("topleft",                    # Add legend to plot
       legend = c("Ca","Fe","Ti","K","K/Ti", "(Ti+K)/Ca"),
       col = c("#fde725","#7ad151","#22a884","#2a788e",
               "#414487","#440154"),pch = 16)


#calculate root mean square of correlations for above scenario:
#1) calculate squared distance between each point for each proxy combination
#2) take mean of all those squared distances
#3) take root of mean squared distance
#do this on the middle 50%, 80% of data to avoid open/closed end effects and to avoid dealing with unequal length
#statistically not ideal to compare Ca with ratio containing Ca but it's an approximation.
#dataset max is 7127, let's take 1000-6000 and 2000-5000

#combinations: 15
  #Ca-Fe, Ca-Ti, Ca-K, Ca-TiKCa, Ca-TiK
  #Fe-Ti, Fe-K, Fe-TiKCa, Fe-TiK
  #Ti-K, Ti-TiKCa, Ti-TiK
  #K-TiKCa, K-TiK
  #TiKCa-TiK

#make dataframes
AW_Ca1<-as.data.frame(cbind(AW_Ca$index1, AW_Ca$index2))
AW_Fe1<-as.data.frame(cbind(AW_Fe$index1, AW_Fe$index2))
AW_Ti1<-as.data.frame(cbind(AW_Ti$index1, AW_Ti$index2))
AW_K1<-as.data.frame(cbind(AW_K$index1, AW_K$index2))
AW_KT1<-as.data.frame(cbind(AW_KT$index1, AW_KT$index2))
AW_tkc1<-as.data.frame(cbind(AW_tkc$index1, AW_tkc$index2))

#remove duplicates in V1; this shouldn't influence RMS too much as there are no hiatuses inferred (ie major jumps)
AW_Ca1<-AW_Ca1[!duplicated(AW_Ca1[,1]), ]
AW_Fe1<-AW_Fe1[!duplicated(AW_Fe1[,1]), ]
AW_Ti1<-AW_Ti1[!duplicated(AW_Ti1[,1]), ]
AW_K1<-AW_K1[!duplicated(AW_K1[,1]), ]
AW_KT1<-AW_KT1[!duplicated(AW_KT1[,1]), ]
AW_tkc1<-AW_tkc1[!duplicated(AW_tkc1[,1]), ]

#now they can be bound together
#in kind of a brute forcing way: structure the input dataframe in such a way 
#that all possible correlations are sequential (because figuring out the looping took too long otherwise)
AW_dat2<-cbind(AW_Ca1$V1, AW_Ca1$V2,AW_Fe1$V2, AW_Ca1$V2,AW_Ti1$V2, AW_Ca1$V2,AW_K1$V2, 
               AW_Ca1$V2,AW_KT1$V2, AW_Ca1$V2,AW_tkc1$V2,
               AW_Fe1$V2,AW_Ti1$V2, AW_Fe1$V2,AW_K1$V2, AW_Fe1$V2,AW_KT1$V2, AW_Fe1$V2,AW_tkc1$V2,
               AW_Ti1$V2,AW_K1$V2, AW_Ti1$V2,AW_KT1$V2, AW_Ti1$V2,AW_tkc1$V2,
               AW_K1$V2,AW_KT1$V2, AW_K1$V2,AW_tkc1$V2,
               AW_KT1$V2,AW_tkc1$V2)

colnames(AW_dat2)<-c("ar", "W_Ca","W_Fe",  "W_Ca","W_Ti", "W_Ca","W_K", "W_Ca","W_KT", "W_Ca","W_tkc",
                     "W_Fe","W_Ti", "W_Fe","W_K", "W_Fe","W_KT", "W_Fe","W_tkc",
                     "W_Ti","W_K", "W_Ti","W_KT", "W_Ti","W_tkc",
                     "W_K","W_KT", "W_K","W_tkc",
                     "W_KT","W_tkc")
view(AW_dat2)

AW_dat2<-as.data.frame(AW_dat2)

AW_dat2a<-subset(AW_dat2, ar>=1000 & ar<=6000)
AW_dat2b<-subset(AW_dat2, ar>=2000 & ar<=5000)





#full data:
m=(ncol(AW_dat2)-1)/2

n=length(AW_dat2$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
   for (i in 1:n){
     
    y[i,j]<-(AW_dat2[i,j*2]-AW_dat2[i,(j*2)+1])^2
    
   }
}

#view(y) #looks good

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #990.7801



#1000-6000 points:
m=(ncol(AW_dat2a)-1)/2

n=length(AW_dat2a$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AW_dat2a[i,j*2]-AW_dat2a[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #1031.336



#2000-5000 points:
m=(ncol(AW_dat2b)-1)/2

n=length(AW_dat2b$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AW_dat2b[i,j*2]-AW_dat2b[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #1063.19










#---------------------------Closed-closed------------------------------------------ 

AW_Ca<-dtw(ULA_Ca_ni$Ca, ULW_Ca_n$Ca,open.begin=F, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
AW_Fe<-dtw(ULA_Fe_ni$Fe, ULW_Fe_n$Fe,open.begin=F, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
AW_Ti<-dtw(ULA_Ti_ni$Ti, ULW_Ti_n$Ti,open.begin=F, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
AW_K<-dtw(ULA_K_ni$K, ULW_K_n$K,open.begin=F, open.end=F,
          step.pattern=asymmetricP05, keep.internals=T)
AW_KT<-dtw(ULA_KT_ni$ratio, ULW_KTn$ratio, open.begin=F, open.end=F, 
           step.pattern=asymmetricP05, keep.internals=T)
AW_tkc<-dtw(ULA_TKC_ni$ratio, ULW_TKCn$ratio, open.begin=F, open.end=F,
            step.pattern=asymmetricP05,keep.internals=T)

plot(AW_Ca, col='#fde725',ylab="Winsenberg", xlab="Arfeld", main="Closed - Closed")
lines(x=AW_Ca$index1, y=AW_Ca$index2, col="#fde725", lw=3)
lines(x=AW_Fe$index1, y=AW_Fe$index2, col="#7ad151", lw=3)
lines(x=AW_Ti$index1, y=AW_Ti$index2, col="#22a884", lw=3)
lines(x=AW_K$index1, y=AW_K$index2, col="#2a788e", lw=3)
lines(x=AW_KT$index1, y=AW_KT$index2, col="#414487", lw=3)
lines(x=AW_tkc$index1, y=AW_tkc$index2, col="#440154", lw=3)
legend("topleft",                    # Add legend to plot
       legend = c("Ca","Fe","Ti","K","K/Ti", "(Ti+K)/Ca"),
       col = c("#fde725","#7ad151","#22a884","#2a788e",
               "#414487","#440154"),pch = 16)

#make dataframes
AW_Ca1<-as.data.frame(cbind(AW_Ca$index1, AW_Ca$index2))
AW_Fe1<-as.data.frame(cbind(AW_Fe$index1, AW_Fe$index2))
AW_Ti1<-as.data.frame(cbind(AW_Ti$index1, AW_Ti$index2))
AW_K1<-as.data.frame(cbind(AW_K$index1, AW_K$index2))
AW_KT1<-as.data.frame(cbind(AW_KT$index1, AW_KT$index2))
AW_tkc1<-as.data.frame(cbind(AW_tkc$index1, AW_tkc$index2))

#remove duplicates in V1; this shouldn't influence RMS too much as there are no hiatuses inferred (ie major jumps)
AW_Ca1<-AW_Ca1[!duplicated(AW_Ca1[,1]), ]
AW_Fe1<-AW_Fe1[!duplicated(AW_Fe1[,1]), ]
AW_Ti1<-AW_Ti1[!duplicated(AW_Ti1[,1]), ]
AW_K1<-AW_K1[!duplicated(AW_K1[,1]), ]
AW_KT1<-AW_KT1[!duplicated(AW_KT1[,1]), ]
AW_tkc1<-AW_tkc1[!duplicated(AW_tkc1[,1]), ]

#now they can be bound together
AW_dat2<-cbind(AW_Ca1$V1, AW_Ca1$V2,AW_Fe1$V2, AW_Ca1$V2,AW_Ti1$V2, AW_Ca1$V2,AW_K1$V2, 
               AW_Ca1$V2,AW_KT1$V2, AW_Ca1$V2,AW_tkc1$V2,
               AW_Fe1$V2,AW_Ti1$V2, AW_Fe1$V2,AW_K1$V2, AW_Fe1$V2,AW_KT1$V2, AW_Fe1$V2,AW_tkc1$V2,
               AW_Ti1$V2,AW_K1$V2, AW_Ti1$V2,AW_KT1$V2, AW_Ti1$V2,AW_tkc1$V2,
               AW_K1$V2,AW_KT1$V2, AW_K1$V2,AW_tkc1$V2,
               AW_KT1$V2,AW_tkc1$V2)

colnames(AW_dat2)<-c("ar", "W_Ca","W_Fe",  "W_Ca","W_Ti", "W_Ca","W_K", "W_Ca","W_KT", "W_Ca","W_tkc",
                     "W_Fe","W_Ti", "W_Fe","W_K", "W_Fe","W_KT", "W_Fe","W_tkc",
                     "W_Ti","W_K", "W_Ti","W_KT", "W_Ti","W_tkc",
                     "W_K","W_KT", "W_K","W_tkc",
                     "W_KT","W_tkc")

AW_dat2<-as.data.frame(AW_dat2)

AW_dat2a<-subset(AW_dat2, ar>=1000 & ar<=6000)
AW_dat2b<-subset(AW_dat2, ar>=2000 & ar<=5000)





#full data:
m=(ncol(AW_dat2)-1)/2

n=length(AW_dat2$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AW_dat2[i,j*2]-AW_dat2[i,(j*2)+1])^2
    
  }
}

view(y)

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #846.4547



#1000-6000 points:
m=(ncol(AW_dat2a)-1)/2

n=length(AW_dat2a$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AW_dat2a[i,j*2]-AW_dat2a[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #975.556



#2000-5000 points:
m=(ncol(AW_dat2b)-1)/2

n=length(AW_dat2b$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AW_dat2b[i,j*2]-AW_dat2b[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #1024.524








#----------------------------Open-closed---------------------------------------------
AW_Ca<-dtw(ULA_Ca_ni$Ca, ULW_Ca_n$Ca,open.begin=T, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
AW_Fe<-dtw(ULA_Fe_ni$Fe, ULW_Fe_n$Fe,open.begin=T, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
AW_Ti<-dtw(ULA_Ti_ni$Ti, ULW_Ti_n$Ti,open.begin=T, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
AW_K<-dtw(ULA_K_ni$K, ULW_K_n$K,open.begin=T, open.end=F,
          step.pattern=asymmetricP05, keep.internals=T)
AW_KT<-dtw(ULA_KT_ni$ratio, ULW_KTn$ratio, open.begin=T, open.end=F, 
           step.pattern=asymmetricP05, keep.internals=T)
AW_tkc<-dtw(ULA_TKC_ni$ratio, ULW_TKCn$ratio, open.begin=T, open.end=F,
            step.pattern=asymmetricP05,keep.internals=T)

plot(AW_Ca, col='#fde725',ylab="Winsenberg", xlab="Arfeld", main="Open - Closed")
lines(x=AW_Ca$index1, y=AW_Ca$index2, col="#fde725", lw=3)
lines(x=AW_Fe$index1, y=AW_Fe$index2, col="#7ad151", lw=3)
lines(x=AW_Ti$index1, y=AW_Ti$index2, col="#22a884", lw=3)
lines(x=AW_K$index1, y=AW_K$index2, col="#2a788e", lw=3)
lines(x=AW_KT$index1, y=AW_KT$index2, col="#414487", lw=3)
lines(x=AW_tkc$index1, y=AW_tkc$index2, col="#440154", lw=3)
legend("topleft",                    # Add legend to plot
       legend = c("Ca","Fe","Ti","K","K/Ti", "(Ti+K)/Ca"),
       col = c("#fde725","#7ad151","#22a884","#2a788e",
               "#414487","#440154"),pch = 16)




#make dataframes
AW_Ca1<-as.data.frame(cbind(AW_Ca$index1, AW_Ca$index2))
AW_Fe1<-as.data.frame(cbind(AW_Fe$index1, AW_Fe$index2))
AW_Ti1<-as.data.frame(cbind(AW_Ti$index1, AW_Ti$index2))
AW_K1<-as.data.frame(cbind(AW_K$index1, AW_K$index2))
AW_KT1<-as.data.frame(cbind(AW_KT$index1, AW_KT$index2))
AW_tkc1<-as.data.frame(cbind(AW_tkc$index1, AW_tkc$index2))

#remove duplicates in V1; this shouldn't influence RMS too much as there are no hiatuses inferred (ie major jumps)
AW_Ca1<-AW_Ca1[!duplicated(AW_Ca1[,1]), ]
AW_Fe1<-AW_Fe1[!duplicated(AW_Fe1[,1]), ]
AW_Ti1<-AW_Ti1[!duplicated(AW_Ti1[,1]), ]
AW_K1<-AW_K1[!duplicated(AW_K1[,1]), ]
AW_KT1<-AW_KT1[!duplicated(AW_KT1[,1]), ]
AW_tkc1<-AW_tkc1[!duplicated(AW_tkc1[,1]), ]

#now they can be bound together
AW_dat2<-cbind(AW_Ca1$V1, AW_Ca1$V2,AW_Fe1$V2, AW_Ca1$V2,AW_Ti1$V2, AW_Ca1$V2,AW_K1$V2, 
               AW_Ca1$V2,AW_KT1$V2, AW_Ca1$V2,AW_tkc1$V2,
               AW_Fe1$V2,AW_Ti1$V2, AW_Fe1$V2,AW_K1$V2, AW_Fe1$V2,AW_KT1$V2, AW_Fe1$V2,AW_tkc1$V2,
               AW_Ti1$V2,AW_K1$V2, AW_Ti1$V2,AW_KT1$V2, AW_Ti1$V2,AW_tkc1$V2,
               AW_K1$V2,AW_KT1$V2, AW_K1$V2,AW_tkc1$V2,
               AW_KT1$V2,AW_tkc1$V2)

colnames(AW_dat2)<-c("ar", "W_Ca","W_Fe",  "W_Ca","W_Ti", "W_Ca","W_K", "W_Ca","W_KT", "W_Ca","W_tkc",
                     "W_Fe","W_Ti", "W_Fe","W_K", "W_Fe","W_KT", "W_Fe","W_tkc",
                     "W_Ti","W_K", "W_Ti","W_KT", "W_Ti","W_tkc",
                     "W_K","W_KT", "W_K","W_tkc",
                     "W_KT","W_tkc")

AW_dat2<-as.data.frame(AW_dat2)

AW_dat2a<-subset(AW_dat2, ar>=1000 & ar<=6000)
AW_dat2b<-subset(AW_dat2, ar>=2000 & ar<=5000)





#full data:
m=(ncol(AW_dat2)-1)/2

n=length(AW_dat2$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AW_dat2[i,j*2]-AW_dat2[i,(j*2)+1])^2
    
  }
}


#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #964.2763



#1000-6000 points:
m=(ncol(AW_dat2a)-1)/2

n=length(AW_dat2a$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AW_dat2a[i,j*2]-AW_dat2a[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #1016.474



#2000-5000 points:
m=(ncol(AW_dat2b)-1)/2

n=length(AW_dat2b$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AW_dat2b[i,j*2]-AW_dat2b[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #1043.608








#--------------------------------Closed-open--------------------------------------
AW_Ca<-dtw(ULA_Ca_ni$Ca, ULW_Ca_n$Ca,open.begin=F, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
AW_Fe<-dtw(ULA_Fe_ni$Fe, ULW_Fe_n$Fe,open.begin=F, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
AW_Ti<-dtw(ULA_Ti_ni$Ti, ULW_Ti_n$Ti,open.begin=F, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
AW_K<-dtw(ULA_K_ni$K, ULW_K_n$K,open.begin=F, open.end=T,
          step.pattern=asymmetricP05, keep.internals=T)
AW_KT<-dtw(ULA_KT_ni$ratio, ULW_KTn$ratio, open.begin=F, open.end=T, 
           step.pattern=asymmetricP05, keep.internals=T)
AW_tkc<-dtw(ULA_TKC_ni$ratio, ULW_TKCn$ratio, open.begin=F, open.end=T,
            step.pattern=asymmetricP05,keep.internals=T)

plot(AW_Ca, col='#fde725',ylab="Winsenberg", xlab="Arfeld", main="Closed - Open")
lines(x=AW_Ca$index1, y=AW_Ca$index2, col="#fde725", lw=3)
lines(x=AW_Fe$index1, y=AW_Fe$index2, col="#7ad151", lw=3)
lines(x=AW_Ti$index1, y=AW_Ti$index2, col="#22a884", lw=3)
lines(x=AW_K$index1, y=AW_K$index2, col="#2a788e", lw=3)
lines(x=AW_KT$index1, y=AW_KT$index2, col="#414487", lw=3)
lines(x=AW_tkc$index1, y=AW_tkc$index2, col="#440154", lw=3)
legend("topleft",                    # Add legend to plot
       legend = c("Ca","Fe","Ti","K","K/Ti", "(Ti+K)/Ca"),
       col = c("#fde725","#7ad151","#22a884","#2a788e",
               "#414487","#440154"),pch = 16)



#make dataframes
AW_Ca1<-as.data.frame(cbind(AW_Ca$index1, AW_Ca$index2))
AW_Fe1<-as.data.frame(cbind(AW_Fe$index1, AW_Fe$index2))
AW_Ti1<-as.data.frame(cbind(AW_Ti$index1, AW_Ti$index2))
AW_K1<-as.data.frame(cbind(AW_K$index1, AW_K$index2))
AW_KT1<-as.data.frame(cbind(AW_KT$index1, AW_KT$index2))
AW_tkc1<-as.data.frame(cbind(AW_tkc$index1, AW_tkc$index2))

#remove duplicates in V1; this shouldn't influence RMS too much as there are no hiatuses inferred (ie major jumps)
AW_Ca1<-AW_Ca1[!duplicated(AW_Ca1[,1]), ]
AW_Fe1<-AW_Fe1[!duplicated(AW_Fe1[,1]), ]
AW_Ti1<-AW_Ti1[!duplicated(AW_Ti1[,1]), ]
AW_K1<-AW_K1[!duplicated(AW_K1[,1]), ]
AW_KT1<-AW_KT1[!duplicated(AW_KT1[,1]), ]
AW_tkc1<-AW_tkc1[!duplicated(AW_tkc1[,1]), ]

#now they can be bound together
AW_dat2<-cbind(AW_Ca1$V1, AW_Ca1$V2,AW_Fe1$V2, AW_Ca1$V2,AW_Ti1$V2, AW_Ca1$V2,AW_K1$V2, 
               AW_Ca1$V2,AW_KT1$V2, AW_Ca1$V2,AW_tkc1$V2,
               AW_Fe1$V2,AW_Ti1$V2, AW_Fe1$V2,AW_K1$V2, AW_Fe1$V2,AW_KT1$V2, AW_Fe1$V2,AW_tkc1$V2,
               AW_Ti1$V2,AW_K1$V2, AW_Ti1$V2,AW_KT1$V2, AW_Ti1$V2,AW_tkc1$V2,
               AW_K1$V2,AW_KT1$V2, AW_K1$V2,AW_tkc1$V2,
               AW_KT1$V2,AW_tkc1$V2)

colnames(AW_dat2)<-c("ar", "W_Ca","W_Fe",  "W_Ca","W_Ti", "W_Ca","W_K", "W_Ca","W_KT", "W_Ca","W_tkc",
                     "W_Fe","W_Ti", "W_Fe","W_K", "W_Fe","W_KT", "W_Fe","W_tkc",
                     "W_Ti","W_K", "W_Ti","W_KT", "W_Ti","W_tkc",
                     "W_K","W_KT", "W_K","W_tkc",
                     "W_KT","W_tkc")

AW_dat2<-as.data.frame(AW_dat2)

AW_dat2a<-subset(AW_dat2, ar>=1000 & ar<=6000)
AW_dat2b<-subset(AW_dat2, ar>=2000 & ar<=5000)





#full data:
m=(ncol(AW_dat2)-1)/2

n=length(AW_dat2$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AW_dat2[i,j*2]-AW_dat2[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #1115.367



#1000-6000 points:
m=(ncol(AW_dat2a)-1)/2

n=length(AW_dat2a$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AW_dat2a[i,j*2]-AW_dat2a[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #1211.944



#2000-5000 points:
m=(ncol(AW_dat2b)-1)/2

n=length(AW_dat2b$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AW_dat2b[i,j*2]-AW_dat2b[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #1260.359










#-------------------------WINSENBERG-STEINBRUCH SCHMIDT---------------------------
#Winsenberg = W, Steinbruch Schmidt = S

#Ca:
SW_Ca<-dtw(ULS_Ca_ni$Ca, ULW_Ca_n$Ca, open.begin=F, open.end=F, step.pattern=asymmetricP05,
             window.type="none",  keep.internals=T)

dtwPlotThreeWay(SW_Ca)
dtwPlotTwoWay(SW_Ca, offset=8)

s<-ULS_Ca_ni[SW_Ca$index1,1]
w<-ULW_Ca_n[SW_Ca$index2,1]
SW<-cbind(s,w)
SW<-SW[!duplicated(SW[,2]), ]
SW<-SW[!duplicated(SW[,1]), ]
ULS_Ca_ni_t<-tune(ULS_Ca_ni, SW)


ggplot()+
  geom_line(data=ULW_Ca_n, aes(x=depth_cm, y=Ca), col="black")+
  geom_line(data=ULS_Ca_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

SW_Ca$normalizedDistance




#Fe:
SW_Fe<-dtw(ULS_Fe_ni$Fe, ULW_Fe_n$Fe, open.begin=F, open.end=F, step.pattern=asymmetricP05,
           window.type="none",  keep.internals=T)

dtwPlotThreeWay(SW_Fe)
dtwPlotTwoWay(SW_Fe, offset=8)

s<-ULS_Fe_ni[SW_Fe$index1,1]
w<-ULW_Fe_n[SW_Fe$index2,1]
SW<-cbind(s,w)
SW<-SW[!duplicated(SW[,2]), ]
SW<-SW[!duplicated(SW[,1]), ]
ULS_Fe_ni_t<-tune(ULS_Fe_ni, SW)


ggplot()+
  geom_line(data=ULW_Fe_n, aes(x=depth_cm, y=Fe), col="black")+
  geom_line(data=ULS_Fe_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

SW_Fe$normalizedDistance
SW_Fe$distance



#K:
SW_K<-dtw(ULS_K_ni$K, ULW_K_n$K, open.begin=F, open.end=F, step.pattern=asymmetricP05,
           window.type="none",  keep.internals=T)

dtwPlotThreeWay(SW_K)
dtwPlotTwoWay(SW_K, offset=8)

s<-ULS_K_ni[SW_K$index1,1]
w<-ULW_K_n[SW_K$index2,1]
SW<-cbind(s,w)
SW<-SW[!duplicated(SW[,2]), ]
SW<-SW[!duplicated(SW[,1]), ]
ULS_K_ni_t<-tune(ULS_K_ni, SW)


ggplot()+
  geom_line(data=ULW_K_n, aes(x=depth_cm, y=K), col="black")+
  geom_line(data=ULS_K_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

SW_K$normalizedDistance
SW_K$distance




#Ti:
SW_Ti<-dtw(ULS_Ti_ni$Ti, ULW_Ti_n$Ti, open.begin=F, open.end=F, step.pattern=asymmetricP05,
           window.type="none",  keep.internals=T)

dtwPlotThreeWay(SW_Ti)
dtwPlotTwoWay(SW_Ti, offset=8)

s<-ULS_Ti_ni[SW_Ti$index1,1]
w<-ULW_Ti_n[SW_Ti$index2,1]
SW<-cbind(s,w)
SW<-SW[!duplicated(SW[,2]), ]
SW<-SW[!duplicated(SW[,1]), ]
ULS_Ti_ni_t<-tune(ULS_Ti_ni, SW)


ggplot()+
  geom_line(data=ULW_Ti_n, aes(x=depth_cm, y=Ti), col="black")+
  geom_line(data=ULS_Ti_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

SW_Ti$normalizedDistance
SW_Ti$distance





#(Ti+K)/Ca:
SW_tkc<-dtw(ULS_TKC_ni$ratio, ULW_TKCn$ratio, open.begin=F, open.end=F, step.pattern=asymmetricP05,
            window.type="none", keep.internals=T)

dtwPlotThreeWay(SW_tkc)
#dtwPlotTwoWay(SW_tkc, offset=8)

s<-ULS_TKC_ni[SW_tkc$index1,1]
w<-ULW_TKCn[SW_tkc$index2,1]
SW<-cbind(s,w)
SW<-SW[!duplicated(SW[,2]), ]
SW<-SW[!duplicated(SW[,1]), ]
ULS_tkc_ni_t<-tune(ULS_TKC_ni, SW)

ggplot()+
  geom_line(data=ULW_TKCn, aes(x=depth_cm, y=ratio), col="black")+
  geom_line(data=ULS_tkc_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

#SW_tkc$distance
SW_tkc$normalizedDistance





#K/Ti:
SW_KT<-dtw(ULS_KT_ni$ratio, ULW_KTn$ratio, open.begin=F, open.end=F, step.pattern=asymmetricP05,
           keep.internals=T)

dtwPlotThreeWay(SW_KT)
#dtwPlotTwoWay(SW_KT, offset=8)

s<-ULS_KT_ni[SW_KT$index1,1]
w<-ULW_KTn[SW_KT$index2,1]
SW<-cbind(s,w)
SW<-SW[!duplicated(SW[,2]), ]
SW<-SW[!duplicated(SW[,1]), ]
ULS_KT_ni_t<-tune(ULS_KT_ni, SW)

ggplot()+
  geom_line(data=ULW_KTn, aes(x=depth_cm, y=ratio), col="black")+
  geom_line(data=ULS_KT_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

#SW_KT$distance
SW_KT$normalizedDistance







#---------------------------alignment plots S-W----------------------------------------
#-----------------------------Open-open-------------------------------------------
SW_Ca<-dtw(ULS_Ca_ni$Ca, ULW_Ca_n$Ca,open.begin=T, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
SW_Fe<-dtw(ULS_Fe_ni$Fe, ULW_Fe_n$Fe,open.begin=T, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
SW_Ti<-dtw(ULS_Ti_ni$Ti, ULW_Ti_n$Ti,open.begin=T, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
SW_K<-dtw(ULS_K_ni$K, ULW_K_n$K,open.begin=T, open.end=T,
          step.pattern=asymmetricP05, keep.internals=T)
SW_KT<-dtw(ULS_KT_ni$ratio, ULW_KTn$ratio, open.begin=T, open.end=T, 
           step.pattern=asymmetricP05, keep.internals=T)
SW_tkc<-dtw(ULS_TKC_ni$ratio, ULW_TKCn$ratio, open.begin=T, open.end=T,
            step.pattern=asymmetricP05,keep.internals=T)

plot(SW_Ca, col='#fde725',ylab="Winsenberg", xlab="Steinbruch Schmidt", main="Open - Open")
lines(x=SW_Ca$index1, y=SW_Ca$index2, col="#fde725", lw=3)
lines(x=SW_Fe$index1, y=SW_Fe$index2, col="#7ad151", lw=3)
lines(x=SW_Ti$index1, y=SW_Ti$index2, col="#22a884", lw=3)
lines(x=SW_K$index1, y=SW_K$index2, col="#2a788e", lw=3)
lines(x=SW_KT$index1, y=SW_KT$index2, col="#414487", lw=3)
lines(x=SW_tkc$index1, y=SW_tkc$index2, col="#440154", lw=3)
legend("topleft",                    # Add legend to plot
       legend = c("Ca","Fe","Ti","K","K/Ti", "(Ti+K)/Ca"),
       col = c("#fde725","#7ad151","#22a884","#2a788e",
               "#414487","#440154"),pch = 16)



#make dataframes
SW_Ca1<-as.data.frame(cbind(SW_Ca$index1, SW_Ca$index2))
SW_Fe1<-as.data.frame(cbind(SW_Fe$index1, SW_Fe$index2))
SW_Ti1<-as.data.frame(cbind(SW_Ti$index1, SW_Ti$index2))
SW_K1<-as.data.frame(cbind(SW_K$index1, SW_K$index2))
SW_KT1<-as.data.frame(cbind(SW_KT$index1, SW_KT$index2))
SW_tkc1<-as.data.frame(cbind(SW_tkc$index1, SW_tkc$index2))

#remove duplicates in V1; this shouldn't influence RMS too much as there are no hiatuses inferred (ie major jumps)
SW_Ca1<-SW_Ca1[!duplicated(SW_Ca1[,1]), ]
SW_Fe1<-SW_Fe1[!duplicated(SW_Fe1[,1]), ]
SW_Ti1<-SW_Ti1[!duplicated(SW_Ti1[,1]), ]
SW_K1<-SW_K1[!duplicated(SW_K1[,1]), ]
SW_KT1<-SW_KT1[!duplicated(SW_KT1[,1]), ]
SW_tkc1<-SW_tkc1[!duplicated(SW_tkc1[,1]), ]

#now they can be bound together
SW_dat2<-cbind(SW_Ca1$V1, SW_Ca1$V2,SW_Fe1$V2, SW_Ca1$V2,SW_Ti1$V2, SW_Ca1$V2,SW_K1$V2, 
               SW_Ca1$V2,SW_KT1$V2, SW_Ca1$V2,SW_tkc1$V2,
               SW_Fe1$V2,SW_Ti1$V2, SW_Fe1$V2,SW_K1$V2, SW_Fe1$V2,SW_KT1$V2, SW_Fe1$V2,SW_tkc1$V2,
               SW_Ti1$V2,SW_K1$V2, SW_Ti1$V2,SW_KT1$V2, SW_Ti1$V2,SW_tkc1$V2,
               SW_K1$V2,SW_KT1$V2, SW_K1$V2,SW_tkc1$V2,
               SW_KT1$V2,SW_tkc1$V2)

colnames(SW_dat2)<-c("ar", "W_Ca","W_Fe",  "W_Ca","W_Ti", "W_Ca","W_K", "W_Ca","W_KT", "W_Ca","W_tkc",
                     "W_Fe","W_Ti", "W_Fe","W_K", "W_Fe","W_KT", "W_Fe","W_tkc",
                     "W_Ti","W_K", "W_Ti","W_KT", "W_Ti","W_tkc",
                     "W_K","W_KT", "W_K","W_tkc",
                     "W_KT","W_tkc")

SW_dat2<-as.data.frame(SW_dat2)

SW_dat2a<-subset(SW_dat2, ar>=1000 & ar<=6000)
SW_dat2b<-subset(SW_dat2, ar>=2000 & ar<=5000)





#full data:
m=(ncol(SW_dat2)-1)/2

n=length(SW_dat2$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SW_dat2[i,j*2]-SW_dat2[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #469.3521



#1000-6000 points:
m=(ncol(SW_dat2a)-1)/2

n=length(SW_dat2a$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SW_dat2a[i,j*2]-SW_dat2a[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #440.6746



#2000-5000 points:
m=(ncol(SW_dat2b)-1)/2

n=length(SW_dat2b$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SW_dat2b[i,j*2]-SW_dat2b[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #421.7167




#------------------------------closed-closed-------------------------------------
SW_Ca<-dtw(ULS_Ca_ni$Ca, ULW_Ca_n$Ca,open.begin=F, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
SW_Fe<-dtw(ULS_Fe_ni$Fe, ULW_Fe_n$Fe,open.begin=F, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
SW_Ti<-dtw(ULS_Ti_ni$Ti, ULW_Ti_n$Ti,open.begin=F, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
SW_K<-dtw(ULS_K_ni$K, ULW_K_n$K,open.begin=F, open.end=F,
          step.pattern=asymmetricP05, keep.internals=T)
SW_KT<-dtw(ULS_KT_ni$ratio, ULW_KTn$ratio, open.begin=F, open.end=F, 
           step.pattern=asymmetricP05, keep.internals=T)
SW_tkc<-dtw(ULS_TKC_ni$ratio, ULW_TKCn$ratio, open.begin=F, open.end=F,
            step.pattern=asymmetricP05,keep.internals=T)

plot(SW_Ca, col='#fde725',ylab="Winsenberg", xlab="Steinbruch Schmidt", main="Closed - Closed")
lines(x=SW_Ca$index1, y=SW_Ca$index2, col="#fde725", lw=3)
lines(x=SW_Fe$index1, y=SW_Fe$index2, col="#7ad151", lw=3)
lines(x=SW_Ti$index1, y=SW_Ti$index2, col="#22a884", lw=3)
lines(x=SW_K$index1, y=SW_K$index2, col="#2a788e", lw=3)
lines(x=SW_KT$index1, y=SW_KT$index2, col="#414487", lw=3)
lines(x=SW_tkc$index1, y=SW_tkc$index2, col="#440154", lw=3)
legend("topleft",                    # Add legend to plot
       legend = c("Ca","Fe","Ti","K","K/Ti", "(Ti+K)/Ca"),
       col = c("#fde725","#7ad151","#22a884","#2a788e",
               "#414487","#440154"),pch = 16)





#make dataframes
SW_Ca1<-as.data.frame(cbind(SW_Ca$index1, SW_Ca$index2))
SW_Fe1<-as.data.frame(cbind(SW_Fe$index1, SW_Fe$index2))
SW_Ti1<-as.data.frame(cbind(SW_Ti$index1, SW_Ti$index2))
SW_K1<-as.data.frame(cbind(SW_K$index1, SW_K$index2))
SW_KT1<-as.data.frame(cbind(SW_KT$index1, SW_KT$index2))
SW_tkc1<-as.data.frame(cbind(SW_tkc$index1, SW_tkc$index2))

#remove duplicates in V1; this shouldn't influence RMS too much as there are no hiatuses inferred (ie major jumps)
SW_Ca1<-SW_Ca1[!duplicated(SW_Ca1[,1]), ]
SW_Fe1<-SW_Fe1[!duplicated(SW_Fe1[,1]), ]
SW_Ti1<-SW_Ti1[!duplicated(SW_Ti1[,1]), ]
SW_K1<-SW_K1[!duplicated(SW_K1[,1]), ]
SW_KT1<-SW_KT1[!duplicated(SW_KT1[,1]), ]
SW_tkc1<-SW_tkc1[!duplicated(SW_tkc1[,1]), ]

#now they can be bound together
SW_dat2<-cbind(SW_Ca1$V1, SW_Ca1$V2,SW_Fe1$V2, SW_Ca1$V2,SW_Ti1$V2, SW_Ca1$V2,SW_K1$V2, 
               SW_Ca1$V2,SW_KT1$V2, SW_Ca1$V2,SW_tkc1$V2,
               SW_Fe1$V2,SW_Ti1$V2, SW_Fe1$V2,SW_K1$V2, SW_Fe1$V2,SW_KT1$V2, SW_Fe1$V2,SW_tkc1$V2,
               SW_Ti1$V2,SW_K1$V2, SW_Ti1$V2,SW_KT1$V2, SW_Ti1$V2,SW_tkc1$V2,
               SW_K1$V2,SW_KT1$V2, SW_K1$V2,SW_tkc1$V2,
               SW_KT1$V2,SW_tkc1$V2)

colnames(SW_dat2)<-c("ar", "W_Ca","W_Fe",  "W_Ca","W_Ti", "W_Ca","W_K", "W_Ca","W_KT", "W_Ca","W_tkc",
                     "W_Fe","W_Ti", "W_Fe","W_K", "W_Fe","W_KT", "W_Fe","W_tkc",
                     "W_Ti","W_K", "W_Ti","W_KT", "W_Ti","W_tkc",
                     "W_K","W_KT", "W_K","W_tkc",
                     "W_KT","W_tkc")

SW_dat2<-as.data.frame(SW_dat2)

SW_dat2a<-subset(SW_dat2, ar>=1000 & ar<=6000)
SW_dat2b<-subset(SW_dat2, ar>=2000 & ar<=5000)





#full data:
m=(ncol(SW_dat2)-1)/2

n=length(SW_dat2$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SW_dat2[i,j*2]-SW_dat2[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #592.3632



#1000-6000 points:
m=(ncol(SW_dat2a)-1)/2

n=length(SW_dat2a$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SW_dat2a[i,j*2]-SW_dat2a[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #672.2157



#2000-5000 points:
m=(ncol(SW_dat2b)-1)/2

n=length(SW_dat2b$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SW_dat2b[i,j*2]-SW_dat2b[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #673.0173




#------------------------------Open-closed----------------------------------------
SW_Ca<-dtw(ULS_Ca_ni$Ca, ULW_Ca_n$Ca,open.begin=T, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
SW_Fe<-dtw(ULS_Fe_ni$Fe, ULW_Fe_n$Fe,open.begin=T, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
SW_Ti<-dtw(ULS_Ti_ni$Ti, ULW_Ti_n$Ti,open.begin=T, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
SW_K<-dtw(ULS_K_ni$K, ULW_K_n$K,open.begin=T, open.end=F,
          step.pattern=asymmetricP05, keep.internals=T)
SW_KT<-dtw(ULS_KT_ni$ratio, ULW_KTn$ratio, open.begin=T, open.end=F, 
           step.pattern=asymmetricP05, keep.internals=T)
SW_tkc<-dtw(ULS_TKC_ni$ratio, ULW_TKCn$ratio, open.begin=T, open.end=F,
            step.pattern=asymmetricP05,keep.internals=T)

plot(SW_Ca, col='#fde725',ylab="Winsenberg", xlab="Steinbruch Schmidt", main="Open - Closed")
lines(x=SW_Ca$index1, y=SW_Ca$index2, col="#fde725", lw=3)
lines(x=SW_Fe$index1, y=SW_Fe$index2, col="#7ad151", lw=3)
lines(x=SW_Ti$index1, y=SW_Ti$index2, col="#22a884", lw=3)
lines(x=SW_K$index1, y=SW_K$index2, col="#2a788e", lw=3)
lines(x=SW_KT$index1, y=SW_KT$index2, col="#414487", lw=3)
lines(x=SW_tkc$index1, y=SW_tkc$index2, col="#440154", lw=3)
legend("topleft",                    # Add legend to plot
       legend = c("Ca","Fe","Ti","K","K/Ti", "(Ti+K)/Ca"),
       col = c("#fde725","#7ad151","#22a884","#2a788e",
               "#414487","#440154"),pch = 16)





#make dataframes
SW_Ca1<-as.data.frame(cbind(SW_Ca$index1, SW_Ca$index2))
SW_Fe1<-as.data.frame(cbind(SW_Fe$index1, SW_Fe$index2))
SW_Ti1<-as.data.frame(cbind(SW_Ti$index1, SW_Ti$index2))
SW_K1<-as.data.frame(cbind(SW_K$index1, SW_K$index2))
SW_KT1<-as.data.frame(cbind(SW_KT$index1, SW_KT$index2))
SW_tkc1<-as.data.frame(cbind(SW_tkc$index1, SW_tkc$index2))

#remove duplicates in V1; this shouldn't influence RMS too much as there are no hiatuses inferred (ie major jumps)
SW_Ca1<-SW_Ca1[!duplicated(SW_Ca1[,1]), ]
SW_Fe1<-SW_Fe1[!duplicated(SW_Fe1[,1]), ]
SW_Ti1<-SW_Ti1[!duplicated(SW_Ti1[,1]), ]
SW_K1<-SW_K1[!duplicated(SW_K1[,1]), ]
SW_KT1<-SW_KT1[!duplicated(SW_KT1[,1]), ]
SW_tkc1<-SW_tkc1[!duplicated(SW_tkc1[,1]), ]

#now they can be bound together
SW_dat2<-cbind(SW_Ca1$V1, SW_Ca1$V2,SW_Fe1$V2, SW_Ca1$V2,SW_Ti1$V2, SW_Ca1$V2,SW_K1$V2, 
               SW_Ca1$V2,SW_KT1$V2, SW_Ca1$V2,SW_tkc1$V2,
               SW_Fe1$V2,SW_Ti1$V2, SW_Fe1$V2,SW_K1$V2, SW_Fe1$V2,SW_KT1$V2, SW_Fe1$V2,SW_tkc1$V2,
               SW_Ti1$V2,SW_K1$V2, SW_Ti1$V2,SW_KT1$V2, SW_Ti1$V2,SW_tkc1$V2,
               SW_K1$V2,SW_KT1$V2, SW_K1$V2,SW_tkc1$V2,
               SW_KT1$V2,SW_tkc1$V2)

colnames(SW_dat2)<-c("ar", "W_Ca","W_Fe",  "W_Ca","W_Ti", "W_Ca","W_K", "W_Ca","W_KT", "W_Ca","W_tkc",
                     "W_Fe","W_Ti", "W_Fe","W_K", "W_Fe","W_KT", "W_Fe","W_tkc",
                     "W_Ti","W_K", "W_Ti","W_KT", "W_Ti","W_tkc",
                     "W_K","W_KT", "W_K","W_tkc",
                     "W_KT","W_tkc")

SW_dat2<-as.data.frame(SW_dat2)

SW_dat2a<-subset(SW_dat2, ar>=1000 & ar<=6000)
SW_dat2b<-subset(SW_dat2, ar>=2000 & ar<=5000)





#full data:
m=(ncol(SW_dat2)-1)/2

n=length(SW_dat2$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SW_dat2[i,j*2]-SW_dat2[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #393.5645



#1000-6000 points:
m=(ncol(SW_dat2a)-1)/2

n=length(SW_dat2a$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SW_dat2a[i,j*2]-SW_dat2a[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #404.7166



#2000-5000 points:
m=(ncol(SW_dat2b)-1)/2

n=length(SW_dat2b$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SW_dat2b[i,j*2]-SW_dat2b[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #412.6714



#-------------------------------Closed-open----------------------------------------
SW_Ca<-dtw(ULS_Ca_ni$Ca, ULW_Ca_n$Ca,open.begin=F, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
SW_Fe<-dtw(ULS_Fe_ni$Fe, ULW_Fe_n$Fe,open.begin=F, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
SW_Ti<-dtw(ULS_Ti_ni$Ti, ULW_Ti_n$Ti,open.begin=F, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
SW_K<-dtw(ULS_K_ni$K, ULW_K_n$K,open.begin=F, open.end=T,
          step.pattern=asymmetricP05, keep.internals=T)
SW_KT<-dtw(ULS_KT_ni$ratio, ULW_KTn$ratio, open.begin=F, open.end=T, 
           step.pattern=asymmetricP05, keep.internals=T)
SW_tkc<-dtw(ULS_TKC_ni$ratio, ULW_TKCn$ratio, open.begin=F, open.end=T,
            step.pattern=asymmetricP05,keep.internals=T)

plot(SW_Ca, col='#fde725',ylab="Winsenberg", xlab="Steinbruch Schmidt", main="Closed - Open")
lines(x=SW_Ca$index1, y=SW_Ca$index2, col="#fde725", lw=3)
lines(x=SW_Fe$index1, y=SW_Fe$index2, col="#7ad151", lw=3)
lines(x=SW_Ti$index1, y=SW_Ti$index2, col="#22a884", lw=3)
lines(x=SW_K$index1, y=SW_K$index2, col="#2a788e", lw=3)
lines(x=SW_KT$index1, y=SW_KT$index2, col="#414487", lw=3)
lines(x=SW_tkc$index1, y=SW_tkc$index2, col="#440154", lw=3)
legend("topleft",                    # Add legend to plot
       legend = c("Ca","Fe","Ti","K","K/Ti", "(Ti+K)/Ca"),
       col = c("#fde725","#7ad151","#22a884","#2a788e",
               "#414487","#440154"),pch = 16)



#make dataframes
SW_Ca1<-as.data.frame(cbind(SW_Ca$index1, SW_Ca$index2))
SW_Fe1<-as.data.frame(cbind(SW_Fe$index1, SW_Fe$index2))
SW_Ti1<-as.data.frame(cbind(SW_Ti$index1, SW_Ti$index2))
SW_K1<-as.data.frame(cbind(SW_K$index1, SW_K$index2))
SW_KT1<-as.data.frame(cbind(SW_KT$index1, SW_KT$index2))
SW_tkc1<-as.data.frame(cbind(SW_tkc$index1, SW_tkc$index2))

#remove duplicates in V1; this shouldn't influence RMS too much as there are no hiatuses inferred (ie major jumps)
SW_Ca1<-SW_Ca1[!duplicated(SW_Ca1[,1]), ]
SW_Fe1<-SW_Fe1[!duplicated(SW_Fe1[,1]), ]
SW_Ti1<-SW_Ti1[!duplicated(SW_Ti1[,1]), ]
SW_K1<-SW_K1[!duplicated(SW_K1[,1]), ]
SW_KT1<-SW_KT1[!duplicated(SW_KT1[,1]), ]
SW_tkc1<-SW_tkc1[!duplicated(SW_tkc1[,1]), ]

#now they can be bound together
SW_dat2<-cbind(SW_Ca1$V1, SW_Ca1$V2,SW_Fe1$V2, SW_Ca1$V2,SW_Ti1$V2, SW_Ca1$V2,SW_K1$V2, 
               SW_Ca1$V2,SW_KT1$V2, SW_Ca1$V2,SW_tkc1$V2,
               SW_Fe1$V2,SW_Ti1$V2, SW_Fe1$V2,SW_K1$V2, SW_Fe1$V2,SW_KT1$V2, SW_Fe1$V2,SW_tkc1$V2,
               SW_Ti1$V2,SW_K1$V2, SW_Ti1$V2,SW_KT1$V2, SW_Ti1$V2,SW_tkc1$V2,
               SW_K1$V2,SW_KT1$V2, SW_K1$V2,SW_tkc1$V2,
               SW_KT1$V2,SW_tkc1$V2)

colnames(SW_dat2)<-c("ar", "W_Ca","W_Fe",  "W_Ca","W_Ti", "W_Ca","W_K", "W_Ca","W_KT", "W_Ca","W_tkc",
                     "W_Fe","W_Ti", "W_Fe","W_K", "W_Fe","W_KT", "W_Fe","W_tkc",
                     "W_Ti","W_K", "W_Ti","W_KT", "W_Ti","W_tkc",
                     "W_K","W_KT", "W_K","W_tkc",
                     "W_KT","W_tkc")

SW_dat2<-as.data.frame(SW_dat2)

SW_dat2a<-subset(SW_dat2, ar>=1000 & ar<=6000)
SW_dat2b<-subset(SW_dat2, ar>=2000 & ar<=5000)





#full data:
m=(ncol(SW_dat2)-1)/2

n=length(SW_dat2$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SW_dat2[i,j*2]-SW_dat2[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #546.2128



#1000-6000 points:
m=(ncol(SW_dat2a)-1)/2

n=length(SW_dat2a$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SW_dat2a[i,j*2]-SW_dat2a[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #572.2006



#2000-5000 points:
m=(ncol(SW_dat2b)-1)/2

n=length(SW_dat2b$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SW_dat2b[i,j*2]-SW_dat2b[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #555.5163






#--------------------ARFELD-STEINBRUCH SCHMIDT--------------------------------------
#Arfeld = A, Steinbruch Schmidt = S

#Ca:
SA_Ca<-dtw(ULS_Ca_ni$Ca, ULA_Ca_ni$Ca, open.begin=F, open.end=F, step.pattern=asymmetricP05,
           window.type="none", keep.internals=T)

dtwPlotThreeWay(SA_Ca)
#dtwPlotTwoWay(SA_Ca, offset=8)

s<-ULS_Ca_ni[SA_Ca$index1,1]
a<-ULA_Ca_ni[SA_Ca$index2,1]
SA<-cbind(s,a)
SA<-SA[!duplicated(SA[,2]), ]
SA<-SA[!duplicated(SA[,1]), ]
ULS_Ca_ni_t<-tune(ULS_Ca_ni, SA)

ggplot()+
  geom_line(data=ULA_Ca_ni, aes(x=depth_cm, y=Ca), col="black")+
  geom_line(data=ULS_Ca_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

SA_Ca$normalizedDistance





#Fe:
SA_Fe<-dtw(ULS_Fe_ni$Fe, ULA_Fe_ni$Fe, open.begin=F, open.end=F, step.pattern=asymmetricP05,
           window.type="none", keep.internals=T)

dtwPlotThreeWay(SA_Fe)
#dtwPlotTwoWay(SA_Fe, offset=8)

s<-ULS_Fe_ni[SA_Fe$index1,1]
a<-ULA_Fe_ni[SA_Fe$index2,1]
SA<-cbind(s,a)
SA<-SA[!duplicated(SA[,2]), ]
SA<-SA[!duplicated(SA[,1]), ]
ULS_Fe_ni_t<-tune(ULS_Fe_ni, SA)

ggplot()+
  geom_line(data=ULA_Fe_ni, aes(x=depth_cm, y=Fe), col="black")+
  geom_line(data=ULS_Fe_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

#SA_Fe$distance
SA_Fe$normalizedDistance




#K:
SA_K<-dtw(ULS_K_ni$K, ULA_K_ni$K, open.begin=F, open.end=F, step.pattern=asymmetricP05,
          window.type="none", keep.internals=T)

dtwPlotThreeWay(SA_K)
#dtwPlotTwoWay(SA_K, offset=8)

s<-ULS_K_ni[SA_K$index1,1]
a<-ULA_K_ni[SA_K$index2,1]
SA<-cbind(s,a)
SA<-SA[!duplicated(SA[,2]), ]
SA<-SA[!duplicated(SA[,1]), ]
ULS_K_ni_t<-tune(ULS_K_ni, SA)

ggplot()+
  geom_line(data=ULA_K_ni, aes(x=depth_cm, y=K), col="black")+
  geom_line(data=ULS_K_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

#SA_K$distance
SA_K$normalizedDistance





#Ti:
SA_Ti<-dtw(ULS_Ti_ni$Ti, ULA_Ti_ni$Ti, open.begin=F, open.end=F, step.pattern=asymmetricP05,
           window.type="none", keep.internals=T)

dtwPlotThreeWay(SA_Ti)
#dtwPlotTwoWay(SA_Ti, offset=8)

s<-ULS_Ti_ni[SA_Ti$index1,1]
a<-ULA_Ti_ni[SA_Ti$index2,1]
SA<-cbind(s,a)
SA<-SA[!duplicated(SA[,2]), ]
SA<-SA[!duplicated(SA[,1]), ]
ULS_Ti_ni_t<-tune(ULS_Ti_ni, SA)

ggplot()+
  geom_line(data=ULA_Ti_ni, aes(x=depth_cm, y=Ti), col="black")+
  geom_line(data=ULS_Ti_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

#SA_Ti$distance
SA_Ti$normalizedDistance





#(Ti+K)/Ca:
SA_tkc<-dtw(ULS_TKC_ni$ratio, ULA_TKC_ni$ratio, open.begin=T, open.end=F, step.pattern=asymmetricP05,
            window.type="none", keep.internals=T)

dtwPlotThreeWay(SA_tkc)
#dtwPlotTwoWay(SA_tkc, offset=8)

s<-ULS_TKC_ni[SA_tkc$index1,1]
a<-ULA_TKC_ni[SA_tkc$index2,1]
SA<-cbind(s,a)
SA<-SA[!duplicated(SA[,2]), ]
SA<-SA[!duplicated(SA[,1]), ]
ULS_tkc_ni_t<-tune(ULS_TKC_ni, SA)

ggplot()+
  geom_line(data=ULA_TKC_ni, aes(x=depth_cm, y=ratio), col="black")+
  geom_line(data=ULS_tkc_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

#SA_tkc$distance
SA_tkc$normalizedDistance





#K/Ti:
SA_KT<-dtw(ULS_KT_ni$ratio, ULA_KT_ni$ratio, open.begin=F, open.end=F, step.pattern=asymmetricP05,
           keep.internals=T)

dtwPlotThreeWay(SA_KT)
#dtwPlotTwoWay(SA_KT, offset=8)

s<-ULS_KT_ni[SA_KT$index1,1]
a<-ULA_KT_ni[SA_KT$index2,1]
SA<-cbind(s,a)
SA<-SA[!duplicated(SA[,2]), ]
SA<-SA[!duplicated(SA[,1]), ]
ULS_KT_ni_t<-tune(ULS_KT_ni, SA)

ggplot()+
  geom_line(data=ULA_KT_ni, aes(x=depth_cm, y=ratio), col="black")+
  geom_line(data=ULS_KT_ni_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

#SA_KT$distance
SA_KT$normalizedDistance




#---------------------------alignment plots S-A----------------------------------------
#----------------------------Open-open---------------------------------------------
SA_Ca<-dtw(ULS_Ca_ni$Ca, ULA_Ca_ni$Ca,open.begin=T, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
SA_Fe<-dtw(ULS_Fe_ni$Fe, ULA_Fe_ni$Fe,open.begin=T, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
SA_Ti<-dtw(ULS_Ti_ni$Ti, ULA_Ti_ni$Ti,open.begin=T, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
SA_K<-dtw(ULS_K_ni$K, ULA_K_ni$K,open.begin=T, open.end=T,
          step.pattern=asymmetricP05, keep.internals=T)
SA_KT<-dtw(ULS_KT_ni$ratio, ULA_KT_ni$ratio, open.begin=T, open.end=T, 
           step.pattern=asymmetricP05, keep.internals=T)
SA_tkc<-dtw(ULS_TKC_ni$ratio, ULA_TKC_ni$ratio, open.begin=T, open.end=T,
            step.pattern=asymmetricP05,keep.internals=T)

plot(SA_Ca, col='#fde725',ylab="Arfeld", xlab="Steinbruch Schmidt", main="Open - Open")
lines(x=SA_Ca$index1, y=SA_Ca$index2, col="#fde725", lw=3)
lines(x=SA_Fe$index1, y=SA_Fe$index2, col="#7ad151", lw=3)
lines(x=SA_Ti$index1, y=SA_Ti$index2, col="#22a884", lw=3)
lines(x=SA_K$index1, y=SA_K$index2, col="#2a788e", lw=3)
lines(x=SA_KT$index1, y=SA_KT$index2, col="#414487", lw=3)
lines(x=SA_tkc$index1, y=SA_tkc$index2, col="#440154", lw=3)
legend("topleft",                    # Add legend to plot
       legend = c("Ca","Fe","Ti","K","K/Ti", "(Ti+K)/Ca"),
       col = c("#fde725","#7ad151","#22a884","#2a788e",
               "#414487","#440154"),pch = 16)





#make dataframes
SA_Ca1<-as.data.frame(cbind(SA_Ca$index1, SA_Ca$index2))
SA_Fe1<-as.data.frame(cbind(SA_Fe$index1, SA_Fe$index2))
SA_Ti1<-as.data.frame(cbind(SA_Ti$index1, SA_Ti$index2))
SA_K1<-as.data.frame(cbind(SA_K$index1, SA_K$index2))
SA_KT1<-as.data.frame(cbind(SA_KT$index1, SA_KT$index2))
SA_tkc1<-as.data.frame(cbind(SA_tkc$index1, SA_tkc$index2))

#remove duplicates in V1; this shouldn't influence RMS too much as there are no hiatuses inferred (ie major jumps)
SA_Ca1<-SA_Ca1[!duplicated(SA_Ca1[,1]), ]
SA_Fe1<-SA_Fe1[!duplicated(SA_Fe1[,1]), ]
SA_Ti1<-SA_Ti1[!duplicated(SA_Ti1[,1]), ]
SA_K1<-SA_K1[!duplicated(SA_K1[,1]), ]
SA_KT1<-SA_KT1[!duplicated(SA_KT1[,1]), ]
SA_tkc1<-SA_tkc1[!duplicated(SA_tkc1[,1]), ]

#now they can be bound together
SA_dat2<-cbind(SA_Ca1$V1, SA_Ca1$V2,SA_Fe1$V2, SA_Ca1$V2,SA_Ti1$V2, SA_Ca1$V2,SA_K1$V2, 
               SA_Ca1$V2,SA_KT1$V2, SA_Ca1$V2,SA_tkc1$V2,
               SA_Fe1$V2,SA_Ti1$V2, SA_Fe1$V2,SA_K1$V2, SA_Fe1$V2,SA_KT1$V2, SA_Fe1$V2,SA_tkc1$V2,
               SA_Ti1$V2,SA_K1$V2, SA_Ti1$V2,SA_KT1$V2, SA_Ti1$V2,SA_tkc1$V2,
               SA_K1$V2,SA_KT1$V2, SA_K1$V2,SA_tkc1$V2,
               SA_KT1$V2,SA_tkc1$V2)

colnames(SA_dat2)<-c("ar", "W_Ca","W_Fe",  "W_Ca","W_Ti", "W_Ca","W_K", "W_Ca","W_KT", "W_Ca","W_tkc",
                     "W_Fe","W_Ti", "W_Fe","W_K", "W_Fe","W_KT", "W_Fe","W_tkc",
                     "W_Ti","W_K", "W_Ti","W_KT", "W_Ti","W_tkc",
                     "W_K","W_KT", "W_K","W_tkc",
                     "W_KT","W_tkc")

SA_dat2<-as.data.frame(SA_dat2)

SA_dat2a<-subset(SA_dat2, ar>=1000 & ar<=6000)
SA_dat2b<-subset(SA_dat2, ar>=2000 & ar<=5000)





#full data:
m=(ncol(SA_dat2)-1)/2

n=length(SA_dat2$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SA_dat2[i,j*2]-SA_dat2[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #392.1655



#1000-6000 points:
m=(ncol(SA_dat2a)-1)/2

n=length(SA_dat2a$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SA_dat2a[i,j*2]-SA_dat2a[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #436.275



#2000-5000 points:
m=(ncol(SA_dat2b)-1)/2

n=length(SA_dat2b$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SA_dat2b[i,j*2]-SA_dat2b[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #374.0009



#---------------------------Closed-closed----------------------------------------
SA_Ca<-dtw(ULS_Ca_ni$Ca, ULA_Ca_ni$Ca,open.begin=F, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
SA_Fe<-dtw(ULS_Fe_ni$Fe, ULA_Fe_ni$Fe,open.begin=F, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
SA_Ti<-dtw(ULS_Ti_ni$Ti, ULA_Ti_ni$Ti,open.begin=F, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
SA_K<-dtw(ULS_K_ni$K, ULA_K_ni$K,open.begin=F, open.end=F,
          step.pattern=asymmetricP05, keep.internals=T)
SA_KT<-dtw(ULS_KT_ni$ratio, ULA_KT_ni$ratio, open.begin=F, open.end=F, 
           step.pattern=asymmetricP05, keep.internals=T)
SA_tkc<-dtw(ULS_TKC_ni$ratio, ULA_TKC_ni$ratio, open.begin=F, open.end=F,
            step.pattern=asymmetricP05,keep.internals=T)

plot(SA_Ca, col='#fde725',ylab="Arfeld", xlab="Steinbruch Schmidt", main="Closed - Closed")
lines(x=SA_Ca$index1, y=SA_Ca$index2, col="#fde725", lw=3)
lines(x=SA_Fe$index1, y=SA_Fe$index2, col="#7ad151", lw=3)
lines(x=SA_Ti$index1, y=SA_Ti$index2, col="#22a884", lw=3)
lines(x=SA_K$index1, y=SA_K$index2, col="#2a788e", lw=3)
lines(x=SA_KT$index1, y=SA_KT$index2, col="#414487", lw=3)
lines(x=SA_tkc$index1, y=SA_tkc$index2, col="#440154", lw=3)
legend("topleft",                    # Add legend to plot
       legend = c("Ca","Fe","Ti","K","K/Ti", "(Ti+K)/Ca"),
       col = c("#fde725","#7ad151","#22a884","#2a788e",
               "#414487","#440154"),pch = 16)



#make dataframes
SA_Ca1<-as.data.frame(cbind(SA_Ca$index1, SA_Ca$index2))
SA_Fe1<-as.data.frame(cbind(SA_Fe$index1, SA_Fe$index2))
SA_Ti1<-as.data.frame(cbind(SA_Ti$index1, SA_Ti$index2))
SA_K1<-as.data.frame(cbind(SA_K$index1, SA_K$index2))
SA_KT1<-as.data.frame(cbind(SA_KT$index1, SA_KT$index2))
SA_tkc1<-as.data.frame(cbind(SA_tkc$index1, SA_tkc$index2))

#remove duplicates in V1; this shouldn't influence RMS too much as there are no hiatuses inferred (ie major jumps)
SA_Ca1<-SA_Ca1[!duplicated(SA_Ca1[,1]), ]
SA_Fe1<-SA_Fe1[!duplicated(SA_Fe1[,1]), ]
SA_Ti1<-SA_Ti1[!duplicated(SA_Ti1[,1]), ]
SA_K1<-SA_K1[!duplicated(SA_K1[,1]), ]
SA_KT1<-SA_KT1[!duplicated(SA_KT1[,1]), ]
SA_tkc1<-SA_tkc1[!duplicated(SA_tkc1[,1]), ]

#now they can be bound together
SA_dat2<-cbind(SA_Ca1$V1, SA_Ca1$V2,SA_Fe1$V2, SA_Ca1$V2,SA_Ti1$V2, SA_Ca1$V2,SA_K1$V2, 
               SA_Ca1$V2,SA_KT1$V2, SA_Ca1$V2,SA_tkc1$V2,
               SA_Fe1$V2,SA_Ti1$V2, SA_Fe1$V2,SA_K1$V2, SA_Fe1$V2,SA_KT1$V2, SA_Fe1$V2,SA_tkc1$V2,
               SA_Ti1$V2,SA_K1$V2, SA_Ti1$V2,SA_KT1$V2, SA_Ti1$V2,SA_tkc1$V2,
               SA_K1$V2,SA_KT1$V2, SA_K1$V2,SA_tkc1$V2,
               SA_KT1$V2,SA_tkc1$V2)

colnames(SA_dat2)<-c("ar", "W_Ca","W_Fe",  "W_Ca","W_Ti", "W_Ca","W_K", "W_Ca","W_KT", "W_Ca","W_tkc",
                     "W_Fe","W_Ti", "W_Fe","W_K", "W_Fe","W_KT", "W_Fe","W_tkc",
                     "W_Ti","W_K", "W_Ti","W_KT", "W_Ti","W_tkc",
                     "W_K","W_KT", "W_K","W_tkc",
                     "W_KT","W_tkc")

SA_dat2<-as.data.frame(SA_dat2)

SA_dat2a<-subset(SA_dat2, ar>=1000 & ar<=6000)
SA_dat2b<-subset(SA_dat2, ar>=2000 & ar<=5000)





#full data:
m=(ncol(SA_dat2)-1)/2

n=length(SA_dat2$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SA_dat2[i,j*2]-SA_dat2[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #382.7797



#1000-6000 points:
m=(ncol(SA_dat2a)-1)/2

n=length(SA_dat2a$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SA_dat2a[i,j*2]-SA_dat2a[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #436.0663



#2000-5000 points:
m=(ncol(SA_dat2b)-1)/2

n=length(SA_dat2b$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SA_dat2b[i,j*2]-SA_dat2b[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #374.0009



#-------------------------------Open-closed--------------------------------------
SA_Ca<-dtw(ULS_Ca_ni$Ca, ULA_Ca_ni$Ca,open.begin=T, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
SA_Fe<-dtw(ULS_Fe_ni$Fe, ULA_Fe_ni$Fe,open.begin=T, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
SA_Ti<-dtw(ULS_Ti_ni$Ti, ULA_Ti_ni$Ti,open.begin=T, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
SA_K<-dtw(ULS_K_ni$K, ULA_K_ni$K,open.begin=T, open.end=F,
          step.pattern=asymmetricP05, keep.internals=T)
SA_KT<-dtw(ULS_KT_ni$ratio, ULA_KT_ni$ratio, open.begin=T, open.end=F, 
           step.pattern=asymmetricP05, keep.internals=T)
SA_tkc<-dtw(ULS_TKC_ni$ratio, ULA_TKC_ni$ratio, open.begin=T, open.end=F,
            step.pattern=asymmetricP05,keep.internals=T)

plot(SA_Ca, col='#fde725',ylab="Arfeld", xlab="Steinbruch Schmidt", main="Open - Closed")
lines(x=SA_Ca$index1, y=SA_Ca$index2, col="#fde725", lw=3)
lines(x=SA_Fe$index1, y=SA_Fe$index2, col="#7ad151", lw=3)
lines(x=SA_Ti$index1, y=SA_Ti$index2, col="#22a884", lw=3)
lines(x=SA_K$index1, y=SA_K$index2, col="#2a788e", lw=3)
lines(x=SA_KT$index1, y=SA_KT$index2, col="#414487", lw=3)
lines(x=SA_tkc$index1, y=SA_tkc$index2, col="#440154", lw=3)
legend("topleft",                    # Add legend to plot
       legend = c("Ca","Fe","Ti","K","K/Ti", "(Ti+K)/Ca"),
       col = c("#fde725","#7ad151","#22a884","#2a788e",
               "#414487","#440154"),pch = 16)



#make dataframes
SA_Ca1<-as.data.frame(cbind(SA_Ca$index1, SA_Ca$index2))
SA_Fe1<-as.data.frame(cbind(SA_Fe$index1, SA_Fe$index2))
SA_Ti1<-as.data.frame(cbind(SA_Ti$index1, SA_Ti$index2))
SA_K1<-as.data.frame(cbind(SA_K$index1, SA_K$index2))
SA_KT1<-as.data.frame(cbind(SA_KT$index1, SA_KT$index2))
SA_tkc1<-as.data.frame(cbind(SA_tkc$index1, SA_tkc$index2))

#remove duplicates in V1; this shouldn't influence RMS too much as there are no hiatuses inferred (ie major jumps)
SA_Ca1<-SA_Ca1[!duplicated(SA_Ca1[,1]), ]
SA_Fe1<-SA_Fe1[!duplicated(SA_Fe1[,1]), ]
SA_Ti1<-SA_Ti1[!duplicated(SA_Ti1[,1]), ]
SA_K1<-SA_K1[!duplicated(SA_K1[,1]), ]
SA_KT1<-SA_KT1[!duplicated(SA_KT1[,1]), ]
SA_tkc1<-SA_tkc1[!duplicated(SA_tkc1[,1]), ]

#now they can be bound together
SA_dat2<-cbind(SA_Ca1$V1, SA_Ca1$V2,SA_Fe1$V2, SA_Ca1$V2,SA_Ti1$V2, SA_Ca1$V2,SA_K1$V2, 
               SA_Ca1$V2,SA_KT1$V2, SA_Ca1$V2,SA_tkc1$V2,
               SA_Fe1$V2,SA_Ti1$V2, SA_Fe1$V2,SA_K1$V2, SA_Fe1$V2,SA_KT1$V2, SA_Fe1$V2,SA_tkc1$V2,
               SA_Ti1$V2,SA_K1$V2, SA_Ti1$V2,SA_KT1$V2, SA_Ti1$V2,SA_tkc1$V2,
               SA_K1$V2,SA_KT1$V2, SA_K1$V2,SA_tkc1$V2,
               SA_KT1$V2,SA_tkc1$V2)

colnames(SA_dat2)<-c("ar", "W_Ca","W_Fe",  "W_Ca","W_Ti", "W_Ca","W_K", "W_Ca","W_KT", "W_Ca","W_tkc",
                     "W_Fe","W_Ti", "W_Fe","W_K", "W_Fe","W_KT", "W_Fe","W_tkc",
                     "W_Ti","W_K", "W_Ti","W_KT", "W_Ti","W_tkc",
                     "W_K","W_KT", "W_K","W_tkc",
                     "W_KT","W_tkc")

SA_dat2<-as.data.frame(SA_dat2)

SA_dat2a<-subset(SA_dat2, ar>=1000 & ar<=6000)
SA_dat2b<-subset(SA_dat2, ar>=2000 & ar<=5000)





#full data:
m=(ncol(SA_dat2)-1)/2

n=length(SA_dat2$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SA_dat2[i,j*2]-SA_dat2[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #391.7716



#1000-6000 points:
m=(ncol(SA_dat2a)-1)/2

n=length(SA_dat2a$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SA_dat2a[i,j*2]-SA_dat2a[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #436.275



#2000-5000 points:
m=(ncol(SA_dat2b)-1)/2

n=length(SA_dat2b$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SA_dat2b[i,j*2]-SA_dat2b[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #374.0009



#-------------------------------Closed-open---------------------------------------
SA_Ca<-dtw(ULS_Ca_ni$Ca, ULA_Ca_ni$Ca,open.begin=F, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
SA_Fe<-dtw(ULS_Fe_ni$Fe, ULA_Fe_ni$Fe,open.begin=F, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
SA_Ti<-dtw(ULS_Ti_ni$Ti, ULA_Ti_ni$Ti,open.begin=F, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
SA_K<-dtw(ULS_K_ni$K, ULA_K_ni$K,open.begin=F, open.end=T,
          step.pattern=asymmetricP05, keep.internals=T)
SA_KT<-dtw(ULS_KT_ni$ratio, ULA_KT_ni$ratio, open.begin=F, open.end=T, 
           step.pattern=asymmetricP05, keep.internals=T)
SA_tkc<-dtw(ULS_TKC_ni$ratio, ULA_TKC_ni$ratio, open.begin=F, open.end=T,
            step.pattern=asymmetricP05,keep.internals=T)

plot(SA_Ca, col='#fde725',ylab="Arfeld", xlab="Steinbruch Schmidt", main="Closed - Open")
lines(x=SA_Ca$index1, y=SA_Ca$index2, col="#fde725", lw=3)
lines(x=SA_Fe$index1, y=SA_Fe$index2, col="#7ad151", lw=3)
lines(x=SA_Ti$index1, y=SA_Ti$index2, col="#22a884", lw=3)
lines(x=SA_K$index1, y=SA_K$index2, col="#2a788e", lw=3)
lines(x=SA_KT$index1, y=SA_KT$index2, col="#414487", lw=3)
lines(x=SA_tkc$index1, y=SA_tkc$index2, col="#440154", lw=3)
legend("topleft",                    # Add legend to plot
       legend = c("Ca","Fe","Ti","K","K/Ti", "(Ti+K)/Ca"),
       col = c("#fde725","#7ad151","#22a884","#2a788e",
               "#414487","#440154"),pch = 16)




#make dataframes
SA_Ca1<-as.data.frame(cbind(SA_Ca$index1, SA_Ca$index2))
SA_Fe1<-as.data.frame(cbind(SA_Fe$index1, SA_Fe$index2))
SA_Ti1<-as.data.frame(cbind(SA_Ti$index1, SA_Ti$index2))
SA_K1<-as.data.frame(cbind(SA_K$index1, SA_K$index2))
SA_KT1<-as.data.frame(cbind(SA_KT$index1, SA_KT$index2))
SA_tkc1<-as.data.frame(cbind(SA_tkc$index1, SA_tkc$index2))

#remove duplicates in V1; this shouldn't influence RMS too much as there are no hiatuses inferred (ie major jumps)
SA_Ca1<-SA_Ca1[!duplicated(SA_Ca1[,1]), ]
SA_Fe1<-SA_Fe1[!duplicated(SA_Fe1[,1]), ]
SA_Ti1<-SA_Ti1[!duplicated(SA_Ti1[,1]), ]
SA_K1<-SA_K1[!duplicated(SA_K1[,1]), ]
SA_KT1<-SA_KT1[!duplicated(SA_KT1[,1]), ]
SA_tkc1<-SA_tkc1[!duplicated(SA_tkc1[,1]), ]

#now they can be bound together
SA_dat2<-cbind(SA_Ca1$V1, SA_Ca1$V2,SA_Fe1$V2, SA_Ca1$V2,SA_Ti1$V2, SA_Ca1$V2,SA_K1$V2, 
               SA_Ca1$V2,SA_KT1$V2, SA_Ca1$V2,SA_tkc1$V2,
               SA_Fe1$V2,SA_Ti1$V2, SA_Fe1$V2,SA_K1$V2, SA_Fe1$V2,SA_KT1$V2, SA_Fe1$V2,SA_tkc1$V2,
               SA_Ti1$V2,SA_K1$V2, SA_Ti1$V2,SA_KT1$V2, SA_Ti1$V2,SA_tkc1$V2,
               SA_K1$V2,SA_KT1$V2, SA_K1$V2,SA_tkc1$V2,
               SA_KT1$V2,SA_tkc1$V2)

colnames(SA_dat2)<-c("ar", "W_Ca","W_Fe",  "W_Ca","W_Ti", "W_Ca","W_K", "W_Ca","W_KT", "W_Ca","W_tkc",
                     "W_Fe","W_Ti", "W_Fe","W_K", "W_Fe","W_KT", "W_Fe","W_tkc",
                     "W_Ti","W_K", "W_Ti","W_KT", "W_Ti","W_tkc",
                     "W_K","W_KT", "W_K","W_tkc",
                     "W_KT","W_tkc")

SA_dat2<-as.data.frame(SA_dat2)

SA_dat2a<-subset(SA_dat2, ar>=1000 & ar<=6000)
SA_dat2b<-subset(SA_dat2, ar>=2000 & ar<=5000)





#full data:
m=(ncol(SA_dat2)-1)/2

n=length(SA_dat2$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SA_dat2[i,j*2]-SA_dat2[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #383.1829



#1000-6000 points:
m=(ncol(SA_dat2a)-1)/2

n=length(SA_dat2a$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SA_dat2a[i,j*2]-SA_dat2a[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #436.0663



#2000-5000 points:
m=(ncol(SA_dat2b)-1)/2

n=length(SA_dat2b$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(SA_dat2b[i,j*2]-SA_dat2b[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #374.0009




#---------------------------alignment plots A-S----------------------------------------
#Arfeld is thicker, but Steinbruch Schmidt (might) have more cycles and thus contain more time.
#Therefore, matching Arfeld to Steinbruch Schmidt instead of the other way around (as above)
#may result in a more realistic fit.

#----------------------------------Open-open----------------------------------------
AS_Ca<-dtw(ULA_Ca_ni$Ca, ULS_Ca_ni$Ca,open.begin=T, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
AS_Fe<-dtw(ULA_Fe_ni$Fe, ULS_Fe_ni$Fe,open.begin=T, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
AS_Ti<-dtw(ULA_Ti_ni$Ti, ULS_Ti_ni$Ti,open.begin=T, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
AS_K<-dtw(ULA_K_ni$K, ULS_K_ni$K,open.begin=T, open.end=T,
          step.pattern=asymmetricP05, keep.internals=T)
AS_KT<-dtw(ULA_KT_ni$ratio, ULS_KT_ni$ratio, open.begin=T, open.end=T, 
           step.pattern=asymmetricP05, keep.internals=T)
AS_tkc<-dtw(ULA_TKC_ni$ratio, ULS_TKC_ni$ratio, open.begin=T, open.end=T,
            step.pattern=asymmetricP05,keep.internals=T)

plot(AS_Ca, col='#fde725',ylab="Steinbruch Schmidt", xlab="Arfeld", main="Open - Open")
lines(x=AS_Ca$index1, y=AS_Ca$index2, col="#fde725", lw=3)
lines(x=AS_Fe$index1, y=AS_Fe$index2, col="#7ad151", lw=3)
lines(x=AS_Ti$index1, y=AS_Ti$index2, col="#22a884", lw=3)
lines(x=AS_K$index1, y=AS_K$index2, col="#2a788e", lw=3)
lines(x=AS_KT$index1, y=AS_KT$index2, col="#414487", lw=3)
lines(x=AS_tkc$index1, y=AS_tkc$index2, col="#440154", lw=3)
legend("topleft",                    # Add legend to plot
       legend = c("Ca","Fe","Ti","K","K/Ti", "(Ti+K)/Ca"),
       col = c("#fde725","#7ad151","#22a884","#2a788e",
               "#414487","#440154"),pch = 16)


#make dataframes
AS_Ca1<-as.data.frame(cbind(AS_Ca$index1, AS_Ca$index2))
AS_Fe1<-as.data.frame(cbind(AS_Fe$index1, AS_Fe$index2))
AS_Ti1<-as.data.frame(cbind(AS_Ti$index1, AS_Ti$index2))
AS_K1<-as.data.frame(cbind(AS_K$index1, AS_K$index2))
AS_KT1<-as.data.frame(cbind(AS_KT$index1, AS_KT$index2))
AS_tkc1<-as.data.frame(cbind(AS_tkc$index1, AS_tkc$index2))

#remove duplicates in V1; this shouldn't influence RMS too much as there are no hiatuses inferred (ie major jumps)
AS_Ca1<-AS_Ca1[!duplicated(AS_Ca1[,1]), ]
AS_Fe1<-AS_Fe1[!duplicated(AS_Fe1[,1]), ]
AS_Ti1<-AS_Ti1[!duplicated(AS_Ti1[,1]), ]
AS_K1<-AS_K1[!duplicated(AS_K1[,1]), ]
AS_KT1<-AS_KT1[!duplicated(AS_KT1[,1]), ]
AS_tkc1<-AS_tkc1[!duplicated(AS_tkc1[,1]), ]

#now they can be bound together
AS_dat2<-cbind(AS_Ca1$V1, AS_Ca1$V2,AS_Fe1$V2, AS_Ca1$V2,AS_Ti1$V2, AS_Ca1$V2,AS_K1$V2, 
               AS_Ca1$V2,AS_KT1$V2, AS_Ca1$V2,AS_tkc1$V2,
               AS_Fe1$V2,AS_Ti1$V2, AS_Fe1$V2,AS_K1$V2, AS_Fe1$V2,AS_KT1$V2, AS_Fe1$V2,AS_tkc1$V2,
               AS_Ti1$V2,AS_K1$V2, AS_Ti1$V2,AS_KT1$V2, AS_Ti1$V2,AS_tkc1$V2,
               AS_K1$V2,AS_KT1$V2, AS_K1$V2,AS_tkc1$V2,
               AS_KT1$V2,AS_tkc1$V2)

colnames(AS_dat2)<-c("ar", "W_Ca","W_Fe",  "W_Ca","W_Ti", "W_Ca","W_K", "W_Ca","W_KT", "W_Ca","W_tkc",
                     "W_Fe","W_Ti", "W_Fe","W_K", "W_Fe","W_KT", "W_Fe","W_tkc",
                     "W_Ti","W_K", "W_Ti","W_KT", "W_Ti","W_tkc",
                     "W_K","W_KT", "W_K","W_tkc",
                     "W_KT","W_tkc")

AS_dat2<-as.data.frame(AS_dat2)

AS_dat2a<-subset(AS_dat2, ar>=1000 & ar<=6000)
AS_dat2b<-subset(AS_dat2, ar>=2000 & ar<=5000)





#full data:
m=(ncol(AS_dat2)-1)/2

n=length(AS_dat2$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AS_dat2[i,j*2]-AS_dat2[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #420.2844



#1000-6000 points:
m=(ncol(AS_dat2a)-1)/2

n=length(AS_dat2a$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AS_dat2a[i,j*2]-AS_dat2a[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #463.8847



#2000-5000 points:
m=(ncol(AS_dat2b)-1)/2

n=length(AS_dat2b$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AS_dat2b[i,j*2]-AS_dat2b[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #409.41




#-----------------------------Closed-closed--------------------------------------
AS_Ca<-dtw(ULA_Ca_ni$Ca, ULS_Ca_ni$Ca,open.begin=F, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
AS_Fe<-dtw(ULA_Fe_ni$Fe, ULS_Fe_ni$Fe,open.begin=F, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
AS_Ti<-dtw(ULA_Ti_ni$Ti, ULS_Ti_ni$Ti,open.begin=F, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
AS_K<-dtw(ULA_K_ni$K, ULS_K_ni$K,open.begin=F, open.end=F,
          step.pattern=asymmetricP05, keep.internals=T)
AS_KT<-dtw(ULA_KT_ni$ratio, ULS_KT_ni$ratio, open.begin=F, open.end=F, 
           step.pattern=asymmetricP05, keep.internals=T)
AS_tkc<-dtw(ULA_TKC_ni$ratio, ULS_TKC_ni$ratio, open.begin=F, open.end=F,
            step.pattern=asymmetricP05,keep.internals=T)

plot(AS_Ca, col='#fde725',ylab="Steinbruch Schmidt", xlab="Arfeld", main="Closed - Closed")
lines(x=AS_Ca$index1, y=AS_Ca$index2, col="#fde725", lw=3)
lines(x=AS_Fe$index1, y=AS_Fe$index2, col="#7ad151", lw=3)
lines(x=AS_Ti$index1, y=AS_Ti$index2, col="#22a884", lw=3)
lines(x=AS_K$index1, y=AS_K$index2, col="#2a788e", lw=3)
lines(x=AS_KT$index1, y=AS_KT$index2, col="#414487", lw=3)
lines(x=AS_tkc$index1, y=AS_tkc$index2, col="#440154", lw=3)
legend("topleft",                    # Add legend to plot
       legend = c("Ca","Fe","Ti","K","K/Ti", "(Ti+K)/Ca"),
       col = c("#fde725","#7ad151","#22a884","#2a788e",
               "#414487","#440154"),pch = 16)



#make dataframes
AS_Ca1<-as.data.frame(cbind(AS_Ca$index1, AS_Ca$index2))
AS_Fe1<-as.data.frame(cbind(AS_Fe$index1, AS_Fe$index2))
AS_Ti1<-as.data.frame(cbind(AS_Ti$index1, AS_Ti$index2))
AS_K1<-as.data.frame(cbind(AS_K$index1, AS_K$index2))
AS_KT1<-as.data.frame(cbind(AS_KT$index1, AS_KT$index2))
AS_tkc1<-as.data.frame(cbind(AS_tkc$index1, AS_tkc$index2))

#remove duplicates in V1; this shouldn't influence RMS too much as there are no hiatuses inferred (ie major jumps)
AS_Ca1<-AS_Ca1[!duplicated(AS_Ca1[,1]), ]
AS_Fe1<-AS_Fe1[!duplicated(AS_Fe1[,1]), ]
AS_Ti1<-AS_Ti1[!duplicated(AS_Ti1[,1]), ]
AS_K1<-AS_K1[!duplicated(AS_K1[,1]), ]
AS_KT1<-AS_KT1[!duplicated(AS_KT1[,1]), ]
AS_tkc1<-AS_tkc1[!duplicated(AS_tkc1[,1]), ]

#now they can be bound together
AS_dat2<-cbind(AS_Ca1$V1, AS_Ca1$V2,AS_Fe1$V2, AS_Ca1$V2,AS_Ti1$V2, AS_Ca1$V2,AS_K1$V2, 
               AS_Ca1$V2,AS_KT1$V2, AS_Ca1$V2,AS_tkc1$V2,
               AS_Fe1$V2,AS_Ti1$V2, AS_Fe1$V2,AS_K1$V2, AS_Fe1$V2,AS_KT1$V2, AS_Fe1$V2,AS_tkc1$V2,
               AS_Ti1$V2,AS_K1$V2, AS_Ti1$V2,AS_KT1$V2, AS_Ti1$V2,AS_tkc1$V2,
               AS_K1$V2,AS_KT1$V2, AS_K1$V2,AS_tkc1$V2,
               AS_KT1$V2,AS_tkc1$V2)

colnames(AS_dat2)<-c("ar", "W_Ca","W_Fe",  "W_Ca","W_Ti", "W_Ca","W_K", "W_Ca","W_KT", "W_Ca","W_tkc",
                     "W_Fe","W_Ti", "W_Fe","W_K", "W_Fe","W_KT", "W_Fe","W_tkc",
                     "W_Ti","W_K", "W_Ti","W_KT", "W_Ti","W_tkc",
                     "W_K","W_KT", "W_K","W_tkc",
                     "W_KT","W_tkc")

AS_dat2<-as.data.frame(AS_dat2)

AS_dat2a<-subset(AS_dat2, ar>=1000 & ar<=6000)
AS_dat2b<-subset(AS_dat2, ar>=2000 & ar<=5000)





#full data:
m=(ncol(AS_dat2)-1)/2

n=length(AS_dat2$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AS_dat2[i,j*2]-AS_dat2[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #409.7838



#1000-6000 points:
m=(ncol(AS_dat2a)-1)/2

n=length(AS_dat2a$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AS_dat2a[i,j*2]-AS_dat2a[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #463.7021



#2000-5000 points:
m=(ncol(AS_dat2b)-1)/2

n=length(AS_dat2b$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AS_dat2b[i,j*2]-AS_dat2b[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #409.4097




#----------------------------Open-closed------------------------------------------
AS_Ca<-dtw(ULA_Ca_ni$Ca, ULS_Ca_ni$Ca,open.begin=T, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
AS_Fe<-dtw(ULA_Fe_ni$Fe, ULS_Fe_ni$Fe,open.begin=T, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
AS_Ti<-dtw(ULA_Ti_ni$Ti, ULS_Ti_ni$Ti,open.begin=T, open.end=F,
           step.pattern=asymmetricP05, keep.internals=T)
AS_K<-dtw(ULA_K_ni$K, ULS_K_ni$K,open.begin=T, open.end=F,
          step.pattern=asymmetricP05, keep.internals=T)
AS_KT<-dtw(ULA_KT_ni$ratio, ULS_KT_ni$ratio, open.begin=T, open.end=F, 
           step.pattern=asymmetricP05, keep.internals=T)
AS_tkc<-dtw(ULA_TKC_ni$ratio, ULS_TKC_ni$ratio, open.begin=T, open.end=F,
            step.pattern=asymmetricP05,keep.internals=T)

plot(AS_Ca, col='#fde725',ylab="Steinbruch Schmidt", xlab="Arfeld", main="Open - Closed")
lines(x=AS_Ca$index1, y=AS_Ca$index2, col="#fde725", lw=3)
lines(x=AS_Fe$index1, y=AS_Fe$index2, col="#7ad151", lw=3)
lines(x=AS_Ti$index1, y=AS_Ti$index2, col="#22a884", lw=3)
lines(x=AS_K$index1, y=AS_K$index2, col="#2a788e", lw=3)
lines(x=AS_KT$index1, y=AS_KT$index2, col="#414487", lw=3)
lines(x=AS_tkc$index1, y=AS_tkc$index2, col="#440154", lw=3)
legend("topleft",                    # Add legend to plot
       legend = c("Ca","Fe","Ti","K","K/Ti", "(Ti+K)/Ca"),
       col = c("#fde725","#7ad151","#22a884","#2a788e",
               "#414487","#440154"),pch = 16)




#make dataframes
AS_Ca1<-as.data.frame(cbind(AS_Ca$index1, AS_Ca$index2))
AS_Fe1<-as.data.frame(cbind(AS_Fe$index1, AS_Fe$index2))
AS_Ti1<-as.data.frame(cbind(AS_Ti$index1, AS_Ti$index2))
AS_K1<-as.data.frame(cbind(AS_K$index1, AS_K$index2))
AS_KT1<-as.data.frame(cbind(AS_KT$index1, AS_KT$index2))
AS_tkc1<-as.data.frame(cbind(AS_tkc$index1, AS_tkc$index2))

#remove duplicates in V1; this shouldn't influence RMS too much as there are no hiatuses inferred (ie major jumps)
AS_Ca1<-AS_Ca1[!duplicated(AS_Ca1[,1]), ]
AS_Fe1<-AS_Fe1[!duplicated(AS_Fe1[,1]), ]
AS_Ti1<-AS_Ti1[!duplicated(AS_Ti1[,1]), ]
AS_K1<-AS_K1[!duplicated(AS_K1[,1]), ]
AS_KT1<-AS_KT1[!duplicated(AS_KT1[,1]), ]
AS_tkc1<-AS_tkc1[!duplicated(AS_tkc1[,1]), ]

#now they can be bound together
AS_dat2<-cbind(AS_Ca1$V1, AS_Ca1$V2,AS_Fe1$V2, AS_Ca1$V2,AS_Ti1$V2, AS_Ca1$V2,AS_K1$V2, 
               AS_Ca1$V2,AS_KT1$V2, AS_Ca1$V2,AS_tkc1$V2,
               AS_Fe1$V2,AS_Ti1$V2, AS_Fe1$V2,AS_K1$V2, AS_Fe1$V2,AS_KT1$V2, AS_Fe1$V2,AS_tkc1$V2,
               AS_Ti1$V2,AS_K1$V2, AS_Ti1$V2,AS_KT1$V2, AS_Ti1$V2,AS_tkc1$V2,
               AS_K1$V2,AS_KT1$V2, AS_K1$V2,AS_tkc1$V2,
               AS_KT1$V2,AS_tkc1$V2)

colnames(AS_dat2)<-c("ar", "W_Ca","W_Fe",  "W_Ca","W_Ti", "W_Ca","W_K", "W_Ca","W_KT", "W_Ca","W_tkc",
                     "W_Fe","W_Ti", "W_Fe","W_K", "W_Fe","W_KT", "W_Fe","W_tkc",
                     "W_Ti","W_K", "W_Ti","W_KT", "W_Ti","W_tkc",
                     "W_K","W_KT", "W_K","W_tkc",
                     "W_KT","W_tkc")

AS_dat2<-as.data.frame(AS_dat2)

AS_dat2a<-subset(AS_dat2, ar>=1000 & ar<=6000)
AS_dat2b<-subset(AS_dat2, ar>=2000 & ar<=5000)





#full data:
m=(ncol(AS_dat2)-1)/2

n=length(AS_dat2$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AS_dat2[i,j*2]-AS_dat2[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #420.0687



#1000-6000 points:
m=(ncol(AS_dat2a)-1)/2

n=length(AS_dat2a$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AS_dat2a[i,j*2]-AS_dat2a[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #463.8847



#2000-5000 points:
m=(ncol(AS_dat2b)-1)/2

n=length(AS_dat2b$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AS_dat2b[i,j*2]-AS_dat2b[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #409.41




#------------------------------Closed-open---------------------------------------
AS_Ca<-dtw(ULA_Ca_ni$Ca, ULS_Ca_ni$Ca,open.begin=F, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
AS_Fe<-dtw(ULA_Fe_ni$Fe, ULS_Fe_ni$Fe,open.begin=F, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
AS_Ti<-dtw(ULA_Ti_ni$Ti, ULS_Ti_ni$Ti,open.begin=F, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)
AS_K<-dtw(ULA_K_ni$K, ULS_K_ni$K,open.begin=F, open.end=T,
          step.pattern=asymmetricP05, keep.internals=T)
AS_KT<-dtw(ULA_KT_ni$ratio, ULS_KT_ni$ratio, open.begin=F, open.end=T, 
           step.pattern=asymmetricP05, keep.internals=T)
AS_tkc<-dtw(ULA_TKC_ni$ratio, ULS_TKC_ni$ratio, open.begin=F, open.end=T,
            step.pattern=asymmetricP05,keep.internals=T)

plot(AS_Ca, col='#fde725',ylab="Steinbruch Schmidt", xlab="Arfeld", main="Closed - Open")
lines(x=AS_Ca$index1, y=AS_Ca$index2, col="#fde725", lw=3)
lines(x=AS_Fe$index1, y=AS_Fe$index2, col="#7ad151", lw=3)
lines(x=AS_Ti$index1, y=AS_Ti$index2, col="#22a884", lw=3)
lines(x=AS_K$index1, y=AS_K$index2, col="#2a788e", lw=3)
lines(x=AS_KT$index1, y=AS_KT$index2, col="#414487", lw=3)
lines(x=AS_tkc$index1, y=AS_tkc$index2, col="#440154", lw=3)
legend("topleft",                    # Add legend to plot
       legend = c("Ca","Fe","Ti","K","K/Ti", "(Ti+K)/Ca"),
       col = c("#fde725","#7ad151","#22a884","#2a788e",
               "#414487","#440154"),pch = 16)




#make dataframes
AS_Ca1<-as.data.frame(cbind(AS_Ca$index1, AS_Ca$index2))
AS_Fe1<-as.data.frame(cbind(AS_Fe$index1, AS_Fe$index2))
AS_Ti1<-as.data.frame(cbind(AS_Ti$index1, AS_Ti$index2))
AS_K1<-as.data.frame(cbind(AS_K$index1, AS_K$index2))
AS_KT1<-as.data.frame(cbind(AS_KT$index1, AS_KT$index2))
AS_tkc1<-as.data.frame(cbind(AS_tkc$index1, AS_tkc$index2))

#remove duplicates in V1; this shouldn't influence RMS too much as there are no hiatuses inferred (ie major jumps)
AS_Ca1<-AS_Ca1[!duplicated(AS_Ca1[,1]), ]
AS_Fe1<-AS_Fe1[!duplicated(AS_Fe1[,1]), ]
AS_Ti1<-AS_Ti1[!duplicated(AS_Ti1[,1]), ]
AS_K1<-AS_K1[!duplicated(AS_K1[,1]), ]
AS_KT1<-AS_KT1[!duplicated(AS_KT1[,1]), ]
AS_tkc1<-AS_tkc1[!duplicated(AS_tkc1[,1]), ]

#now they can be bound together
AS_dat2<-cbind(AS_Ca1$V1, AS_Ca1$V2,AS_Fe1$V2, AS_Ca1$V2,AS_Ti1$V2, AS_Ca1$V2,AS_K1$V2, 
               AS_Ca1$V2,AS_KT1$V2, AS_Ca1$V2,AS_tkc1$V2,
               AS_Fe1$V2,AS_Ti1$V2, AS_Fe1$V2,AS_K1$V2, AS_Fe1$V2,AS_KT1$V2, AS_Fe1$V2,AS_tkc1$V2,
               AS_Ti1$V2,AS_K1$V2, AS_Ti1$V2,AS_KT1$V2, AS_Ti1$V2,AS_tkc1$V2,
               AS_K1$V2,AS_KT1$V2, AS_K1$V2,AS_tkc1$V2,
               AS_KT1$V2,AS_tkc1$V2)

colnames(AS_dat2)<-c("ar", "W_Ca","W_Fe",  "W_Ca","W_Ti", "W_Ca","W_K", "W_Ca","W_KT", "W_Ca","W_tkc",
                     "W_Fe","W_Ti", "W_Fe","W_K", "W_Fe","W_KT", "W_Fe","W_tkc",
                     "W_Ti","W_K", "W_Ti","W_KT", "W_Ti","W_tkc",
                     "W_K","W_KT", "W_K","W_tkc",
                     "W_KT","W_tkc")

AS_dat2<-as.data.frame(AS_dat2)

AS_dat2a<-subset(AS_dat2, ar>=1000 & ar<=6000)
AS_dat2b<-subset(AS_dat2, ar>=2000 & ar<=5000)





#full data:
m=(ncol(AS_dat2)-1)/2

n=length(AS_dat2$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AS_dat2[i,j*2]-AS_dat2[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #410.0049



#1000-6000 points:
m=(ncol(AS_dat2a)-1)/2

n=length(AS_dat2a$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AS_dat2a[i,j*2]-AS_dat2a[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #463.7021



#2000-5000 points:
m=(ncol(AS_dat2b)-1)/2

n=length(AS_dat2b$ar) 

y<-matrix(data=NA, nrow=n, ncol=m) #in separate matrix in case you want to look at individual correlations

for (j in 1:m){
  
  for (i in 1:n){
    
    y[i,j]<-(AS_dat2b[i,j*2]-AS_dat2b[i,(j*2)+1])^2
    
  }
}

#to calculate rms: sum all values, then *1/n and sqrt
rms<-sqrt(sum(y)*(1/(ncol(y)*nrow(y))))
rms #409.4097







#Ca:
SA_Ca<-dtw(ULA_Ca_n$Ca, ULS_Ca_n$Ca, open.begin=T, open.end=T, step.pattern=asymmetricP05,
           window.type="none", keep.internals=T)

dtwPlotThreeWay(SA_Ca)
#dtwPlotTwoWay(SA_Ca, offset=8)

a<-ULA_Ca_n[SA_Ca$index1,1]
s<-ULS_Ca_n[SA_Ca$index2,1]
SA<-cbind(a,s)
SA<-SA[!duplicated(SA[,2]), ]
SA<-SA[!duplicated(SA[,1]), ]
ULA_Ca_n_t<-tune(ULA_Ca_n, SA)

ggplot()+
  geom_line(data=ULS_Ca_n, aes(x=depth_cm, y=Ca), col="black")+
  geom_line(data=ULA_Ca_n_t, aes(x=X1, y=X2+3), col="blue")+
  theme_classic()

SA_Ca$normalizedDistance
SA_Ca$distance


#Closed-open
AS_Ca<-dtw(ULA_Ca_n$Ca, ULS_Ca_n$Ca,open.begin=F, open.end=T,
           step.pattern=asymmetricP05, keep.internals=T)

plot(AS_Ca, col='#fde725',ylab="Steinbruch Schmidt", xlab="Arfeld", main="Closed - Open")
lines(x=AS_Ca$index1, y=AS_Ca$index2, col="#fde725", lw=3)




#-----------------------------(Ti+K)/Ca plot----------------------------------------
#Arfeld and Steinbruch Schmidt fitted to Winsenberg:

#(Ti+K)/Ca:
AW_tkc<-dtw(ULA_TKC_ni$ratio, ULW_TKCn$ratio, open.begin=T, open.end=F, step.pattern=asymmetricP05,
            window.type="none", keep.internals=T)

a<-ULA_TKC_ni[AW_tkc$index1,1]
w<-ULW_TKCn[AW_tkc$index2,1]
aw<-cbind(a,w)
aw<-aw[!duplicated(aw[,2]), ]
aw<-aw[!duplicated(aw[,1]), ]
ULA_tkc_ni_t<-tune(ULA_TKC_ni, aw)

#(Ti+K)/Ca:
SW_tkc<-dtw(ULS_TKC_ni$ratio, ULW_TKCn$ratio, open.begin=T, open.end=F, step.pattern=asymmetricP05,
            window.type="none", keep.internals=T)

s<-ULS_TKC_ni[SW_tkc$index1,1]
w<-ULW_TKCn[SW_tkc$index2,1]
SW<-cbind(s,w)
SW<-SW[!duplicated(SW[,2]), ]
SW<-SW[!duplicated(SW[,1]), ]
ULS_tkc_ni_t<-tune(ULS_TKC_ni, SW)

w<-ggplot()+
  geom_line(data=ULW_TKCn, aes(x=depth_cm, y=ratio), col="black")+
  labs(x="depth (cm)", y="(Ti+K)/Ca [-]", title="Winsenberg (Ti+K)/Ca")+
  xlim(0,75)+
  coord_flip()+
  theme_classic()

s<-ggplot()+ #fitted to winsenberg
  geom_line(data=ULS_tkc_ni_t, aes(x=X1, y=X2), col="black")+
  labs(x="depth (cm)", y="(Ti+K)/Ca [-]", title="Steinbruch Schmidt (Ti+K)/Ca")+
  xlim(0,75)+
  coord_flip()+
  theme_classic()

a<-ggplot()+
  geom_line(data=ULA_tkc_ni_t, aes(x=X1, y=X2), col="black")+
  labs(x="depth (cm)", y="(Ti+K)/Ca [-]", title="Arfeld (Ti+K)/Ca")+
  xlim(0,75)+
  coord_flip()+
  theme_classic()

cowplot::plot_grid(w,a,s, ncol=3)




#Arfeld fitted to Winsenberg, Steinbruch Schmidt fitted to Arfeld:
#in order to get plot right
#(Ti+K)/Ca:
AW_tkc<-dtw(ULA_TKC_ni$ratio, ULW_TKCn$ratio, open.begin=T, open.end=F, step.pattern=asymmetricP05,
            window.type="none", keep.internals=T)

a<-ULA_TKC_ni[AW_tkc$index1,1]
w<-ULW_TKCn[AW_tkc$index2,1]
aw<-cbind(a,w)
aw<-aw[!duplicated(aw[,2]), ]
aw<-aw[!duplicated(aw[,1]), ]
ULA_tkc_ni_t<-tune(ULA_TKC_ni, aw)

dtwPlotThreeWay(AW_tkc)


#(Ti+K)/Ca:
SA_tkc<-dtw(ULS_TKC_ni$ratio, ULA_tkc_ni_t$X2, open.begin=T, open.end=F, step.pattern=asymmetricP05,
            window.type="none", keep.internals=T)

s<-ULS_TKC_ni[SA_tkc$index1,1]
a<-ULA_tkc_ni_t[SA_tkc$index2,1]
SA<-cbind(s,a)
SA<-SA[!duplicated(SA[,2]), ]
SA<-SA[!duplicated(SA[,1]), ]
ULS_tkc_ni_t<-tune(ULS_TKC_ni, SA) 

dtwPlotThreeWay(SA_tkc)


w<-ggplot()+
  geom_line(data=ULW_TKCn, aes(x=depth_cm, y=ratio), col="black")+
  labs(x="depth (cm)", y="(Ti+K)/Ca [-]", title="Winsenberg (Ti+K)/Ca")+
  xlim(0,75)+
  coord_flip()+
  theme_classic()

s<-ggplot()+ #fitted to Arfeld
  geom_line(data=ULS_tkc_ni_t, aes(x=X1, y=X2), col="black")+
  labs(x="depth (cm)", y="(Ti+K)/Ca [-]", title="Steinbruch Schmidt (Ti+K)/Ca")+
  xlim(0,75)+
  coord_flip()+
  theme_classic()

a<-ggplot()+
  geom_line(data=ULA_tkc_ni_t, aes(x=X1, y=X2), col="black")+
  labs(x="depth (cm)", y="(Ti+K)/Ca [-]", title="Arfeld (Ti+K)/Ca")+
  xlim(0,75)+
  coord_flip()+
  theme_classic()

cowplot::plot_grid(w,a,s, ncol=3)

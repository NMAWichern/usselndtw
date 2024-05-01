#code for generating elemental depth records from micro-XRF generated maps
#Dynamic Time Warping (DTW) has been applied on these depth records
#This code is for the Winsenberg Usseln Limestone (ULW)
#ULW consists of nine sections, ULW a1, a2, a3, b1, c1, c2, d1, d2, and f1.

#steps within this code (no non-rock area removal required):
#1) formatting data into excel files similar to ULA and ULS
#2) plot all individual segments together as a composite map
#3) determine overlap between segments and combine them
#4) average all rows in order to generate a composite depth record

#Note that this dataset was measured on a different Bruker M4 Tornado at MARUM, Bremen.
#Data format and required steps differ slightly to ULA and ULS.

library(tidyverse)
library(astrochron)
library(devtools)
library(ptw)
library(RColorBrewer)

setwd("insert wd")
#----------------Initial data checking and formatting------------------------------

#Load in data:
ULWa1_Ca<-read.table("ULWa1_Ca.txt", sep=";")
ULWa1_Fe<-read.table("ULWa1_Fe.txt", sep=";")
ULWa1_Ti<-read.table("ULWa1_Ti.txt", sep=";")
ULWa1_K<-read.table("ULWa1_K.txt", sep=";")
ULWa1_Mn<-read.table("ULWa1_Mn.txt", sep=";")
ULWa1_Si<-read.table("ULWa1_Si.txt", sep=";")

ULWa2_Ca<-read.table("ULWa2_Ca.txt", sep=";")
ULWa2_Fe<-read.table("ULWa2_Fe.txt", sep=";")
ULWa2_Ti<-read.table("ULWa2_Ti.txt", sep=";")
ULWa2_K<-read.table("ULWa2_K.txt", sep=";")
ULWa2_Mn<-read.table("ULWa2_Mn.txt", sep=";")
ULWa2_Si<-read.table("ULWa2_Si.txt", sep=";")

ULWa3_Ca<-read.table("ULWa3_Ca.txt", sep=";")
ULWa3_Fe<-read.table("ULWa3_Fe.txt", sep=";")
ULWa3_Ti<-read.table("ULWa3_Ti.txt", sep=";")
ULWa3_K<-read.table("ULWa3_K.txt", sep=";")
ULWa3_Mn<-read.table("ULWa3_Mn.txt", sep=";")
ULWa3_Si<-read.table("ULWa3_Si.txt", sep=";")

ULWb1_Ca<-read.table("ULWb1_Ca.txt", sep=";")
ULWb1_Fe<-read.table("ULWb1_Fe.txt", sep=";")
ULWb1_Ti<-read.table("ULWb1_Ti.txt", sep=";")
ULWb1_K<-read.table("ULWb1_K.txt", sep=";")
ULWb1_Mn<-read.table("ULWb1_Mn.txt", sep=";")
ULWb1_Si<-read.table("ULWb1_Si.txt", sep=";")

ULWc1_Ca<-read.table("ULWc1_Ca.txt", sep=";")
ULWc1_Fe<-read.table("ULWc1_Fe.txt", sep=";")
ULWc1_Ti<-read.table("ULWc1_Ti.txt", sep=";")
ULWc1_K<-read.table("ULWc1_K.txt", sep=";")
ULWc1_Mn<-read.table("ULWc1_Mn.txt", sep=";")
ULWc1_Si<-read.table("ULWc1_Si.txt", sep=";")

ULWc2_Ca<-read.table("ULWc2_Ca.txt", sep=";")
ULWc2_Fe<-read.table("ULWc2_Fe.txt", sep=";")
ULWc2_Ti<-read.table("ULWc2_Ti.txt", sep=";")
ULWc2_K<-read.table("ULWc2_K.txt", sep=";")
ULWc2_Mn<-read.table("ULWc2_Mn.txt", sep=";")
ULWc2_Si<-read.table("ULWc2_Si.txt", sep=";")

ULWd1_Ca<-read.table("ULWd1_Ca.txt", sep=";")
ULWd1_Fe<-read.table("ULWd1_Fe.txt", sep=";")
ULWd1_Ti<-read.table("ULWd1_Ti.txt", sep=";")
ULWd1_K<-read.table("ULWd1_K.txt", sep=";")
ULWd1_Mn<-read.table("ULWd1_Mn.txt", sep=";")
ULWd1_Si<-read.table("ULWd1_Si.txt", sep=";")

ULWd2_Ca<-read.table("ULWd2_Ca.txt", sep=";")
ULWd2_Fe<-read.table("ULWd2_Fe.txt", sep=";")
ULWd2_Ti<-read.table("ULWd2_Ti.txt", sep=";")
ULWd2_K<-read.table("ULWd2_K.txt", sep=";")
ULWd2_Mn<-read.table("ULWd2_Mn.txt", sep=";")
ULWd2_Si<-read.table("ULWd2_Si.txt", sep=";")

ULWf1_Ca<-read.table("ULWf1_Ca.txt", sep=";")
ULWf1_Fe<-read.table("ULWf1_Fe.txt", sep=";")
ULWf1_Ti<-read.table("ULWf1_Ti.txt", sep=";")
ULWf1_K<-read.table("ULWf1_K.txt", sep=";")
ULWf1_Mn<-read.table("ULWf1_Mn.txt", sep=";")
ULWf1_Si<-read.table("ULWf1_Si.txt", sep=";")


#reverse all datasets as measuring was done from top to bottom:
#(and rename colnames, as it's confusing if they're in reverse order)

ULWa1_Ca<-ULWa1_Ca[157:1,526:1]
ULWa1_Fe<-ULWa1_Fe[157:1,526:1]
ULWa1_K<-ULWa1_K[157:1,526:1]
ULWa1_Ti<-ULWa1_Ti[157:1,526:1]
ULWa1_Si<-ULWa1_Si[157:1,526:1]
ULWa1_Mn<-ULWa1_Mn[157:1,526:1]

colnames(ULWa1_Ca)<-c(1:526) 
rownames(ULWa1_Ca)<-c(1:157)
colnames(ULWa1_Fe)<-c(1:526) 
rownames(ULWa1_Fe)<-c(1:157)
colnames(ULWa1_K)<-c(1:526) 
rownames(ULWa1_K)<-c(1:157)
colnames(ULWa1_Ti)<-c(1:526) 
rownames(ULWa1_Ti)<-c(1:157)
colnames(ULWa1_Si)<-c(1:526) 
rownames(ULWa1_Si)<-c(1:157)
colnames(ULWa1_Mn)<-c(1:526) 
rownames(ULWa1_Mn)<-c(1:157)

ULWa2_Ca<-ULWa2_Ca[113:1,1400:1]
ULWa2_Fe<-ULWa2_Fe[113:1,1400:1]
ULWa2_K<-ULWa2_K[113:1,1400:1]
ULWa2_Ti<-ULWa2_Ti[113:1,1400:1]
ULWa2_Si<-ULWa2_Si[113:1,1400:1]
ULWa2_Mn<-ULWa2_Mn[113:1,1400:1]

colnames(ULWa2_Ca)<-c(1:1400) 
rownames(ULWa2_Ca)<-c(1:113)
colnames(ULWa2_Fe)<-c(1:1400) 
rownames(ULWa2_Fe)<-c(1:113)
colnames(ULWa2_K)<-c(1:1400) 
rownames(ULWa2_K)<-c(1:113)
colnames(ULWa2_Ti)<-c(1:1400) 
rownames(ULWa2_Ti)<-c(1:113)
colnames(ULWa2_Si)<-c(1:1400) 
rownames(ULWa2_Si)<-c(1:113)
colnames(ULWa2_Mn)<-c(1:1400) 
rownames(ULWa2_Mn)<-c(1:113)

ULWa3_Ca<-ULWa3_Ca[115:1,1288:1]
ULWa3_Fe<-ULWa3_Fe[115:1,1288:1]
ULWa3_K<-ULWa3_K[115:1,1288:1]
ULWa3_Ti<-ULWa3_Ti[115:1,1288:1]
ULWa3_Si<-ULWa3_Si[115:1,1288:1]
ULWa3_Mn<-ULWa3_Mn[115:1,1288:1]

colnames(ULWa3_Ca)<-c(1:1288) 
rownames(ULWa3_Ca)<-c(1:115)
colnames(ULWa3_Fe)<-c(1:1288) 
rownames(ULWa3_Fe)<-c(1:115)
colnames(ULWa3_K)<-c(1:1288) 
rownames(ULWa3_K)<-c(1:115)
colnames(ULWa3_Ti)<-c(1:1288) 
rownames(ULWa3_Ti)<-c(1:115)
colnames(ULWa3_Si)<-c(1:1288) 
rownames(ULWa3_Si)<-c(1:115)
colnames(ULWa3_Mn)<-c(1:1288) 
rownames(ULWa3_Mn)<-c(1:115)

ULWb1_Ca<-ULWb1_Ca[90:1,1440:1]
ULWb1_Fe<-ULWb1_Fe[90:1,1440:1]
ULWb1_K<-ULWb1_K[90:1,1440:1]
ULWb1_Ti<-ULWb1_Ti[90:1,1440:1]
ULWb1_Si<-ULWb1_Si[90:1,1440:1]
ULWb1_Mn<-ULWb1_Mn[90:1,1440:1]

colnames(ULWb1_Ca)<-c(1:1440) 
rownames(ULWb1_Ca)<-c(1:90)
colnames(ULWb1_Fe)<-c(1:1440) 
rownames(ULWb1_Fe)<-c(1:90)
colnames(ULWb1_K)<-c(1:1440) 
rownames(ULWb1_K)<-c(1:90)
colnames(ULWb1_Ti)<-c(1:1440) 
rownames(ULWb1_Ti)<-c(1:90)
colnames(ULWb1_Si)<-c(1:1440) 
rownames(ULWb1_Si)<-c(1:90)
colnames(ULWb1_Mn)<-c(1:1440) 
rownames(ULWb1_Mn)<-c(1:90)

ULWc1_Ca<-ULWc1_Ca[120:1,1260:1]
ULWc1_Fe<-ULWc1_Fe[120:1,1260:1]
ULWc1_K<-ULWc1_K[120:1,1260:1]
ULWc1_Ti<-ULWc1_Ti[120:1,1260:1]
ULWc1_Si<-ULWc1_Si[120:1,1260:1]
ULWc1_Mn<-ULWc1_Mn[120:1,1260:1]

colnames(ULWc1_Ca)<-c(1:1260) 
rownames(ULWc1_Ca)<-c(1:120)
colnames(ULWc1_Fe)<-c(1:1260) 
rownames(ULWc1_Fe)<-c(1:120)
colnames(ULWc1_K)<-c(1:1260) 
rownames(ULWc1_K)<-c(1:120)
colnames(ULWc1_Ti)<-c(1:1260) 
rownames(ULWc1_Ti)<-c(1:120)
colnames(ULWc1_Si)<-c(1:1260) 
rownames(ULWc1_Si)<-c(1:120)
colnames(ULWc1_Mn)<-c(1:1260) 
rownames(ULWc1_Mn)<-c(1:120)

ULWc2_Ca<-ULWc2_Ca[106:1,1338:1]
ULWc2_Fe<-ULWc2_Fe[106:1,1338:1]
ULWc2_K<-ULWc2_K[106:1,1338:1]
ULWc2_Ti<-ULWc2_Ti[106:1,1338:1]
ULWc2_Si<-ULWc2_Si[106:1,1338:1]
ULWc2_Mn<-ULWc2_Mn[106:1,1338:1]

colnames(ULWc2_Ca)<-c(1:1338) 
rownames(ULWc2_Ca)<-c(1:106)
colnames(ULWc2_Fe)<-c(1:1338) 
rownames(ULWc2_Fe)<-c(1:106)
colnames(ULWc2_K)<-c(1:1338) 
rownames(ULWc2_K)<-c(1:106)
colnames(ULWc2_Ti)<-c(1:1338) 
rownames(ULWc2_Ti)<-c(1:106)
colnames(ULWc2_Si)<-c(1:1338) 
rownames(ULWc2_Si)<-c(1:106)
colnames(ULWc2_Mn)<-c(1:1338) 
rownames(ULWc2_Mn)<-c(1:106)

ULWd1_Ca<-ULWd1_Ca[92:1,1158:1]
ULWd1_Fe<-ULWd1_Fe[92:1,1158:1]
ULWd1_K<-ULWd1_K[92:1,1158:1]
ULWd1_Ti<-ULWd1_Ti[92:1,1158:1]
ULWd1_Si<-ULWd1_Si[92:1,1158:1]
ULWd1_Mn<-ULWd1_Mn[92:1,1158:1]

colnames(ULWd1_Ca)<-c(1:1158) 
rownames(ULWd1_Ca)<-c(1:92)
colnames(ULWd1_Fe)<-c(1:1158) 
rownames(ULWd1_Fe)<-c(1:92)
colnames(ULWd1_K)<-c(1:1158) 
rownames(ULWd1_K)<-c(1:92)
colnames(ULWd1_Ti)<-c(1:1158) 
rownames(ULWd1_Ti)<-c(1:92)
colnames(ULWd1_Si)<-c(1:1158) 
rownames(ULWd1_Si)<-c(1:92)
colnames(ULWd1_Mn)<-c(1:1158) 
rownames(ULWd1_Mn)<-c(1:92)

ULWd2_Ca<-ULWd2_Ca[49:1,493:1]
ULWd2_Fe<-ULWd2_Fe[49:1,493:1]
ULWd2_K<-ULWd2_K[49:1,493:1]
ULWd2_Ti<-ULWd2_Ti[49:1,493:1]
ULWd2_Si<-ULWd2_Si[49:1,493:1]
ULWd2_Mn<-ULWd2_Mn[49:1,493:1]

colnames(ULWd2_Ca)<-c(1:493) 
rownames(ULWd2_Ca)<-c(1:49)
colnames(ULWd2_Fe)<-c(1:493) 
rownames(ULWd2_Fe)<-c(1:49)
colnames(ULWd2_K)<-c(1:493) 
rownames(ULWd2_K)<-c(1:49)
colnames(ULWd2_Ti)<-c(1:493) 
rownames(ULWd2_Ti)<-c(1:49)
colnames(ULWd2_Si)<-c(1:493) 
rownames(ULWd2_Si)<-c(1:49)
colnames(ULWd2_Mn)<-c(1:493) 
rownames(ULWd2_Mn)<-c(1:49)

ULWf1_Ca<-ULWf1_Ca[107:1,1586:1]
ULWf1_Fe<-ULWf1_Fe[107:1,1586:1]
ULWf1_K<-ULWf1_K[107:1,1586:1]
ULWf1_Ti<-ULWf1_Ti[107:1,1586:1]
ULWf1_Si<-ULWf1_Si[107:1,1586:1]
ULWf1_Mn<-ULWf1_Mn[107:1,1586:1]

colnames(ULWf1_Ca)<-c(1:1586) 
rownames(ULWf1_Ca)<-c(1:107)
colnames(ULWf1_Fe)<-c(1:1586) 
rownames(ULWf1_Fe)<-c(1:107)
colnames(ULWf1_K)<-c(1:1586) 
rownames(ULWf1_K)<-c(1:107)
colnames(ULWf1_Ti)<-c(1:1586) 
rownames(ULWf1_Ti)<-c(1:107)
colnames(ULWf1_Si)<-c(1:1586) 
rownames(ULWf1_Si)<-c(1:107)
colnames(ULWf1_Mn)<-c(1:1586) 
rownames(ULWf1_Mn)<-c(1:107)



#Plot as heatmaps:
heatmap(as.matrix(ULWa1_Ca), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWa1_Fe), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWa1_K), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWa1_Ti), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWa1_Si), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWa1_Mn), Rowv=NA, Colv=NA, asp=1)

heatmap(as.matrix(ULWa2_Ca), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWa2_Fe), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWa2_K), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWa2_Ti), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWa2_Si), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWa2_Mn), Rowv=NA, Colv=NA, asp=1)

heatmap(as.matrix(ULWa3_Ca), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWa3_Fe), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWa3_K), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWa3_Ti), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWa3_Si), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWa3_Mn), Rowv=NA, Colv=NA, asp=1)

heatmap(as.matrix(ULWb1_Ca), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWb1_Fe), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWb1_K), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWb1_Ti), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWb1_Si), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWb1_Mn), Rowv=NA, Colv=NA, asp=1)

heatmap(as.matrix(ULWc1_Ca), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWc1_Fe), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWc1_K), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWc1_Ti), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWc1_Si), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWc1_Mn), Rowv=NA, Colv=NA, asp=1)

heatmap(as.matrix(ULWc2_Ca), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWc2_Fe), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWc2_K), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWc2_Ti), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWc2_Si), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWc2_Mn), Rowv=NA, Colv=NA, asp=1)

heatmap(as.matrix(ULWd1_Ca), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWd1_Fe), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWd1_K), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWd1_Ti), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWd1_Si), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWd1_Mn), Rowv=NA, Colv=NA, asp=1)

heatmap(as.matrix(ULWd2_Ca), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWd2_Fe), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWd2_K), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWd2_Ti), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWd2_Si), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWd2_Mn), Rowv=NA, Colv=NA, asp=1)

heatmap(as.matrix(ULWf1_Ca), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWf1_Fe), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWf1_K), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWf1_Ti), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWf1_Si), Rowv=NA, Colv=NA, asp=1)
heatmap(as.matrix(ULWf1_Mn), Rowv=NA, Colv=NA, asp=1)


#save as csvs:
write.csv(ULWa1_Ca, "ULWa1_Ca.csv")
write.csv(ULWa1_Fe, "ULWa1_Fe.csv")
write.csv(ULWa1_K, "ULWa1_K.csv")
write.csv(ULWa1_Ti, "ULWa1_Ti.csv")
write.csv(ULWa1_Si, "ULWa1_Si.csv")
write.csv(ULWa1_Mn, "ULWa1_Mn.csv")

write.csv(ULWa2_Ca, "ULWa2_Ca.csv")
write.csv(ULWa2_Fe, "ULWa2_Fe.csv")
write.csv(ULWa2_K, "ULWa2_K.csv")
write.csv(ULWa2_Ti, "ULWa2_Ti.csv")
write.csv(ULWa2_Si, "ULWa2_Si.csv")
write.csv(ULWa2_Mn, "ULWa2_Mn.csv")

write.csv(ULWa3_Ca, "ULWa3_Ca.csv")
write.csv(ULWa3_Fe, "ULWa3_Fe.csv")
write.csv(ULWa3_K, "ULWa3_K.csv")
write.csv(ULWa3_Ti, "ULWa3_Ti.csv")
write.csv(ULWa3_Si, "ULWa3_Si.csv")
write.csv(ULWa3_Mn, "ULWa3_Mn.csv")

write.csv(ULWb1_Ca, "ULWb1_Ca.csv")
write.csv(ULWb1_Fe, "ULWb1_Fe.csv")
write.csv(ULWb1_K, "ULWb1_K.csv")
write.csv(ULWb1_Ti, "ULWb1_Ti.csv")
write.csv(ULWb1_Si, "ULWb1_Si.csv")
write.csv(ULWb1_Mn, "ULWb1_Mn.csv")

write.csv(ULWc1_Ca, "ULWc1_Ca.csv")
write.csv(ULWc1_Fe, "ULWc1_Fe.csv")
write.csv(ULWc1_K, "ULWc1_K.csv")
write.csv(ULWc1_Ti, "ULWc1_Ti.csv")
write.csv(ULWc1_Si, "ULWc1_Si.csv")
write.csv(ULWc1_Mn, "ULWc1_Mn.csv")

write.csv(ULWc2_Ca, "ULWc2_Ca.csv")
write.csv(ULWc2_Fe, "ULWc2_Fe.csv")
write.csv(ULWc2_K, "ULWc2_K.csv")
write.csv(ULWc2_Ti, "ULWc2_Ti.csv")
write.csv(ULWc2_Si, "ULWc2_Si.csv")
write.csv(ULWc2_Mn, "ULWc2_Mn.csv")

write.csv(ULWd1_Ca, "ULWd1_Ca.csv")
write.csv(ULWd1_Fe, "ULWd1_Fe.csv")
write.csv(ULWd1_K, "ULWd1_K.csv")
write.csv(ULWd1_Ti, "ULWd1_Ti.csv")
write.csv(ULWd1_Si, "ULWd1_Si.csv")
write.csv(ULWd1_Mn, "ULWd1_Mn.csv")

write.csv(ULWd2_Ca, "ULWd2_Ca.csv")
write.csv(ULWd2_Fe, "ULWd2_Fe.csv")
write.csv(ULWd2_K, "ULWd2_K.csv")
write.csv(ULWd2_Ti, "ULWd2_Ti.csv")
write.csv(ULWd2_Si, "ULWd2_Si.csv")
write.csv(ULWd2_Mn, "ULWd2_Mn.csv")

write.csv(ULWf1_Ca, "ULWf1_Ca.csv")
write.csv(ULWf1_Fe, "ULWf1_Fe.csv")
write.csv(ULWf1_K, "ULWf1_K.csv")
write.csv(ULWf1_Ti, "ULWf1_Ti.csv")
write.csv(ULWf1_Si, "ULWf1_Si.csv")
write.csv(ULWf1_Mn, "ULWf1_Mn.csv")




#No need for removing non-rock areas.

#-------------------------Composite maps-------------------------------------------
#Read in data; don't have to re-run above code
ULWa1_Ca<-read.csv("ULWa1_Ca.csv")
ULWa1_Fe<-read.csv("ULWa1_Fe.csv")
ULWa1_K<-read.csv("ULWa1_K.csv")
ULWa1_Ti<-read.csv("ULWa1_Ti.csv")
ULWa1_Si<-read.csv("ULWa1_Si.csv")
ULWa1_Mn<-read.csv("ULWa1_Mn.csv")

ULWa2_Ca<-read.csv("ULWa2_Ca.csv")
ULWa2_Fe<-read.csv("ULWa2_Fe.csv")
ULWa2_K<-read.csv("ULWa2_K.csv")
ULWa2_Ti<-read.csv("ULWa2_Ti.csv")
ULWa2_Si<-read.csv("ULWa2_Si.csv")
ULWa2_Mn<-read.csv("ULWa2_Mn.csv")

ULWa3_Ca<-read.csv("ULWa3_Ca.csv")
ULWa3_Fe<-read.csv("ULWa3_Fe.csv")
ULWa3_K<-read.csv("ULWa3_K.csv")
ULWa3_Ti<-read.csv("ULWa3_Ti.csv")
ULWa3_Si<-read.csv("ULWa3_Si.csv")
ULWa3_Mn<-read.csv("ULWa3_Mn.csv")

ULWb1_Ca<-read.csv("ULWb1_Ca.csv")
ULWb1_Fe<-read.csv("ULWb1_Fe.csv")
ULWb1_K<-read.csv("ULWb1_K.csv")
ULWb1_Ti<-read.csv("ULWb1_Ti.csv")
ULWb1_Si<-read.csv("ULWb1_Si.csv")
ULWb1_Mn<-read.csv("ULWb1_Mn.csv")

ULWc1_Ca<-read.csv("ULWc1_Ca.csv")
ULWc1_Fe<-read.csv("ULWc1_Fe.csv")
ULWc1_K<-read.csv("ULWc1_K.csv")
ULWc1_Ti<-read.csv("ULWc1_Ti.csv")
ULWc1_Si<-read.csv("ULWc1_Si.csv")
ULWc1_Mn<-read.csv("ULWc1_Mn.csv")

ULWc2_Ca<-read.csv("ULWc2_Ca.csv")
ULWc2_Fe<-read.csv("ULWc2_Fe.csv")
ULWc2_K<-read.csv("ULWc2_K.csv")
ULWc2_Ti<-read.csv("ULWc2_Ti.csv")
ULWc2_Si<-read.csv("ULWc2_Si.csv")
ULWc2_Mn<-read.csv("ULWc2_Mn.csv")

ULWd1_Ca<-read.csv("ULWd1_Ca.csv")
ULWd1_Fe<-read.csv("ULWd1_Fe.csv")
ULWd1_K<-read.csv("ULWd1_K.csv")
ULWd1_Ti<-read.csv("ULWd1_Ti.csv")
ULWd1_Si<-read.csv("ULWd1_Si.csv")
ULWd1_Mn<-read.csv("ULWd1_Mn.csv")

ULWd2_Ca<-read.csv("ULWd2_Ca.csv")
ULWd2_Fe<-read.csv("ULWd2_Fe.csv")
ULWd2_K<-read.csv("ULWd2_K.csv")
ULWd2_Ti<-read.csv("ULWd2_Ti.csv")
ULWd2_Si<-read.csv("ULWd2_Si.csv")
ULWd2_Mn<-read.csv("ULWd2_Mn.csv")

ULWf1_Ca<-read.csv("ULWf1_Ca.csv")
ULWf1_Fe<-read.csv("ULWf1_Fe.csv")
ULWf1_K<-read.csv("ULWf1_K.csv")
ULWf1_Ti<-read.csv("ULWf1_Ti.csv")
ULWf1_Si<-read.csv("ULWf1_Si.csv")
ULWf1_Mn<-read.csv("ULWf1_Mn.csv")

#---------------------------Formatting---------------------------------------------
#Remove row index:
ULWa1_Ca<-ULWa1_Ca[,2:ncol(ULWa1_Ca)]
ULWa1_Fe<-ULWa1_Fe[,2:ncol(ULWa1_Fe)]
ULWa1_K<-ULWa1_K[,2:ncol(ULWa1_K)]
ULWa1_Ti<-ULWa1_Ti[,2:ncol(ULWa1_Ti)]
ULWa1_Si<-ULWa1_Si[,2:ncol(ULWa1_Si)]
ULWa1_Mn<-ULWa1_Mn[,2:ncol(ULWa1_Mn)]

ULWa2_Ca<-ULWa2_Ca[,2:ncol(ULWa2_Ca)]
ULWa2_Fe<-ULWa2_Fe[,2:ncol(ULWa2_Fe)]
ULWa2_K<-ULWa2_K[,2:ncol(ULWa2_K)]
ULWa2_Ti<-ULWa2_Ti[,2:ncol(ULWa2_Ti)]
ULWa2_Si<-ULWa2_Si[,2:ncol(ULWa2_Si)]
ULWa2_Mn<-ULWa2_Mn[,2:ncol(ULWa2_Mn)]

ULWa3_Ca<-ULWa3_Ca[,2:ncol(ULWa3_Ca)]
ULWa3_Fe<-ULWa3_Fe[,2:ncol(ULWa3_Fe)]
ULWa3_K<-ULWa3_K[,2:ncol(ULWa3_K)]
ULWa3_Ti<-ULWa3_Ti[,2:ncol(ULWa3_Ti)]
ULWa3_Si<-ULWa3_Si[,2:ncol(ULWa3_Si)]
ULWa3_Mn<-ULWa3_Mn[,2:ncol(ULWa3_Mn)]

ULWb1_Ca<-ULWb1_Ca[,2:ncol(ULWb1_Ca)]
ULWb1_Fe<-ULWb1_Fe[,2:ncol(ULWb1_Fe)]
ULWb1_K<-ULWb1_K[,2:ncol(ULWb1_K)]
ULWb1_Ti<-ULWb1_Ti[,2:ncol(ULWb1_Ti)]
ULWb1_Si<-ULWb1_Si[,2:ncol(ULWb1_Si)]
ULWb1_Mn<-ULWb1_Mn[,2:ncol(ULWb1_Mn)]

ULWc1_Ca<-ULWc1_Ca[,2:ncol(ULWc1_Ca)]
ULWc1_Fe<-ULWc1_Fe[,2:ncol(ULWc1_Fe)]
ULWc1_K<-ULWc1_K[,2:ncol(ULWc1_K)]
ULWc1_Ti<-ULWc1_Ti[,2:ncol(ULWc1_Ti)]
ULWc1_Si<-ULWc1_Si[,2:ncol(ULWc1_Si)]
ULWc1_Mn<-ULWc1_Mn[,2:ncol(ULWc1_Mn)]

ULWc2_Ca<-ULWc2_Ca[,2:ncol(ULWc2_Ca)]
ULWc2_Fe<-ULWc2_Fe[,2:ncol(ULWc2_Fe)]
ULWc2_K<-ULWc2_K[,2:ncol(ULWc2_K)]
ULWc2_Ti<-ULWc2_Ti[,2:ncol(ULWc2_Ti)]
ULWc2_Si<-ULWc2_Si[,2:ncol(ULWc2_Si)]
ULWc2_Mn<-ULWc2_Mn[,2:ncol(ULWc2_Mn)]

ULWd1_Ca<-ULWd1_Ca[,2:ncol(ULWd1_Ca)]
ULWd1_Fe<-ULWd1_Fe[,2:ncol(ULWd1_Fe)]
ULWd1_K<-ULWd1_K[,2:ncol(ULWd1_K)]
ULWd1_Ti<-ULWd1_Ti[,2:ncol(ULWd1_Ti)]
ULWd1_Si<-ULWd1_Si[,2:ncol(ULWd1_Si)]
ULWd1_Mn<-ULWd1_Mn[,2:ncol(ULWd1_Mn)]

ULWd2_Ca<-ULWd2_Ca[,2:ncol(ULWd2_Ca)]
ULWd2_Fe<-ULWd2_Fe[,2:ncol(ULWd2_Fe)]
ULWd2_K<-ULWd2_K[,2:ncol(ULWd2_K)]
ULWd2_Ti<-ULWd2_Ti[,2:ncol(ULWd2_Ti)]
ULWd2_Si<-ULWd2_Si[,2:ncol(ULWd2_Si)]
ULWd2_Mn<-ULWd2_Mn[,2:ncol(ULWd2_Mn)]

ULWf1_Ca<-ULWf1_Ca[,2:ncol(ULWf1_Ca)]
ULWf1_Fe<-ULWf1_Fe[,2:ncol(ULWf1_Fe)]
ULWf1_K<-ULWf1_K[,2:ncol(ULWf1_K)]
ULWf1_Ti<-ULWf1_Ti[,2:ncol(ULWf1_Ti)]
ULWf1_Si<-ULWf1_Si[,2:ncol(ULWf1_Si)]
ULWf1_Mn<-ULWf1_Mn[,2:ncol(ULWf1_Mn)]


#----------------------------determine overlap-------------------------------------
#--------------------------Ca----------------------------------------------------
#ULWa1-ULWa2:
#directly on top of each other, no overlap.
#526+1400=1926 total length
ULWa1_Ca_padded<-padzeros(data=ULWa1_Ca, nzeros=(1926-526), side="right")
ULWa2_Ca_padded<-padzeros(data=ULWa2_Ca, nzeros=(1926-1400), side="left")
colnames(ULWa1_Ca_padded)<-c(1:1926)
colnames(ULWa2_Ca_padded)<-c(1:1926)
ULWa12_Ca_padded<-rbind(ULWa2_Ca_padded, ULWa1_Ca_padded)
#heatmap(as.matrix(ULWa12_Ca_padded), Rowv=NA, Colv=NA, asp=1) 

#ULWa2-ULWb1:
#complete overlap, ends at same spot as a2. 
#total length = 1926
ULWb1_Ca_padded<-padzeros(data=ULWb1_Ca, nzeros=(1926-1440), side="left")
colnames(ULWb1_Ca_padded)<-c(1:1926)
ULWa12b_Ca_padded<-rbind(ULWa12_Ca_padded, ULWb1_Ca_padded)
#heatmap(as.matrix(ULWa12b_Ca_padded), Rowv=NA, Colv=NA, asp=1) 

#UlWa2-ULWa3:
#directly on top of each other, no overlap.
#1926+1288=3214 total length
ULWa12b_Ca_padded<-padzeros(data=ULWa12b_Ca_padded, nzeros=(3214-1926), side="right")
ULWa3_Ca_padded<-padzeros(data=ULWa3_Ca, nzeros=(3214-1288), side="left")
colnames(ULWa12b_Ca_padded)<-c(1:3214)
colnames(ULWa3_Ca_padded)<-c(1:3214)
ULWa12b3_Ca_padded<-rbind(ULWa3_Ca_padded, ULWa12b_Ca_padded)
#heatmap(as.matrix(ULWa12b3_Ca_padded), Rowv=NA, Colv=NA, asp=1)

#ULWa3-ULWc1:
#372 overlap
#(3214+1260)-372=4102 total length
ULWa12b3_Ca_padded<-padzeros(data=ULWa12b3_Ca_padded, nzeros=(4102-3214), side="right")
ULWc1_Ca_padded<-padzeros(data=ULWc1_Ca, nzeros=(4102-1260), side="left")
colnames(ULWa12b3_Ca_padded)<-c(1:4102)
colnames(ULWc1_Ca_padded)<-c(1:4102)
ULWa12b3c1_Ca_padded<-rbind(ULWc1_Ca_padded, ULWa12b3_Ca_padded)
#heatmap(as.matrix(ULWa12b3c1_Ca_padded), Rowv=NA, Colv=NA, asp=1)

#ULWc1-ULWc2:
#directly on top of each other, no overlap.
#1338+4102=5440 total length
ULWa12b3c1_Ca_padded<-padzeros(data=ULWa12b3c1_Ca_padded, nzeros=(5440-4102), side="right")
ULWc2_Ca_padded<-padzeros(data=ULWc2_Ca, nzeros=(5440-1338), side="left")
colnames(ULWa12b3c1_Ca_padded)<-c(1:5440)
colnames(ULWc2_Ca_padded)<-c(1:5440)
ULWa12b3c12_Ca_padded<-rbind(ULWc2_Ca_padded, ULWa12b3c1_Ca_padded)
#heatmap(as.matrix(ULWa12b3c12_Ca_padded), Rowv=NA, Colv=NA, asp=1)

#ULWc2-ULWd1:
#complete overlap. 5440 total length
#d1 ends 318 before c2, total d1 length is 1158
#so 318 padding on right, 5440-1158-318=3964 padding on left
ULWd1_Ca_padded<-padzeros(ULWd1_Ca, nzeros=318, side="right")
ULWd1_Ca_padded<-padzeros(ULWd1_Ca_padded, nzeros=3964, side="left")
colnames(ULWd1_Ca_padded)<-c(1:5440)
ULWa12b3c12d1_Ca_padded<-rbind(ULWa12b3c12_Ca_padded, ULWd1_Ca_padded)
#heatmap(as.matrix(ULWa12b3c12d1_Ca_padded), Rowv=NA, Colv=NA, asp=1) 

#ULWd1-ULWd2:
#263 overlap between c2 and d2.
#(5440+493)-263=5670 total lenght.
ULWa12b3c12d1_Ca_padded<-padzeros(data=ULWa12b3c12d1_Ca_padded, nzeros=(5670-5440), side="right")
ULWd2_Ca_padded<-padzeros(data=ULWd2_Ca, nzeros=(5670-493), side="left")
colnames(ULWa12b3c12d1_Ca_padded)<-c(1:5670)
colnames(ULWd2_Ca_padded)<-c(1:5670)
ULWa12b3c12d12_Ca_padded<-rbind(ULWd2_Ca_padded, ULWa12b3c12d1_Ca_padded)
#heatmap(as.matrix(ULWa12b3c12d12_Ca_padded), Rowv=NA, Colv=NA, asp=1)

#ULWd2-ULWf1:
#130 overlap. 
#(5670+1586)-130=7126 total length
ULWa12b3c12d12_Ca_padded<-padzeros(data=ULWa12b3c12d12_Ca_padded, nzeros=(7126-5670), side="right")
ULWf1_Ca_padded<-padzeros(data=ULWf1_Ca, nzeros=(7126-1586), side="left")
colnames(ULWa12b3c12d12_Ca_padded)<-c(1:7126)
colnames(ULWf1_Ca_padded)<-c(1:7126)
ULWa12b3c12d12f1_Ca_padded<-rbind(ULWf1_Ca_padded, ULWa12b3c12d12_Ca_padded)
ULWa12b3c12d12f1_Ca_padded[ULWa12b3c12d12f1_Ca_padded==0]<-NA
heatmap(as.matrix(ULWa12b3c12d12f1_Ca_padded), Rowv=NA, Colv=NA, asp=1)





#--------------------------Fe----------------------------------------------------
#ULWa1-ULWa2:
#directly on top of each other, no overlap.
#526+1400=1926 total length
ULWa1_Fe_padded<-padzeros(data=ULWa1_Fe, nzeros=(1926-526), side="right")
ULWa2_Fe_padded<-padzeros(data=ULWa2_Fe, nzeros=(1926-1400), side="left")
colnames(ULWa1_Fe_padded)<-c(1:1926)
colnames(ULWa2_Fe_padded)<-c(1:1926)
ULWa12_Fe_padded<-rbind(ULWa2_Fe_padded, ULWa1_Fe_padded)
#heatmap(as.matrix(ULWa12_Fe_padded), Rowv=NA, Colv=NA, asp=1) 

#ULWa2-ULWb1:
#complete overlap, ends at same spot as a2. 
#total length = 1926
ULWb1_Fe_padded<-padzeros(data=ULWb1_Fe, nzeros=(1926-1440), side="left")
colnames(ULWb1_Fe_padded)<-c(1:1926)
ULWa12b_Fe_padded<-rbind(ULWa12_Fe_padded, ULWb1_Fe_padded)
#heatmap(as.matrix(ULWa12b_Fe_padded), Rowv=NA, Colv=NA, asp=1) 

#UlWa2-ULWa3:
#directly on top of each other, no overlap.
#1926+1288=3214 total length
ULWa12b_Fe_padded<-padzeros(data=ULWa12b_Fe_padded, nzeros=(3214-1926), side="right")
ULWa3_Fe_padded<-padzeros(data=ULWa3_Fe, nzeros=(3214-1288), side="left")
colnames(ULWa12b_Fe_padded)<-c(1:3214)
colnames(ULWa3_Fe_padded)<-c(1:3214)
ULWa12b3_Fe_padded<-rbind(ULWa3_Fe_padded, ULWa12b_Fe_padded)
#heatmap(as.matrix(ULWa12b3_Fe_padded), Rowv=NA, Colv=NA, asp=1)

#ULWa3-ULWc1:
#372 overlap
#(3214+1260)-372=4102 total length
ULWa12b3_Fe_padded<-padzeros(data=ULWa12b3_Fe_padded, nzeros=(4102-3214), side="right")
ULWc1_Fe_padded<-padzeros(data=ULWc1_Fe, nzeros=(4102-1260), side="left")
colnames(ULWa12b3_Fe_padded)<-c(1:4102)
colnames(ULWc1_Fe_padded)<-c(1:4102)
ULWa12b3c1_Fe_padded<-rbind(ULWc1_Fe_padded, ULWa12b3_Fe_padded)
#heatmap(as.matrix(ULWa12b3c1_Fe_padded), Rowv=NA, Colv=NA, asp=1)

#ULWc1-ULWc2:
#directly on top of each other, no overlap.
#1338+4102=5440 total length
ULWa12b3c1_Fe_padded<-padzeros(data=ULWa12b3c1_Fe_padded, nzeros=(5440-4102), side="right")
ULWc2_Fe_padded<-padzeros(data=ULWc2_Fe, nzeros=(5440-1338), side="left")
colnames(ULWa12b3c1_Fe_padded)<-c(1:5440)
colnames(ULWc2_Fe_padded)<-c(1:5440)
ULWa12b3c12_Fe_padded<-rbind(ULWc2_Fe_padded, ULWa12b3c1_Fe_padded)
#heatmap(as.matrix(ULWa12b3c12_Fe_padded), Rowv=NA, Colv=NA, asp=1)

#ULWc2-ULWd1:
#complete overlap. 5440 total length
#d1 ends 318 before c2, total d1 length is 1158
#so 318 padding on right, 5440-1158-318=3964 padding on left
ULWd1_Fe_padded<-padzeros(ULWd1_Fe, nzeros=318, side="right")
ULWd1_Fe_padded<-padzeros(ULWd1_Fe_padded, nzeros=3964, side="left")
colnames(ULWd1_Fe_padded)<-c(1:5440)
ULWa12b3c12d1_Fe_padded<-rbind(ULWa12b3c12_Fe_padded, ULWd1_Fe_padded)
#heatmap(as.matrix(ULWa12b3c12d1_Fe_padded), Rowv=NA, Colv=NA, asp=1) 

#ULWd1-ULWd2:
#263 overlap between c2 and d2.
#(5440+493)-263=5670 total lenght.
ULWa12b3c12d1_Fe_padded<-padzeros(data=ULWa12b3c12d1_Fe_padded, nzeros=(5670-5440), side="right")
ULWd2_Fe_padded<-padzeros(data=ULWd2_Fe, nzeros=(5670-493), side="left")
colnames(ULWa12b3c12d1_Fe_padded)<-c(1:5670)
colnames(ULWd2_Fe_padded)<-c(1:5670)
ULWa12b3c12d12_Fe_padded<-rbind(ULWd2_Fe_padded, ULWa12b3c12d1_Fe_padded)
#heatmap(as.matrix(ULWa12b3c12d12_Fe_padded), Rowv=NA, Colv=NA, asp=1)

#ULWd2-ULWf1:
#130 overlap. 
#(5670+1586)-130=7126 total length
ULWa12b3c12d12_Fe_padded<-padzeros(data=ULWa12b3c12d12_Fe_padded, nzeros=(7126-5670), side="right")
ULWf1_Fe_padded<-padzeros(data=ULWf1_Fe, nzeros=(7126-1586), side="left")
colnames(ULWa12b3c12d12_Fe_padded)<-c(1:7126)
colnames(ULWf1_Fe_padded)<-c(1:7126)
ULWa12b3c12d12f1_Fe_padded<-rbind(ULWf1_Fe_padded, ULWa12b3c12d12_Fe_padded)
ULWa12b3c12d12f1_Fe_padded[ULWa12b3c12d12f1_Fe_padded==0]<-NA
heatmap(as.matrix(ULWa12b3c12d12f1_Fe_padded), Rowv=NA, Colv=NA, asp=1)




#--------------------------K----------------------------------------------------
#ULWa1-ULWa2:
#directly on top of each other, no overlap.
#526+1400=1926 total length
ULWa1_K_padded<-padzeros(data=ULWa1_K, nzeros=(1926-526), side="right")
ULWa2_K_padded<-padzeros(data=ULWa2_K, nzeros=(1926-1400), side="left")
colnames(ULWa1_K_padded)<-c(1:1926)
colnames(ULWa2_K_padded)<-c(1:1926)
ULWa12_K_padded<-rbind(ULWa2_K_padded, ULWa1_K_padded)
#heatmap(as.matrix(ULWa12_K_padded), Rowv=NA, Colv=NA, asp=1) 

#ULWa2-ULWb1:
#complete overlap, ends at same spot as a2. 
#total length = 1926
ULWb1_K_padded<-padzeros(data=ULWb1_K, nzeros=(1926-1440), side="left")
colnames(ULWb1_K_padded)<-c(1:1926)
ULWa12b_K_padded<-rbind(ULWa12_K_padded, ULWb1_K_padded)
#heatmap(as.matrix(ULWa12b_K_padded), Rowv=NA, Colv=NA, asp=1) 

#UlWa2-ULWa3:
#directly on top of each other, no overlap.
#1926+1288=3214 total length
ULWa12b_K_padded<-padzeros(data=ULWa12b_K_padded, nzeros=(3214-1926), side="right")
ULWa3_K_padded<-padzeros(data=ULWa3_K, nzeros=(3214-1288), side="left")
colnames(ULWa12b_K_padded)<-c(1:3214)
colnames(ULWa3_K_padded)<-c(1:3214)
ULWa12b3_K_padded<-rbind(ULWa3_K_padded, ULWa12b_K_padded)
#heatmap(as.matrix(ULWa12b3_K_padded), Rowv=NA, Colv=NA, asp=1)

#ULWa3-ULWc1:
#372 overlap
#(3214+1260)-372=4102 total length
ULWa12b3_K_padded<-padzeros(data=ULWa12b3_K_padded, nzeros=(4102-3214), side="right")
ULWc1_K_padded<-padzeros(data=ULWc1_K, nzeros=(4102-1260), side="left")
colnames(ULWa12b3_K_padded)<-c(1:4102)
colnames(ULWc1_K_padded)<-c(1:4102)
ULWa12b3c1_K_padded<-rbind(ULWc1_K_padded, ULWa12b3_K_padded)
#heatmap(as.matrix(ULWa12b3c1_K_padded), Rowv=NA, Colv=NA, asp=1)

#ULWc1-ULWc2:
#directly on top of each other, no overlap.
#1338+4102=5440 total length
ULWa12b3c1_K_padded<-padzeros(data=ULWa12b3c1_K_padded, nzeros=(5440-4102), side="right")
ULWc2_K_padded<-padzeros(data=ULWc2_K, nzeros=(5440-1338), side="left")
colnames(ULWa12b3c1_K_padded)<-c(1:5440)
colnames(ULWc2_K_padded)<-c(1:5440)
ULWa12b3c12_K_padded<-rbind(ULWc2_K_padded, ULWa12b3c1_K_padded)
#heatmap(as.matrix(ULWa12b3c12_K_padded), Rowv=NA, Colv=NA, asp=1)

#ULWc2-ULWd1:
#complete overlap. 5440 total length
#d1 ends 318 before c2, total d1 length is 1158
#so 318 padding on right, 5440-1158-318=3964 padding on left
ULWd1_K_padded<-padzeros(ULWd1_K, nzeros=318, side="right")
ULWd1_K_padded<-padzeros(ULWd1_K_padded, nzeros=3964, side="left")
colnames(ULWd1_K_padded)<-c(1:5440)
ULWa12b3c12d1_K_padded<-rbind(ULWa12b3c12_K_padded, ULWd1_K_padded)
#heatmap(as.matrix(ULWa12b3c12d1_K_padded), Rowv=NA, Colv=NA, asp=1) 

#ULWd1-ULWd2:
#263 overlap between c2 and d2.
#(5440+493)-263=5670 total lenght.
ULWa12b3c12d1_K_padded<-padzeros(data=ULWa12b3c12d1_K_padded, nzeros=(5670-5440), side="right")
ULWd2_K_padded<-padzeros(data=ULWd2_K, nzeros=(5670-493), side="left")
colnames(ULWa12b3c12d1_K_padded)<-c(1:5670)
colnames(ULWd2_K_padded)<-c(1:5670)
ULWa12b3c12d12_K_padded<-rbind(ULWd2_K_padded, ULWa12b3c12d1_K_padded)
#heatmap(as.matrix(ULWa12b3c12d12_K_padded), Rowv=NA, Colv=NA, asp=1)

#ULWd2-ULWf1:
#130 overlap. 
#(5670+1586)-130=7126 total length
ULWa12b3c12d12_K_padded<-padzeros(data=ULWa12b3c12d12_K_padded, nzeros=(7126-5670), side="right")
ULWf1_K_padded<-padzeros(data=ULWf1_K, nzeros=(7126-1586), side="left")
colnames(ULWa12b3c12d12_K_padded)<-c(1:7126)
colnames(ULWf1_K_padded)<-c(1:7126)
ULWa12b3c12d12f1_K_padded<-rbind(ULWf1_K_padded, ULWa12b3c12d12_K_padded)
ULWa12b3c12d12f1_K_padded[ULWa12b3c12d12f1_K_padded==0]<-NA
heatmap(as.matrix(ULWa12b3c12d12f1_K_padded), Rowv=NA, Colv=NA, asp=1)





#--------------------------Ti----------------------------------------------------
#ULWa1-ULWa2:
#directly on top of each other, no overlap.
#526+1400=1926 total length
ULWa1_Ti_padded<-padzeros(data=ULWa1_Ti, nzeros=(1926-526), side="right")
ULWa2_Ti_padded<-padzeros(data=ULWa2_Ti, nzeros=(1926-1400), side="left")
colnames(ULWa1_Ti_padded)<-c(1:1926)
colnames(ULWa2_Ti_padded)<-c(1:1926)
ULWa12_Ti_padded<-rbind(ULWa2_Ti_padded, ULWa1_Ti_padded)
#heatmap(as.matrix(ULWa12_Ti_padded), Rowv=NA, Colv=NA, asp=1) 

#ULWa2-ULWb1:
#complete overlap, ends at same spot as a2. 
#total length = 1926
ULWb1_Ti_padded<-padzeros(data=ULWb1_Ti, nzeros=(1926-1440), side="left")
colnames(ULWb1_Ti_padded)<-c(1:1926)
ULWa12b_Ti_padded<-rbind(ULWa12_Ti_padded, ULWb1_Ti_padded)
#heatmap(as.matrix(ULWa12b_Ti_padded), Rowv=NA, Colv=NA, asp=1) 

#UlWa2-ULWa3:
#directly on top of each other, no overlap.
#1926+1288=3214 total length
ULWa12b_Ti_padded<-padzeros(data=ULWa12b_Ti_padded, nzeros=(3214-1926), side="right")
ULWa3_Ti_padded<-padzeros(data=ULWa3_Ti, nzeros=(3214-1288), side="left")
colnames(ULWa12b_Ti_padded)<-c(1:3214)
colnames(ULWa3_Ti_padded)<-c(1:3214)
ULWa12b3_Ti_padded<-rbind(ULWa3_Ti_padded, ULWa12b_Ti_padded)
#heatmap(as.matrix(ULWa12b3_Ti_padded), Rowv=NA, Colv=NA, asp=1)

#ULWa3-ULWc1:
#372 overlap
#(3214+1260)-372=4102 total length
ULWa12b3_Ti_padded<-padzeros(data=ULWa12b3_Ti_padded, nzeros=(4102-3214), side="right")
ULWc1_Ti_padded<-padzeros(data=ULWc1_Ti, nzeros=(4102-1260), side="left")
colnames(ULWa12b3_Ti_padded)<-c(1:4102)
colnames(ULWc1_Ti_padded)<-c(1:4102)
ULWa12b3c1_Ti_padded<-rbind(ULWc1_Ti_padded, ULWa12b3_Ti_padded)
#heatmap(as.matrix(ULWa12b3c1_Ti_padded), Rowv=NA, Colv=NA, asp=1)

#ULWc1-ULWc2:
#directly on top of each other, no overlap.
#1338+4102=5440 total length
ULWa12b3c1_Ti_padded<-padzeros(data=ULWa12b3c1_Ti_padded, nzeros=(5440-4102), side="right")
ULWc2_Ti_padded<-padzeros(data=ULWc2_Ti, nzeros=(5440-1338), side="left")
colnames(ULWa12b3c1_Ti_padded)<-c(1:5440)
colnames(ULWc2_Ti_padded)<-c(1:5440)
ULWa12b3c12_Ti_padded<-rbind(ULWc2_Ti_padded, ULWa12b3c1_Ti_padded)
#heatmap(as.matrix(ULWa12b3c12_Ti_padded), Rowv=NA, Colv=NA, asp=1)

#ULWc2-ULWd1:
#complete overlap. 5440 total length
#d1 ends 318 before c2, total d1 length is 1158
#so 318 padding on right, 5440-1158-318=3964 padding on left
ULWd1_Ti_padded<-padzeros(ULWd1_Ti, nzeros=318, side="right")
ULWd1_Ti_padded<-padzeros(ULWd1_Ti_padded, nzeros=3964, side="left")
colnames(ULWd1_Ti_padded)<-c(1:5440)
ULWa12b3c12d1_Ti_padded<-rbind(ULWa12b3c12_Ti_padded, ULWd1_Ti_padded)
#heatmap(as.matrix(ULWa12b3c12d1_Ti_padded), Rowv=NA, Colv=NA, asp=1) 

#ULWd1-ULWd2:
#263 overlap between c2 and d2.
#(5440+493)-263=5670 total lenght.
ULWa12b3c12d1_Ti_padded<-padzeros(data=ULWa12b3c12d1_Ti_padded, nzeros=(5670-5440), side="right")
ULWd2_Ti_padded<-padzeros(data=ULWd2_Ti, nzeros=(5670-493), side="left")
colnames(ULWa12b3c12d1_Ti_padded)<-c(1:5670)
colnames(ULWd2_Ti_padded)<-c(1:5670)
ULWa12b3c12d12_Ti_padded<-rbind(ULWd2_Ti_padded, ULWa12b3c12d1_Ti_padded)
#heatmap(as.matrix(ULWa12b3c12d12_Ti_padded), Rowv=NA, Colv=NA, asp=1)

#ULWd2-ULWf1:
#130 overlap. 
#(5670+1586)-130=7126 total length
ULWa12b3c12d12_Ti_padded<-padzeros(data=ULWa12b3c12d12_Ti_padded, nzeros=(7126-5670), side="right")
ULWf1_Ti_padded<-padzeros(data=ULWf1_Ti, nzeros=(7126-1586), side="left")
colnames(ULWa12b3c12d12_Ti_padded)<-c(1:7126)
colnames(ULWf1_Ti_padded)<-c(1:7126)
ULWa12b3c12d12f1_Ti_padded<-rbind(ULWf1_Ti_padded, ULWa12b3c12d12_Ti_padded)
ULWa12b3c12d12f1_Ti_padded[ULWa12b3c12d12f1_Ti_padded==0]<-NA
heatmap(as.matrix(ULWa12b3c12d12f1_Ti_padded), Rowv=NA, Colv=NA, asp=1)





#--------------------------Si----------------------------------------------------
#ULWa1-ULWa2:
#directly on top of each other, no overlap.
#526+1400=1926 total length
ULWa1_Si_padded<-padzeros(data=ULWa1_Si, nzeros=(1926-526), side="right")
ULWa2_Si_padded<-padzeros(data=ULWa2_Si, nzeros=(1926-1400), side="left")
colnames(ULWa1_Si_padded)<-c(1:1926)
colnames(ULWa2_Si_padded)<-c(1:1926)
ULWa12_Si_padded<-rbind(ULWa2_Si_padded, ULWa1_Si_padded)
#heatmap(as.matrix(ULWa12_Si_padded), Rowv=NA, Colv=NA, asp=1) 

#ULWa2-ULWb1:
#complete overlap, ends at same spot as a2. 
#total length = 1926
ULWb1_Si_padded<-padzeros(data=ULWb1_Si, nzeros=(1926-1440), side="left")
colnames(ULWb1_Si_padded)<-c(1:1926)
ULWa12b_Si_padded<-rbind(ULWa12_Si_padded, ULWb1_Si_padded)
#heatmap(as.matrix(ULWa12b_Si_padded), Rowv=NA, Colv=NA, asp=1) 

#UlWa2-ULWa3:
#directly on top of each other, no overlap.
#1926+1288=3214 total length
ULWa12b_Si_padded<-padzeros(data=ULWa12b_Si_padded, nzeros=(3214-1926), side="right")
ULWa3_Si_padded<-padzeros(data=ULWa3_Si, nzeros=(3214-1288), side="left")
colnames(ULWa12b_Si_padded)<-c(1:3214)
colnames(ULWa3_Si_padded)<-c(1:3214)
ULWa12b3_Si_padded<-rbind(ULWa3_Si_padded, ULWa12b_Si_padded)
#heatmap(as.matrix(ULWa12b3_Si_padded), Rowv=NA, Colv=NA, asp=1)

#ULWa3-ULWc1:
#372 overlap
#(3214+1260)-372=4102 total length
ULWa12b3_Si_padded<-padzeros(data=ULWa12b3_Si_padded, nzeros=(4102-3214), side="right")
ULWc1_Si_padded<-padzeros(data=ULWc1_Si, nzeros=(4102-1260), side="left")
colnames(ULWa12b3_Si_padded)<-c(1:4102)
colnames(ULWc1_Si_padded)<-c(1:4102)
ULWa12b3c1_Si_padded<-rbind(ULWc1_Si_padded, ULWa12b3_Si_padded)
#heatmap(as.matrix(ULWa12b3c1_Si_padded), Rowv=NA, Colv=NA, asp=1)

#ULWc1-ULWc2:
#directly on top of each other, no overlap.
#1338+4102=5440 total length
ULWa12b3c1_Si_padded<-padzeros(data=ULWa12b3c1_Si_padded, nzeros=(5440-4102), side="right")
ULWc2_Si_padded<-padzeros(data=ULWc2_Si, nzeros=(5440-1338), side="left")
colnames(ULWa12b3c1_Si_padded)<-c(1:5440)
colnames(ULWc2_Si_padded)<-c(1:5440)
ULWa12b3c12_Si_padded<-rbind(ULWc2_Si_padded, ULWa12b3c1_Si_padded)
#heatmap(as.matrix(ULWa12b3c12_Si_padded), Rowv=NA, Colv=NA, asp=1)

#ULWc2-ULWd1:
#complete overlap. 5440 total length
#d1 ends 318 before c2, total d1 length is 1158
#so 318 padding on right, 5440-1158-318=3964 padding on left
ULWd1_Si_padded<-padzeros(ULWd1_Si, nzeros=318, side="right")
ULWd1_Si_padded<-padzeros(ULWd1_Si_padded, nzeros=3964, side="left")
colnames(ULWd1_Si_padded)<-c(1:5440)
ULWa12b3c12d1_Si_padded<-rbind(ULWa12b3c12_Si_padded, ULWd1_Si_padded)
#heatmap(as.matrix(ULWa12b3c12d1_Si_padded), Rowv=NA, Colv=NA, asp=1) 

#ULWd1-ULWd2:
#263 overlap between c2 and d2.
#(5440+493)-263=5670 total lenght.
ULWa12b3c12d1_Si_padded<-padzeros(data=ULWa12b3c12d1_Si_padded, nzeros=(5670-5440), side="right")
ULWd2_Si_padded<-padzeros(data=ULWd2_Si, nzeros=(5670-493), side="left")
colnames(ULWa12b3c12d1_Si_padded)<-c(1:5670)
colnames(ULWd2_Si_padded)<-c(1:5670)
ULWa12b3c12d12_Si_padded<-rbind(ULWd2_Si_padded, ULWa12b3c12d1_Si_padded)
#heatmap(as.matrix(ULWa12b3c12d12_Si_padded), Rowv=NA, Colv=NA, asp=1)

#ULWd2-ULWf1:
#130 overlap. 
#(5670+1586)-130=7126 total length
ULWa12b3c12d12_Si_padded<-padzeros(data=ULWa12b3c12d12_Si_padded, nzeros=(7126-5670), side="right")
ULWf1_Si_padded<-padzeros(data=ULWf1_Si, nzeros=(7126-1586), side="left")
colnames(ULWa12b3c12d12_Si_padded)<-c(1:7126)
colnames(ULWf1_Si_padded)<-c(1:7126)
ULWa12b3c12d12f1_Si_padded<-rbind(ULWf1_Si_padded, ULWa12b3c12d12_Si_padded)
ULWa12b3c12d12f1_Si_padded[ULWa12b3c12d12f1_Si_padded==0]<-NA
heatmap(as.matrix(ULWa12b3c12d12f1_Si_padded), Rowv=NA, Colv=NA, asp=1)




#--------------------------Mn----------------------------------------------------
#ULWa1-ULWa2:
#directly on top of each other, no overlap.
#526+1400=1926 total length
ULWa1_Mn_padded<-padzeros(data=ULWa1_Mn, nzeros=(1926-526), side="right")
ULWa2_Mn_padded<-padzeros(data=ULWa2_Mn, nzeros=(1926-1400), side="left")
colnames(ULWa1_Mn_padded)<-c(1:1926)
colnames(ULWa2_Mn_padded)<-c(1:1926)
ULWa12_Mn_padded<-rbind(ULWa2_Mn_padded, ULWa1_Mn_padded)
#heatmap(as.matrix(ULWa12_Mn_padded), Rowv=NA, Colv=NA, asp=1) 

#ULWa2-ULWb1:
#complete overlap, ends at same spot as a2. 
#total length = 1926
ULWb1_Mn_padded<-padzeros(data=ULWb1_Mn, nzeros=(1926-1440), side="left")
colnames(ULWb1_Mn_padded)<-c(1:1926)
ULWa12b_Mn_padded<-rbind(ULWa12_Mn_padded, ULWb1_Mn_padded)
#heatmap(as.matrix(ULWa12b_Mn_padded), Rowv=NA, Colv=NA, asp=1) 

#UlWa2-ULWa3:
#directly on top of each other, no overlap.
#1926+1288=3214 total length
ULWa12b_Mn_padded<-padzeros(data=ULWa12b_Mn_padded, nzeros=(3214-1926), side="right")
ULWa3_Mn_padded<-padzeros(data=ULWa3_Mn, nzeros=(3214-1288), side="left")
colnames(ULWa12b_Mn_padded)<-c(1:3214)
colnames(ULWa3_Mn_padded)<-c(1:3214)
ULWa12b3_Mn_padded<-rbind(ULWa3_Mn_padded, ULWa12b_Mn_padded)
#heatmap(as.matrix(ULWa12b3_Mn_padded), Rowv=NA, Colv=NA, asp=1)

#ULWa3-ULWc1:
#372 overlap
#(3214+1260)-372=4102 total length
ULWa12b3_Mn_padded<-padzeros(data=ULWa12b3_Mn_padded, nzeros=(4102-3214), side="right")
ULWc1_Mn_padded<-padzeros(data=ULWc1_Mn, nzeros=(4102-1260), side="left")
colnames(ULWa12b3_Mn_padded)<-c(1:4102)
colnames(ULWc1_Mn_padded)<-c(1:4102)
ULWa12b3c1_Mn_padded<-rbind(ULWc1_Mn_padded, ULWa12b3_Mn_padded)
#heatmap(as.matrix(ULWa12b3c1_Mn_padded), Rowv=NA, Colv=NA, asp=1)

#ULWc1-ULWc2:
#directly on top of each other, no overlap.
#1338+4102=5440 total length
ULWa12b3c1_Mn_padded<-padzeros(data=ULWa12b3c1_Mn_padded, nzeros=(5440-4102), side="right")
ULWc2_Mn_padded<-padzeros(data=ULWc2_Mn, nzeros=(5440-1338), side="left")
colnames(ULWa12b3c1_Mn_padded)<-c(1:5440)
colnames(ULWc2_Mn_padded)<-c(1:5440)
ULWa12b3c12_Mn_padded<-rbind(ULWc2_Mn_padded, ULWa12b3c1_Mn_padded)
#heatmap(as.matrix(ULWa12b3c12_Mn_padded), Rowv=NA, Colv=NA, asp=1)

#ULWc2-ULWd1:
#complete overlap. 5440 total length
#d1 ends 318 before c2, total d1 length is 1158
#so 318 padding on right, 5440-1158-318=3964 padding on left
ULWd1_Mn_padded<-padzeros(ULWd1_Mn, nzeros=318, side="right")
ULWd1_Mn_padded<-padzeros(ULWd1_Mn_padded, nzeros=3964, side="left")
colnames(ULWd1_Mn_padded)<-c(1:5440)
ULWa12b3c12d1_Mn_padded<-rbind(ULWa12b3c12_Mn_padded, ULWd1_Mn_padded)
#heatmap(as.matrix(ULWa12b3c12d1_Mn_padded), Rowv=NA, Colv=NA, asp=1) 

#ULWd1-ULWd2:
#263 overlap between c2 and d2.
#(5440+493)-263=5670 total lenght.
ULWa12b3c12d1_Mn_padded<-padzeros(data=ULWa12b3c12d1_Mn_padded, nzeros=(5670-5440), side="right")
ULWd2_Mn_padded<-padzeros(data=ULWd2_Mn, nzeros=(5670-493), side="left")
colnames(ULWa12b3c12d1_Mn_padded)<-c(1:5670)
colnames(ULWd2_Mn_padded)<-c(1:5670)
ULWa12b3c12d12_Mn_padded<-rbind(ULWd2_Mn_padded, ULWa12b3c12d1_Mn_padded)
#heatmap(as.matrix(ULWa12b3c12d12_Mn_padded), Rowv=NA, Colv=NA, asp=1)

#ULWd2-ULWf1:
#130 overlap. 
#(5670+1586)-130=7126 total length
ULWa12b3c12d12_Mn_padded<-padzeros(data=ULWa12b3c12d12_Mn_padded, nzeros=(7126-5670), side="right")
ULWf1_Mn_padded<-padzeros(data=ULWf1_Mn, nzeros=(7126-1586), side="left")
colnames(ULWa12b3c12d12_Mn_padded)<-c(1:7126)
colnames(ULWf1_Mn_padded)<-c(1:7126)
ULWa12b3c12d12f1_Mn_padded<-rbind(ULWf1_Mn_padded, ULWa12b3c12d12_Mn_padded)
ULWa12b3c12d12f1_Mn_padded[ULWa12b3c12d12f1_Mn_padded==0]<-NA
heatmap(as.matrix(ULWa12b3c12d12f1_Mn_padded), Rowv=NA, Colv=NA, asp=1)



#write all to .csv:
write.csv(ULWa12b3c12d12f1_Ca_padded, "ULW_comp_map_Ca.csv")
write.csv(ULWa12b3c12d12f1_Fe_padded, "ULW_comp_map_Fe.csv")
write.csv(ULWa12b3c12d12f1_K_padded, "ULW_comp_map_K.csv")
write.csv(ULWa12b3c12d12f1_Ti_padded, "ULW_comp_map_Ti.csv")
write.csv(ULWa12b3c12d12f1_Si_padded, "ULW_comp_map_Si.csv")
write.csv(ULWa12b3c12d12f1_Mn_padded, "ULW_comp_map_Mn.csv")




#------------------------Composite records-----------------------------------------
#read in data so the whole script doesn't have to be run
ULW_Ca<-read.csv("ULW_comp_map_Ca.csv")
ULW_Fe<-read.csv("ULW_comp_map_Fe.csv")
ULW_K<-read.csv("ULW_comp_map_K.csv")
ULW_Ti<-read.csv("ULW_comp_map_Ti.csv")
ULW_Si<-read.csv("ULW_comp_map_Si.csv")
ULW_Mn<-read.csv("ULW_comp_map_Mn.csv")

#Once again remove row index
ULW_Ca<-ULW_Ca[,2:ncol(ULW_Ca)]
ULW_Fe<-ULW_Fe[,2:ncol(ULW_Fe)]
ULW_K<-ULW_K[,2:ncol(ULW_K)]
ULW_Ti<-ULW_Ti[,2:ncol(ULW_Ti)]
ULW_Si<-ULW_Si[,2:ncol(ULW_Si)]
ULW_Mn<-ULW_Mn[,2:ncol(ULW_Mn)]


#------------------------Ca-------------------------------------------------------
#Average data:
ULW_Ca_mean<-colMeans(ULW_Ca, na.rm=T)

#add depth column, step size was 100 um = 0.1 mm
ULW_Ca_d<-seq(0, length(ULW_Ca_mean)-1, 1)*0.1
ULW_Ca_d<-ULW_Ca_d+7.5 #starts at 7.4 mm
ULW_Ca_depth<-cbind(ULW_Ca_d, ULW_Ca_mean)

#colnames
colnames(ULW_Ca_depth)<-c("depth_mm", "Ca")

#to dataframe
ULW_Ca_depth<-as.data.frame(ULW_Ca_depth)

#plot
ggplot()+
  geom_line(data=ULW_Ca_depth, aes(x=depth_mm, y=Ca), col="navy")+
  theme_classic()


#------------------------Fe-------------------------------------------------------
#Average data:
ULW_Fe_mean<-colMeans(ULW_Fe, na.rm=T)

#add depth column, step size was 100 um = 0.1 mm
ULW_Fe_d<-seq(0, length(ULW_Fe_mean)-1, 1)*0.1
ULW_Fe_d<-ULW_Fe_d+7.5 #starts at 7.4 mm
ULW_Fe_depth<-cbind(ULW_Fe_d, ULW_Fe_mean)

#colnames
colnames(ULW_Fe_depth)<-c("depth_mm", "Fe")

#to dataframe
ULW_Fe_depth<-as.data.frame(ULW_Fe_depth)

#plot
ggplot()+
  geom_line(data=ULW_Fe_depth, aes(x=depth_mm, y=Fe), col="navy")+
  theme_classic()


#------------------------K-------------------------------------------------------
#Average data:
ULW_K_mean<-colMeans(ULW_K, na.rm=T)

#add depth column, step size was 100 um = 0.1 mm
ULW_K_d<-seq(0, length(ULW_K_mean)-1, 1)*0.1
ULW_K_d<-ULW_K_d+7.5 #starts at 7.4 mm
ULW_K_depth<-cbind(ULW_K_d, ULW_K_mean)

#colnames
colnames(ULW_K_depth)<-c("depth_mm", "K")

#to dataframe
ULW_K_depth<-as.data.frame(ULW_K_depth)

#plot
ggplot()+
  geom_line(data=ULW_K_depth, aes(x=depth_mm, y=K), col="navy")+
  theme_classic()


#------------------------Ti-------------------------------------------------------
#Average data:
ULW_Ti_mean<-colMeans(ULW_Ti, na.rm=T)

#add depth column, step size was 100 um = 0.1 mm
ULW_Ti_d<-seq(0, length(ULW_Ti_mean)-1, 1)*0.1
ULW_Ti_d<-ULW_Ti_d+7.5 #starts at 7.4 mm
ULW_Ti_depth<-cbind(ULW_Ti_d, ULW_Ti_mean)

#colnames
colnames(ULW_Ti_depth)<-c("depth_mm", "Ti")

#to dataframe
ULW_Ti_depth<-as.data.frame(ULW_Ti_depth)

#plot
ggplot()+
  geom_line(data=ULW_Ti_depth, aes(x=depth_mm, y=Ti), col="navy")+
  theme_classic()


#------------------------Si-------------------------------------------------------
#Average data:
ULW_Si_mean<-colMeans(ULW_Si, na.rm=T)

#add depth column, step size was 100 um = 0.1 mm
ULW_Si_d<-seq(0, length(ULW_Si_mean)-1, 1)*0.1
ULW_Si_d<-ULW_Si_d+7.5 #starts at 7.4 mm
ULW_Si_depth<-cbind(ULW_Si_d, ULW_Si_mean)

#colnames
colnames(ULW_Si_depth)<-c("depth_mm", "Si")

#to dataframe
ULW_Si_depth<-as.data.frame(ULW_Si_depth)

#plot
ggplot()+
  geom_line(data=ULW_Si_depth, aes(x=depth_mm, y=Si), col="navy")+
  theme_classic()


#------------------------Mn-------------------------------------------------------
#Average data:
ULW_Mn_mean<-colMeans(ULW_Mn, na.rm=T)

#add depth column, step size was 100 um = 0.1 mm
ULW_Mn_d<-seq(0, length(ULW_Mn_mean)-1, 1)*0.1
ULW_Mn_d<-ULW_Mn_d+7.5 #starts at 7.4 mm
ULW_Mn_depth<-cbind(ULW_Mn_d, ULW_Mn_mean)

#colnames
colnames(ULW_Mn_depth)<-c("depth_mm", "Mn")

#to dataframe
ULW_Mn_depth<-as.data.frame(ULW_Mn_depth)

#plot
ggplot()+
  geom_line(data=ULW_Mn_depth, aes(x=depth_mm, y=Mn), col="navy")+
  theme_classic()




#write all to csv
write.csv(ULW_Ca_depth, "ULW_Ca_depth.csv")
write.csv(ULW_Fe_depth, "ULW_Fe_depth.csv")
write.csv(ULW_K_depth, "ULW_K_depth.csv")
write.csv(ULW_Ti_depth, "ULW_Ti_depth.csv")
write.csv(ULW_Si_depth, "ULW_Si_depth.csv")
write.csv(ULW_Mn_depth, "ULW_Mn_depth.csv")



#------------------------Quick plotting of records-------------------------------
#read in data:
ULW_Ca_depth<-read.csv("ULW_Ca_depth.csv")
ULW_Fe_depth<-read.csv("ULW_Fe_depth.csv")
ULW_K_depth<-read.csv("ULW_K_depth.csv")
ULW_Ti_depth<-read.csv("ULW_Ti_depth.csv")
ULW_Mn_depth<-read.csv("ULW_Mn_depth.csv")
ULW_Si_depth<-read.csv("ULW_Si_depth.csv")

#formatting: remove first column
ULW_Ca_depth<-ULW_Ca_depth[,2:3]
ULW_Fe_depth<-ULW_Fe_depth[,2:3]
ULW_K_depth<-ULW_K_depth[,2:3]
ULW_Ti_depth<-ULW_Ti_depth[,2:3]
ULW_Mn_depth<-ULW_Mn_depth[,2:3]
ULW_Si_depth<-ULW_Si_depth[,2:3]

#Plot individual records:
ggplot()+
  geom_line(data=ULW_Ca_depth, aes(x=depth_mm, y=Ca), col="navy")+
  theme_classic()

ggplot()+
  geom_line(data=ULW_Fe_depth, aes(x=depth_mm, y=Fe), col="navy")+
  theme_classic()

ggplot()+
  geom_line(data=ULW_K_depth, aes(x=depth_mm, y=K), col="navy")+
  theme_classic()

ggplot()+
  geom_line(data=ULW_Ti_depth, aes(x=depth_mm, y=Ti), col="navy")+
  theme_classic()

ggplot()+
  geom_line(data=ULW_Si_depth, aes(x=depth_mm, y=Si), col="navy")+
  theme_classic()

ggplot()+
  geom_line(data=ULW_Mn_depth, aes(x=depth_mm, y=Mn), col="navy")+
  theme_classic()



#Plot together:
Ca_p<-ggplot()+
  geom_line(data=ULW_Ca_depth, aes(x=depth_mm, y=Ca), col="navy")+
  coord_flip()+
  theme_classic()

Fe_p<-ggplot()+
  geom_line(data=ULW_Fe_depth, aes(x=depth_mm, y=Fe), col="navy")+
  coord_flip()+
  theme_classic()

K_p<-ggplot()+
  geom_line(data=ULW_K_depth, aes(x=depth_mm, y=K), col="navy")+
  coord_flip()+
  theme_classic()

Ti_p<-ggplot()+
  geom_line(data=ULW_Ti_depth, aes(x=depth_mm, y=Ti), col="navy")+
  coord_flip()+
  theme_classic()

Si_p<-ggplot()+
  geom_line(data=ULW_Si_depth, aes(x=depth_mm, y=Si), col="navy")+
  coord_flip()+
  theme_classic()

Mn_p<-ggplot()+
  geom_line(data=ULW_Mn_depth, aes(x=depth_mm, y=Mn), col="navy")+
  coord_flip()+
  theme_classic()


cowplot::plot_grid(Ca_p, Fe_p, K_p, Ti_p, Si_p, Mn_p, ncol=6)
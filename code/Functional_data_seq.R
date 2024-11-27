library(openxlsx)
library(ggplot2)
library(Rmisc)

date<-"20241121"

#rename experiments based on exp, file of origin, by
setwd("~/MBC Dropbox/Lab Poli PhD/Aurora/Projects_wd/BC networks")
colony_num<-read.xlsx("data/Functional data/Functional data table - colony.xlsx", 1)
colony_num$treatment<-factor(colony_num$treatment, levels=c("WT", "EV", "TFDP1", "E2F3"))
colony_num$exp_name<-paste(colony_num$by, colony_num$file.of.origin, colony_num$exp)

colony_size<-read.xlsx("data/Functional data/Functional data table - colony.xlsx", 2)
colony_size$treatment<-factor(colony_size$treatment, levels=c("WT", "EV", "TFDP1", "E2F3"))
colony_size$exp_name<-paste(colony_size$by, colony_size$file.of.origin, colony_size$exp)

prolif<-read.xlsx("data/Functional data/Functional data table - prolif.xlsx", 1)
prolif$treatment<-factor(prolif$treatment, levels=c("WT", "EV", "TFDP1", "E2F3"))
prolif$exp_name<-paste(prolif$by, prolif$file.of.origin, prolif$exp)

###########
## only sequenced clones
##########

p<-matrix(ncol=5, nrow=3)
d<-matrix(ncol=5, nrow=3)

rownames(p)<-c("Colony_num", "Colony_size", "Prolif")
rownames(d)<-c("Colony_num", "Colony_size", "Prolif")
colnames(d)<-c("Hs TFDP1", "231 TFDP1", "231 E2F3", "468 TFDP1", "468 E2F3")
colnames(p)<-c("Hs TFDP1", "231 TFDP1", "231 E2F3", "468 TFDP1", "468 E2F3")

pc<-matrix(ncol=10, nrow=3)
dc<-matrix(ncol=10, nrow=3)

rownames(pc)<-c("Colony_num", "Colony_size", "Prolif")
rownames(dc)<-c("Colony_num", "Colony_size", "Prolif")
colnames(dc)<-c("Hs TFDP1 14", "Hs TFDP1 18", "231 TFDP1 21","231 TFDP1 6", 
                "231 E2F3 18", "231 E2F3 20", "468 TFDP1 4", "468 TFDP1 9", 
                "468 E2F3 1", "468 E2F3 3")
colnames(pc)<-c("Hs TFDP1 14", "Hs TFDP1 18", "231 TFDP1 21","231 TFDP1 6", 
                "231 E2F3 18", "231 E2F3 20", "468 TFDP1 4", "468 TFDP1 9", 
                "468 E2F3 1", "468 E2F3 3")

colony_num_sel<-colony_num[c(which(colony_num$cell.line %in% c("HS-578T", "MDA-468", "MDA-231") & colony_num$treatment=="EV" & colony_num$clone.ID=="bulk" ),
                             which(colony_num$cell.line %in% c("HS-578T") & colony_num$treatment=="TFDP1" & colony_num$clone.ID %in% c("14", "18")),
                             which(colony_num$cell.line %in% c("MDA-468") & colony_num$treatment=="TFDP1" & colony_num$clone.ID %in% c("4", "9")),
                             which(colony_num$cell.line %in% c("MDA-468") & colony_num$treatment=="E2F3" & colony_num$clone.ID %in% c("1", "3")),
                             which(colony_num$cell.line %in% c("MDA-231") & colony_num$treatment=="TFDP1" & colony_num$clone.ID %in% c("6", "21")),
                             which(colony_num$cell.line %in% c("MDA-231") & colony_num$treatment=="E2F3" & colony_num$clone.ID %in% c("18", "20"))),]


colony_num_sel$treatment<-factor(colony_num_sel$treatment, levels=c("WT", "EV", "TFDP1", "E2F3"))
ggplot(colony_num_sel, aes(x=treatment, y=n.of.colonies))+geom_boxplot()+facet_wrap(~cell.line)
ggplot(colony_num_sel, aes(x=treatment, y=n.of.colonies))+geom_boxplot()+facet_wrap(~cell.line+exp_name)

png(paste("results/",date, "/colony_num_seq_exp.png", sep=""), res=300, 2500, 2500)
ggplot(colony_num_sel, aes(x=treatment, y=n.of.colonies, fill=clone.ID))+geom_boxplot()+facet_wrap(~cell.line+exp_name)+                                                                # Change font size
  theme(strip.text.x = element_text(size = 4))
dev.off()

ggplot(colony_num_sel, aes(x=treatment, y=n.of.colonies, fill=clone.ID))+geom_boxplot()+facet_wrap(~cell.line)

colony_num_sel$treatment2<-paste(colony_num_sel$treatment, colony_num_sel$clone.ID)

av<-aov(n.of.colonies~cell.line+treatment+clone.ID+exp_name, colony_num_sel)
summary(av)

#as clone is nested within treatment, we need to test treatment and clone separately
avHs<-aov(n.of.colonies~treatment+exp_name, subset(colony_num_sel, cell.line=="HS-578T"))
t<-TukeyHSD(avHs, which=c("treatment"))
d["Colony_num",1]<-t[[1]][,"diff"]
p["Colony_num",1]<-t[[1]][,"p adj"]

av231<-aov(n.of.colonies~treatment+exp_name, subset(colony_num_sel, cell.line=="MDA-231"))
t<-TukeyHSD(av231, which=c("treatment"))
d["Colony_num",c(2,3)]<-t[[1]][,"diff"][c(1,2)]
p["Colony_num",c(2,3)]<-t[[1]][,"p adj"][c(1,2)]

av468<-aov(n.of.colonies~treatment+exp_name, subset(colony_num_sel, cell.line=="MDA-468"))
t<-TukeyHSD(av468, which=c("treatment"))
d["Colony_num",c(4,5)]<-t[[1]][,"diff"][c(1,2)]
p["Colony_num",c(4,5)]<-t[[1]][,"p adj"][c(1,2)]



av<-aov(n.of.colonies~cell.line+treatment2+exp_name, colony_num_sel)
summary(av)

avHs<-aov(n.of.colonies~treatment2+exp_name, subset(colony_num_sel, cell.line=="HS-578T"))
t<-TukeyHSD(avHs, which=c("treatment2"))
dc["Colony_num",c(1,2)]<-t[[1]][,"diff"][c(1,2)]
pc["Colony_num",c(1,2)]<-t[[1]][,"p adj"][c(1,2)]

av231<-aov(n.of.colonies~treatment2+exp_name, subset(colony_num_sel, cell.line=="MDA-231"))
t<-TukeyHSD(av231, which=c("treatment2"))
dc["Colony_num",c(3,4)]<-t[[1]][,"diff"][c(8,9)]
dc["Colony_num",c(5,6)]<- -t[[1]][,"diff"][c(2,5)]
pc["Colony_num",c(3:6)]<-t[[1]][,"p adj"][c(8,9,2,5)]

av468<-aov(n.of.colonies~treatment2+exp_name, subset(colony_num_sel, cell.line=="MDA-468"))
t<-TukeyHSD(av468, which=c("treatment2"))
dc["Colony_num",c(7,8)]<-t[[1]][,"diff"][c(8,9)]
dc["Colony_num",c(9,10)]<- -t[[1]][,"diff"][c(2,5)]
pc["Colony_num",c(7:10)]<-t[[1]][,"p adj"][c(8,9,2,5)]

colony_size_sel<-colony_size[c(which(colony_size$cell.line %in% c("HS-578T", "MDA-468", "MDA-231") & colony_size$treatment=="EV" & colony_size$clone.ID=="bulk" ),
                               which(colony_size$cell.line %in% c("HS-578T") & colony_size$treatment=="TFDP1" & colony_size$clone.ID %in% c("14", "18")),
                               which(colony_size$cell.line %in% c("MDA-468") & colony_size$treatment=="TFDP1" & colony_size$clone.ID %in% c("4", "9")),
                               which(colony_size$cell.line %in% c("MDA-468") & colony_size$treatment=="E2F3" & colony_size$clone.ID %in% c("1", "3")),
                               which(colony_size$cell.line %in% c("MDA-231") & colony_size$treatment=="TFDP1" & colony_size$clone.ID %in% c("6", "21")),
                               which(colony_size$cell.line %in% c("MDA-231") & colony_size$treatment=="E2F3" & colony_size$clone.ID %in% c("18", "20"))),]

colony_size_sel$treatment<-factor(colony_size_sel$treatment, levels=c("WT", "EV", "TFDP1", "E2F3"))

ggplot(colony_size_sel, aes(x=treatment, y=avg.size))+geom_boxplot()+facet_wrap(~cell.line)
ggplot(colony_size_sel, aes(x=treatment, y=avg.size))+geom_boxplot()+facet_wrap(~cell.line+exp_name)

png(paste("results/",date, "/colony_size_seq_exp.png", sep=""), res=300, 2500, 2500)
ggplot(colony_size_sel, aes(x=treatment, y=avg.size, fill=clone.ID))+geom_boxplot()+facet_wrap(~cell.line+exp_name)+                                                                # Change font size
  theme(strip.text.x = element_text(size = 4))
dev.off()


ggplot(colony_size_sel, aes(x=treatment, y=avg.size, fill=clone.ID))+geom_boxplot()+facet_wrap(~cell.line)

colony_size_sel$treatment2<-paste(colony_size_sel$treatment, colony_size_sel$clone.ID)
av<-aov(avg.size~cell.line+treatment+clone.ID+exp_name, colony_size_sel)
summary(av)

avHs<-aov(avg.size~treatment+exp_name, subset(colony_size_sel, cell.line=="HS-578T"))
t<-TukeyHSD(avHs, which=c("treatment"))
d["Colony_size",1]<-t[[1]][,"diff"]
p["Colony_size",1]<-t[[1]][,"p adj"]

av231<-aov(avg.size~treatment+exp_name, subset(colony_size_sel, cell.line=="MDA-231"))
t<-TukeyHSD(av231, which=c("treatment"))
d["Colony_size",c(2,3)]<-t[[1]][,"diff"][c(1,2)]
p["Colony_size",c(2,3)]<-t[[1]][,"p adj"][c(1,2)]

av468<-aov(avg.size~treatment+exp_name, subset(colony_size_sel, cell.line=="MDA-468"))
t<-TukeyHSD(av468, which=c("treatment"))
d["Colony_size",c(4,5)]<-t[[1]][,"diff"][c(1,2)]
p["Colony_size",c(4,5)]<-t[[1]][,"p adj"][c(1,2)]

av<-aov(avg.size~cell.line+treatment2+exp_name, colony_size_sel)
summary(av)

avHs<-aov(avg.size~treatment2+exp_name, subset(colony_size_sel, cell.line=="HS-578T"))
t<-TukeyHSD(avHs, which=c("treatment2"))
dc["Colony_size",c(1,2)]<-t[[1]][,"diff"][c(1,2)]
pc["Colony_size",c(1,2)]<-t[[1]][,"p adj"][c(1,2)]

av231<-aov(avg.size~treatment2+exp_name, subset(colony_size_sel, cell.line=="MDA-231"))
t<-TukeyHSD(av231, which=c("treatment2"))
dc["Colony_size",c(3,4)]<-t[[1]][,"diff"][c(8,9)]
dc["Colony_size",c(5,6)]<- -t[[1]][,"diff"][c(2,5)]
pc["Colony_size",c(3:6)]<-t[[1]][,"p adj"][c(8,9,2,5)]

av468<-aov(avg.size~treatment2+exp_name, subset(colony_size_sel, cell.line=="MDA-468"))
t<-TukeyHSD(av468, which=c("treatment2"))
dc["Colony_size",c(7,8)]<-t[[1]][,"diff"][c(1,2)]
pc["Colony_size",c(7,8)]<-t[[1]][,"p adj"][c(1,2)]


prolif_sel<-prolif[c(which(prolif$cell.line %in% c("HS-578T", "MDA-468", "MDA-231") & prolif$treatment=="EV" & prolif$clone.ID=="bulk" ),
                     which(prolif$cell.line %in% c("HS-578T") & prolif$treatment=="TFDP1" & prolif$clone.ID %in% c("14", "18")),
                     which(prolif$cell.line %in% c("MDA-468") & prolif$treatment=="TFDP1" & prolif$clone.ID %in% c("4", "9")),
                     which(prolif$cell.line %in% c("MDA-468") & prolif$treatment=="E2F3" & prolif$clone.ID %in% c("1", "3")),
                     which(prolif$cell.line %in% c("MDA-231") & prolif$treatment=="TFDP1" & prolif$clone.ID %in% c("6", "21")),
                     which(prolif$cell.line %in% c("MDA-231") & prolif$treatment=="E2F3" & prolif$clone.ID %in% c("18", "20"))),]

prolif_se <- summarySE(prolif_sel, na.rm=T, measurevar="value", groupvars=c("time","cell.line","treatment"))

png(paste("results/",date, "/prolif_seq_summary.png", sep=""), res=300, 2500, 2000)
ggplot(prolif_se, aes(x = time, y = value, color = treatment)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, size=1) +
  geom_line(linewidth=1) +
  geom_point()+theme_bw()+facet_grid(.~cell.line)+ylab("Norm OD")
dev.off()



prolif_se <- summarySE(prolif_sel, na.rm=T, measurevar="value", groupvars=c("time","cell.line","treatment", "clone.ID","exp_name"))

png(paste("results/",date, "/prolif_seq_exp.png", sep=""), res=300, 7000, 2000)
ggplot(prolif_se, aes(x = time, y = value, color = treatment, shape=clone.ID)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, size=1) +
  geom_line(linewidth=1) +
  geom_point()+theme_bw()+facet_grid(.~cell.line+exp_name)+ylab("Norm OD")+theme(strip.text.x = element_text(size = 5))
dev.off()

prolif_se <- summarySE(prolif_sel, na.rm=T, measurevar="value", groupvars=c("time","cell.line","treatment", "clone.ID"))

png(paste("results/",date, "/prolif_seq_summary_clones.png", sep=""), res=300, 2500, 2000)
ggplot(prolif_se, aes(x = time, y = value, color = treatment, shape=clone.ID)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, size=1) +
  geom_line(linewidth=1) +
  geom_point()+theme_bw()+facet_grid(.~cell.line)+ylab("Norm OD")
dev.off()



avHs<-aov(value~treatment+exp_name+as.factor(time), subset(prolif_sel, cell.line=="HS-578T"))
t<-TukeyHSD(avHs, which=c("treatment"))
d["Prolif",1]<-t[[1]][,"diff"]
p["Prolif",1]<-t[[1]][,"p adj"]

av231<-aov(value~treatment+exp_name+as.factor(time), subset(prolif_sel, cell.line=="MDA-231"))
t<-TukeyHSD(av231, which=c("treatment"))
d["Prolif",c(2,3)]<-t[[1]][,"diff"][c(1,2)]
p["Prolif",c(2,3)]<-t[[1]][,"p adj"][c(1,2)]

av468<-aov(value~treatment+exp_name+as.factor(time), subset(prolif_sel, cell.line=="MDA-468"))
t<-TukeyHSD(av468, which=c("treatment"))
d["Prolif",c(4,5)]<-t[[1]][,"diff"][c(1,2)]
p["Prolif",c(4,5)]<-t[[1]][,"p adj"][c(1,2)]


prolif_sel$treatment2<-paste(prolif_sel$treatment, prolif_sel$clone.ID)

av<-aov(value~cell.line+treatment2+exp_name, prolif_sel)
summary(av)

avHs<-aov(value~treatment2+exp_name+as.factor(time), subset(prolif_sel, cell.line=="HS-578T"))
t<-TukeyHSD(avHs, which=c("treatment2"))
dc["Prolif",c(1,2)]<-t[[1]][,"diff"][c(1,2)]
pc["Prolif",c(1,2)]<-t[[1]][,"p adj"][c(1,2)]

av231<-aov(value~treatment2+exp_name+as.factor(time), subset(prolif_sel, cell.line=="MDA-231"))
t<-TukeyHSD(av231, which=c("treatment2"))
dc["Prolif",c(3,4)]<- t[[1]][,"diff"][c(4,5)]
dc["Prolif",c(6)]<- -t[[1]][,"diff"][c(1)]
pc["Prolif",c(3,4)]<-t[[1]][,"p adj"][c(4,5)]
pc["Prolif",c(6)]<-t[[1]][,"p adj"][c(1)]


av468<-aov(value~treatment2+exp_name+as.factor(time), subset(prolif_sel, cell.line=="MDA-468"))
t<-TukeyHSD(av468, which=c("treatment2"))
dc["Prolif",c(7,8)]<-t[[1]][,"diff"][c(8,9)]
dc["Prolif",c(9,10)]<- -t[[1]][,"diff"][c(2,5)]
pc["Prolif",c(7:10)]<-t[[1]][,"p adj"][c(8,9,2,5)]

save(dc, file=paste("results/",date, "/dc.RData", sep=""))
save(pc, file=paste("results/",date, "/pc.RData", sep=""))
save(d, file=paste("results/",date, "/d.RData", sep=""))
save(p, file=paste("results/",date, "/p.RData", sep=""))

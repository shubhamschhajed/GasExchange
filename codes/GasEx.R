#### Import library(s) -----------
library(tidyverse)
library(ggplot2)
library(emmeans)
library(lsmeans)
library(xtable)
library(car)
library(ape)
library(nlme)


#### Clearing the environment -----------
rm(list=ls())

#### Compiling data from different LICORs -----------
LI5 <- read.csv("D://R//RCodesCNR1//Input//LI5.csv")
LI6 <- read.csv("D://R//RCodesCNR1//Input//LI6.csv")
LI <- bind_rows(LI5,LI6)
write.csv(LI, "D:\\R\\RCodesCNR1\\Output\\GasEx.csv")  # exporting

#### Importing data & converting it to correct format -----------

dsA <- openxlsx::read.xlsx("D://R//RCodesCNR1//Input//GasEx_Castlereagh.xlsx", 1)
dsR <- openxlsx::read.xlsx("D://R//RCodesCNR1//Input//GasEx_Castlereagh.xlsx", 2)

str(dsA)
dsA$Spp <- as.factor(dsA$Spp)
dsA$Rep <- as.factor(dsA$Rep)

str(dsR)
dsR$Spp <- as.factor(dsR$Spp)
dsR$Rep <- as.factor(dsR$Rep)



#### Step 1: Filter for those quality parameters from the LiCOR -----------

plot(dsA$gsw, dsA$A)
points(dsA$gsw[dsA$Humidifier_.< 0], dsA$A[dsA$Humidifier_.< 0], col="red")

# dsA <- dsA[dsA$Humidifier_.>= 0,]
# dsA <- dsA[dsA$Ci >= 50,]  # keep only positive leaf-internal CO2 concentration values
# dsA <- dsA[dsA$Ci <= 380,]  # keep data with leaf-internal CO2 concentration values less than ambient CO2 concentration values
# dsA <- dsA[dsA$gsw >= 0,]  # keep only positive stomatal conductance values
# dsA <- dsA[dsA$CiCa >= 0.1,]

dsA$Ci[dsA$Ci <= 50] <- NA
dsA$Ci[dsA$Ci >= 380] <- NA
dsA$gsw[dsA$gsw <= 0] <- NA
dsA$CiCa[dsA$CiCa <= 0.1] <- NA
dsA$A[dsA$A <= 0] <- NA
dsA$A[dsA$A <= 3 & dsA$Spp == "CG"] <- NA
dsA$A[dsA$A <= 3 & dsA$Spp == "ER"] <- NA
dsA$A[dsA$CiCa < 3.5 & dsA$Spp == "PL"] <- NA

dsA <- subset(dsA, !(Ci <= 50))
dsA <- subset(dsA, !(Ci >=380))
dsA <- subset(dsA, !(CiCa <= 0.1))
dsA <- subset(dsA, !(gsw <= 0))
dsA <- subset(dsA, !(A <= 0))
dsA <- subset(dsA, !(Spp == "CG" & A <= 3))
dsA <- subset(dsA, !(Spp == "ER" & A <= 3))
dsA <- subset(dsA, !(Spp == "PL" & CiCa < 0.35))
#### Means -----------
dsAm <- dsA %>% group_by(Spp, Rep) %>% slice(which.max(A))
# summarise_if(is.numeric, .fun = mean, na.rm = TRUE)  # Removing pseudo-replicates

plyr::count(dsAm, "Spp")

dsA1 <- dsAm %>% group_by(Spp) %>%
  summarise_if(is.numeric, .fun = mean, na.rm = TRUE)  # Calculating species means

dsRm <- dsR %>% group_by(Spp, Rep) %>%
  summarise_if(is.numeric, .fun = mean, na.rm = TRUE)  # Removing pseudo-replicates

dsR1 <- dsR %>% group_by(Spp) %>%
  summarise_if(is.numeric, .fun = mean, na.rm = TRUE)  # Calculating species means




quantile(dsA1$A, na.rm = TRUE)

quantile(dsA1$A, c(0.25, 0.8), na.rm = TRUE)

myquant <- function(x){
  quantile(x, c(0.8), na.rm = TRUE)}

dsA2 <- dsA1 %>% group_by(Spp) %>%
  summarise_if(is.numeric, .fun = myquant)  # Calculating species means

dsA3 <- dsA1 %>% group_by(Spp) %>%
  summarise_if(is.numeric, .fun = median)  # Calculating species means

# dsA2 <- aggregate(dsA1[,6:16], by=list(Spp), FUN=myquant)

# dsA3 <- aggregate(dsA1, by=list(Spp), FUN=median)


# cor.test(dsAm$A, dsA3$A)
# cor.test(dsA2$A, dsAm$A)

plot(dsA3$gsw, dsA2$A, col="red", type = "p")
points(dsA3$gsw, dsA3$A, col="blue", type = "p")

ggplot(dsAm, aes(x=gsw, y=CiCa)) +
  geom_point() +
  geom_smooth(method="lm", se=F) + ggpubr::stat_cor() +
  # geom_abline(slope=0.85, col="blue") +
  # labs(y="A (median)", x="A (mean)") +
  directlabels::geom_dl(aes(label = Spp), method = list("last.points", cex = 0.8)) +
  theme_bw(base_size=11)+
  theme(axis.title=element_text(size=25),
        axis.text=element_text(size=15),
        legend.position="None",
        legend.title=element_blank(),
        legend.text=element_blank(),
        # legend.title=element_text(size=25),
        # legend.text=element_text(size=25),
        # panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent"),
        panel.border=element_rect(size=1))

ggsave("Output/GasExPlots/median-mean.tiff", compression = "lzw", width=220, height=200, units="mm")


write.csv(dsA1, "D:\\R\\RCodesCNR1\\Output\\ACleaned2.csv")  # exporting

# aa<-dsA2$A
# write.csv(aa, "D:\\R\\RCodesCNR1\\Output\\A80.csv")  # exporting
# aa<-dsA3$A
# write.csv(aa, "D:\\R\\RCodesCNR1\\Output\\Amedian.csv")  # exporting

#### Step 2: variation and outliers -----------

# identify outliers using univariate analysis 

ggplot(dsAm, aes(x = reorder(Spp, A, FUN = median), y = A)) +
  geom_boxplot() +
  geom_jitter() +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent"))

ggplot(dsA, aes(x=gsw, y=A)) +
  geom_point(aes(colour = VPDleaf)) +
  theme_classic()

ggplot(dsA, aes(x=gsw, y=A)) +
  geom_point(aes(colour = Tleaf)) +
  theme_classic()

ggplot(dsA, aes(x=Tleaf, y=A)) +
  geom_point(aes(colour = RHcham)) +
  theme_classic()

ggplot(dsA, aes(x=Ci/Ca, y=A)) +
  geom_point(aes(colour = Tleaf)) +
  theme_classic()

ggplot(dsA, aes(x=Ci/Ca, y=A)) +
  geom_point(aes(colour = RHcham)) +
  theme_classic()

#### Step 3: Variance Component Analysis -----------

lme1 <- lme(log10(A)~1, random=~1|Spp, data=dsA, na.action=na.exclude)
l1 <- varcomp(lme1, cum=FALSE, scale=TRUE)

lme2 <- lme(log10(gsw)~1, random=~1|Spp, data=dsA, na.action=na.exclude)
l2 <- varcomp(lme2, cum=FALSE, scale=TRUE)

lme3 <- lme(log10(Ci/Ca)~1, random=~1|Spp, data=dsA, na.action=na.exclude)
l3 <- varcomp(lme3, cum=FALSE, scale=TRUE)

lme4 <- lme(log10(E)~1, random=~1|Spp, data=dsA, na.action=na.exclude)
l4 <- varcomp(lme4, cum=FALSE, scale=TRUE)

lme5 <- lme(log10(gsw)~1, random=~1|Spp, data=dsA, na.action=na.exclude)
l5 <- varcomp(lme5, cum=FALSE, scale=TRUE)

lme6 <- lme(log10(gsw)~1, random=~1|Spp, data=dsA, na.action=na.exclude)
l6 <- varcomp(lme6, cum=FALSE, scale=TRUE)

lme7 <- lme(log10(gsw)~1, random=~1|Spp, data=dsA, na.action=na.exclude)
l7 <- varcomp(lme7, cum=FALSE, scale=TRUE)

# lme8<-lme(log10(trunk_wood_Rmass20_nmolg1s1) ~1, random=~1|site/Species,data=dat,na.action=na.exclude)
# l8<-varcomp(lme8, cum = FALSE,scale=TRUE)

ll1<-as.data.frame(as.matrix(l1[1:3]))
ll2<-as.data.frame(as.matrix(l2[1:3]))
ll3<-as.data.frame(as.matrix(l3[1:3]))
ll4<-as.data.frame(as.matrix(l4[1:3]))
ll5<-as.data.frame(as.matrix(l5[1:3]))
ll6<-as.data.frame(as.matrix(l6[1:3]))
ll7<-as.data.frame(as.matrix(l7[1:3]))

ll1$Tissue<-"Leaf"
ll2$Tissue<-"Stem bark"
ll3$Tissue<-"Stem wood"
ll4$Tissue<-"Trunk bark"
ll5$Tissue<-"Trunk wood"
ll6$Tissue<-"Trunk wood"
ll7$Tissue<-"Trunk wood"

bigms2<-rbind(ll1, ll2, ll3, ll4, ll5, ll6, ll7)
bigms2<-rbind(ll1, ll2, ll3, ll5)
r1<-rownames(bigms2)
bigms2$Source<-r1
bigms2$Source<-as.factor(bigms2$Source)

levels(bigms2$Source)[levels(bigms2$Source)=="site"] <- "Site"
levels(bigms2$Source)[levels(bigms2$Source)=="Within"] <- "Residual"
levels(bigms2$Source)[levels(bigms2$Source)=="site1"] <- "Site"
levels(bigms2$Source)[levels(bigms2$Source)=="Within1"] <- "Residual"
levels(bigms2$Source)[levels(bigms2$Source)=="site2"] <- "Site"
levels(bigms2$Source)[levels(bigms2$Source)=="Within2"] <- "Residual"
levels(bigms2$Source)[levels(bigms2$Source)=="site3"] <- "Site"
levels(bigms2$Source)[levels(bigms2$Source)=="Within3"] <- "Residual"
levels(bigms2$Source)[levels(bigms2$Source)=="site4"] <- "Site"
levels(bigms2$Source)[levels(bigms2$Source)=="Within4"] <- "Residual"
levels(bigms2$Source)[levels(bigms2$Source)=="Species1"] <- "Species"
levels(bigms2$Source)[levels(bigms2$Source)=="Species2"] <- "Species"
levels(bigms2$Source)[levels(bigms2$Source)=="Species3"] <- "Species"
levels(bigms2$Source)[levels(bigms2$Source)=="Species4"] <- "Species"

bigms2$Variance<-bigms2$V1

g1<-ggplot(bigms2, aes(fill=Source,y=Variance,x=Tissue))+
  geom_bar(position="stack",stat="identity",color="grey",width=0.85)+
  theme_classic()+theme(text = element_text(size=10),
                        axis.text.x = element_text(angle=45, hjust=1,size=10))+
  scale_fill_manual(values=c("midnightblue", "dodgerblue", "skyblue"))+
  labs(y="Variance components", x = " ")+
  scale_x_discrete(labels=c("Leaf","Stem bark", "Stem wood","Trunk bark","Trunk wood"))+
  guides(fill=guide_legend(title=" "))
#+facet_wrap(~trait,nrow=5,scales="free")

library(cowplot)
png("Appendix2.tiff",res=600,height=11,width=14,units="cm",pointsize = 1,bg="white",type=c("cairo"))
plot_grid(g1,
          labels=c(""),axis=c("r"),nrow=1,ncol=1,label_size = 12)
dev.off()

#Venables, W. N., & Ripley, B. D. (2002). Random and mixed effects. In Modern applied statistics with S (pp. 271-300). Springer, New York, NY.



#### Step 3: One-way ANOVA -----------

model <- lm(data = dsA, A ~ Spp)

model2 <- lm(data = dsA, log10(A) ~ Spp)
#test
res <- Anova(model,type=3)
res
#drop1(model,test="Chisq")
ano <- summary(aov(data = dsA, A ~ Spp))
ano

shapiro.test(residuals(model2))
hist(dsA$A)
hist(log10(dsA$A)) #do for all variables

leveneTest(model)

plot(model)


#refgrid<-ref.grid(model)
#summary(refgrid)


#### Step 4: Merge Gas Exchange data with Water Potential and Cs data -----------
dsA4 <- merge(dsAm, dsA2, by="Spp")
#### Step 5: Regression analysis -----------
regr <- read.csv("D://R//RCodesCNR1//Input//CNRSppGasExPsiCsMeans.csv")
regr <- regr[regr$Spp != "EL",]

M <- cor(regr)
res1 <- cor.mtest(regr, conf.level = .95)
corrplot(M, type = "upper", tl.pos = "td", method = "number", p.mat = res1$p, number.cex=0.6, cl.pos="r")

######### plot the graph
ggplot(regr, aes(x= CsPD, y= gsw))+
  #stat_summary(fun.y=mean,geom="point",aes(group=Geno),size=2)+
  geom_point()+
  geom_smooth(method="lm", se=F)+
  # xlim(0, 1500)+
  # geom_point()+
  # geom_smooth(method="lm", se=F, colour="Black")+
  # stat_poly_eq(formula = y~x, aes(label = paste(..rr.label.., sep = "~~~")),
  #              # stat_poly_eq(formula = y~x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE,
  #              label.x.npc ="left",
  #              label.y.npc = 0.9, size=5)+
  # stat_fit_glance(method = 'lm',method.args = list(formula = y~x), geom = 'text',
  #                 aes(label = paste("p = ", signif(..p.value.., digits = 4), sep = "")),
  #                 label.x.npc ="left",
#                 label.y.npc = 0.8, size=5)+
# geom_abline(intercept=0,linetype="longdash",color="Maroon",size=1.5)+
# geom_hline(yintercept=0, linetype="dashed",color="Maroon",size=1.5)+
#labs(x=expression("Photsynthetic assimilation ("*mu*"mol/cm"^2*"/sec)"),y="RGR (per day)")+#, colour = "Treatment")+
#labs(x=expression("log"[10]*"(Root dry mass)"), y=expression("log"[10]*"(Leaf dry mass)"))+#~Delta*"RGR (in g/g/day)")+
# labs(y=expression("RGR"[stress]*" (g/g/day)"), x=expression("RGR"[max]*" (g/g/day)"))+#,y=~Delta*"RGR (g/g/day)")+
# labs(y="Capacitance at pre-dawn / kg/m3", x="Predawn water potential / MPa")+#,y=~Delta*"RGR (g/g/day)")+
theme_bw(base_size=35)+
  # scale_x_continuous(limits=c(0.17,NA))+
  # scale_y_continuous(limits=c(NA,0.266))+
  theme(axis.title=element_text(size=27),
        #axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        # axis.text.x = element_text(size=25),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent"),
        # legend.position = "None",
        #legend.title=element_text(size=25),
        #legend.text=element_text(size=18),
        #legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        panel.border=element_rect(size=3))

ggsave("D:\\R\\RCodesCNR1\\Output\\CNRCsPDPsiPDPlot.pdf", width = 20, height = 20, units = "cm")


pval.lm = lm(data$RGR ~ data$W1, data=data)
#pval.lm = lmList(data$TotalH2 ~ data$AT_RGRls | data$Trt, data=data)
summary(pval.lm)
#dev.off()


# Vcmax function (one-pt method) ====

#The equation for calculation Vcmax using one-point method
# Asat : net photosynthesis under ambient CO2 and saturated irradiance condition, µmol CO2 m-² s-1
# Tl: leaf temperature at Asast, in degree
# Ci :internal leaf CO2 at Asat, in ppm
# Rd/Ra: leaf mitochondrial respiration in the light, µmol CO2 m-² s-1
# Tref: reference temperature, it can be 25 or other wanted temperature

Vcmax_calc <- function(data,Tref=25,Rd=FALSE){
  Asat <- data$A
  Ci <- data$Ci
  Tl <- data$Tleaf
  Ra <- data$Rabs
  R<-8.314
  O2<-210000 #umol/mol
  T0<-273.15
  #### Ha from in Bernacchi et al., 2001 in unit J·mol-1
  H_kc <-79430
  H_ko <-36380
  H_gstar <- 37830
  ####k25 from Bernacchi et al., 2001 in μmol·mol-1
  k25_kc <- 404.9
  k25_ko <- 278.4*1000
  k25_gstar <-42.75
  ####Kinetic constants at 25 °C (Table 2.3, von Caemmerer (2000)):
  #####################################################
  Kc<-k25_kc*exp(H_kc*(Tl-Tref)/(R*(Tl+T0)*298.15))
  Ko<-k25_ko*exp(H_ko*(Tl-Tref)/(R*(Tl+T0)*298.15))
  gStar<-k25_gstar*exp(H_gstar*(Tl-Tref)/(R*(Tl+T0)*298.15))
  
  Km<-Kc*(1+O2/Ko)
  
  if(Rd == TRUE){
    Vcmax<-(Asat+Ra)*((Ci+Km)/(Ci-gStar))   
  }
  if(Rd == FALSE){
    Vcmax<-(Asat)*((Ci+Km)/(Ci-gStar)-0.015)
  }
  
  # bernacchi equation and paramaters
  Vcmax_25<-Vcmax/exp(65330*(Tl-Tref)/(R*(Tl+T0)*298.15))
  
  # mdelyn equation and paramaters
  #kt<-exp(65330*(Tl-Tref)/(R*(Tl+TKR)*(Tref+TKR)))
  #kpeak<-(1+exp(((Tref+TKR)*656-200000)/(R*(Tref+TKR))))/(1+exp(((Tl+TKR)*656-200000)/(R*(Tl+TKR))))
  # Vcmax_std_med<-Vcmax/(kt*kpeak)
  
  Vcmaxfile <- data.frame(data$Spp, data$Rep, Vcmax, Vcmax_25)
}

dsAm$Vcmax <- NA
dsAm$Vcmax <- dsAm2$Vcmax
dsAm$Vcmax25 <- NA
dsAm$Vcmax25 <- dsAm2$Vcmax_25
# ====

dsAm2 <- Vcmax_calc(data=dsAm)
dsAm2$Vcmax[dsAm2$Vcmax > 500] <- NA
dsAm2$Vcmax_25[dsAm2$Vcmax_25 > 500] <- NA

plot(dsAm$gsw, dsAm$Vcmax)
plot(dsAm$gsw, dsAm$Vcmax25)

dsA1 <- dsAm %>% group_by(Spp) %>%
  summarise_if(is.numeric, .fun = mean, na.rm = TRUE)  # Calculating species means

plot(dsA1$gsw, dsA1$Vcmax)
plot(dsA1$gsw, dsA1$Vcmax25)


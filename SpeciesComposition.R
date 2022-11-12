library(tidyverse)
library(picante)


getwd()
setwd("D:/ECNU/Lilab/AFEC-X/PPT form Professors/Aki/Practice")
tree=read.tree("dummy_tree.newick")
tree
plot(tree)
trait=read.csv("dummy_trait.csv")
trait
summary(trait)
soil=read.csv("soil_env_use.csv")
soil
summary(soil)

samp=read.csv("sp data.csv",header = T,row.names = 1)
summary(samp)
samp

env=read.csv("env data.csv")
summary(env)
head(env)

res_mds <- metaMDS(samp)
plot(res_mds)
ordiplot(res_mds, type = "n")
orditorp(res_mds, display = "species", col = "red", air = 0.01)
orditorp(res_mds, display = "sites", cex = 1, air = 0.01)

NMDS.inv<-metaMDS(samp, distance="bray",autotransform =F)
NMDS.inv
str(NMDS.inv)
plot(NMDS.inv)

#using ggord to generate the ordination instead of the above
#extract NMDS scores (x and y coordinates)
data.scores <- as.data.frame(scores(NMDS.inv))
summary(data.scores)

#add columns to data frame 
data.scores$site <- env$ForestType
data.scores$treat = as.factor(env$Elevation)

#ggplot
(plot.2 <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
    geom_point(size = 4, aes( shape = site, colour = treat))+ 
    labs(x = "NMDS1", colour = "sample", y = "NMDS2", shape = "Type")  + 
    scale_colour_manual(values = c("red", "blue","green","yellow"))+ 
    theme(axis.text.y = element_text(colour = "black", size = 12), 
          axis.text.x = element_text(colour = "black", size = 12), 
          axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA,size = 1.2)
    )
)

#install.packages("ggfortify")
#installed.packages("ggvegan")
library(ggfortify)
library(ggvegan)
library(devtools)
#devtools::install_github("gavinsimpson/ggvegan")


spsite=samp
spsite$Elevation=env$Elevation
spsite$forest=env$ForestType
head(spsite)

sumNMDS<-NMDS.inv$points
sumNMDS<-as.data.frame(sumNMDS)
sumNMDS$Elevation <- env$Elevation
sumNMDS$forest <- env$ForestType

##autoplot can plot more beautiful plots####
autoplot(NMDS.inv, data = spsite, layers = c("biplot"))+ 
  geom_point(data=sumNMDS,aes(x=MDS1,y=MDS2,color=Elevation,shape=forest),size=4)+
  theme_light()+ 
  theme(
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA))

##and can also use ggplot to add the points of the species####
sumspecies <-as.data.frame (NMDS.inv$species)
ggplot(data=sumspecies,aes(x=MDS1,y=MDS2))+geom_point(size=3,color="pink") +
  geom_point(data=sumNMDS,aes(x=MDS1,y=MDS2,color=Elevation,shape=forest),size=4)+theme_light()+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))

#ggsave("D:/ECNU/Lilab/AFEC-X/PPT form Professors/Aki/Practice/NMDS.pdf",width=6,height=4.5, dpi=600,units="in")


#Statistical tests PERMANOVA
#set permutation
permAll <- how(nperm = 9999)
permAll

treat=env[,c(7:11)]
treat
samp
inv.ad<-adonis2(samp~., data=treat, permutations = permAll) 
inv.ad


### Ordistep function to chose the best model. We will first use "backward" to be conservative
summary(samp)
env.data=env[,-c(1:6)]
summary(env.data)

env.data2=env.data[,c(3,5)]


inv.rda<-capscale(samp~., data = env.data, distance="bray", scale=T)
summary(inv.rda)
plot(inv.rda, display = c("wa","cn"))
anova(inv.rda) #model not sig explain the variation in the multivariate community data


inv.rda1<-capscale(samp~1, data=env.data, distance="bray", scale=T) #need this if stepwise selection is used
inv.rda1

ordistep.inv.rda<-ordistep(inv.rda, permutations = permAll ,direction = "backward") 
ordistep.inv.rda
plot(ordistep.inv.rda) 

sumrda <- as.data.frame(unlist(summary(inv.rda)$sites))
sumrda$Elevation <- env$Elevation
sumrda$forest <- env$ForestType
autoplot(inv.rda, data = spsite, layers = c("biplot"))+ geom_point(data=sumrda,aes(x=CAP1,y=CAP2,color=Elevation,shape=forest),size=4)+theme_light()+ theme(
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA))
#ggsave("dbRDA.all.factors.pdf",width=6,height=4.5,dpi=600,units="in")

#delete un-significant factors####
##RDA2####
inv.rda2<-capscale(samp~., data = env.data2, distance="bray", scale=T)
summary(inv.rda2)
plot(inv.rda2, display = c("wa","cn"))
anova(inv.rda2,by="term",permutations = 2000) #model not sig explain the variation in the multivariate community data
anova(inv.rda,by="term",permutations = 2000)

sumrda2 <- as.data.frame(unlist(summary(inv.rda2)$sites))
sumrda2$Elevation <- env$Elevation
sumrda2$forest <- env$ForestType
autoplot(inv.rda2, data = spsite, layers = c("biplot"))+ 
  geom_point(data=sumrda2,aes(x=CAP1,y=CAP2,color=Elevation,shape=forest),size=4)+
  theme_light()+ 
  theme(
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA))

ggsave("dbRDA.pdf",width=6,height=4.5,dpi=600,units="in")

#######################################################################################

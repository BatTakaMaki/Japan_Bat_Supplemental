library(tidyverse)
library(betapart)
library(vegan)
library(ade4)
library(adespatial)
library(adegraphics)
library(spdep)
library(rgdal)
library(sf)
library(gridExtra)


setwd("C:/Users/uni23/Documents/Japanese_island_bat")

#Read distribution and factor data
RDA_dis1<- read.csv("bat_dist_Entire.csv", row.names = 1)
env<- read.csv("env_Entire.csv", row.names = 1)
his<- read.csv("his_Entire.csv", row.names = 1)
coo<- st_read("coo_Entire.shp")


#Calculate beta diversity
beta.dist <- beta.pair(RDA_dis1, index.family = "sorensen")
sor <- beta.dist$beta.sor #betasor
sim <- beta.dist$beta.sim #betasim
sne <- beta.dist$beta.sne #betasne
matsor<- as.matrix(sor) #matrix data of betasor

#############################################
#Calculate MEM based on Monadjem et al. 2022#
#############################################

## Useful function from monadjem et al. 2022
# *****************
nb2ggplot <- function(nb, coord) {
  # 'coord' must be a matrix/dataframe with two columns (called "long" and "lat")
  # 'nb' is an object of class nb
  # take out the connections from the nb object and assign them the lat and long in a dataframe
  n <- length(attributes(nb$neighbours)$region.id)
  DA <- data.frame(
    from = rep(1:n, sapply(nb$neighbours, length)),
    to = unlist(nb$neighbours),
    weight = unlist(nb$weights)
  )
  DA <- cbind(DA, coord[DA$from, 1:2], coord[DA$to, 1:2])
  colnames(DA)[4:7] = c("long", "lat", "long_to", "lat_to")
  return(DA)
}
set.seed(2020)

bat_xy<- as.data.frame(st_coordinates(x=coo)) %>% 
  rename(long = X, lat = Y)
row.names(bat_xy)<- coo$id


#Construct graph-based connectivity schemes
bat.nbtri <- tri2nb(bat_xy) # Delaunay triangulation
bat.nbgab <- graph2nb(gabrielneigh(bat_xy), sym = TRUE) # Gabriel
bat.nbrel <- graph2nb(relativeneigh(bat_xy), sym = TRUE) # Relative neighbourhood
bat.nbmst <- mst.nb(dist(bat_xy)) # minimum spanning tree

#Transform the nb object into listw objects for spdep packages
bat.nbtri_listw <- nb2listw(bat.nbtri) 
bat.nbgab_listw <- nb2listw(bat.nbgab)
bat.nbrel_listw <- nb2listw(bat.nbrel)
bat.nbmst_listw <- nb2listw(bat.nbmst)

#Create the MEMs based on the network
bat.DA_tri <- nb2ggplot(bat.nbtri_listw, bat_xy)
bat.tri_g <- bat_xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) +
  geom_point() +
  coord_fixed() +
  geom_segment(
    data = bat.DA_tri,
    aes(xend = long_to, yend = lat_to),
    size = 0.3,
    alpha = 0.5,
    colour = "darkred"
  ) +
  labs(title = "Delaunay triangulation") +
  theme_bw()

bat.DA_gab <- nb2ggplot(bat.nbgab_listw, bat_xy)
bat.gab_g <- bat_xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) +
  geom_point() +
  coord_fixed() +
  geom_segment(
    data = bat.DA_gab,
    aes(xend = long_to, yend = lat_to),
    size = 0.3,
    alpha = 0.5,
    colour = "darkred"
  ) +
  labs(title = "Gabriel") +
  theme_bw()

bat.DA_rel <- nb2ggplot(bat.nbrel_listw, bat_xy)
bat.rel_g <- bat_xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) +
  geom_point() +
  coord_fixed() +
  geom_segment(
    data = bat.DA_rel,
    aes(xend = long_to, yend = lat_to),
    size = 0.3,
    alpha = 0.5,
    colour = "chocolate"
  ) +
  labs(title = "Relative neighbourhood") +
  theme_bw()

bat.DA_mst <- nb2ggplot(bat.nbmst_listw, bat_xy)
bat.mst_g <- bat_xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) +
  geom_point() +
  coord_fixed() +
  geom_segment(
    data = bat.DA_mst,
    aes(xend = long_to, yend = lat_to),
    size = 0.3,
    alpha = 0.5,
    colour = "goldenrod"
  ) +
  labs(title = "Min. span tree") +
  theme_bw()

grid.arrange(bat.tri_g, bat.gab_g, bat.rel_g, bat.mst_g, ncol = 2, nrow = 2)

# Binary forms (no weighths added to the connections):
bat.nbtri_edited_b <- nb2listw(bat.nbtri) # changed names here
bat.nbgab_edited_b <- nb2listw(bat.nbgab) # changed names here
bat.nbrel_edited_b <- nb2listw(bat.nbrel) # changed names here
bat.nbmst_edited_b <- nb2listw(bat.nbmst) # changed names here

# Construction of the list of candidate spatial weighting matrices:
# *****************************************************************
# Four candidates
bat.candidates <- list(bat.tri_b = bat.nbtri_edited_b,
                       bat.gab_b = bat.nbgab_edited_b,
                       bat.rel_b = bat.nbrel_edited_b,
                       bat.nbm_b = bat.nbmst_edited_b)

# Optimisation of a subset of spatial predictors
# optimise the selection of a subset of spatial predictors from the
# best-suited SWM.
bat.select <-
  listw.select(
    RDA_dis1,
    bat.candidates,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )

# Optimised selected SWM:
bat.select$best.id

# Summarises the candidate SWMs:
bat.select$candidates

# Summarises MEM var. selected and Moran's I and p-val in residuals after inclusion of the MEM var.
bat.select$best$summary

#best-suited SWM based on beta sor
bat.select.sor <-
  listw.select(
    matsor,
    bat.candidates,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )

# MEM variables to use in models and simulations:
# ***********************************************
spa<-bat.select.sor$best$MEM.select

#Variation partitioning based on dbRDA
varsor<-varpart(sor,spa, env, his)
varsim<-varpart(sim,spa, env, his)
varsne<-varpart(sne,spa, env, his)

plot(varsor)
plot(varsim)
plot(varsne)

#Extract adjusted R square
sorfract<-varsor$part$fract$Adj.R.square
sorindfract<-varsor$part$indfract$Adj.R.square
simfract<-varsim$part$fract$Adj.R.square
simindfract<-varsim$part$indfract$Adj.R.square
snefract<-varsne$part$fract$Adj.R.square
sneindfract<-varsne$part$indfract$Adj.R.square


#Permutation with variation partitioning
RDA_ran<-RDA_dis1[sample(nrow(RDA_dis1), replace = FALSE),]
beta.dist2 <- beta.pair(RDA_ran, index.family = "sorensen")
sor2 <- beta.dist2$beta.sor
sim2 <- beta.dist2$beta.sim
sne2 <- beta.dist2$beta.sne

#Variation partitining analysis on permuted distribution data
varran<- varpart(sor2,spa, env, his)
varransim<- varpart(sim2,spa, env, his)
varransne<- varpart(sne2,spa, env, his)

#Extract adjusted R square
rfract<-varran$part$fract %>% 
  dplyr::select(Adj.R.square)
rfractsim<-varransim$part$fract %>% 
  dplyr::select(Adj.R.square)
rfractsne<-varransne$part$fract %>% 
  dplyr::select(Adj.R.square)

rindfract<-varran$part$indfract %>% 
  dplyr::select(Adj.R.square)
rindfractsim<-varransim$part$indfract %>% 
  dplyr::select(Adj.R.square)
rindfractsne<-varransne$part$indfract %>% 
  dplyr::select(Adj.R.square)

#Additional permutation with variation partitioning with 9998 iteration
for (i in 1:9998){
  RDA_ran<-RDA_dis1[sample(nrow(RDA_dis1), replace = FALSE),]
  beta.dist2 <- beta.pair(RDA_ran, index.family = "sorensen")
  sor2 <- beta.dist2$beta.sor
  sim2 <- beta.dist2$beta.sim
  sne2 <- beta.dist2$beta.sne
  
  varran<- varpart(sor2,spa, env, his)
  varransim<- varpart(sim2,spa, env, his)
  varransne<- varpart(sne2,spa, env, his)
  
  rfract<-cbind(rfract, varran$part$fract$Adj.R.square)
  rfractsim<-cbind(rfractsim, varransim$part$fract$Adj.R.square)
  rfractsne<-cbind(rfractsne, varransne$part$fract$Adj.R.square)
  
  rindfract<-cbind(rindfract, varran$part$indfract$Adj.R.square)
  rindfractsim<-cbind(rindfractsim, varransim$part$indfract$Adj.R.square)
  rindfractsne<-cbind(rindfractsne, varransne$part$indfract$Adj.R.square)
}  

#Arrange the data of R square
rfractsor_t<-as.data.frame(t(rfract))
rfractsim_t<-as.data.frame(t(rfractsim))
rfractsne_t<-as.data.frame(t(rfractsne))

colnames(rfractsor_t)<- c("S", "E", "H", "S_E", "S_H", "E_H", "All")
colnames(rfractsim_t)<- c("S", "E", "H", "S_E", "S_H", "E_H", "All")
colnames(rfractsne_t)<- c("S", "E", "H", "S_E", "S_H", "E_H", "All")

rindfractsor_t<-as.data.frame(t(rindfract))
rindfractsim_t<-as.data.frame(t(rindfractsim))
rindfractsne_t<-as.data.frame(t(rindfractsne))

colnames(rindfractsor_t)<- c("indS", "indE", "indH", "shareS_E",
                             "shareE_H", "shareS_H", "shareS_E_H", "residuals")
colnames(rindfractsim_t)<- c("indS", "indE", "indH", "shareS_E",
                             "shareE_H", "shareS_H", "shareS_E_H", "residuals")
colnames(rindfractsne_t)<- c("indS", "indE", "indH", "shareS_E",
                             "shareE_H", "shareS_H", "shareS_E_H", "residuals")

#Test observed R square greater than the value from permutation
sorS<- (10000-(sum(rfractsor_t[,1]<(sorfract[1]))))/10000
sorE<- (10000-(sum(rfractsor_t[,2]<(sorfract[2]))))/10000
sorH<- (10000-(sum(rfractsor_t[,3]<(sorfract[3]))))/10000
sorS_E<- (10000-(sum(rfractsor_t[,4]<(sorfract[4]))))/10000
sorS_H<- (10000-(sum(rfractsor_t[,5]<(sorfract[5]))))/10000
sorE_H<- (10000-(sum(rfractsor_t[,6]<(sorfract[6]))))/10000
sorAll<- (10000-(sum(rfractsor_t[,7]<(sorfract[7]))))/10000

sorindS<- (10000-(sum(rindfractsor_t[,1]<(sorindfract[1]))))/10000
sorindE<- (10000-(sum(rindfractsor_t[,2]<(sorindfract[2]))))/10000
sorindH<- (10000-(sum(rindfractsor_t[,3]<(sorindfract[3]))))/10000
sorshareS_E<- (10000-(sum(rindfractsor_t[,4]<(sorindfract[4]))))/10000
sorshareE_H<- (10000-(sum(rindfractsor_t[,5]<(sorindfract[5]))))/10000
sorshareS_H<- (10000-(sum(rindfractsor_t[,6]<(sorindfract[6]))))/10000
sorshareS_E_H<- (10000-(sum(rindfractsor_t[,7]<(sorindfract[7]))))/10000
sorresiduals<- (10000-(sum(rindfractsor_t[,8]<(sorindfract[8]))))/10000


simS<- (10000-(sum(rfractsim_t[,1]<(simfract[1]))))/10000
simE<- (10000-(sum(rfractsim_t[,2]<(simfract[2]))))/10000
simH<- (10000-(sum(rfractsim_t[,3]<(simfract[3]))))/10000
simS_E<- (10000-(sum(rfractsim_t[,4]<(simfract[4]))))/10000
simS_H<- (10000-(sum(rfractsim_t[,5]<(simfract[5]))))/10000
simE_H<- (10000-(sum(rfractsim_t[,6]<(simfract[6]))))/10000
simAll<- (10000-(sum(rfractsim_t[,7]<(simfract[7]))))/10000

simindS<- (10000-(sum(rindfractsim_t[,1]<(simindfract[1]))))/10000
simindE<- (10000-(sum(rindfractsim_t[,2]<(simindfract[2]))))/10000
simindH<- (10000-(sum(rindfractsim_t[,3]<(simindfract[3]))))/10000
simshareS_E<- (10000-(sum(rindfractsim_t[,4]<(simindfract[4]))))/10000
simshareE_H<- (10000-(sum(rindfractsim_t[,5]<(simindfract[5]))))/10000
simshareS_H<- (10000-(sum(rindfractsim_t[,6]<(simindfract[6]))))/10000
simshareS_E_H<- (10000-(sum(rindfractsim_t[,7]<(simindfract[7]))))/10000
simresiduals<- (10000-(sum(rindfractsim_t[,8]<(simindfract[8]))))/10000


sneS<- (10000-(sum(rfractsne_t[,1]<(snefract[1]))))/10000
sneE<- (10000-(sum(rfractsne_t[,2]<(snefract[2]))))/10000
sneH<- (10000-(sum(rfractsne_t[,3]<(snefract[3]))))/10000
sneS_E<- (10000-(sum(rfractsne_t[,4]<(snefract[4]))))/10000
sneS_H<- (10000-(sum(rfractsne_t[,5]<(snefract[5]))))/10000
sneE_H<- (10000-(sum(rfractsne_t[,6]<(snefract[6]))))/10000
sneAll<- (10000-(sum(rfractsne_t[,7]<(snefract[7]))))/10000

sneindS<- (10000-(sum(rindfractsne_t[,1]<(sneindfract[1]))))/10000
sneindE<- (10000-(sum(rindfractsne_t[,2]<(sneindfract[2]))))/10000
sneindH<- (10000-(sum(rindfractsne_t[,3]<(sneindfract[3]))))/10000
sneshareS_E<- (10000-(sum(rindfractsne_t[,4]<(sneindfract[4]))))/10000
sneshareE_H<- (10000-(sum(rindfractsne_t[,5]<(sneindfract[5]))))/10000
sneshareS_H<- (10000-(sum(rindfractsne_t[,6]<(sneindfract[6]))))/10000
sneshareS_E_H<- (10000-(sum(rindfractsne_t[,7]<(sneindfract[7]))))/10000
sneresiduals<- (10000-(sum(rindfractsne_t[,8]<(sneindfract[8]))))/10000

sorfract_RP<- c(sorS, sorE, sorH, sorS_E, sorS_H, sorE_H, sorAll)
sorindfract_RP<- c(sorindS, sorindE, sorindH, sorshareS_E, sorshareE_H, sorshareS_H, sorshareS_E_H, sorresiduals)

simfract_RP<- c(simS, simE, simH, simS_E, simS_H, simE_H, simAll)
simindfract_RP<- c(simindS, simindE, simindH, simshareS_E, simshareE_H, simshareS_H, simshareS_E_H, simresiduals)

snefract_RP<- c(sneS, sneE, sneH, sneS_E, sneS_H, sneE_H, sneAll)
sneindfract_RP<- c(sneindS, sneindE, sneindH, sneshareS_E, sneshareE_H, sneshareS_H, sneshareS_E_H, sneresiduals)
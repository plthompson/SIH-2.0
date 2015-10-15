#Metacommunity parameters####
regions<-10 #number of regions
patches<-8 #number of local patches in each region
tot_patches<-regions*patches

species<-80 #number of species

dispV<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1)
landscape_dispersal_V<-c(0,0.25,0.5,0.75,1)

#Environmental parameters####
ePeriod<-1000
eAMP<-1
Tmax<-ePeriod*3
env_gradient<-5 #magintude of the environmental gradient
env_increase<-2

sync<-seq(0,1,length=5)


loc_env<-array(NA,dim=c(Tmax,tot_patches,length(sync)))
for(s in 1:length(sync)){
  for(TS in 1:Tmax){
    loc_env[TS,,s]<-0.5*eAMP*(sin((2*pi/ePeriod)*TS+1+(seq(0,(1-sync[s])*(patches-1),length=patches)+1)*2*pi/patches)+1)
  }}

par(mfrow=c(3,2))
for(s in 1:length(sync)){
matplot(loc_env[1:ePeriod*2,,s], type='l',lty=1, col=1:patches)
}
par(mfrow=c(1,1))

reg_env<-matrix(NA, Tmax,tot_patches)
for(TS in 1:Tmax){
  reg_env[TS,]<-0.5*eAMP*(sin((2*pi/ePeriod*1.5)*TS+1+(rep(1:regions,each=patches))*2*pi/(regions))+1)
}

matplot(reg_env[1:ePeriod*2,], type='l',lty=1, col=1:10)

#environmental increase
landscape_increase<-matrix(seq(0,env_increase, length=Tmax),Tmax,tot_patches)
#environmental gradient
gradient<-matrix(rep(seq(0,(env_gradient-1),length=regions),each=Tmax*patches),Tmax,tot_patches)
#Overall environment
#Env<-gradient+landscape_increase+loc_env[,,s]+reg_env
Env<-loc_env[,,4]+gradient

matplot(Env[1:3000,],type='l', lty=1, col=1:patches)



#community####
Com_types<-c("NoInt","Comp","Mixed","Foodweb")
BB<-array(0,dim=c(species, species, length(Com_types)),dimnames=list(1:species,1:species,Com_types))

#no interactions
diag(BB[,,"NoInt"])<--.2

#competition
b11=-.15
BB[,,"Comp"]=b11*matrix(runif(species*species),species,species)
diag(BB[,,"Comp"])<--.2

#mixed
hold<-matrix(-1,species,species)
int.n<-sum(hold[upper.tri(hold)])*-1
hold[upper.tri(hold)][sample(int.n, replace=F)<=(0.35*int.n)]<-0.5
hold[lower.tri(hold)][t(hold)[lower.tri(hold)]>0][sample(0.35*int.n, replace=F)<(0.10*int.n)]<-0.5
BB[,,"Mixed"]<-hold*-BB[,,"Comp"]

#no foodweb at the moment

#weight by species number
weight=1/species*3
BB<-BB*weight

#intrinsic rate of increase
C<-0.1

#environmental response
T_Opt<-seq(min(Env[1,]),max(Env[1,]),length=species)


#dispersal####
#local
loc_dispersal_matrix<-matrix(0, tot_patches,tot_patches)
for(reg in 1:regions){
  loc_dispersal_matrix[(((1:patches)+(reg-1)*patches)),(((1:patches)+(reg-1)*patches))]<-1/(patches-1)
}
diag(loc_dispersal_matrix)<-0

#regional
reg_dispersal_matrix<-matrix(0, tot_patches,tot_patches)
for(reg in 2:(regions-1)){
  reg_dispersal_matrix[c(((1:patches)+(patches*reg)-patches*2),((1:patches)+(patches*reg))),(1:patches)+(patches*(reg-1))]<-0.5/patches
}
reg_dispersal_matrix[1:(2*patches),1:patches]<-0.5/patches
reg_dispersal_matrix[(tot_patches-(2*patches-1)):tot_patches,(tot_patches-(patches-1)):tot_patches]<-0.5/patches


#data storage####
Spatial.Insurance<-data.frame(SR=NA,Biomass=NA,Biomass_CV=NA,Dispersal=dispV,Local_sync=rep(sync,each=length(dispV)),Scale=rep(c("Local","Regional","Landscape"),each=length(dispV)*length(sync)),Landscape_dispersal=rep(landscape_dispersal_V,each=length(dispV)*length(sync)*3))

#abundance matrix
X<-array(10,dim=c(species,Tmax,tot_patches))

for(land_disp in 1:length(landscape_dispersal_V)){
for(s in 1:length(sync)){
  Env<-loc_env[,,s]+gradient
  matplot(Env[1:3000,],type='l', lty=1, col=1:patches)
for(d in 1:length(dispV)){
  dispersal_matrix<-loc_dispersal_matrix*(1-landscape_dispersal_V[land_disp])+reg_dispersal_matrix*landscape_dispersal_V[land_disp]
  for(l in 1:(Tmax-1)){
  for(com in 2){
  Env_resp<--abs(matrix(rep(Env[l,],each=species),species,tot_patches)-T_Opt)*0.2
  X[,l+1,]<-X[,l,]*exp(C+BB[,,com]%*%X[,l,]+Env_resp)+(X[,l,]%*%dispersal_matrix)*dispV[d]-X[,l,]*dispV[d]
  X[,l+1,][(X[,l+1,]<0.01)]<-0
  }}
  matplot(t(X[,-c(1:500),1]),type='l',col=1)
  
  hold3<-hold2<-hold<-rep(NA,regions)
  for(reg in 1:regions){
    hold[reg]<-mean(colSums(apply(X[,-c(1:2000),(1:patches)+patches*(reg-1)],2,rowSums)>0))
    #hold2[reg]<-mean(rowSums(apply((Prod_sample[,,(1:patches)+patches*(reg-1)]),3,rowSums)),na.rm=T)
    #hold3[reg]<-sd(rowSums(apply((Prod_sample[,,(1:patches)+patches*(reg-1)]),3,rowSums)),na.rm=T)/mean(rowSums(apply((Prod_sample[,,(1:patches)+patches*(reg-1)]),3,rowSums)))
  }
  
  Spatial.Insurance$SR[Spatial.Insurance$Dispersal==dispV[d] & Spatial.Insurance$Local_sync==sync[s] & Spatial.Insurance$Scale=="Local" & Spatial.Insurance$Landscape_dispersal == landscape_dispersal_V[land_disp]] <-mean(apply(X[,-c(1:2000),]>0,3,colSums))
  Spatial.Insurance$SR[Spatial.Insurance$Dispersal==dispV[d] & Spatial.Insurance$Local_sync==sync[s] & Spatial.Insurance$Scale=="Regional" & Spatial.Insurance$Landscape_dispersal == landscape_dispersal_V[land_disp]] <-mean(hold)
  Spatial.Insurance$SR[Spatial.Insurance$Dispersal==dispV[d] & Spatial.Insurance$Local_sync==sync[s] & Spatial.Insurance$Scale=="Landscape" & Spatial.Insurance$Landscape_dispersal == landscape_dispersal_V[land_disp]] <-mean(colSums(apply(X[,-c(1:2000),],2,rowSums)>0))
  }}}

require(ggplot2)

ggplot(Spatial.Insurance,aes(x=Dispersal,y=SR,color=Scale))+
  geom_line(size=2)+
  facet_grid(Local_sync~Landscape_dispersal)+
  scale_x_log10()+
  theme_bw(base_size = 15)

require(dplyr)
group_by(Spatial.Insurance,Scale,Landscape_dispersal)
filter(Spatial.Insurance,Dispersal==0.0001)
max_SI<-summarize(group_by(Spatial.Insurance,Scale,Landscape_dispersal,Local_sync),max_SR=max(SR))
null_SR<-arrange(select(filter(group_by(Spatial.Insurance,Scale,Landscape_dispersal,Local_sync),Dispersal == 0.001),SR),Scale,Landscape_dispersal,Local_sync)

max_SI$max_SI<-max_SI$max_SR/null_SR$SR

ggplot(max_SI,aes(x=Local_sync,y=Landscape_dispersal,z=max_SI))+
  stat_contour(binwidth=0.25,aes(colour=..level..))+
  scale_color_continuous(low="dodgerblue",high="red")+
  theme_bw(base_size = 15)+
  facet_wrap("Scale")

ggplot(filter(max_SI, Scale=="Landscape"),aes(x=Local_sync,y=Landscape_dispersal,z=max_SI))+
  stat_contour(binwidth=0.05,aes(colour=..level..))+
  scale_color_continuous(low="dodgerblue",high="red")+
  theme_bw(base_size = 15)
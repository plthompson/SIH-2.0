patches<-30 #number of patches
species<-80 #number of species
dispV<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1)


#environment####
ePeriod<-1000
eAMP<-1
regional_gradient<-3
Tmax<-ePeriod*3
env_gradient<-5 #magintude of the environmental gradient
env_increase<-2

sync<-seq(0,1,length=5)

loc_env<-array(NA,dim=c(Tmax,patches,length(sync),2))
for(s in 1:length(sync)){
  for(TS in 1:Tmax){
    loc_env[TS,,s,1]<-0.5*eAMP*(sin((2*pi/ePeriod)*TS+1+(seq(0,(1-sync[s])*(patches-1),length=patches)+1)*2*pi/patches)+1)
  }}
loc_env[,,,2]<-loc_env[,,,1]+rep(seq(0,regional_gradient,length=patches),each=Tmax)

matplot(loc_env[,,5,2], type='l', lty=1)

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
T_Opt<-seq(0,max(loc_env),length=species)


#dispersal####
dispersal_mat<-matrix(1/(patches-1), patches,patches)
diag(dispersal_mat)<-0

#data storage####
Spatial.Insurance<-data.frame(SR=NA,Biomass=NA,Biomass_CV=NA,Dispersal=dispV,Local_sync=rep(sync,each=length(dispV)),Scale=rep(c("Local","Regional"),each=length(dispV)*length(sync)),Vary_mean=rep(c("No","Yes"),each=length(dispV)*length(sync)*2))

#abundance matrix
X<-array(10,dim=c(species,Tmax,patches))

for(vm in 1:2){
  for(s in 1:length(sync)){
    Env<-loc_env[,,s,vm]
    matplot(Env[1:3000,],type='l', lty=1, col=1:patches)
    for(d in 1:length(dispV)){
      for(l in 1:(Tmax-1)){
        for(com in 2){
          Env_resp<--abs(matrix(rep(Env[l,],each=species),species,patches)-T_Opt)*0.2
          X[,l+1,]<-X[,l,]*exp(C+BB[,,com]%*%X[,l,]+Env_resp)+(X[,l,]%*%dispersal_mat)*dispV[d]-X[,l,]*dispV[d]
          X[,l+1,][(X[,l+1,]<0.01)]<-0
        }}
      matplot(t(X[,-c(1:500),2]),type='l',col=1, main=paste(s,d,sep="_"))
      
      Spatial.Insurance$SR[Spatial.Insurance$Dispersal==dispV[d] & Spatial.Insurance$Local_sync==sync[s] & Spatial.Insurance$Scale=="Local" & Spatial.Insurance$Vary_mean == c("No","Yes")[vm]] <-mean(apply(X[,-c(1:2000),]>0,3,colSums))
      Spatial.Insurance$SR[Spatial.Insurance$Dispersal==dispV[d] & Spatial.Insurance$Local_sync==sync[s] & Spatial.Insurance$Scale=="Regional" & Spatial.Insurance$Vary_mean == c("No","Yes")[vm]] <-mean(colSums(apply(X[,-c(1:2000),],2,rowSums)>0))
    }}}

require(ggplot2)

ggplot(Spatial.Insurance,aes(x=Dispersal,y=SR,color=Scale))+
  geom_line(size=2)+
  facet_grid(Local_sync~Vary_mean)+
  scale_x_log10()+
  theme_bw(base_size = 15)
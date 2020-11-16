
require(tidyverse)
require(reshape2)
require(vegan) 
require(sp)
require(viridis)
require(sf)
require(geosphere)
proRD<-"+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.417,50.3319,465.552,-0.398957,0.343988,-1.8774,4.0725 +units=m +no_defs"
proWG<-CRS("+proj=longlat +datum=WGS84")
proUTM<-CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")
NL<-readRDS("../gadm36_NLD_0_sp.rds")

dataDir<-"./data/derived_data"
outDir <-"./output"

# data inladen
load(file.path(dataDir,"dats.Rdata"))
samps <- samps %>%
  mutate(area=substr(sttcd,4,5))

# The complete table has been transformed in readdate.R into a cross-table, 
# and then split into a stations part and a community matrix. 
# Summary statistics, such as total density and number of species, are added here to the stations part.

samps$totdens<-apply(abund,1,sum) 
samps$totldens<-apply(abund,1,FUN=function(x)sum(log(x+1)))
samps$nspec<-apply(abund,1,FUN = function(x) length(x[x>0]))

# we perform an exploratory PCA to compare the composition of the differ-
#  ent samples. We use different symbols and colors
# for the five areas. 

#remove infrequent species
freq<-apply(abund,2,FUN=function(x)length(x[x>0]))
abfre<-abund[,which(freq>10)]

dpca<-rda(log(abfre+1)) 
samps$ax1<-scores(dpca)$sites[,1]
samps$ax2<-scores(dpca)$sites[,2] 

# preparation for plots
depcol<-rainbow(5)
arpch<-tibble(area= as.character(c("ZH","NH","TX","AM","SC")),
              num = as.numeric(c(1,2,3,4,5)),
              pch = as.numeric(c(15,16,17,18,19)),
              cc  = as.character(depcol[c(1,2,3,4,5)]))
samps<- samps %>% left_join(arpch[,c("area","num")],by="area")

#some summary statistics
samps$areaf<-as.factor(samps$area)
summarie<-samps %>% group_by(areaf) %>% summarize(Mean_nspec=mean(nspec),Sd_nspec=sd(nspec),
                                       Mean_totdens=mean(totdens),Sd_totdens=sd(totdens),
                                       Mean_totldens=mean(totldens),Sd_totldens=sd(totldens),
                                       num=mean(num)) %>%
  arrange(by_group=num)
print(summarie)

## plotting section

pdf(file.path(outDir,"basplots.pdf")) 

# preparation
areas<-c(unique(samps$area),"All")
dums<-samps[1:2,]
dums[,]<-NA
rcls<-viridis_pal()(100)
plf<-function(title,ar,var,scal){
  if(scal<(max(var,na.rm=T)*1.05))scal<-(max(var,na.rm=T)*1.05)
  if(min(var,na.rm=T)<0)var<-var-min(var,na.rm=T)+0.01
  if(ar=="All")cexp<-1 else cexp<-2
  plot(arstats,cex=cexp, 
       main=paste(ar,title),
       pch=19,col=rcls[100*var/scal+1]) 
  #points(arstats,cex=1)
  plot(NL,add=T)
}
for (ar in areas){ 
  if(ar=="All")arstats<-samps else arstats<-samps[samps$area==ar,]
  dums[1,9:10]<-c(min(arstats$geometriepunt_x)-0.01,min(arstats$geometriepunt_y)-0.01)
  dums[2,9:10]<-c(max(arstats$geometriepunt_x)+0.01,max(arstats$geometriepunt_y)+0.01)
  arstats<-rbind(arstats,dums)
  coordinates(arstats)<- ~geometriepunt_x+geometriepunt_y
  proj4string(arstats)<-proWG

  plf("Tot_log_dens",ar,arstats$totldens,30)
  plf("Tot_dens",ar,arstats$totdens,5500)
  plf("Nspec",ar,arstats$nspec,5)
  plf("Axis1",ar,arstats$ax1+3,1)
  plf("Axis2",ar,arstats$ax2+3,1)
} 
dev.off() 

pdf(file.path(outDir,"specplots.pdf")) 

for (i in 1:ncol(abfre)){
    ar<-""
    arstats<-samps
    arstats$ab<-ifelse(abfre[,i]>0,log(abfre[,i]),NA)
    dums<-arstats[1:2,]
    dums[,]<-NA

    dums[1,9:10]<-c(min(arstats$geometriepunt_x)-0.01,min(arstats$geometriepunt_y)-0.01)
    dums[2,9:10]<-c(max(arstats$geometriepunt_x)+0.01,max(arstats$geometriepunt_y)+0.01)
    arstats<-rbind(arstats,dums)
    coordinates(arstats)<- ~geometriepunt_x+geometriepunt_y
    proj4string(arstats)<-proWG
    
    plf(names(abfre)[i],ar,arstats$ab,1)
}

dev.off()


# find the distance to the coast for all sample points
sfsamps<-samps
coordinates(sfsamps)<- ~geometriepunt_x+geometriepunt_y
proj4string(sfsamps)<-proWG
sfsamps<-st_as_sf(sfsamps)

sfNL<-st_as_sf(NL)

areas<-unique(samps$area)
samps$dist2Cst<-rep(NA,nrow(samps))
for(ar in areas){
  slst<-which(sfsamps$area==ar)
  xmin<-min(st_coordinates(sfsamps[sfsamps$area==ar,])[,1])-0.1
  xmax<-max(st_coordinates(sfsamps[sfsamps$area==ar,])[,1])+0.1
  ymin<-min(st_coordinates(sfsamps[sfsamps$area==ar,])[,2])-0.1
  ymax<-max(st_coordinates(sfsamps[sfsamps$area==ar,])[,2])+0.1
  cst<-st_crop(sfNL,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
  dd<-dist2Line(p=st_coordinates(sfsamps[slst,]),line=st_coordinates(cst)[,1:2])
  samps$dist2Cst[slst]<-dd[,1]
}

pdf(file.path(outDir,"factorplots.pdf"))

# plots of the PCA
plot(samps$ax1,samps$ax2,pch=arpch$pch[samps$num],
     col=depcol[samps$num]) 
legend("bottomleft",legend=arpch$area,pch=arpch$pch,col=arpch$cc)

freq<-apply(abfre,2,FUN=function(x)length(x[x>0]))
tt<-which(freq>20)
plot(dpca)
text(scores(dpca)$species[tt,1],scores(dpca)$species[tt,2],names(abfre)[tt],cex=0.5,col="red")


# make some plots versus distance from coast
plDst<-function(labx,laby,varx,vary,logx=FALSE){
  plot(varx,vary,main=laby,pch=arpch$pch[samps$num],
       col=depcol[samps$num],xlab=labx,ylab=laby,log=ifelse(logx,"x",""))
  legend("bottomright",legend=arpch$area,pch=arpch$pch,col=arpch$cc)
}
plDst("Distance to Coast (m)","Number of species",samps$dist2Cst,samps$nspec)
plDst("Distance to Coast (m)","SumLogDens",samps$dist2Cst,samps$totldens)
plDst("Distance to Coast (m)","Log10 Total Density",samps$dist2Cst,log(samps$totdens+1)/log(10))
plDst("Distance to Coast (m)","Axis 1 scores",samps$dist2Cst,samps$ax1)
plDst("Distance to Coast (m)","Axis 2 scores",samps$dist2Cst,samps$ax2)

# make some plots versus depth
plDst("Depth (m)","Number of species",samps$begindiepte_m,samps$nspec)
plDst("Depth (m)","SumLogDens",samps$begindiepte_m,samps$totldens)
plDst("Depth (m)","Log10 Total Density",samps$begindiepte_m,log(samps$totdens+1)/log(10))
plDst("Depth (m)","Axis 1 scores",samps$begindiepte_m,samps$ax1)
plDst("Depth (m)","Axis 2 scores",samps$begindiepte_m,samps$ax2)

plDst("D50 (um)","Number of species",samps$D50,samps$nspec,logx=TRUE)
plDst("D50 (um)","SumLogDens",samps$D50,samps$totldens,logx=TRUE)
plDst("D50 (um)","Log10 Total Density",samps$D50,log(samps$totdens+1)/log(10),logx=TRUE)
plDst("D50 (um)","Axis 1 scores",samps$D50,samps$ax1,logx=TRUE)
plDst("D50 (um)","Axis 2 scores",samps$D50,samps$ax2,logx=TRUE)

plDst("Perc <63 um","Number of species",samps$perc0_2+samps$perc2_63,samps$nspec)
plDst("Perc <63 um","SumLogDens",samps$perc0_2+samps$perc2_63,samps$totldens)
plDst("Perc <63 um","Log10 Total Density",samps$perc0_2+samps$perc2_63,log(samps$totdens+1)/log(10))
plDst("Perc <63 um","Axis 1 scores",samps$perc0_2+samps$perc2_63,samps$ax1)
plDst("Perc <63 um","Axis 2 scores",samps$perc0_2+samps$perc2_63,samps$ax2)

samps$lD50<-log(samps$D50)

plot(exp(samps$lD50),samps$begindiepte_m,cex=samps$ax1+2.5,col=depcol[samps$num],pch=19,log="x",
     main="Scores as 1 (grootte symbolen)",xlab="Mediaan (um)",ylab="Diepte(m)")
lines(c(exp(9.08/1.82),exp((15+9.08/0.144)*0.144/1.82)),c(0,15),lwd=2)
legend("bottomright",legend=arpch$area,pch=19,col=arpch$cc)
plot(exp(samps$lD50),samps$begindiepte_m,cex=samps$totldens/50,col=depcol[samps$num],pch=19,log="x",
     main="SumLogDens (grootte symbolen)",xlab="Mediaan (um)",ylab="Diepte(m)")
legend("bottomright",legend=arpch$area,pch=19,col=arpch$cc)

# abline(-9.08/0.144,1.82/0.144)
# for(A in seq(-10-5*1.82/0.144,+10-5*1.82/0.144,by=1)){
#   a=A
#   b=1.82/0.144
#   abline(a,b)
# }

samps <- samps %>% mutate(trnct=as.numeric(substr(sttcd,7,8))) %>%
  mutate(pos=as.numeric(substr(sttcd,9,9)))
trncts<-unique(samps$trnct)
trncts<-trncts[order(trncts)]
plot(samps$ax1,samps$ax2,asp=1,pch=NA)
for(t in trncts){
  pnts<-which(samps$trnct==t)
  pnts<-pnts[order(samps$pos[pnts])]
  lines(samps$ax1[pnts],samps$ax2[pnts],col=arpch$cc[unique(samps$num[pnts])])
  text(samps$ax1[pnts],samps$ax2[pnts],labels = samps$pos[pnts],cex=0.7)
}
legend("bottomright",legend=arpch$area,col=arpch$cc,lty=1)


splist<-specsacc$ScientificName_accepted[specsacc$AphiaID_accepted>0]
for(sp in splist){
  plot(as.factor(samps$area),log10(abund[,sp]+1),ylab="Log10(abundance+1)",main=sp)
}
dev.off()

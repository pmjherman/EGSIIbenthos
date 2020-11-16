require(tidyverse)
require(reshape2)
require(vegan) 
require(sp)
require(viridis)
require(sf)
require(geosphere)
proRD  <- "+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.417,50.3319,465.552,-0.398957,0.343988,-1.8774,4.0725 +units=m +no_defs"
proWG  <- CRS("+proj=longlat +datum=WGS84")
proUTM <- CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")
NL     <- readRDS("./data/gadm36_NLD_0_sp.rds")

dataDir <- "./data/derived_data"
outDir  <- "./output"

# data inladen
load(file.path(dataDir,"dats.Rdata"))


depcol <- rainbow(5,alpha = 1)
allsamps <- samps %>% select(sttcd)
allsamps <- allsamps %>% 
  mutate(stno = as.numeric(substr(sttcd,7,9)),
         area = ifelse(stno<200,1,
                    ifelse(stno<300,2,
                           ifelse(stno<400,3,
                                  ifelse(stno<500,4,5))))
         ) %>% arrange(by = stno)


specsweg <- c("Crangon crangon","Spisula solida","Spisula",
            "Chamelea striatula","Mactra stultorum","Bivalvia")
datsL <- datsL %>% filter(!ScientificName_accepted %in% specsweg)

pdf(file.path(outDir,"histograms.pdf"))

for(sp in unique(datsL$ScientificName_accepted)){
  tt <- datsL %>% filter(ScientificName_accepted == sp) %>%
    select(sttcd,numeriekewaarde) %>%
    left_join(allsamps,by = "sttcd") %>%
    mutate(st = paste0(area,substr(sttcd,4,5))) %>%
    mutate(`lengte (mm)` = numeriekewaarde)
  arlst <- unique(tt$area)
  arlst <- arlst[order(arlst)]
   pl <- ggplot(tt, aes(x = numeriekewaarde, fill = st,col = st)) + 
    scale_fill_manual(values = depcol[arlst]) +
    scale_color_manual(values = depcol[arlst]) +
    geom_histogram(alpha = 0.3, position = "identity")  +
     ggtitle(sp)
    
  plot(pl)
    #plot(as.factor(tt$stno),tt$numeriekewaarde,main=sp,col=allsamps$col)
}
dev.off()


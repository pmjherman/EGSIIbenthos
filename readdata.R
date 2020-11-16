require(tidyverse)
require(reshape2)

downloadDir<- "./data/raw_data"
dataDir    <- "./data/derived_data"

#read in downloaded data file
dats<-read_csv(file.path(downloadDir,"NV_B_All_OD_Geo_view.csv"))

# remove non-informative columns that have only one value over all rows
nuni<-vector(length=ncol(dats))
for(i in 1:ncol(dats)) nuni[i]<-nrow(unique(dats[,i]))
dats <- dats %>% select(-which(nuni==1))

#make separate grootheid_code for dg and adg
dats <- dats %>% 
  mutate(grootheid_code=ifelse(grootheid_code=="MASSA",
                        paste(grootheid_code,hoedanigheid_code),
                        grootheid_code)) %>%
  select(-hoedanigheid_code)

# determine for each sample the surface sampled and add it to the sample characteristics
surf<- dats %>% 
  filter(grootheid_code=="BEMSROPVK") %>% 
  group_by(ident_m) %>%
  summarise(surf=mean(numeriekewaarde))

# merge surf with rest of file and remove records with this information
dats <- dats %>% 
  left_join(surf, by="ident_m") %>%
  filter (grootheid_code != "BEMSROPVK")

# Make a list of all species and submit them to Worms. Why isn't that done already?
specsreported<- dats %>% group_by (biotaxon_naam)  %>% tally() %>% select(biotaxon_naam)
write_csv(specsreported,file.path(dataDir,"specsreported.csv"))
specsacc<-read_csv(file.path(dataDir,"specsreported_matched.csv"))
dats<- dats %>% left_join(specsacc,by="biotaxon_naam")

# derive WMR's sample code from the remarks column

dats <- dats %>%
  mutate(sttcd=substr(meetwaardeopmerking,1,9))

# now dats is ready.
#neglect length measurements for now
datsAB <- dats %>% filter (grootheid_code!="LENGTE")

# and cast abundance and biomass, after calculating per m2

dats2<-dcast(datsAB,ident_m+sttcd+meetpunt_identificatie+bemonsteringsmethode_code+
                  bemonsteringsapparaat_code+begindiepte_m+referentievlak_code+monsterophaaldatum+
                  monsterophaaltijd+geometriepunt_x+geometriepunt_y+wmr_gisid+
                  wmr_giscode+wmr_year_m+wmr_projnr_m+waardebepalingsmethode_code+
                  begindatum+begintijd+resultaatdatum+resultaattijd+
                  wmr_year_t+wmr_point+surf+AphiaID_accepted+ScientificName_accepted+
                  Kingdom+Phylum+Class+Order+Family+Genus+Subgenus+Species ~ 
               grootheid_code,value.var="numeriekewaarde",fun.aggregate = function(x)sum(x,na.rm=T))

dats3<- dats2 %>%
  mutate(abund=AANTAL/surf,biomdg=`MASSA dg`/surf,biomadg=`MASSA adg`/surf)

dats4abund<-dcast(dats3,ident_m+sttcd+meetpunt_identificatie+bemonsteringsmethode_code+
               bemonsteringsapparaat_code+begindiepte_m+referentievlak_code+monsterophaaldatum+
               monsterophaaltijd+geometriepunt_x+geometriepunt_y+wmr_gisid+
               wmr_giscode+wmr_year_m+wmr_projnr_m+
               begindatum+begintijd+resultaatdatum+
               resultaattijd+wmr_year_t+wmr_point+surf~
               ScientificName_accepted,value.var="abund",fun.aggregate = function(x)sum(x,na.rm=T))

dats4biomadg<-dcast(dats3,ident_m+sttcd+meetpunt_identificatie+bemonsteringsmethode_code+
                    bemonsteringsapparaat_code+begindiepte_m+referentievlak_code+monsterophaaldatum+
                    monsterophaaltijd+geometriepunt_x+geometriepunt_y+wmr_gisid+
                    wmr_giscode+wmr_year_m+wmr_projnr_m+
                    begindatum+begintijd+resultaatdatum+
                    resultaattijd+wmr_year_t+wmr_point+surf~
                    ScientificName_accepted,value.var="biomadg",fun.aggregate = function(x)sum(x,na.rm=T))

dats4biomdg<-dcast(dats3,ident_m+sttcd+meetpunt_identificatie+bemonsteringsmethode_code+
                      bemonsteringsapparaat_code+begindiepte_m+referentievlak_code+monsterophaaldatum+
                      monsterophaaltijd+geometriepunt_x+geometriepunt_y+wmr_gisid+
                      wmr_giscode+wmr_year_m+wmr_projnr_m+
                      begindatum+begintijd+resultaatdatum+
                      resultaattijd+wmr_year_t+wmr_point+surf~
                      ScientificName_accepted,value.var="biomdg",fun.aggregate = function(x)sum(x,na.rm=T))

# prepare the necessary data frames
# 1. the species

specsacc<-specsacc %>% select (-biotaxon_naam)
specsacc <- unique(specsacc)

# 2. the samples

samps<-dats4abund[,1:22]

# 3. the abundances

abund<-dats4abund[,-(1:22)]

# 4. the biomasses, for what it's worth

biom<-dats4biomadg[,-(1:22)]

# samps is adapted by adding sediment grain size information and completing all depths
# it is read in again here, and merged to what we already had

sampsed <- read_delim(file.path(dataDir,"SampsSed.csv"),delim=",")
sampsed <- sampsed %>%
  select(ident_m,sttcd,begindiepte_m,D10,D50,D90,perc0_2,
         perc2_63,perc63_125,perc125_250,perc250_500,
         perc500_1000,perc1000_2000)
samps <- samps %>%
  select(-begindiepte_m) %>%
  left_join(sampsed,by=c("ident_m","sttcd"))

# separate data with lengths
datsL<- dats %>% filter(grootheid_code=="LENGTE")

# save all these files
save(specsacc,samps,abund,biom,datsL,file=file.path(dataDir,"dats.Rdata"))
# and write them into csv's
write_csv(specsacc,path=file.path(dataDir,"specsacc.csv"))
write_csv(samps   ,path=file.path(dataDir,"samps.csv"))
write_csv(abund,path=file.path(dataDir,"abund.csv"))
write_csv(biom,path=file.path(dataDir,"biom.csv"))
write_csv(datsL,path=file.path(dataDir,"datsL.csv"))





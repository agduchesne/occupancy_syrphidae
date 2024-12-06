# Script used to run variant models as part of the sensitivity analysis.
# Recommended to run in a separate session from the main script, as this script uses many of the same object names.
# This script differs from the original model script in 3 ways: (1) either CellSize, Erasize, or IntervalSize is altered each time. (2) species must have >=500 records to be analyzed (to reduce computational load). (3) MCMC sampling is reduced (only 30k iterations, to reduce computational load)
# Results of each variant model are saved in a folder that has the name of the variant model (specified below)
# For this analysis, the variant models were named: "BaseModel", "50x50km", "200x200km", "6eras", "2eras", "10intervals", and "3intervals"
# Please set working directory to syrphid_occupancy folder


# Name the variant model. This will be the name of the folder in which these results are saved
NameOfVariantModel <- "BaseModel"
dir.create("SensitivityAnalysis")
dir.create(paste("SensitivityAnalysis/",NameOfVariantModel,sep=""))

# Variant model specifications
CellSize <- 100 #Set Gridcell size (in km). default 100
EraSize <- 30       #Size of eras, in years. default 30
IntervalSize <- 5   #Size of intervals (the repeated 'surveys' within an era), in years. default 5





########################
#### Load Libraries ####
########################

library(ggplot2) #plotting tools
library(tidyr) #drop_na function
library(abind) #abind function
library(dplyr) #Creating tables and other simple functions
options(dplyr.summarise.inform = FALSE) #Suppress the info messages when using dyplr::summarise()
library(scales) #Displaying percentages on ggplot

# Libraries for Spatial Analysis:
library(maps) #acquires maps
library(rnaturalearth) #maps with Canadian province borders (note this requires rnaturalearthhires; you can install it by running the two lines below)
#install.packages("devtools")
#devtools::install_github("ropensci/rnaturalearthhires")
library(sf)

# Libraries for JAGS model
library(parallel) #parallel.seeds function
library(rjags)
library(R2jags)
library(runjags)
library(coda) #For visualizing posteriors




##############
#### Maps ####
##############

#Create projection
na_laea <- st_crs("+proj=laea +lat_0=50 +lon_0=-100 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") #Use the Lambert Azimuthal Equal Area projection. Use NAD83 datum (ESPG:4269, more accurate than WGS84 in North America).

# Create map to make the grid
m <- map(regions = c("Canada","USA(?!:Hawaii)"), xlim = c(-180,-50), ylim = c(20,85), plot = F, fill = T) #Acquire maps and set coordinate limits
m <- m %>% st_as_sf(coords = c('x', 'y'), crs = 4269) %>% st_transform(crs = na_laea) #Project map
m.buffered <- st_buffer(m, dist=50 * 1000) #create 50km buffer zone (in meters) around edges

# Create grid
grid <- m.buffered %>%
  st_make_grid(cellsize = CellSize * 1000) %>% st_sf %>% st_cast("POLYGON") %>% #Set cell size (in meters) and make sure grid is treated as polygons (not lines, points, etc.)
  st_intersection(st_union(m.buffered)) %>%  #Clip grid to outline of buffered map (map is first turned into a single geometry)
  st_collection_extract("POLYGON") %>% st_cast() %>% #extract only the polygons (because sometimes extraneous geometries are created)
  mutate(grid.id = row_number()) %>% #Add grid ID's to each cell
  st_transform(crs=na_laea) #reproject grid (just in case)

# Create map with state/province borders, for making nice visual maps
m.states <- rnaturalearth::ne_states(c("united states of america", "canada"), returnclass="sf")
m.states <- m.states %>% st_as_sf(coords = c('x', 'y'), crs = 4269) %>% st_transform(crs = na_laea) %>% dplyr::select(name,geometry) %>% st_intersection(m) #Project map and make sure it has same shape as the map used to make grid
#plot(grid$geometry, border="lightblue")
#plot(m.states$geometry, add=T)




###################################
#### Load and Combine Datasets ####
###################################

#Function to return the year of identification from a string of ID label text. Used for the CNC Dataset.
findIDyear <- function(text) {
  if (grepl("'\\d{2}|19\\d{2}|20\\d{2}",text)==TRUE){
    if(grepl("19\\d{2}|20\\d{2}",text)==TRUE){
      ifelse(regmatches(text, regexpr("19\\d{2}|20\\d{2}", text))<=2024,return(regmatches(text, regexpr("19\\d{2}|20\\d{2}", text))),return("NA"))
    }else{
      vagueYear <- substr(regmatches(text, regexpr("'\\d{2}", text)),2,3)
      ifelse(vagueYear>20,return(paste("19",vagueYear,sep="")),return("Error: Unknown Century"))}
  }else{return("NA")}}

# CNC Dataset (Ottawa, Ontario). Obtained via CNC Online Database. https://www.cnc.agr.gc.ca/taxonomy/SpecSearchD15.php
CNCDataset <- read.csv("CNC_Syrphidae.csv") %>% mutate(Start.year = ifelse(Start.year > 0, Start.year, End.year)) # If start year has no value but end year does, move end year value to start year
CNCDataset <- mutate(CNCDataset,Collection="CNC", Year=Start.year, Month=Start.month, IdentificationYear= sapply(CNCDataset$identification,findIDyear)) #Rename columns and create identification year column

# UGuelph Dataset (Guelph, Ontario). Obtained via GBIF. https://doi.org/10.15468/dl.cs23h6
GuelphDataset <- read.csv("Guelph_Syrphidae.csv") %>% mutate(Collection="Guelph", Location=paste(locality,county,sep=", "), Province=stateProvince, Country=countryCode, Latitude=decimalLatitude, Longitude=decimalLongitude, Year=year, Month=month, Order=order, Family=family, Subfamily=subfamily, Genus=genus, Species=specificEpithet, Collector=recordedBy, IdentificationYear=substr(dateIdentified,1,4))

# E.H. Strickland Dataset (Edmonton, Alberta). Obtained via GBIF. https://doi.org/10.15468/dl.fbtax9
StricklandDataset <- read.csv("Strickland_Syrphidae.csv") %>% mutate(Collection="Strickland", Location=locality, Province=stateProvince, Country=countryCode, Latitude=decimalLatitude, Longitude=decimalLongitude, Year=year, Month=month, Order=order, Family=family, Subfamily="", Genus=genus, Species=sub(".* ","",species), Collector=recordedBy, IdentificationYear=substr(dateIdentified,1,4))

# UBC Beaty Biodiversity Museum Dataset (Vancouver, British Columbia). Obtained via Karen Needham.
SpencerDataset <- read.csv("Spencer_Syrphidae.csv") %>% mutate(Collection="Spencer", Location=paste(Collection.Location,Region,sep=", "), Year=Year.1, Month=Month.1, IdentificationYear=Determiner.Date)

# UAlaska Museum Dataset (Fairbanks, Alaska). Obtained via GBIF. https://doi.org/10.15468/dl.jnrxky
AlaskaDataset <- read.csv("Alaska_Syrphidae.csv") %>% mutate(Collection="Alaska", Location=locality, Province=stateProvince, Country=countryCode, Latitude=decimalLatitude, Longitude=decimalLongitude, Year=year, Month=month, Order=order, Family=family, Subfamily="", Genus=genus, Species=sub(".* ","",species), Collector=recordedBy, IdentificationYear=substr(dateIdentified,1,4))

# UTexas Biodiversity Center Dataset (Austin, Texas). Obtained via GBIF. https://doi.org/10.15468/dl.2u4ma4
TexasDataset <- read.csv("Texas_Syrphidae.csv") %>% mutate(Collection="Texas", Location=locality, Province=stateProvince, Country=countryCode, Latitude=decimalLatitude, Longitude=decimalLongitude, Year=year, Month=month, Order=order, Family=family, Subfamily="", Genus=genus, Species=sub(".* ","",species), Collector=recordedBy, IdentificationYear=substr(dateIdentified,1,4))

# Los Angeles County Natural History Museum Dataset (Los Angeles, California). Obtained via GBIF. https://doi.org/10.15468/dl.y2ksjz
LACMDataset <- read.csv("LACM_Syrphidae.csv") %>% mutate(Collection="LACM", Location=locality, Province=stateProvince, Country=countryCode, Latitude=decimalLatitude, Longitude=decimalLongitude, Year=year, Month=month, Order=order, Family=family, Subfamily="", Genus=genus, Species=sub(".* ","",species), Collector=recordedBy, IdentificationYear=substr(dateIdentified,1,4))

# Function to keep only desired columns
select_desired <- function(dataframe){
  return(dplyr::select(dataframe,Collection, Location, Province, Country, Latitude, Longitude, Year, Month, Order, Family, Subfamily, Genus, Species, Collector, IdentificationYear))}
# Combine the desired columns of all datasets
SyrphidDataset <- do.call(rbind,lapply(list(CNCDataset,GuelphDataset,StricklandDataset,SpencerDataset,AlaskaDataset,TexasDataset,LACMDataset), select_desired))
# Remove original datasets (to save memory)
rm(CNCDataset,GuelphDataset,StricklandDataset,SpencerDataset,AlaskaDataset,TexasDataset,LACMDataset)

# Make sure numbers are treated as such
SyrphidDataset$Latitude <- as.numeric(SyrphidDataset$Latitude)
SyrphidDataset$Longitude <- as.numeric(SyrphidDataset$Longitude)
SyrphidDataset$Year <- as.numeric(SyrphidDataset$Year)
SyrphidDataset$IdentificationYear <- as.numeric(SyrphidDataset$IdentificationYear)



############################
#### Data Cleanup Pt. 1 ####
############################

Begin_year <- 1900  #Set time period boundaries
Final_year <- 2020
records_threshold <- 500  #Minimum number of records/specimens a species must have to be included in analysis.
eras_threshold <- 4      #Minimum number of eras a species must be recorded in to be included in analysis


# Remove records outside continental CAN/USA:
SyrphidDataset <- subset(SyrphidDataset, Country %in% c("Canada","CAN","CA","U.S.A.","USA","U.S.","US","United States","United States of America")) %>% filter(Province != "Hawaii" & Province != "Puerto Rico")

# Remove records with undesired/no date and assign Era/Interval to each record:
message("Removing ", nrow(subset(SyrphidDataset, Year < Begin_year | Year > Final_year)), " records (",round(100*nrow(subset(SyrphidDataset, Year < Begin_year | Year > Final_year)) / nrow(SyrphidDataset), digits=2)  ,"%) with undesired/no date.") #Print message
EraLabels <- c(paste(Begin_year,"-",Begin_year+EraSize,sep="")) # Create the labels for eras, based on era size
for(i in seq(Begin_year+EraSize,Final_year-EraSize,EraSize)){
  EraLabels <- append(EraLabels, paste(i+1,"-",i+EraSize,sep=""))}
IntervalLabels <- c(paste(Begin_year,"-",Begin_year+IntervalSize,sep="")) # Create the labels for intervals, based on interval size
for(i in seq(Begin_year+IntervalSize,Final_year-IntervalSize,IntervalSize)){
  IntervalLabels <- append(IntervalLabels, paste(i+1,"-",i+IntervalSize,sep=""))}
SyrphidDataset <- SyrphidDataset %>% subset(Year >= Begin_year & Year <= Final_year) %>% # Remove records with no date or date outside of desired range
  mutate(Era = cut(Year, breaks = c(-Inf, seq(Begin_year+EraSize,Final_year,EraSize)), labels = EraLabels, right = TRUE)) %>% # Create column for era
  mutate(Interval = cut(Year, breaks = c(-Inf, seq(Begin_year+IntervalSize,Final_year,IntervalSize)), labels = IntervalLabels, right = TRUE)) # Create column for interval

# Remove impossible/missing coordinates:
# Southern most point in continental U.S. is Key West, FL. Western and Eastern boundaries are the Aleutian Islands and Newfoundland.
BadCoords <- filter(SyrphidDataset, is.na(Latitude) | is.na(Longitude) | Latitude==0 | Longitude==0 | Latitude < 24.5 | Latitude > 85 | Longitude > -52.6 | Longitude < -180) #List of records with missing or impossible coordinates. Out of curiosity.
message(nrow(SyrphidDataset)-nrow(BadCoords), " records (", round(100*(nrow(SyrphidDataset)-nrow(BadCoords))/nrow(SyrphidDataset), digits=2),"%) are geocoded. ", nrow(BadCoords), " records (", round(100*(nrow(BadCoords)/nrow(SyrphidDataset)), digits=2),"%) are not geocoded.") #Print message stating percent geocoded
rm(BadCoords)
SyrphidDataset <- SyrphidDataset %>% drop_na(Latitude, Longitude) %>% filter(Latitude > 24.5 & Latitude < 85 & Longitude < -52.6 & Longitude > -180.0) #Remove impossible/missing coordinates

# Remove records with no species ID
message("Removing ", nrow(subset(SyrphidDataset, Genus == "" | Species %in% c("","sp","sp."))), " records (", round(100*nrow(subset(SyrphidDataset, Genus == "" | Species %in% c("","sp","sp."))) / nrow(SyrphidDataset), digits=2), "%) without species ID.")
SyrphidDataset <- filter(SyrphidDataset, Genus != "" & !(Species %in% c("","sp","sp.")))

# Create index column
SyrphidDataset <- mutate(SyrphidDataset, index = 1:nrow(SyrphidDataset), .before="Collection")

# Create a column of Genus + Species:
SyrphidDataset$genusspecies <- paste(SyrphidDataset$Genus,SyrphidDataset$Species)

# Table showing number of records and number of eras for each species
Species_table.0 <- SyrphidDataset %>%
  group_by(genusspecies) %>% #group by species. While grouped, any given operations (ex. "summarise") will be performed on each grouping.
  summarise(sprecords = n(), num.eras = n_distinct(Era)) %>% #For each group, display the total records and the number of eras it was found in
  arrange(.by_group=TRUE) %>%
  ungroup() #Ungroup




######################################
#### Combining Synonymous Species ####
######################################

# Create dataframe with the synonym names you want to replace/combine (synonyms will be named under a single "Operational Name"). Each set of synonyms is a different row. This dataframe also corrects some spelling errors.
# Synonym information was taken from GBIF.org, the CNC Online Database Taxon Search (https://www.cnc.agr.gc.ca/taxonomy/Taxonomy.php?id=18567), the Guelph Collection Dataset on GBIF.org (which marked their records containing synonyms; https://doi.org/10.15468/dl.cs23h6), and Skevington, Jeffrey H. (2019). Field Guide to the Flower Flies of Northeastern North America. ISBN 9780691189406, and BugGuide.net.
Synonyms_table <- matrix(ncol=2, byrow=TRUE, dimnames = list(NULL,c("Synonyms","OperationalName")),
                   c("Anasimyia chrysostomus", "Anasimyia chrysostoma",
                     "Anasimyia bilinearis, Lejops billinearis, Lejops bilinearis", "Anasimyia (Lejops) bilinearis",
                     "Anasimyia distinctus", "Anasimyia distincta",
                     "Anasimyia perfidiosus, Lejops perfidiosus", "Anasimyia (Lejops) perfidiosus",
                     "Arctosyrphus willingii, Lejops willingii", "Arctosyrphus (Lejops) willingii",
                     "Brachypalpus femorata", "Brachypalpus femoratus",
                     "Brachypalpus margaritus, Brachypalpus oarus", "Brachypalpus oarus (margaritus)",
                     "Callicera auripila, Callicera erratica", "Callicera erratica (auripila)",
                     "Ceriana sayi, Polybiomyia sayi", "Ceriana (Polybiomyia) sayi",
                     "Ceriana townsendi, Polybiomyia townsendi", "Ceriana (Polybiomyia) townsendi",
                     "Chalcosyrphus libo, Xylota libo", "Chalcosyrphus (Xylota) libo",
                     "Chalcosyrphus nemorum, Xylota nemorum", "Chalcosyrphus (Xylota) nemorum",
                     "Chalcosyrphus pigra, Chalcosyrphus piger, Xylota pigra", "Chalcosyrphus (Xylota) piger (pigra)",
                     "Chalcosyrphus satanica, Xylota satanica", "Chalcosyrphus (Xylota) satanica",
                     "Chalcosyrphus vecors, Xylota vecors", "Chalcosyrphus (Xylota) vecors",
                     "Cheilosia nigroseta, Cheilosia borealis", "Cheilosia borealis (nigroseta)",
                     "Cheilosia browni, Cheilosia nigroapicata, Cheilosia lasiophthalmus", "Cheilosia lasiophthalmus (browni,nigroapicata)",
                     "Cheilosia caltha, Cheilosia sensua, Cheilosia comosa", "Cheilosia comosa (caltha, sensua)",
                     "Cheilosia consentiens, Cheilosia orilliaensis", "Cheilosia orilliaensis (consentiens)",
                     "Cheilosia latrans, Cheilosia tristis", "Cheilosia latrans (tristis)",
                     "Cheilosia cineralis, Cheilosia cynoprosopa", "Cheilosia cynoprosopa (cineralis)",
                     "Cheilosia nigrofasciata, Cheilosia hunteri", "Cheilosia hunteri (nigrofasciata)",
                     "Cheilosia alpinensis, Cheilosia rita", "Cheilosia rita (alpinensis)",
                     "Cheilosia varipila, Cheilosia variseta, Cheilosia yukonensis", "Cheilosia yukonensis (varipila,variseta)",
                     "Chrysogaster ontario, Chrysogaster antithea, Chrysogaster antitheus", "Chrysogaster antitheus (ontario)",
                     "Chrysosyrphus versipellis, Chrysosyrphus latus", "Chrysosyrphus latus (versipellis)",
                     "Chrysotoxum coloradense, Chrysotoxum fasciata, Chrysotoxum fasciatum", "Chrysotoxum coloradense (fasciatum)",
                     "Chrysotoxum ventricosum, Chrysotoxum integre, Chrysotoxum columbianum, Chrysotoxum minor, Chrysotoxum derivatum","Chrysotoxum derivatum (columbianum,integre,minor,ventricosum)",
                     "Chrysotoxum perplexum, Chrysotoxum plumeum", "Chrysotoxum plumeum (perplexum)",
                     "Chrysotoxum luteopilosum, Chrysotoxum pubescens", "Chrysotoxum pubescens (luteopilosum)",
                     "Copestylum satur, Volucella satur", "Copestylum (Volucella) satur",
                     "Copestylum vittatum, Volucella americana", "Copestylum (Volucella) vittatum (americana)",
                     "Criorhina luna, Criorhina bubulcus", "Criorhina bubulcus (luna)",
                     "Criorhina maritima, Criorhina mystaceae, Criorhina nigriventris", "Criorhina nigriventris (maritima,mystaceae)",
                     "Criorhina occidentalis, Sphecomyia occidentalis", "Criorhina (Sphecomyia) occidentalis",
                     "Cynorhinella canadensis, Cynorhinella bella", "Cynorhinella bella (canadensis)",
                     "Dasysyrphus venustus, Syrphus arcuatus", "Dasysyrphus (Syrphus) venustus (arcuatus)",
                     "Dioprosopa clavatus, Pseudodoros clavatus", "Dioprosopa (Pseudodoros) clavatus",
                     "Epistrophe emarginata, Epistrophella emarginata", "Epistrophella (Epistrophe) emarginata",
                     "Eristalis bastardii, Eristalis anthophorinus, Eristalis anthophorina", "Eristalis anthophorina (bastardii)",
                     "Eristalis compactus, Eristalis cryptarum", "Eristalis cryptarum (compactus)",
                     "Eristalis dimidiatus", "Eristalis dimidiata",
                     "Eristalis hirtus", "Eristalis hirta",
                     "Eristalis interrupta, Eristalis nemorum", "Eristalis nemorum (interrupta)",
                     "Eristalis obscurus", "Eristalis obscura",
                     "Eupeodes vockerothi, Eupeodes lungifer, Eupeodes luniger", "Eupeodes luniger (vockerothi)",
                     "Eurimyia stipatus", "Eurimyia stipata",
                     "Eurimyia lineata, Eurimyia lineatus, Lejops lineatus", "Eurimyia (Lejops) lineatus",
                     "Fagisyrphus cincta, Melangyna cincta, Meligramma cincta", "Fagisyrphus (Melangyna, Meligramma) cincta",
                     "Fazia micrura, Allograpta micrura", "Fazia (Allograpta) micrura",
                     "Ferdinandea dives, Ferdinandea buccata", "Ferdinandea buccata (dives)",
                     "Hadromyia aldrichi, Caliprobola aldrichi", "Hadromyia (Caliprobola) aldrichi",
                     "Hadromyia crawfordi, Caliprobola crawfordi", "Hadromyia (Caliprobola) crawfordi",
                     "Hadromyia pulcher", "Hadromyia pulchra",
                     "Hammerschmidtia ferruginea, Hammerschmidtia rufa, Brachyopa ferruginea", "Hammerschmidtia (Brachyopa) rufa (ferruginea)",
                     "Helophilus stricklandi, Helophilus bottnicus", "Helophilus bottnicus (stricklandi)",
                     "Hybobathus lineatus, Ocyptamus lineatus", "Hybobathus (Ocyptamus) lineatus",
                     "Hypocritanus fascipennis, Ocyptamus fascipennis", "Hypocritanus (Ocyptamus) fascipennis",
                     "Hypocritanus lemur, Ocyptamus lemur", "Hypocritanus (Ocyptamus) lemur",
                     "Lapposyrphus aberrantis, Lapposyrphus abberrantis, Eupeodes abberrantis, Eupeodes aberrantis", "Lapposyrphus (Eupeodes) aberrantis",
                     "Lapposyrphus lapponicus, Eupeodes lapponicus, Metasyrphus lapponicus", "Lapposyrphus (Eupeodes,Metasyrphus) lapponicus",
                     "Leucozona velutinus", "Leucozona velutina",
                     "Myathropa florea, Mallota florea", "Myathropa (Mallota) florea",
                     "Mallota palmerae, Mallota illinoensis", "Mallota illinoensis (palmerae)",
                     "Mallota posticata, Mallota separata, Imatisma posticata, Imatisma bardus", "Mallota (Imatisma) posticata (bardus,separata)",
                     "Mallota bautias, Imatisma bautias", "Mallota (Imatisma) bautias",
                     "Mallota columbiae, Mallota sackeni", "Mallota sackeni (columbiae)",
                     "Megasyrphus laxus, Eriozona laxa", "Megasyrphus (Eriozona) laxus",
                     "Melanostoma mellina, Melanostoma angustatum, Melanostoma fallax, Melanostoma melanderi, Melanostoma pallitarsis, Melanostoma mellinum", "Melanostoma mellinum (angustatum,fallax,melanderi,pallitarsis)",
                     "Melangyna triangulifera, Meligramma triangulifera", "Meligramma (Melangyna) triangulifera",
                     "Melangyna guttata, Meligramma guttata", "Meligramma (Melangyna) guttata",
                     "Microdon hutchingsi, Microdon globosus", "Microdon globosus (hutchingsi)",
                     "Microdon basicornis, Microdon ruficrus", "Microdon ruficrus (basicornis)",
                     "Neoascia conica, Neoascia geniculata", "Neoascia geniculata (conica)",
                     "Neoascia distincta, Neoascia globosa", "Neoascia globosa (distincta)",
                     "Neocnemodon calcaratus", "Neocnemodon calcarata",
                     "Neocnemodon carinata, Heringia carinata", "Neocnemodon (Heringia) carinata",
                     "Neocnemodon coxalis, Heringia coxalis", "Neocnemodon (Heringia) coxalis",
                     "Orthonevra parva, Chrysogaster parva", "Orthonevra (Chrysogaster) parva",
                     "Paragus hemorrhous", "Paragus haemorrhous",
                     "Parasyrphus groenlandicus, Parasyrphus interruptus, Parasyrphus bulbosus, Parasyrphus groenlandica", "Parasyrphus groenlandica (bulbosus,interruptus)",
                     "Parasyrphus lineola", "Parasyrphus lineolus",
                     "Parhelophilus obsoletus, Helophilus obsoletus", "Parhelophilus (Helophilus) obsoletus",
                     "Pelecinobaccha costata, Ocyptamus costatus", "Pelecinobaccha (Ocyptamus) costata",
                     "Pipiza oregona, Pipiza crassipes", "Pipiza crassipes (oregona)",
                     "Pipiza nigrotibiata, Pipiza puella", "Pipiza puella (nigrotibiata)",
                     "Platycheirus pauper, Platycheirus aeratus", "Platycheirus aeratus (pauper)",
                     "Platycheirus erraticus, Platycheirus hyperboreus", "Platycheirus hyperboreus (erraticus)",
                     "Platycheirus felix, Platycheirus immarginatus", "Platycheirus immarginatus (felix)",
                     "Platycheirus lata", "Platycheirus latus",
                     "Platycheirus holarcticus, Platycheirus naso", "Platycheirus naso (holarcticus)",
                     "Platycheirus bigelowi, Platycheirus parmatus", "Platycheirus parmatus (bigelowi)",
                     "Platycheirus rufimaculatus, Platycheirus concinnus, Platycheirus pictipes", "Platycheirus pictipes (concinnus,rufimaculatus)",
                     "Platycheirus podagrata", "Platycheirus podagratus",
                     "Platycheirus quadtarus", "Platycheirus quadratus",
                     "Polydontomyia curvipes, Lejops curvipes", "Polydontomyia (Lejops) curvipes",
                     "Pseudoscaeva diversifasciatus, Pseudoscaeva diversifasciata, Ocyptamus diversifasciatus", "Pseudoscaeva (Ocyptamus) diversifasciata",
                     "Pyritis montigena, Pyritis kincaidii", "Pyritis kincaidii (montigena)",
                     "Pyrophaena rosarum, Platycheirus rosarum", "Pyrophaena (Platycheirus) rosarum",
                     "Pyrophaena granditarsis, Pyrophaena granditarsus, Platycheirus granditarsis, Platycheirus granditarsus", "Pyrophaena (Platycheirus) granditarsis",
                     "Rhopalosyrphus carolae, Rhopalosyrphus guentherii", "Rhopalosyrphus guentherii (carolae)",
                     "Scaeva affinis, Scaeva pyrastri", "Scaeva affinis (pyrastri)",
                     "Serichlamys scutifer, Microdon scutifer", "Serichlamys (Microdon) scutifer",
                     "Somula marivirginiae, Somula mississippiensis", "Somula mississippiensis (marivirginiae)",
                     "Sphaerophoria contiqua, Sphaerophoria cylindricus, Sphaerophoria contigua", "Sphaerophoria contigua (cylindricus)",
                     "Sphaerophoria menthastri, Sphaerophoria scripta", "Sphaerophoria scripta (menthastri)",
                     "Sphaerophoria philantha, Sphaerophoria dubia, Sphaerophoria nigritarsi, Sphaerophoria robusta, Sphaerophoria philanthus", "Sphaerophoria philanthus (dubia,nigritarsi,robusta)",
                     "Sphaerophoria guttulata, Sphaerophoria pyrrhina", "Sphaerophoria pyrrhina (guttulata)",
                     "Sphecomyia vittata, Chrysotoxum vittatum", "Sphecomyia (Chrysotoxum) vittata",
                     "Sphegina notata, Sphegina flavomaculata", "Sphegina flavomaculata (notata)",
                     "Sphiximorpha willistoni, Ceriana willistoni", "Sphiximorpha (Ceriana) willistoni",
                     "Syrphus transversalis, Syrphus rectus", "Syrphus rectus (transversalis)",
                     "Syrphus bigelowi, Syrphus ribesii", "Syrphus ribesii (bigelowi)",
                     "Temnostoma daochum", "Temnostoma daochus",
                     "Temnostoma nipigonensis, Temnostoma venustum", "Temnostoma venustum (nipigonensis)",
                     "Temnostoma excentrica, Temnostoma excentricum, Temnostoma vespiforme", "Temnostoma excentrica (vespiforme)",
                     "Toxomerus marginata", "Toxomerus marginatus",
                     "Volucella bombylans-complex", "Volucella bombylans",
                     "Volucella fascialis, Volucella lateralis, Volucella rufomaculata, Volucella facialis", "Volucella facialis (lateralis,rufomaculata)",
                     "Xanthogramma flavipes, Philhelius flavipes", "Xanthogramma (Philhelius) flavipes",
                     "Xylota oregona, Xylota lovetti", "Xylota lovetti (oregona)",
                     "Xylota atlantica, Xylota naknek", "Xylota naknek (atlantica)",
                     "Xylota artemita, Xylota quadrimaculata", "Xylota quadrimaculata (artemita)",
                     "Xylota bicolor, Brachypalpoides bicolor", "Xylota (Brachypalpoides) bicolor",
                     "Xylota tuberculatus", "Xylota tuberculata"
                     )) %>% data.frame()

# For each record in the Dataset, if the species name matches any of the synonyms in the synonym table, substitute the name in the genusspecies column with the operational name.
for(rownum in 1:nrow(SyrphidDataset)){
if(any(grepl(SyrphidDataset$genusspecies[rownum], Synonyms_table$Synonyms))){
  SyrphidDataset$genusspecies[rownum] <- Synonyms_table$OperationalName[grep(SyrphidDataset$genusspecies[rownum], Synonyms_table$Synonyms)]}}





############################
#### Data Cleanup Pt. 2 ####
############################

# Recreate the Species Table showing number of records for each species, now that synonyms have been collapsed
# Species_table.0 shows all species in the Dataset before any are excluded
Species_table.0 <- SyrphidDataset %>% group_by(genusspecies) %>% summarise(sprecords = n(), num.eras = n_distinct(Era)) %>% arrange(.by_group=TRUE) %>% ungroup() #Ungroup

# Remove the species with less than the minimum desired number of records
message("Removing ", length(which(Species_table.0$sprecords < records_threshold)), " species (", round(100*length(which(Species_table.0$sprecords < records_threshold)) / nrow(Species_table.0), digits=2), "%) with less than ", records_threshold, " records." )
SyrphidDataset <- filter(SyrphidDataset, genusspecies %in% Species_table.0$genusspecies[which(Species_table.0$sprecords >= records_threshold)])

# Remove the species with less than the minimum desired number of eras
message("Removing ", length(which(Species_table.0$sprecords >= records_threshold & Species_table.0$num.eras < eras_threshold)), " species (", round(100*length(which(Species_table.0$sprecords >= records_threshold & Species_table.0$num.eras < eras_threshold)) / length(which(Species_table.0$sprecords >= records_threshold)), digits=2), "%) recorded in fewer than ", eras_threshold ," eras.")
SyrphidDataset <- filter(SyrphidDataset, genusspecies %in% Species_table.0$genusspecies[which(Species_table.0$num.eras >= eras_threshold)])

# Remove Volucella bombylans, which is a species of debatable validity in North America.
# Remove Criorhina willistoni, which is likely an unpublished species from Kevin Moran's work.
SyrphidDataset <- filter(SyrphidDataset, !(genusspecies %in% c("Volucella bombylans","Criorhina willistoni")))

# Remove Chrysotoxum specimens from outside the CNC (revision in progress)
SyrphidDataset <- filter(SyrphidDataset, !(Genus=="Chrysotoxum" & Collection!="CNC"))

# Remove Dasysyrphus specimens from outside the CNC that were identified before 2013 (the revision by Locke & Skevington)
SyrphidDataset <- filter(SyrphidDataset,!(Genus=="Dasysyrphus" & Collection!="CNC" & (IdentificationYear<2013|(IdentificationYear=="NA"&Year<2013))))

# Filter the Species Table to only the species being analyzed
Species_table <- SyrphidDataset %>% group_by(genusspecies) %>% summarise(sprecords = n(), num.eras = n_distinct(Era)) %>% arrange(.by_group=TRUE) %>% ungroup()

# Histogram of records per year
#records_hist <- ggplot(data=SyrphidDataset, aes(x=Year)) + geom_histogram(aes(y=..density..), binwidth=1, closed="left", color="black", alpha=0.3) + geom_density() + geom_rug() + scale_x_continuous(breaks = seq(1900,2020,10), labels = seq(1900,2020,10)) + theme(axis.text.x=element_text(angle=30))
#records_hist

# Bar graph showing number of records of each species, grouped by number of eras a species was found in
#Spgraph <- arrange(Species_table.0, num.eras,sprecords) #All species
#Spgraph <- arrange(filter(Species_table.0,sprecords<100), num.eras,sprecords) #Species with less than 100 records
#barplot(Spgraph$sprecords)
#abline(a=30,b=0,col="red")






############################################
#### Assigning coordinates to gridcells ####
############################################

# Turn data coordinates into points on map
SyrphidDataset <- SyrphidDataset %>% st_as_sf(coords = c('Longitude', 'Latitude'), remove=FALSE, crs = 4269) %>%
  st_transform(crs = na_laea)

# Remove records from SyrphidDataset that are outside the gridcells (as these coordinates are either erroneous or on faraway minuscule islands):
outside_gridcells <- which(lengths(st_intersects(SyrphidDataset$geometry,grid$geometry)) == 0) #Find which records lie outside the gridcells
#SyrphidDataset$Location[outside_gridcells] #See where these records are from
# Remove these records from the dataset:
if (length(outside_gridcells) != 0) {
  SyrphidDataset <- SyrphidDataset[-c(outside_gridcells),]}
message("Removed ",length(outside_gridcells)," records outside of the map.") #Print message

# Plot map. Save as .png file in working directory
#png(filename="RecordsMap.png",width=2500,height=2500, res=200)
#plot(grid$geometry, border="lightblue"); plot(m.states$geometry, add=T); 
#plot(SyrphidDataset$geometry[which(SyrphidDataset$Collection=="CNC")], add=T, cex=0.4, pch=20, col=rgb(red=1, green=0, blue=0, alpha=0.3)); #CNC is red
#plot(SyrphidDataset$geometry[which(SyrphidDataset$Collection=="Guelph")], add=T, cex=0.4, pch=20, col=rgb(red=0, green=1, blue=0, alpha=0.3)); #Guelph is green
#plot(SyrphidDataset$geometry[which(SyrphidDataset$Collection=="Strickland")], add=T, cex=0.4, pch=20, col=rgb(red=1, green=1, blue=0, alpha=0.3)); #Strickland is yellow
#plot(SyrphidDataset$geometry[which(SyrphidDataset$Collection=="Spencer")], add=T, cex=0.4, pch=20, col=rgb(red=0, green=0, blue=1, alpha=0.3)); #Spencer is blue
#plot(SyrphidDataset$geometry[which(SyrphidDataset$Collection=="Alaska")], add=T, cex=0.4, pch=20, col=rgb(red=1, green=0, blue=1, alpha=0.3)); #Alaska is pink
#plot(SyrphidDataset$geometry[which(SyrphidDataset$Collection=="Texas")], add=T, cex=0.4, pch=20, col=rgb(red=0, green=1, blue=1, alpha=0.3)); #Texas is cyan
#plot(SyrphidDataset$geometry[which(SyrphidDataset$Collection=="LACM")], add=T, cex=0.4, pch=20, col=rgb(red=0.5, green=0.3, blue=0.6, alpha=0.3)); #LACM is purple
#title(paste("Syrphidae (n = ",nrow(SyrphidDataset),")"),cex.main=3)
#dev.off()

# Assign gridcell id numbers to each point:
SyrphidDataset <- st_intersection(SyrphidDataset, grid) 





##################################
#### Preparing Occupancy Data ####
##################################

# Create list of species to be analyzed
splist <- sort(unique(SyrphidDataset$genusspecies))
# Create array denoting the intervals within each era
EraArray <- t(array(data=IntervalLabels,dim=c(EraSize/IntervalSize,(Final_year-Begin_year)/EraSize),dimnames=list(1:(EraSize/IntervalSize), EraLabels)))

# Function to get the gridcells within each species' assumed range (range is assumed to be all the cells within a convex hull of the cells where the species has ever been found)
get.range <- function(sp){ 
 
  # Find all gridcells within a convex hull of the cells where the species has ever been found
  IDs.observed <- unique(filter(SyrphidDataset,genusspecies==splist[sp])$grid.id) #Obtain the ID's of the gridcells where a species was found.
  geometry.observed <- grid[IDs.observed,] #Obtain the geometry of the cells where a species was found
  centroids.observed <- st_multipoint(st_coordinates(st_centroid(geometry.observed$geometry))) #Create multipoint object that contains all the centroid points of these cells
  convex.hull <- st_convex_hull(centroids.observed) #Create a convex hull polygon around those points
  assumed.range <- grid[st_intersects(grid, convex.hull, sparse = FALSE),]  #Find all gridcells within the convex hull
  
  # In case you wish to plot the range:
  #plot(m.states$geometry)
  #plot(convex.hull,border=rgb(1,0,0, alpha=0.6),lwd=2,add=T)
  #plot(geometry.observed$geometry, col=rgb(1,0,0, alpha=0.7), border=rgb(0,0,0, alpha=0.05), add=T)
  #plot(assumed.range$geometry,col=rgb(1,0,0, alpha=0.1),border=rgb(0,0,0, alpha=0.05), add=T)
  #mtext(paste(splist[sp],", n = ",Species_table$sprecords[sp],sep=""), cex=2)
  
  # Return a logical vector of all sites denoting whether they are within the assumed range
  logical.range <- ifelse(1:nrow(grid) %in% assumed.range$grid.id, TRUE, FALSE)
  return(logical.range)
}


# sp.range is a matrix of [species x site] with logical data (using function above). Shows, for each species, which sites are within their assumed range.
sp.range <- do.call(rbind, lapply(1:length(splist),get.range))
colnames(sp.range) <- paste("grid_",1:nrow(grid),sep='')

# vis.arr is an empty array of each "survey" or "visit" [site x era x interval]. Each interval in a given era is assumed to be a different "visit" to a site. In this analysis, we confine the model to only the relevant visits for each species (more on that below).
vis.arr <-array(data=NA, dim=c(nrow(grid), n_distinct(SyrphidDataset$Era), EraSize/IntervalSize), dimnames=list(paste("grid_",1:nrow(grid),sep=''),sort(as.vector(unique(SyrphidDataset$Era))),paste("interval_",1:(EraSize/IntervalSize),sep='')))

# X is [species x site x era x interval]. The observation data. 0 for undetected, 1 for detected.
get.obs <- function(sp,era,visit){ #function to get the observation data of one species during a specific interval.
  interval <- EraArray[era,visit] #obtain the interval label of this era and interval number (ex. interval 1 of era 1961-1990 is 1961-1965)
  ifelse(1:nrow(grid) %in% filter(SyrphidDataset,genusspecies==sp,Interval==interval)$grid.id, 1, 0)} 
# Create and populate the X array.
for (j in 1:(EraSize/IntervalSize)){ 
  for (i in 1:length(EraLabels)){
    Era <- EraLabels[i]
    Xij <- t(sapply(splist,era=Era,visit=j,get.obs)) #Obtain the species x site observations from each era i, given this interval number j (e.g., interval 1 from each era)
    if(i == 1){
      Xj <- Xij 
    }else{
      Xj <- abind(Xj,Xij, along=3)} #Bind together the sp x site observations from each era i, given this interval number j
  }
  if(j == 1){
    X <- Xj
  }else{
    X <- abind(X,Xj, along=4)} #Bind together the sp x site x era observations from each interval number
}
rm(Xij,Xj) #remove objects no longer needed
dimnames(X) <- list(sp=splist, site=paste("grid_",1:nrow(grid),sep=''), era=sort(as.vector(unique(SyrphidDataset$Era))), visit=paste("interval_",1:(EraSize/IntervalSize),sep=''))

# Combine data objects into a single list object
OccData <- list(sp.range=sp.range,vis.arr=vis.arr,X=X, nsp=length(splist), nsite=nrow(grid), nera=120/EraSize, nvisit=EraSize/IntervalSize)

# Keep only sites that yielded a detection of at least one species. Note gridcell ID numbers will not match site index numbers after removing gridcells that had no species (e.g. the 100th site in the matrix may actually be grid_241).
site.keep <- which(apply(OccData$X, 'site', sum)>0)
OccData$X        <- OccData$X[,site.keep,,,drop=FALSE] #Do not drop dimensions (ex. 4D object must stay 4D)
OccData$sp.range <- OccData$sp.range[,site.keep,drop=FALSE]
OccData$vis.arr  <- OccData$vis.arr[site.keep,,,drop=FALSE]
OccData$nsite    <- length(site.keep)

# Now we will filter data to only the relevant visits (visit = combination of site, era, interval) 
# (this prevents unnecessary iterating through all irrelevant sites and visits, improving model efficiency). 
# A relevant visit is one that occurs in a site within the species' range 
# and in an interval when at least 1 syrphid of ANY species was detected 
# (meaning we know someone was at that site at that time sampling for syrphids, and 
# had the potential to detect the species in question).

# Function to determine indexes of which visits are relevant for a given species
get.indices <- function(sp) {
  vis.arr <- OccData$vis.arr #Begin with a copy of vis.arr
  vis.arr[TRUE] <- 1 #Set all visits to 1
  nsp.detected <- apply(OccData$X, 2:4, sum) #object with same dimensions as vis.arr, recording the number of species detected during each visit (visit = a certain site, era and interval... aka dimensions 2:4 of X)
  vis.arr[nsp.detected==0] <- 0 #Set a visit to 0 ("irrelevant") for this species if no species were detected during that visit/interval
  vis.arr[!OccData$sp.range[sp,],,] <- 0 #Set a visit to 0 if the site is not in species range (i.e. sp.range==FALSE)
  
  tmp <- which(vis.arr==1, arr.ind=TRUE) #Find which cells in vis.arr are 1 (i.e., "relevant") for this species. arr.ind=TRUE means array indices are returned, because vis.arr is an array. Each row will direct you to a 1 within vis.arr[site,era,interval]. dim1 is site, dim2 is era and dim3 is interval.

  cbind(rep(sp,nrow(tmp)),tmp) #return the array indices found above, but add a column denoting the species number that these relevant visits are associated with.
}

master.index <- do.call(rbind, lapply(1:OccData$nsp, get.indices)) #create matrix showing indexes of relevant visits for each species.
colnames(master.index) <- c('sp','site','era','visit') #rename columns

# Create the dataset, containing observations for relevant visits only
relevant.dat <- list(X=OccData$X[master.index],                  #X becomes the observations (0 or 1) at relevant visits only, going through each species in series (one after the other). Sorted by species, then visit, then era, then site
                  master.index=master.index,                     #the master index, allowing the model to know the sp, site, era and interval of each observation in X
                  nsp=length(unique(master.index[,'sp'])),       #number of species
                  nsite=length(unique(master.index[,'site'])),   #number of relevant sites (across all spp.)
                  nera=length(unique(master.index[,'era'])),     #number of eras
                  nvisit=length(unique(master.index[,'visit'])), #number of intervals in an era (may be fewer than expected if some intervals are never relevant)
                  nind=nrow(master.index))                       #number of relevant visits/observations total across all spp.




###########################
#### Running the Model ####
###########################

# Specify the parameters to be monitored
get.params <- function() {
  c('mu.psi',         #Baseline occupancy of syrphids across species and eras (includes default effect of era 1)
    'mu.p',           #Baseline/mean detection probability of syrphids (includes default effect of era 1)
    'psi.sp',         #Random Effect of a certain species on occupancy
    'sigma.psi.sp',   #Hyperprior sd for effects of species on occupancy. Represents variation amongst the effects of different species on occupancy.
    'p.sp',           #Random Effect of a certain species on detection
    'sigma.p.sp',     #Hyperprior sd for effects of species on detection. Represents variation amongst the effects of different species on detection.
    'p.site',         #Random effect of site*era on detection. Includes effect of sampling effort between sites and eras.
    'sigma.p.site',   #Hyperprior sd for effects of site*era on detection. Represents variation amongst the effects of different site*eras on detection.
    'psi.era',        #Species-specific effect of era on occupancy
    'mu.psi.era',     #Hyperprior mean for effects of era on occupancy. Represents the mean value amongst the different species' effects on occupancy in an era. In other words, the mean across-species effect of an era on occupancy.
    'sigma.psi.era',  #Hyperprior sd for effects of era on occupancy. Represents variation amongst the different species' effects on occupancy in an era. In other words, the variation in the across-species effect of an era on occupancy.
    'p.era') }        #Effect of era on detection.

# Function to run the model
run.model <- function(relevant.dat, n.iter, n.burnin, n.adapt, n.thin) { 

  model.txt <- sprintf('SyrphidModel.txt') #Bring in the BUGS-syntax model from the external source file

  # Create data object that will be passed to the model
  jags.data <- list(X=relevant.dat$X,                            #Observations during all relevant visits
                      era=relevant.dat$master.index[,'era'],     #Eras associated with each observation
                      site=relevant.dat$master.index[,'site'],   #Sites associated with each observation
                      sp=relevant.dat$master.index[,'sp'],       #Species number associated with each observation
                      nsp=relevant.dat$nsp,                      #Number of sp
                      nsite=relevant.dat$nsite,                  #Number of sites
                      nera=relevant.dat$nera,                    #Number of eras
                      nind=relevant.dat$nind)                    #Number of total observations

  # Prepare initial values for MCMC
  Zst <- array(1,dim=c(relevant.dat$nsp, relevant.dat$nsite, relevant.dat$nera)) #Creates array of [sp x site x era] filled with 1's
  make.inits <- function() {  #Function that helps initialize the random number generator of each MCMC chain when running parallel chains (on separate cores).
    RNG <- parallel.seeds("base::BaseRNG", 1)
    c(list(Z=Zst), RNG[[1]])}
  inits1 <- make.inits() #Initialize random number generator for each parallel chain you are prepared to run
  inits2 <- make.inits()
  inits3 <- make.inits()
  inits4 <- make.inits()
    
  # Then, run the model
  jags.out <- run.jags(model=model.txt,
                        monitor=get.params(),
                        data=jags.data,
                        inits=list(inits1,inits2,inits3,inits4),
                        n.chains=4, #Four parallel chains being processed at a time
                        burnin=n.burnin,
                        sample=floor(n.iter/n.thin),
                        thin=n.thin,
                        adapt=n.adapt,
                        method='rjags')
  
  # Return the data object used and the model output
  list(jags.data=jags.data, jags.out=jags.out)
}

# Run the model! (using the run.model function defined above)
model.out <- run.model(relevant.dat,
                       n.iter=30e3,   #Number of iterations of each MCMC chain
                       n.burnin=4e3, #Number of burn-in/warmup iterations
                       n.adapt=2e3,  #These iterations are used to tune the behaviour of the MCMC sampler, making it more efficient
                       n.thin=1e1)   #Thinning interval to increase independence of samples within chain. Every 10th sample is kept.









#######################
#### Model Outputs ####
#######################

# Save/load the model output on your computer
#save(model.out,file=paste("SensitivityAnalysis/",NameOfVariantModel,"/model.out.RData",sep=""))     #SAVE
load(paste("SensitivityAnalysis/",NameOfVariantModel,"/model.out.RData",sep=""))                    #LOAD

# Create the Model Summary
#ModelSummary <- summary(model.out$jags.out)
# Save/load the model summary on your computer
#save(ModelSummary,file=paste("SensitivityAnalysis/",NameOfVariantModel,"/ModelSummary.RData",sep=""))  #SAVE
load(paste("SensitivityAnalysis/",NameOfVariantModel,"/ModelSummary.RData",sep=""))                    #LOAD

# Create column in ModelSummary that lists the species name associated with each parameter.
ModelSummary <- mutate(data.frame(ModelSummary), Species = NA, .before = "Lower95")
for(i in 1:nrow(ModelSummary)){
  if(grepl("psi.sp|p.sp|psi.era", rownames(ModelSummary[i,]))){ #If the parameter is any of the following parameters, find the associated species name: psi.sp[sp], p.sp[sp], psi.era[sp,era]
    if(grepl("psi.sp|p.sp", rownames(ModelSummary[i,]))){ #if it is psi.sp or p.sp, the species number will be within square brackets
      ModelSummary$Species[i] <- splist[as.numeric(regmatches(rownames(ModelSummary[i,]), regexec("\\[(\\d+)\\]", rownames(ModelSummary[i,])))[[1]][2])]} #regexec finds the position of the parts of the string that satisfy the condition (i.e., are within square brackets). regmatches then obtains the substrings at those positions. We only want the first match (there should only be one), so we use [[1]]. We want the number only, without the square brackets, so we use [2].
    if(grepl("psi.era", rownames(ModelSummary[i,]))){ #if it is psi.era, the species number will be between a square bracket and a comma
      ModelSummary$Species[i] <- splist[as.numeric(regmatches(rownames(ModelSummary[i,]), regexec("\\[(\\d+)\\,", rownames(ModelSummary[i,])))[[1]][2])]}}} #regexec finds the position of the parts of the string that satisfy the condition (i.e., are within square brackets). regmatches then obtains the substrings at those positions. We only want the first match (there should only be one), so we use [[1]]. We want the number only, without the square brackets, so we use [2].

MCMC <- model.out$jags.out$mcmc #from here on, use MCMC as shorthand for quick/easy access to the mcmc samples in the model output

# Plotting a single specified posterior
#parameter.to.plot <- "mu.psi" #View the traceplot and posterior density of this parameter
#plot(MCMC[,parameter.to.plot], main = parameter.to.plot) #plot the mcmc output (all samples/chains) for the desired parameter. The mcmc.list object is a list of 4 chains (each an mcmc object) with the dimensions [#samples,#parameters]. main is the title




################################
#### Estimates of Occupancy ####
################################

# Function to compute occupancy (as a mcmc.list object containing each sampled chain) for a given species and era
compute_occupancy <- function(sp,era){
  #Occupancy is baseline occupancy + effect of species + effect of era
  mu.psi <- matrix(as.matrix(MCMC[,"mu.psi"]), nrow=model.out$jags.out$sample, ncol=length(MCMC)) #Obtain the necessary parameter distributions
  psi.sp <- matrix(as.matrix(MCMC[,paste("psi.sp[",sp,"]",sep="")]), nrow=model.out$jags.out$sample, ncol=length(MCMC))
  psi.era <- matrix(as.matrix(MCMC[,paste("psi.era[",sp,",",era,"]",sep="")]), nrow=model.out$jags.out$sample, ncol=length(MCMC))
  occupancy <- as.matrix(plogis(as.matrix(mu.psi+psi.sp+psi.era))) #Compute occupancy
  occupancy <- as.mcmc.list(lapply(seq_len(ncol(occupancy)), function(i) as.mcmc(occupancy[,i]))) #convert matrix back into mcmc object (each chain i is converted into an mcmc object by as.mcmc(), which are combined into an mcmc.list object)
  return(occupancy)}

# Function to compute across-species mean occupancy in each era
compute_cross.sp_occupancy <- function(era){
  mu.psi <- matrix(as.matrix(MCMC[,"mu.psi"]), nrow=model.out$jags.out$sample, ncol=length(MCMC)) #Obtain mu.psi distribution
  if (era==1){ #If era is 1, across species occupancy is simply the expit of mu.psi. The effect of era 1 (the "default" era) is already contained in mu.psi.
    mu.psi.era <- 0
  }else{ #If era is not 1, must obtain the mu.psi.era parameter (which represents the mean effect of this era across species)
    mu.psi.era <- matrix(as.matrix(MCMC[,paste("mu.psi.era[",era,"]",sep="")]), nrow=model.out$jags.out$sample, ncol=length(MCMC))}
  occupancy <- as.matrix(plogis(as.matrix(mu.psi+mu.psi.era))) #Compute occupancy
  occupancy <- as.mcmc.list(lapply(seq_len(ncol(occupancy)), function(i) as.mcmc(occupancy[,i]))) #convert matrix back into mcmc object
  return(occupancy)}

# Create a list object containing all occupancy posteriors
Occupancy_list <- list()
meanoccupancies <- list()
for (era in 1:model.out$jags.data$nera){ #Compute cross-species occupancies for each era
  meanoccupancies[[paste(era)]] <- compute_cross.sp_occupancy(era)}
Occupancy_list[["cross.sp.mean"]] <- meanoccupancies #add the cross-species occupancy posteriors to the master list
rm(meanoccupancies)
for (sp in 1:model.out$jags.data$nsp) { #Compute species-specific occupancies
  spoccupancies <- list()
  for (era in 1:model.out$jags.data$nera) {
    spoccupancies[[paste(era)]] <- compute_occupancy(sp,era)}
  Occupancy_list[[splist[sp]]] <- spoccupancies #add this species' occupancy posteriors to the master list
  rm(spoccupancies)}

# Create a dataframe of summary statistics (mean, sd, 95% intervals) for all occupancy estimates.
get_occ_stats <- function(row){ #Function to get the summary stats for each era's occupancy for a given species
  dataframe <- matrix(nrow=model.out$jags.data$nera,ncol=6) #create matrix to hold data.
  for (era in 1:model.out$jags.data$nera){ #Get the species name, era name, lower95, mean, upper95, and SD.
    dataframe[era,] <- c(names(Occupancy_list[row]), names(Occupancy_list[[row]][era]), HPDinterval(as.mcmc(unlist(Occupancy_list[[row]][[era]])), prob=0.95)[[1]] , summary(Occupancy_list[[row]][[era]])$statistics[[1]], HPDinterval(as.mcmc(unlist(Occupancy_list[[row]][[era]])), prob=0.95)[[2]] , summary(Occupancy_list[[row]][[era]])$statistics[[2]])}
  return(dataframe)}
OccupancySummary <- as.data.frame(do.call(rbind,lapply(1:length(Occupancy_list), get_occ_stats))) #Create the dataframe
colnames(OccupancySummary) <- c("Species","Era","Lower95","Mean","Upper95","SD")

# Make sure columns are treated as numeric values
OccupancySummary$Lower95 <- as.numeric(OccupancySummary$Lower95)
OccupancySummary$Mean <- as.numeric(OccupancySummary$Mean)
OccupancySummary$Upper95 <- as.numeric(OccupancySummary$Upper95)
OccupancySummary$SD <- as.numeric(OccupancySummary$SD)

# Save Occupancy Summary
#save(OccupancySummary,file=paste("SensitivityAnalysis/",NameOfVariantModel,"/OccupancySummary.RData",sep=""))
# Load Occupancy Summary
load(paste("SensitivityAnalysis/",NameOfVariantModel,"/OccupancySummary.RData",sep=""))





#############################################
#### Occupancy Changes (delta occupancy) ####
#############################################

# Function to compute the percent occupancy change between two eras
compute_delta.occ <- function(row,era1,era2) { 
  Occ1 <- matrix(as.matrix(Occupancy_list[[row]][[era1]]), nrow=model.out$jags.out$sample, ncol=length(MCMC)) #Occupancy estimate of era1
  Occ2 <- matrix(as.matrix(Occupancy_list[[row]][[era2]]), nrow=model.out$jags.out$sample, ncol=length(MCMC)) #Occupancy estimate of era2
  delta.occ <- as.matrix((Occ1-Occ2)/Occ2) #subtract the distributions and divide by earlier era
  delta.occ <- as.mcmc.list(lapply(seq_len(ncol(delta.occ)), function(i) as.mcmc(delta.occ[,i]))) 
  return(delta.occ)}

# Create a list object containing all occupancy changes between eras (delta occupancy)
delta.occ_list <- list()
for (row in 1:length(Occupancy_list)) {
  occ.changes <- list()
  for (era in model.out$jags.data$nera:2) {
    for (nextera in (era-1):1) {
      occ.changes[[paste(era,"-",nextera, sep="")]] <- compute_delta.occ(row,era,nextera)}} #Occupancy change is computed between the era in question and every era below it. Example with 4 eras: 4-3,4-2,4-1,3-2,3-1,2-1
  delta.occ_list[[names(Occupancy_list[row])]] <- occ.changes #add this species' contrasts to the master list
  rm(occ.changes,row,era,nextera)}

# Create a dataframe of summary statistics (mean, sd, 95% intervals) for all occupancy change estimates.
get_delta_occ_stats <- function(row){ #Function to get the summary stats for all occupancy changes for a given species
  dataframe <- matrix(nrow=model.out$jags.data$nera*(model.out$jags.data$nera-1)/2,ncol=6) #create matrix to hold data. The number of pairwise contrasts is always nera*(nera-1)/2
  for (contrast in 1:(model.out$jags.data$nera*(model.out$jags.data$nera-1)/2)){ #For each contrast, get the species name, contrast name, lower95, mean, upper95, and SD.
    dataframe[contrast,] <- c(names(Occupancy_list[row]), names(delta.occ_list[[row]][contrast]), HPDinterval(as.mcmc(unlist(delta.occ_list[[row]][[contrast]])), prob=0.95)[[1]] , summary(delta.occ_list[[row]][[contrast]])$statistics[[1]], HPDinterval(as.mcmc(unlist(delta.occ_list[[row]][[contrast]])), prob=0.95)[[2]] , summary(delta.occ_list[[row]][[contrast]])$statistics[[2]])}
  return(dataframe)}
delta.occSummary <- as.data.frame(do.call(rbind,lapply(1:length(Occupancy_list), get_delta_occ_stats))) #Create the dataframe
colnames(delta.occSummary) <- c("Species","EraContrast","Lower95","Mean","Upper95","SD")

# Make sure columns are treated as numeric values
delta.occSummary$Lower95 <- as.numeric(delta.occSummary$Lower95)
delta.occSummary$Mean <- as.numeric(delta.occSummary$Mean)
delta.occSummary$Upper95 <- as.numeric(delta.occSummary$Upper95)
delta.occSummary$SD <- as.numeric(delta.occSummary$SD)

# Save delta Occupancy Summary
#save(delta.occSummary,file=paste("SensitivityAnalysis/",NameOfVariantModel,"/delta.occSummary.RData",sep=""))
# Load delta Occupancy Summary
load(paste("SensitivityAnalysis/",NameOfVariantModel,"/delta.occSummary.RData",sep=""))







###############################################
#### Lists of Declining/Increasing Species ####
###############################################


# List of negative delta occupancy posteriors whose 95% credible intervals do not include zero (i.e., declines)
Negative_delta.occ <- filter(delta.occSummary, Lower95 < 0 & Upper95 <= 0)
Negative_sp <- unique(filter(Negative_delta.occ, Species!="cross.sp.mean")$Species) # List of species (do not include the cross.sp.mean)

# List of positive delta occupancy posteriors whose 95% credible intervals do not include zero (i.e., increases)
Positive_delta.occ <- filter(delta.occSummary, Lower95 >= 0 & Upper95 > 0)
Positive_sp <- unique(filter(Positive_delta.occ, Species!="cross.sp.mean")$Species) # List of species

# List of species in recent decline (declined between the 2 most recent eras)
RecentDecline_delta.occ <- filter(delta.occSummary, Lower95 < 0 & Upper95 <= 0 & EraContrast==paste(relevant.dat$nera,"-",relevant.dat$nera-1,sep=""))
RecentDecline_sp <- unique(filter(RecentDecline_delta.occ, Species!="cross.sp.mean")$Species)

# List of species in recent growth (increased between the 2 most recent eras)
RecentGrowth_delta.occ <- filter(delta.occSummary, Lower95 >= 0 & Upper95 > 0 & EraContrast==paste(relevant.dat$nera,"-",relevant.dat$nera-1,sep=""))
RecentGrowth_sp <- unique(filter(RecentGrowth_delta.occ, Species!="cross.sp.mean")$Species)

# List of species with exceptional recent decline -- a lower occupancy now than all historic time periods (a subset of RecentDecline)
ExceptionalRecentDeclines.run <- function(){
  ExceptionalRecentDecline_sp <<- as.character() #Use global assignment operator to make sure this object isn't temporary.
  for (sp in RecentDecline_sp){
    for (historic.era in (relevant.dat$nera-2):1){ #No need to check if it declined between 2 most recent periods, because using RecentDecline_sp ensures that.
      if(!any(Negative_delta.occ$Species==sp & Negative_delta.occ$EraContrast==paste(relevant.dat$nera,"-",historic.era,sep=""))){ #If a contrast between the most recent era and any other historical era is NOT in Negative_delta.occ, skip/ignore this species.
        break}
      if(historic.era==1){ #If all historic eras have been checked and no breaks occurred (i.e., they are all in Negative_delta.occ), add the species to the list
        ExceptionalRecentDecline_sp <<- append(ExceptionalRecentDecline_sp, sp)}}}
  ExceptionalRecentDecline_delta.occ <<- filter(delta.occSummary, Species %in% ExceptionalRecentDecline_sp)}
ExceptionalRecentDeclines.run()

# List of species with exceptional recent growth -- a higher occupancy now than all historic time periods (a subset of RecentGrowth)
ExceptionalRecentGrowth.run <- function(){
  ExceptionalRecentGrowth_sp <<- as.character() #Use global assignment operator to make sure this object isn't temporary.
  for (sp in RecentGrowth_sp){
    for (historic.era in (relevant.dat$nera-2):1){ #No need to check if it increased between 2 most recent periods, because using RecentGrowth_sp ensures that.
      if(!any(Positive_delta.occ$Species==sp & Positive_delta.occ$EraContrast==paste(relevant.dat$nera,"-",historic.era,sep=""))){ #If a contrast between the most recent era and any other historical era is NOT in Positive_delta.occ, skip/ignore this species.
        break}
      if(historic.era==1){ #If all historic eras have been checked and no breaks occurred (i.e., they are all in Positive_delta.occ), add the species to the list
        ExceptionalRecentGrowth_sp <<- append(ExceptionalRecentGrowth_sp, sp)}}}
  ExceptionalRecentGrowth_delta.occ <<- filter(delta.occSummary, Species %in% ExceptionalRecentGrowth_sp)}
ExceptionalRecentGrowth.run()

# List of species with historic decline (declined between the 2 terminal eras)
HistoricDecline_delta.occ <- filter(delta.occSummary, Lower95 < 0 & Upper95 <= 0 & EraContrast==paste(relevant.dat$nera,"-1",sep=""))
HistoricDecline_sp <- unique(filter(HistoricDecline_delta.occ, Species!="cross.sp.mean")$Species)

# List of species with historic growth (increased between the 2 terminal areas)
HistoricGrowth_delta.occ <- filter(delta.occSummary, Lower95 >= 0 & Upper95 > 0 & EraContrast==paste(relevant.dat$nera,"-1",sep=""))
HistoricGrowth_sp <- unique(filter(HistoricGrowth_delta.occ, Species!="cross.sp.mean")$Species)

# List of species with no recent or historic change in occupancy
NoDiscernibleTrend_sp <- splist[!(splist %in% c(RecentDecline_sp,HistoricDecline_sp,RecentGrowth_sp,HistoricGrowth_sp))]

# Create dataframe with all of these lists
max_length <- max(length(RecentDecline_sp),length(ExceptionalRecentDecline_sp),length(HistoricDecline_sp),length(RecentGrowth_sp),length(ExceptionalRecentGrowth_sp),length(HistoricGrowth_sp),length(NoDiscernibleTrend_sp))
Results.table <- data.frame(RecentDecline = c(paste("Declined between", EraLabels[relevant.dat$nera-1] ,"and", EraLabels[relevant.dat$nera]),RecentDecline_sp,rep("",max_length-length(RecentDecline_sp))),
                            ExceptionalRecentDecline = c(paste("Occupancy in",EraLabels[relevant.dat$nera],"is lower than all historical time periods. (A subset of Recent Decline and Historical Decline)"),ExceptionalRecentDecline_sp,rep("",max_length-length(ExceptionalRecentDecline_sp))),
                            HistoricDecline = c(paste("Declined between", EraLabels[1] ,"and", EraLabels[relevant.dat$nera]),HistoricDecline_sp,rep("",max_length-length(HistoricDecline_sp))),
                            RecentGrowth = c(paste("Increased between", EraLabels[relevant.dat$nera-1] ,"and", EraLabels[relevant.dat$nera]),RecentGrowth_sp,rep("",max_length-length(RecentGrowth_sp))),
                            ExceptionalRecentGrowth = c(paste("Occupancy in",EraLabels[relevant.dat$nera],"is greater than all historical time periods (A subset of Recent Growth and Historical Growth)"),ExceptionalRecentGrowth_sp,rep("",max_length-length(ExceptionalRecentGrowth_sp))),
                            HistoricGrowth = c(paste("Increased between", EraLabels[1] ,"and", EraLabels[relevant.dat$nera]),HistoricGrowth_sp,rep("",max_length-length(HistoricGrowth_sp))),
                            NoDiscernibleTrend = c("No significant change over recent or historical time frames", NoDiscernibleTrend_sp,rep("",max_length-length(NoDiscernibleTrend_sp))))

# Save or load results table
#write.csv(Results.table,paste("SensitivityAnalysis/",NameOfVariantModel,"/Results.table.csv",sep=""),row.names=F) #Save results
Results.table <- read.csv(paste("SensitivityAnalysis/",NameOfVariantModel,"/Results.table.csv",sep="")) #Load results














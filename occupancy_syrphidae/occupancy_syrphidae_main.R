# MAIN SCRIPT
# Please set working directory to occupancy_syrphidae folder


########################
#### Load Libraries ####
########################

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

# Libraries for Association Analyses
library(multcomp) #glht function
library(FSA) #Dunn's post-hoc test
library(vcd) #mosaic plotting

# Libraries for plotting
library(ggplot2) #plotting tools
library(scales) #Displaying percentages on ggplot
library(grid)
library(gridExtra) #grid.arrange()

# Misc. Libraries
library(tidyr) #drop_na function
library(abind) #abind function
library(dplyr) #Creating tables and other simple functions
options(dplyr.summarise.inform = FALSE) #Suppress unnecessary info messages when using dyplr::summarise()



##############
#### Maps ####
##############

CellSize <- 100 #Set Gridcell size (in km) 

#Create projection
na_laea <- st_crs("+proj=laea +lat_0=50 +lon_0=-100 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") #Use the Lambert Azimuthal Equal Area projection. Use NAD83 datum (ESPG:4269, more accurate than WGS84 in North America).

# Create map to make the grid
m <- map(regions = c("Canada","USA(?!:Hawaii)"), xlim = c(-180,-50), ylim = c(20,85), plot = F, fill = T) #Acquire maps and set coordinate limits
m <- m %>% st_as_sf(coords = c('x', 'y'), crs = 4269) %>% st_transform(crs = na_laea) #Project map
m.buffered <- st_buffer(m, dist=50 * 1000) #create 50km buffer zone (in meters) around edges of map

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
# Plotting map
#plot(grid$geometry, border="lightblue")
#plot(m.states$geometry, add=T)




###################################
#### Load and Combine Datasets ####
###################################

#Function to return the year of identification from a string of transcribed label text. Used for the CNC Dataset.
findIDyear <- function(text) {
  if (grepl("'\\d{2}|19\\d{2}|20\\d{2}",text)==TRUE){ #\\d specifies any digit {2} exactly 2 times. Search for any 2 digits that follow an apostrophe (e.g., '82, '11,), a 19, or a 20.
    if(grepl("19\\d{2}|20\\d{2}",text)==TRUE){
      ifelse(regmatches(text, regexpr("19\\d{2}|20\\d{2}", text))<=2024,return(regmatches(text, regexpr("19\\d{2}|20\\d{2}", text))),return("NA"))
    }else{ #if year was represented with an apostrophe, century must be inferred. Must be 20th century if the year is greater than 20 (because only data up to 2020 is included)
      vagueYear <- substr(regmatches(text, regexpr("'\\d{2}", text)),2,3) 
      ifelse(vagueYear>20,return(paste("19",vagueYear,sep="")),return("Error: Unknown Century"))}
  }else{return("NA")}}

# Load CNC Dataset (Ottawa, Ontario). Obtained via CNC Online Database. https://www.cnc.agr.gc.ca/taxonomy/SpecSearchD15.php
CNCDataset <- read.csv("Datasets/CNC_Syrphidae.csv") %>% mutate(Start.year = ifelse(Start.year > 0, Start.year, End.year)) # If start year has no value but end year does, move end year value to start year
CNCDataset <- mutate(CNCDataset,Collection="CNC", Year=Start.year, Month=Start.month, IdentificationYear= sapply(CNCDataset$identification,findIDyear)) #Rename columns and create identification year column

# Load UGuelph Dataset (Guelph, Ontario). Obtained via GBIF. https://doi.org/10.15468/dl.cs23h6
GuelphDataset <- read.csv("Datasets/Guelph_Syrphidae.csv") %>% mutate(Collection="Guelph", Location=paste(locality,county,sep=", "), Province=stateProvince, Country=countryCode, Latitude=decimalLatitude, Longitude=decimalLongitude, Year=year, Month=month, Order=order, Family=family, Subfamily=subfamily, Genus=genus, Species=specificEpithet, Collector=recordedBy, IdentificationYear=substr(dateIdentified,1,4))

# Load E.H. Strickland Entomological Museum Dataset (Edmonton, Alberta). Obtained via GBIF. https://doi.org/10.15468/dl.fbtax9
StricklandDataset <- read.csv("Datasets/Strickland_Syrphidae.csv") %>% mutate(Collection="Strickland", Location=locality, Province=stateProvince, Country=countryCode, Latitude=decimalLatitude, Longitude=decimalLongitude, Year=year, Month=month, Order=order, Family=family, Subfamily="", Genus=genus, Species=sub(".* ","",species), Collector=recordedBy, IdentificationYear=substr(dateIdentified,1,4))

# Load UBC Spencer Entomological Collection Dataset (Vancouver, British Columbia). Obtained directly from curator Karen Needham.
SpencerDataset <- read.csv("Datasets/Spencer_Syrphidae.csv") %>% mutate(Collection="Spencer", Location=paste(Collection.Location,Region,sep=", "), Year=Year.1, Month=Month.1, IdentificationYear=Determiner.Date)

# Load UAlaska Museum Dataset (Fairbanks, Alaska). Obtained via GBIF. https://doi.org/10.15468/dl.jnrxky
AlaskaDataset <- read.csv("Datasets/Alaska_Syrphidae.csv") %>% mutate(Collection="Alaska", Location=locality, Province=stateProvince, Country=countryCode, Latitude=decimalLatitude, Longitude=decimalLongitude, Year=year, Month=month, Order=order, Family=family, Subfamily="", Genus=genus, Species=sub(".* ","",species), Collector=recordedBy, IdentificationYear=substr(dateIdentified,1,4))

# Load UTexas Biodiversity Center Dataset (Austin, Texas). Obtained via GBIF. https://doi.org/10.15468/dl.2u4ma4
TexasDataset <- read.csv("Datasets/Texas_Syrphidae.csv") %>% mutate(Collection="Texas", Location=locality, Province=stateProvince, Country=countryCode, Latitude=decimalLatitude, Longitude=decimalLongitude, Year=year, Month=month, Order=order, Family=family, Subfamily="", Genus=genus, Species=sub(".* ","",species), Collector=recordedBy, IdentificationYear=substr(dateIdentified,1,4))

# Load Los Angeles County Natural History Museum Dataset (Los Angeles, California). Obtained via GBIF. https://doi.org/10.15468/dl.y2ksjz
LACMDataset <- read.csv("Datasets/LACM_Syrphidae.csv") %>% mutate(Collection="LACM", Location=locality, Province=stateProvince, Country=countryCode, Latitude=decimalLatitude, Longitude=decimalLongitude, Year=year, Month=month, Order=order, Family=family, Subfamily="", Genus=genus, Species=sub(".* ","",species), Collector=recordedBy, IdentificationYear=substr(dateIdentified,1,4))

# Function to keep only desired columns
select_desired <- function(dataframe){
  return(dplyr::select(dataframe,Collection, Location, Province, Country, Latitude, Longitude, Year, Month, Order, Family, Subfamily, Genus, Species, Collector, IdentificationYear))}
# Combine the desired columns of all datasets
SyrphidDataset <- do.call(rbind,lapply(list(CNCDataset,GuelphDataset,StricklandDataset,SpencerDataset,AlaskaDataset,TexasDataset,LACMDataset), select_desired))
# Remove original datasets (to save memory)
rm(CNCDataset,GuelphDataset,StricklandDataset,SpencerDataset,AlaskaDataset,TexasDataset,LACMDataset)

# Make sure numbers are treated as numbers
SyrphidDataset$Latitude <- as.numeric(SyrphidDataset$Latitude)
SyrphidDataset$Longitude <- as.numeric(SyrphidDataset$Longitude)
SyrphidDataset$Year <- as.numeric(SyrphidDataset$Year)
SyrphidDataset$IdentificationYear <- as.numeric(SyrphidDataset$IdentificationYear)




############################
#### Data Cleanup Pt. 1 ####
############################

Begin_year <- 1900  #Set time period boundaries
Final_year <- 2020
EraSize <- 30       #Size of eras, in years
IntervalSize <- 5   #Size of intervals (the repeated 'surveys' within an era), in years
records_threshold <- 50  #Minimum number of records/specimens a species must have to be included in analysis
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
  group_by(genusspecies) %>% #group records by species name. While grouped, any given operations (ex. "summarise") will be performed on each grouping.
  summarise(sprecords = n(), num_eras = n_distinct(Era)) %>% #For each group, display the total records and the number of eras it was found in
  arrange(.by_group=TRUE) %>%
  ungroup() #Ungroup




######################################
#### Combining Synonymous Species ####
######################################

# Create dataframe with the synonym names you want to replace/combine (synonyms will be named under a single "Operational Name"). Each set of synonyms is a different row. This dataframe also corrects some spelling errors.
# Synonym information was taken from GBIF.org, the CNC Online Database Taxon Search (https://www.cnc.agr.gc.ca/taxonomy/Taxonomy.php?id=18567), the Guelph Collection Dataset on GBIF.org (which marked their records containing synonyms; https://doi.org/10.15468/dl.cs23h6), and Skevington, Jeffrey H. (2019). Field Guide to the Flower Flies of Northeastern North America. ISBN 9780691189406.
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

Synonyms_table <- mutate(Synonyms_table, OperationalName_NoBrackets = gsub("\\s\\(.*?\\)", "", Synonyms_table$OperationalName)) #create column displaying operational name without synonyms in brackets
Synonyms_table$Synonyms_NoOperationalName <- NA
for (i in 1:nrow(Synonyms_table)){ #create column showing synonyms, but excluding the operational name
  Synonyms_table$Synonyms_NoOperationalName[i] <- gsub(paste(Synonyms_table$OperationalName_NoBrackets[i],"\\,\\s|\\,\\s",Synonyms_table$OperationalName_NoBrackets[i],sep=""), "", Synonyms_table$Synonyms[i])}

# Save synonyms table
#write.csv(Synonyms_table[,c(2,1,3,4)], "FinalModel/SynonymsTable.csv") 

# For each record in the SyrphidDataset, if the species name matches any of the synonyms in the synonym table, substitute the name in the genusspecies column with the operational name.
for(rownum in 1:nrow(SyrphidDataset)){
if(any(grepl(SyrphidDataset$genusspecies[rownum], Synonyms_table$Synonyms))){
  SyrphidDataset$genusspecies[rownum] <- Synonyms_table$OperationalName[grep(SyrphidDataset$genusspecies[rownum], Synonyms_table$Synonyms)]}}




############################################
#### Assigning coordinates to gridcells ####
############################################

# Turn data coordinates into points on map
SyrphidDataset <- SyrphidDataset %>% st_as_sf(coords = c('Longitude', 'Latitude'), remove=FALSE, crs = 4269) %>%
  st_transform(crs = na_laea)

# Remove records from SyrphidDataset that are outside the grid cells (as these coordinates are either erroneous or on far-away minuscule islands):
outside_gridcells <- which(lengths(st_intersects(SyrphidDataset$geometry,grid$geometry)) == 0) #Find which records lie outside the grid cells
#SyrphidDataset$Location[outside_gridcells] #See where these erroneous records are from
# Remove these erroneous records from the dataset:
if (length(outside_gridcells) != 0) {
  SyrphidDataset <- SyrphidDataset[-c(outside_gridcells),]}
message("Removed ",length(outside_gridcells)," records outside of the map.") #Print message

# Assign gridcell id numbers to each point:
SyrphidDataset <- st_intersection(SyrphidDataset, grid)





############################
#### Data Cleanup Pt. 2 ####
############################
# Recreate the Species Table displaying each species in the analysis, now that synonyms have been collapsed and grid cells have been assigned

# Species_table.0 shows all species in the Dataset before any are excluded
Species_table.0 <- SyrphidDataset %>% group_by(genusspecies) %>% summarise(sprecords = n(), num_eras = n_distinct(Era), num_grids = n_distinct(grid.id), distinct_site.x.interval = n_distinct(grid.id,Interval)) %>% arrange(.by_group=TRUE) %>% ungroup() %>% st_drop_geometry()

# Remove the species with less than the minimum desired number of records
message("Removing ", length(which(Species_table.0$sprecords < records_threshold)), " species (", round(100*length(which(Species_table.0$sprecords < records_threshold)) / nrow(Species_table.0), digits=2), "%) with less than ", records_threshold, " records." )
SyrphidDataset <- filter(SyrphidDataset, genusspecies %in% Species_table.0$genusspecies[which(Species_table.0$sprecords >= records_threshold)])

# Remove the species with less than the minimum desired number of eras
message("Removing ", length(which(Species_table.0$sprecords >= records_threshold & Species_table.0$num_eras < eras_threshold)), " species (", round(100*length(which(Species_table.0$sprecords >= records_threshold & Species_table.0$num_eras < eras_threshold)) / length(which(Species_table.0$sprecords >= records_threshold)), digits=2), "%) recorded in fewer than ", eras_threshold ," eras.")
SyrphidDataset <- filter(SyrphidDataset, genusspecies %in% Species_table.0$genusspecies[which(Species_table.0$num_eras >= eras_threshold)])

# Remove Volucella bombylans, which is a species of debatable validity in North America (see Jeff Skevington's field guide).
# Remove Criorhina willistoni, which is likely an unpublished species from Kevin Moran's work at the CNC (records only exist in the CNC dataset).
SyrphidDataset <- filter(SyrphidDataset, !(genusspecies %in% c("Volucella bombylans","Criorhina willistoni")))

# Remove Chrysotoxum specimens from outside the CNC (revision in progress)
SyrphidDataset <- filter(SyrphidDataset, !(Genus=="Chrysotoxum" & Collection!="CNC"))

# Remove Dasysyrphus specimens from outside the CNC that were identified before 2013 (the year of the revision by Locke & Skevington)
SyrphidDataset <- filter(SyrphidDataset,!(Genus=="Dasysyrphus" & Collection!="CNC" & (IdentificationYear<2013|(IdentificationYear=="NA"&Year<2013))))

# Filter the Species Table to only the species being analyzed
Species_table <- SyrphidDataset %>% group_by(genusspecies) %>% summarise(sprecords = n(), num_eras = n_distinct(Era), num_grids = n_distinct(grid.id), distinct_site.x.interval = n_distinct(grid.id,Interval)) %>% arrange(.by_group=TRUE) %>% ungroup() %>% st_drop_geometry()

# Add column denoting genus
Species_table$Genus <- regmatches(Species_table$genusspecies, regexpr("\\b\\w+\\b", Species_table$genusspecies))

# Add column denoting subfamily
Eristalinae_genera <- c("Anasimyia","Arctosyrphus","Blera","Brachyopa","Brachypalpus","Caliprobola","Callicera","Ceriana","Chalcosyrphus","Cheilosia","Chrysogaster","Chrysosyrphus","Copestylum","Criorhina","Eristalinus","Eristalis","Eumerus","Eurimyia","Ferdinandea","Hadromyia","Hammerschmidtia","Helophilus","Hiatomyia","Lejops","Lejota","Mallota","Merodon","Meromacrus","Milesia","Myathropa","Myolepta","Nausigaster","Neoascia","Orthonevra","Ornidia","Palpada","Parhelophilus","Pelecocera","Pocota","Polydontomyia","Psilota","Pterallastes","Pyritis","Rhingia","Sericomyia","Somula","Sphecomyia","Sphegina","Sphiximorpha","Spilomyia","Syritta","Temnostoma","Teuchocnemis","Tropidia","Volucella","Xylota")
Syrphinae_genera <- c("Allograpta","Baccha","Chrysotoxum","Dasysyrphus","Didea","Dioprosopa","Doros","Epistrophe","Epistrophella","Eriozona","Eupeodes","Fazia","Hypocritanus","Lapposyrphus","Leucozona","Megasyrphus","Melangyna","Melanostoma","Meligramma","Meliscaeva","Metasyrphus","Ocyptamus","Paragus","Parasyrphus","Pelecinobaccha","Philhelius","Platycheirus","Pseudodoros","Pseudoscaeva","Pyrophaena","Scaeva","Sphaerophoria","Syrphus","Toxomerus","Xanthogramma")
Pipizinae_genera <- c("Heringia", "Neocnemodon","Pipiza","Trichopsomyia")
Microdontinae_genera <- c("Laetodon","Microdon","Mixogaster","Omegasyrphus","Serichlamys")
Subfamily_table <- rbind(data.frame(Genus=Eristalinae_genera,Subfamily="Eristalinae"),data.frame(Genus=Syrphinae_genera,Subfamily="Syrphinae"),data.frame(Genus=Pipizinae_genera,Subfamily="Pipizinae"),data.frame(Genus=Microdontinae_genera,Subfamily="Microdontinae"))
Species_table <- left_join(Species_table, Subfamily_table, by="Genus")
rm(Eristalinae_genera,Syrphinae_genera,Pipizinae_genera,Microdontinae_genera)

# Add column denoting larval trophic guild based on genus (guild info is taken from Klymko et al. 2023, Skevington et al. 2019, and various online sources such as BugGuide.net)
Saprophagous_genera <- c("Anasimyia","Arctosyrphus","Blera","Brachyopa","Brachypalpus","Caliprobola","Callicera","Ceriana","Chalcosyrphus","Chrysogaster","Chrysosyrphus","Copestylum","Criorhina","Eristalinus","Eristalis","Eumerus","Eurimyia","Ferdinandea","Hammerschmidtia","Helophilus","Lejops","Lejota","Mallota","Meromacrus","Milesia","Myathropa","Myolepta","Nausigaster","Neoascia","Orthonevra","Ornidia","Palpada","Parhelophilus","Pocota","Polydontomyia","Psilota","Pyritis","Rhingia","Sericomyia","Somula","Sphecomyia","Sphegina","Sphiximorpha","Spilomyia","Syritta","Temnostoma","Tropidia","Xylota")
Predator_genera <- c("Allograpta","Baccha","Chrysotoxum","Dasysyrphus","Didea","Dioprosopa","Doros","Epistrophe","Epistrophella","Eriozona","Eupeodes","Heringia","Hypocritanus","Lapposyrphus","Leucozona","Megasyrphus","Melangyna","Melanostoma","Meligramma","Meliscaeva","Neocnemodon","Ocyptamus","Paragus","Parasyrphus","Pelecinobaccha","Pipiza","Platycheirus","Pseudoscaeva","Pyrophaena","Scaeva","Sphaerophoria","Syrphus","Toxomerus","Trichopsomyia","Xanthogramma")
Phytophagous_genera <- c("Cheilosia", "Fazia", "Merodon")
BroodParasite_genera <- c("Laetodon", "Microdon", "Mixogaster", "Omegasyrphus", "Serichlamys", "Volucella")
Unknown_genera <- c("Cynorhinella","Hadromyia","Hiatomyia","Pelecocera","Pterallastes","Teuchocnemis")
Larva_table <- rbind(data.frame(Genus=Saprophagous_genera,LarvalGuild="Saprophagous"),data.frame(Genus=Predator_genera,LarvalGuild="Predator"),data.frame(Genus=Phytophagous_genera,LarvalGuild="Phytophagous"),data.frame(Genus=BroodParasite_genera,LarvalGuild="BroodParasite"),data.frame(Genus=Unknown_genera,LarvalGuild="Unknown"))
Species_table <- left_join(Species_table, Larva_table, by="Genus")
rm(Saprophagous_genera,Predator_genera,Phytophagous_genera,BroodParasite_genera,Unknown_genera)

# Save species table
#write.csv(Species_table, file="FinalModel/Species_table.csv")






##################################
#### Preparing Occupancy Data ####
##################################

# Create list of species to be analyzed
splist <- sort(unique(SyrphidDataset$genusspecies))
save(splist, file="FinalModel/splist.RData")
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
  
  # In case you wish to plot the assumed range:
  #plot(m.states$geometry)
  #plot(convex.hull,border=rgb(1,0,0, alpha=0.6),lwd=2,add=T)
  #plot(geometry.observed$geometry, col=rgb(1,0,0, alpha=0.7), border=rgb(0,0,0, alpha=0.05), add=T)
  #plot(assumed.range$geometry,col=rgb(1,0,0, alpha=0.1),border=rgb(0,0,0, alpha=0.05), add=T)
  #mtext(paste(splist[sp],", n = ",Species_table$sprecords[sp],sep=""), cex=2)
  
  # Return a logical vector of all sites denoting whether they are within the assumed range
  logical.range <- ifelse(1:nrow(grid) %in% assumed.range$grid.id, TRUE, FALSE)
  return(logical.range)
}

# sp.range is a matrix of [species x site] with logical data (using function above). It shows, for each species, which sites are within their assumed range.
sp.range <- do.call(rbind, lapply(1:length(splist),get.range))
colnames(sp.range) <- paste("grid_",1:nrow(grid),sep='')
# Add range size column to species table
Species_table$range.size <- NA
for (row in 1:nrow(Species_table)){
Species_table$range.size[row] <- length(which(sp.range[row,])) 
}

# vis.arr is an empty array of each "survey" or "visit" [site x era x interval]. Each interval in a given era is assumed to be a different "visit" to a site. In this analysis, we will confine the model to only the relevant visits for each species (more on that below).
vis.arr <-array(data=NA, dim=c(nrow(grid), n_distinct(SyrphidDataset$Era), EraSize/IntervalSize), dimnames=list(paste("grid_",1:nrow(grid),sep=''),sort(as.vector(unique(SyrphidDataset$Era))),paste("interval_",1:(EraSize/IntervalSize),sep='')))

# X is [species x site x era x interval]. This is the observation data gathered from the presence/absence of specimen records. 0 for undetected, 1 for detected.
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

# Combine these data objects into a single list object
OccData <- list(sp.range=sp.range,vis.arr=vis.arr,X=X, nsp=length(splist), nsite=nrow(grid), nera=120/EraSize, nvisit=EraSize/IntervalSize)

# Keep only sites that yielded a detection of at least one species. Note from this point forward gridcell ID numbers will not match site index numbers, due to the removal of gridcells that had no species (e.g. the 100th site kept in the matrix may actually be grid_241).
site.keep <- which(apply(OccData$X, 'site', sum)>0)
OccData$X        <- OccData$X[,site.keep,,,drop=FALSE] #Do not drop dimensions (ex. 4D object must stay 4D)
OccData$sp.range <- OccData$sp.range[,site.keep,drop=FALSE]
OccData$vis.arr  <- OccData$vis.arr[site.keep,,,drop=FALSE]
OccData$nsite    <- length(site.keep)

# Now we will filter data to only the relevant visits (visit = combination of site, era, interval) 
# A relevant visit is one that occurs in a site within the species' range 
# and in an interval when at least 1 syrphid of ANY species was detected 
# (meaning we know someone was at that site at that time sampling for syrphids, and 
# had the potential to detect the species in question. This justifies our non-detections).
# (this also prevents unnecessary iterating through all irrelevant sites and visits, improving model efficiency). 

# Function to determine indexes of which visits are relevant for a given species
get.indices <- function(sp) {
  vis.arr <- OccData$vis.arr #Begin with a copy of vis.arr
  vis.arr[TRUE] <- 1 #Set all visits to 1 (i.e., assume all visits are relevant to begin with)
  nsp.detected <- apply(OccData$X, 2:4, sum) #object with same dimensions as vis.arr, recording the number of species detected during each visit (visit = site x era x interval... aka dimensions 2:4 of X)
  vis.arr[nsp.detected==0] <- 0 #Set a visit to 0 ("irrelevant") for this species if no species were detected during that visit/interval
  vis.arr[!OccData$sp.range[sp,],,] <- 0 #Set a visit to 0 if the site is not in species range (i.e. sp.range==FALSE)
  
  tmp <- which(vis.arr==1, arr.ind=TRUE) #Find which cells in vis.arr are 1 (i.e., "relevant") for this species. arr.ind=TRUE means array indices are returned, because vis.arr is an array. Each index will direct you to a 1 within vis.arr[site,era,interval]. dim1 is site, dim2 is era and dim3 is interval.

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
















################################################################################
################### ANNOTATED COPY of the Occupancy Model ######################
################################################################################
# Written in BUGS language to be used by JAGS. Not written in R; THIS CODE IS NOT USED.
# This is presented for viewing/annotating only. This section of code is not meant to be run!
# Rjags pulls the model from an external file: "SyrphidModel.txt" (which is identical to the code below, but without the annotations). 

#model{
  
  ### Priors ###
  
  # Mean (i.e., intercept) occupancy and detection probability across species, mu.psi and mu.p
#  mu.psi ~ dnorm(0,0.01) # Sets prior distribution for mean occupancy. JAGS uses precision instead of sd to describe normal distributions. In this case, mean 0 and precision 0.01. SD = sqrt(1/precision) = 10. Low precision means high SD.
#  mu.p   ~ dnorm(0,0.01) # prior for mean detection

  # Random effect of species identity on occupancy and detection, psi.sp and p.sp
  # This is a hierarchical prior, where the prior distribution of each related parameter (i.e. each "effect of species on y") comes from a hyperprior distribution, which represents our knowledge about the entire population of these parameters.
  # The shape of the hyperprior distribution gets updated as each related parameter is estimated. Thus, the model is essentially learning the priors from the sample itself.
#  sigma.psi.sp ~ dunif(0,10) #hyperprior sd in the effects of species on occupancy. Uniform distribution between 0 and 10 on the x axis.
#  sigma.p.sp   ~ dunif(0,10) #hyperprior sd in the effects of species on detection.
#  tau.psi.sp  <- 1/(sigma.psi.sp*sigma.psi.sp) #precision (1/sd^2), because we will be using a normal distribution that needs precision instead of sd
#  tau.p.sp    <- 1/(sigma.p.sp*sigma.p.sp)
#  for(sp in 1:nsp) { #Prior for the random effect. The effect of each species gets its own parameter whose posterior will update based on the observations and on the hyperprior, which also updates as we learn more about each parameter. Random effects are set to mean 0 because they are expressed as deviations from the global mean.
#    psi.sp[sp] ~ dnorm(0, tau.psi.sp)
#    p.sp[sp]   ~ dnorm(0, tau.p.sp)
#  }

  # Random effect of site*era on detection, p.site (includes the effect of sampling effort between sites and eras)
#  sigma.p.site   ~ dunif(0,10) #hyperprior sd in detection across sites
#  tau.p.site    <- 1/(sigma.p.site*sigma.p.site)
#  for(site in 1:nsite) { #Each site*era has its own effect on detection (its own parameter).
#    for(era in 1:nera) {
#      p.site[site,era]   ~ dnorm(0, tau.p.site) 
#    }
#  }

  # Effect of species*era on occupancy, psi.era
#  for (era in 2:nera) {
#    mu.psi.era[era]    ~ dnorm(0,0.01) #hyperprior mean (across-species) effect of an era on occupancy
#    sigma.psi.era[era] ~ dunif(0,10) #hyperprior sd in species-specific effects of an era on occupancy
#    tau.psi.era[era]  <- 1/(sigma.psi.era[era]*sigma.psi.era[era])
#  }
#  for(sp in 1:nsp) {
#    psi.era[sp,1] <- 0  #Effect of era 1 is set to 0, so that it becomes the default effect around which to compare the other eras.
#    for(era in 2:nera) {
#      psi.era[sp,era] ~ dnorm(mu.psi.era[era], tau.psi.era[era]) #Each species is differently affected by each era on its occupancy
#    }
#  }
  
  # Effect of era on detection, p.era (further controls for sampling effort between eras)
#  p.era[1] <- 0 #Effect of era 1 is set to 0, so that it becomes the default effect around which to compare the other eras.
#  for(era in 2:nera) {
#    p.era[era] ~ dnorm(0,0.01)
#  }
  
  ### Occupancy and Detection probabilities ###
  
#  for(sp in 1:nsp) {
#    for(era in 1:nera) {
      # occupancy
#      logit(psi[sp,era]) <- #The occupancy (proportion of sites) of a given species in a given era is a function of the baseline occupancy plus the effect of this species on occupancy plus the effect of this era on occupancy (for this sp). Logit link is used to constrain psi[sp,era] to a probability value (between 0 and 1)
#        mu.psi +
#        psi.sp[sp] +
#        psi.era[sp,era]
      # detection
#      for(site in 1:nsite) { 
#        logit(p[sp,site,era]) <- #The detection probability of a species in a given site and era is a function of the baseline detection plus the effect of species on detection plus the effect of site(per era) on detection plus the effect of era on detection. Logit link is used to constrain to a probability value (between 0 and 1)
#          mu.p +
#          p.sp[sp] +
#          p.site[site,era] +
#          p.era[era]
#      }
#    }
#  }
  
  ### Latent Z state for each site ###
  # Z is the true presence of the species. Either a 0 or 1 and can only change between eras. Bernoulli trial with prob. of success equal to psi[sp,era] (which was defined above)
  
#  for(sp in 1:nsp) {
#    for(site in 1:nsite) {
#      for(era in 1:nera) {
#        Z[sp,site,era] ~ dbern(psi[sp,era])
#      }
#    }
#  }

  ### Model the observations ###
  # Z is true (latent) occupancy state, p is prob of detection, and X is the observations
  
#  for(ind in 1:nind) {
#    p.eff[ind] <-
#      Z[sp[ind],site[ind],era[ind]] *
#      p[sp[ind],site[ind],era[ind]]
#    X[ind] ~ dbern(p.eff[ind]) #The observations, X, are bernoulli trials (0 or 1) whose outcome depends on latent state Z (zeroes are guaranteed non-detections) and detection probability p for that sp*site*era)
#  }
#}
#############################################################################
#############################################################################
#############################################################################














###########################
#### Running the Model ####
###########################

# Specify the parameters to be monitored
get.params <- function() {
  c('mu.psi',         #Baseline occupancy of syrphids across species and eras (includes default effect of era 1)
    'mu.p',           #Baseline/mean detection probability of syrphids across species, sites, and eras (includes default effect of era 1)
    'psi.sp',         #Random Effect of a certain species on occupancy
    'sigma.psi.sp',   #Hyperprior sd for effects of species on occupancy. Represents variation amongst the effects of different species on occupancy.
    'p.sp',           #Random Effect of a certain species on detection
    'sigma.p.sp',     #Hyperprior sd for effects of species on detection. Represents variation amongst the effects of different species on detection.
    'p.site',         #Random effect of site*era on detection. Includes effect of sampling effort between sites and eras.
    'sigma.p.site',   #Hyperprior sd for effects of site*era on detection. Represents variation amongst the effects of different site*eras on detection.
    'psi.era',        #Random species-specific effect of era on occupancy
    'mu.psi.era',     #Hyperprior mean for effects of species*era on occupancy. Represents the mean value amongst the different species' effects on occupancy in an era. In other words, the mean across-species effect of an era on occupancy.
    'sigma.psi.era',  #Hyperprior sd for effects of species*era on occupancy. Represents variation amongst the different species' effects on occupancy in an era. In other words, the variation in the species-specific effects of an era on occupancy.
    'p.era') }        #Effect of era on detection.

# Function to run the model
run.model <- function(relevant.dat, n.iter, n.burnin, n.adapt, n.thin) { 

  model.txt <- sprintf('SyrphidModel.txt') #Bring in the BUGS-syntax model from the external source file

  # Create data object that will be passed to the model
  jags.data <- list(X=relevant.dat$X,                            #Observations for all species for their relevant visits
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
#model.out <- run.model(relevant.dat,
#                       n.iter=50e3,   #Number of iterations of each MCMC chain. How many samples to take.
#                       n.burnin=4e3, #Number of burn-in/warmup iterations. Reduces dependence of the chain on the starting values of MCMC sampling.
#                       n.adapt=2e3,  #These iterations are used to tune the behaviour of the MCMC sampler, making it more efficient
#                       n.thin=1e1)   #Thinning interval to increase independence of samples within each chain. Every 10th sample is kept. In other words, reduces autocorrelation and increases effective sample size.






#######################
#### Model Outputs ####
#######################

# Save/load the model output on your computer
#save(model.out,file="FinalModel/model.out.RData")     #SAVE
#load("FinalModel/model.out.RData")                    #LOAD

# Create Model Summary
#ModelSummary <- summary(model.out$jags.out)
# Save/load the model summary on your computer
#save(ModelSummary,file="FinalModel/ModelSummary.RData")  #SAVE
load("FinalModel/ModelSummary.RData")                    #LOAD

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
#parameter.to.plot <- "mu.psi" #Identify the parameter for which you wish to view the traceplot and posterior density
#plot(MCMC[,parameter.to.plot], main = parameter.to.plot) #plot the mcmc output (all samples/chains) for the desired parameter. The mcmc.list object within model.out is a list of 4 chains (each an mcmc object), each with the dimensions [#samples,#parameters]. main is the title



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
for (era in 1:relevant.dat$nera){ #Compute cross-species occupancies for each era
  meanoccupancies[[paste(era)]] <- compute_cross.sp_occupancy(era)}
Occupancy_list[["cross.sp.mean"]] <- meanoccupancies #add the cross-species occupancy posteriors to the master list
rm(meanoccupancies)
for (sp in 1:relevant.dat$nsp) { #Compute species-specific occupancies
  spoccupancies <- list()
  for (era in 1:relevant.dat$nera) {
    spoccupancies[[paste(era)]] <- compute_occupancy(sp,era)}
  Occupancy_list[[splist[sp]]] <- spoccupancies #add this species' occupancy posteriors to the master list
  rm(spoccupancies)}

# Create a dataframe of summary statistics (mean, sd, 95% intervals) for all occupancy estimates.
get_occ_stats <- function(row){ #Function to get the summary stats for each era's occupancy for a given species
  dataframe <- matrix(nrow=relevant.dat$nera,ncol=6) #create matrix to hold data.
  for (era in 1:relevant.dat$nera){ #Get the species name, era name, lower95, mean, upper95, and SD.
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
#save(OccupancySummary,file="FinalModel/OccupancySummary.RData")
# Load Occupancy Summary
load("FinalModel/OccupancySummary.RData")






#########################################################
#### Estimates of Occupancy Change (delta occupancy) ####
#########################################################

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
  occ.changes <- list() #create empty list to hold the estimates of occupancy change for this species (each row in Occupancy_list is a different species)
  for (era in relevant.dat$nera:2) {
      for (nextera in (era-1):1) {
        occ.changes[[paste(era,"-",nextera, sep="")]] <- compute_delta.occ(row,era,nextera)}} #Occupancy change is computed between the era in question and every era below it. Example with 4 eras: 4-3,4-2,4-1,3-2,3-1,2-1
  delta.occ_list[[names(Occupancy_list[row])]] <- occ.changes #add this species' contrasts to the master list
  rm(occ.changes,row,era,nextera)}

# Create a dataframe of summary statistics (mean, sd, 95% intervals) for all occupancy change estimates.
get_delta_occ_stats <- function(row){ #Function to get the summary stats for all occupancy changes for a given species
  dataframe <- matrix(nrow=relevant.dat$nera*(relevant.dat$nera-1)/2,ncol=6) #create matrix to hold data. The number of pairwise contrasts is always nera*(nera-1)/2
  for (contrast in 1:(relevant.dat$nera*(relevant.dat$nera-1)/2)){ #For each contrast, get the species name, contrast name, lower95, mean, upper95, and SD.
    dataframe[contrast,] <- c(names(Occupancy_list[row]), names(delta.occ_list[[row]][contrast]), HPDinterval(as.mcmc(unlist(delta.occ_list[[row]][[contrast]])), prob=0.95)[[1]] , summary(delta.occ_list[[row]][[contrast]])$statistics[[1]], HPDinterval(as.mcmc(unlist(delta.occ_list[[row]][[contrast]])), prob=0.95)[[2]] , summary(delta.occ_list[[row]][[contrast]])$statistics[[2]])}
  return(dataframe)}
delta.occSummary <- as.data.frame(do.call(rbind,lapply(1:length(Occupancy_list), get_delta_occ_stats))) #Create the dataframe
colnames(delta.occSummary) <- c("Species","EraContrast","Lower95","Mean","Upper95","SD")

# Make sure columns are treated as numeric values
delta.occSummary$Lower95 <- as.numeric(delta.occSummary$Lower95)
delta.occSummary$Mean <- as.numeric(delta.occSummary$Mean)
delta.occSummary$Upper95 <- as.numeric(delta.occSummary$Upper95)
delta.occSummary$SD <- as.numeric(delta.occSummary$SD)

# Save Delta Occupancy Summary
#save(delta.occSummary,file="FinalModel/delta.occSummary.RData")
# Load Delta Occupancy Summary
load("FinalModel/delta.occSummary.RData")







###############################################
#### Lists of Declining/Increasing Species ####
###############################################


# List of negative delta occupancy posteriors whose 95% credible intervals do not include zero (i.e., declines)
Negative_delta.occ <- filter(delta.occSummary, Lower95 < 0 & Upper95 <= 0)
Negative_sp <- unique(filter(Negative_delta.occ, Species!="cross.sp.mean")$Species) # List of species with strongly supported declines (do not include the cross.sp.mean)

# List of positive delta occupancy posteriors whose 95% credible intervals do not include zero (i.e., increases)
Positive_delta.occ <- filter(delta.occSummary, Lower95 >= 0 & Upper95 > 0)
Positive_sp <- unique(filter(Positive_delta.occ, Species!="cross.sp.mean")$Species) # List of species with strongly supported increases (do not include the cross.sp.mean)

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

# List of species with historical decline (declined between the 2 terminal eras)
HistoricDecline_delta.occ <- filter(delta.occSummary, Lower95 < 0 & Upper95 <= 0 & EraContrast==paste(relevant.dat$nera,"-1",sep=""))
HistoricDecline_sp <- unique(filter(HistoricDecline_delta.occ, Species!="cross.sp.mean")$Species)

# List of species with historical growth (increased between the 2 terminal areas)
HistoricGrowth_delta.occ <- filter(delta.occSummary, Lower95 >= 0 & Upper95 > 0 & EraContrast==paste(relevant.dat$nera,"-1",sep=""))
HistoricGrowth_sp <- unique(filter(HistoricGrowth_delta.occ, Species!="cross.sp.mean")$Species)

# List of species with no strongly supported recent or historic change in occupancy
Stable_sp <- splist[!(splist %in% c(RecentDecline_sp,HistoricDecline_sp,RecentGrowth_sp,HistoricGrowth_sp))]

# Create Results dataframe with all of these lists
max_length <- max(length(RecentDecline_sp),length(ExceptionalRecentDecline_sp),length(HistoricDecline_sp),length(RecentGrowth_sp),length(ExceptionalRecentGrowth_sp),length(HistoricGrowth_sp),length(Stable_sp)) #determine length of columns needed
Results.table <- data.frame(RecentDecline = c(paste("Declined between", EraLabels[relevant.dat$nera-1] ,"and", EraLabels[relevant.dat$nera]),RecentDecline_sp,rep("",max_length-length(RecentDecline_sp))),
                            ExceptionalRecentDecline = c(paste("Occupancy in",EraLabels[relevant.dat$nera],"is lower than all historical time periods. (A subset of Recent Decline and Historical Decline)"),ExceptionalRecentDecline_sp,rep("",max_length-length(ExceptionalRecentDecline_sp))),
                            HistoricDecline = c(paste("Declined between", EraLabels[1] ,"and", EraLabels[relevant.dat$nera]),HistoricDecline_sp,rep("",max_length-length(HistoricDecline_sp))),
                            RecentGrowth = c(paste("Increased between", EraLabels[relevant.dat$nera-1] ,"and", EraLabels[relevant.dat$nera]),RecentGrowth_sp,rep("",max_length-length(RecentGrowth_sp))),
                            ExceptionalRecentGrowth = c(paste("Occupancy in",EraLabels[relevant.dat$nera],"is greater than all historical time periods (A subset of Recent Growth and Historical Growth)"),ExceptionalRecentGrowth_sp,rep("",max_length-length(ExceptionalRecentGrowth_sp))),
                            HistoricGrowth = c(paste("Increased between", EraLabels[1] ,"and", EraLabels[relevant.dat$nera]),HistoricGrowth_sp,rep("",max_length-length(HistoricGrowth_sp))),
                            Stable = c("No significant change over recent or historical time frames", Stable_sp,rep("",max_length-length(Stable_sp))))

# Save or load Results dataframe
#write.csv(Results.table,"FinalModel/Results.table.csv",row.names=F) #Save results
#Results.table <- read.csv("FinalModel/Results.table.csv") #Load results








##############################
### Prior Predictive Check ###
##############################
# Run model on a dataset consisting of only NA values, in order to only sample from prior.
relevant.dat.NAs <- relevant.dat
relevant.dat.NAs$X <- rep(NA,length(relevant.dat$X))
PriorCheck <- run.model(relevant.dat.NAs,
                        n.iter=50e3,
                        n.burnin=4e3,
                        n.adapt=2e3, 
                        n.thin=1e1)
# Save output for analysis in another script. See PriorPredictiveCheck.R
save(PriorCheck,file="PriorPredCheck/PriorCheck.RData")





###################
### Simulations ###
###################
# Simulations are run in a different script. See Simulation.R





##################################
### Posterior Predictive Check ###
##################################
# For each observation in X [species x site x era x interval], use the posteriors to generate a simulated observation.
# If pattern of 0's and 1's are often similar between simulated vs. real data, the model is accurate!
run.PostPredictiveCheck <- function(repeatsimulations){ #Function that takes the number of replicate post predictive checks to run, and returns a list object of the summary for each run
  X.simulated <- vector(mode="numeric", length=nrow(relevant.dat$master.index))
  for (n in repeatsimulations) { #Number of repeat simulation datasets to create (i.e., replicates). Results from each one will be listed in PostPredictiveCheck.
    
    # Step 1: Simulate data
    # Step 1a: Simulate latent Z states
    Z <- array(NA, dim=c(relevant.dat$nsp,relevant.dat$nsite,relevant.dat$nera)) #Create empty Z array
    for (spp in 1:relevant.dat$nsp){ #new index variable cannot be the same as a column name in master.index, so alter the index name slightly (i.e., sp becomes spp, era becomes eraa, etc.)
      for (eraa in 1:relevant.dat$nera){
        for (sitee in unique(filter(as.data.frame(relevant.dat$master.index), sp==spp & era==eraa)$site)){ #For each relevant site within each era for each species, simulate Z state. (all species will have at least 1 observation in each era)
          # Simulate an Occupancy Probability for this species*era
          mu.psi <- sample(as.vector(as.matrix(MCMC[,"mu.psi"])),1) #Draw a single sample from the posterior (better than using the means, because it incorporates the uncertainty of the posterior)
          psi.sp <- sample(as.vector(as.matrix(MCMC[,paste("psi.sp[",spp,"]",sep="")])),1)
          psi.era <- sample(as.vector(as.matrix(MCMC[,paste("psi.era[",spp,",",eraa,"]",sep="")])),1)
          psi <- plogis(mu.psi+psi.sp+psi.era)
          # Simulate a Latent State Z (whether the site occupied) for this species*site*era
          Z[spp,sitee,eraa] <- rbinom(1,1,psi)
        }}}
    
    # Step 1b: Simulate observations for each observation X in master.index (only relevant visits are simulated)
    for (x in 1:nrow(relevant.dat$master.index)) { 
      # Simulate a detection probability for this species*site*era
      if(x==1 | any(c(relevant.dat$master.index[x-1,1],relevant.dat$master.index[x-1,2],relevant.dat$master.index[x-1,3])!=c(relevant.dat$master.index[x,1],relevant.dat$master.index[x,2],relevant.dat$master.index[x,3]))){ #if previous row was the same [sp,site,era] (i.e., row x differs from x-1 in interval only), do not sample a new p.eff value; use the previous p.eff value
        mu.p <- sample(as.vector(as.matrix(MCMC[,"mu.p"])),1)
        p.sp <- sample(as.vector(as.matrix(MCMC[,paste("p.sp[",relevant.dat$master.index[x,1],"]",sep="")])),1)
        p.era <- sample(as.vector(as.matrix(MCMC[,paste("p.era[",relevant.dat$master.index[x,3],"]",sep="")])),1)
        p.site <- sample(as.vector(as.matrix(MCMC[,paste("p.site[",relevant.dat$master.index[x,2],",",relevant.dat$master.index[x,3],"]",sep="")])),1)
        p <- plogis(mu.p+p.sp+p.era+p.site)
        # Combine Z state and detection probability into the effective detection probability
        p.eff <- Z[relevant.dat$master.index[x,1],relevant.dat$master.index[x,2],relevant.dat$master.index[x,3]]*p #if Z is 0, detection is impossible
      }
      # Simulate observation X
      X.simulated[x] <- rbinom(1,1,p.eff)
    }
    
    # Step 2: Compare real vs. simulated datasets. The following statistics will be recorded for each species as well as overall:
    # Record percent similarity (what percent of visits (site*era*interval) are identical between real vs. simulated)
    # Record the number of detections (i.e. the number of 1's) in both simulated and real datasets. Also record the difference in number of detections between simulated vs. actual. Simulated - Actual.
    # Record the number of simulated non-detections that were detections in reality ("Wrong Zeroes"), and the number of simulated detections that were non-detections in reality ("Wrong Ones").
    # Record the number of correct (i.e., matching) ones and zeroes between the simulated and real datasets.
    PostPredictiveCheck <- data.frame(SpeciesNum="Total (across all species)",SpeciesName="NA", NumObservations=length(relevant.dat$X),
                                      PercentSimilarity=100*as.numeric(summary(relevant.dat$X==X.simulated)[3])/length(relevant.dat$X),
                                      PositiveObs.Simulated=sum(X.simulated), PositiveObs.real= sum(relevant.dat$X),
                                      DifferenceInPositiveObs=sum(X.simulated)-sum(relevant.dat$X),
                                      WrongZeroes=length(which((X.simulated-relevant.dat$X)==-1)),
                                      WrongOnes=length(which((X.simulated-relevant.dat$X)==1)),
                                      RightZeroes=length(which(X.simulated==0 & relevant.dat$X==0)),
                                      RightOnes=length(which(X.simulated==1 & relevant.dat$X==1)))
    for(i in 1:relevant.dat$nsp){ #Calculate the results for each species
      PostPredictiveCheck[i+1,] <- data.frame(SpeciesNum=i,SpeciesName=splist[i], NumObservations=length(relevant.dat$X[which(relevant.dat$master.index[,1]==i)]),
                                              PercentSimilarity= 100*as.numeric(summary(relevant.dat$X[which(relevant.dat$master.index[,1]==i)]==X.simulated[which(relevant.dat$master.index[,1]==i)])[3])/length(which(relevant.dat$master.index[,1]==i)),
                                              PositiveObs.Simulated= sum(X.simulated[which(relevant.dat$master.index[,1]==i)]), PositiveObs.real= sum(relevant.dat$X[which(relevant.dat$master.index[,1]==i)]),
                                              DifferenceInPositiveObs=sum(X.simulated[which(relevant.dat$master.index[,1]==i)])-sum(relevant.dat$X[which(relevant.dat$master.index[,1]==i)]),
                                              WrongZeroes=length(which((X.simulated[which(relevant.dat$master.index[,1]==i)]-relevant.dat$X[which(relevant.dat$master.index[,1]==i)])==-1)),
                                              WrongOnes=length(which((X.simulated[which(relevant.dat$master.index[,1]==i)]-relevant.dat$X[which(relevant.dat$master.index[,1]==i)])==1)),
                                              RightZeroes=length(which(X.simulated[which(relevant.dat$master.index[,1]==i)]==0 & relevant.dat$X[which(relevant.dat$master.index[,1]==i)]==0)),
                                              RightOnes=length(which(X.simulated[which(relevant.dat$master.index[,1]==i)]==1 & relevant.dat$X[which(relevant.dat$master.index[,1]==i)]==1)))}
    if (n==1){ #Create list of result objects
      PostPredictiveChecks <<- list(PostPredictiveCheck) 
    }else{
      PostPredictiveChecks <<- append(PostPredictiveChecks,list(PostPredictiveCheck))}}
}

run.PostPredictiveCheck(repeatsimulations=1:3) #Run the Posterior Predictive Check. Create 3 simulated datasets.

#save(PostPredictiveChecks, file="PosteriorPredictiveChecks.RData") #save results
load("PosteriorPredictiveChecks.RData") #load results
#View(PostPredictiveChecks[[1]]) #View a Posterior Predictive Check
#mean(PostPredictiveChecks[[3]][-1,]$PercentSimilarity) #View mean percent similarity of a Posterior Predictive Check

# Function for Stacked Barplot of Percent Similarity of observations for each species, for a given simulated dataset
plot.PercentSimilarity.PPC <- function(PPCnum){
  PostPredictiveCheck.plot <- arrange(PostPredictiveChecks[[PPCnum]][-1,],PercentSimilarity) #sort the species in the PPC by Percent Similarity, but remove the first row that displays totals.
  StackedBars <- rbind(100*PostPredictiveCheck.plot$RightZeroes/PostPredictiveCheck.plot$NumObservations,
                       100*PostPredictiveCheck.plot$RightOnes/PostPredictiveCheck.plot$NumObservations,
                       100*PostPredictiveCheck.plot$WrongZeroes/PostPredictiveCheck.plot$NumObservations,
                       100*PostPredictiveCheck.plot$WrongOnes/PostPredictiveCheck.plot$NumObservations)
  barplot(StackedBars, ylim=c(0,100), col=c("#348FA4","#74C0D2","#FFE479","#F6FEAA"), ylab="Percentage of the species' dataset (%)", xlab=paste("Species Identity (n=",length(splist),")",sep=""), space=0, border="transparent", cex.axis=1.5, cex.lab=2, las=1) #Create stacked bar plot
  #axis(side = 2, at = seq(0, 100, by = 10)) #y-axis ticks go up by 10's
  abline(a=mean(PostPredictiveChecks[[PPCnum]][-1,]$PercentSimilarity), b=0, lty=1, col="red") #average percent similarity
  legend(c(170,315),c(36,1),c("Simulated 1, Real 0", "Simulated 0, Real 1","Matching 1's","Matching 0's"), fill = c("#F6FEAA","#FFE479","#74C0D2","#348FA4"), cex = 2, inset = c(0.1, 0.35), y.intersp = 0.7, x.intersp =0.2)}
# Plotting
#par(mar = c(5, 5, 4, 2) + 0.1, mfrow = c(1,1)) #Increase left margin so y axis title doesnt get cut off
#plot.PercentSimilarity.PPC(1)
#par(mar = c(5, 4, 4, 2) + 0.1)  # Default margin values

# Function for scatterplot of Number of simulated vs. real detections for each spp, for a given simulated dataset
plot.numDetections.PPC <- function(PPCnum){
  PostPredictiveCheck.plot <- arrange(PostPredictiveChecks[[PPCnum]][-1,],PositiveObs.real) #sort the PPC by number of real detections, and remove the first row that displays totals.
  plot(PostPredictiveCheck.plot$PositiveObs.real, pch=19, col="blue", cex=0.7, ylim=c(0,1000), ylab="", xlab=paste("Species Identity (n=",length(splist),")",sep=""), xaxt="n", cex.lab=2, las=1, yaxt="n", mgp=c(2,1,0))
  points(PostPredictiveCheck.plot$PositiveObs.Simulated, pch=19, col="red", cex=0.7)
  axis(2, at = c(0,200,400,600,800,1000), cex.axis=2, las=1)}
#Plotting
#par(mar = c(5, 7, 4, 2) + 0.1) #Increase left margin so y axis title doesnt get cut off
#par(bty = "l")
#plot.numDetections.PPC(1)
#mtext("Number of Detections", cex=2, side = 2, line = 5, las = 0)
#par(mar = c(5, 4, 4, 2) + 0.1)  # Default margin values






##############################
#### Sensitivity Analysis ####
##############################
# Comparing results of models with differing spatial/temporal resolutions
# Sensitivity analysis model variants will be run in a separate, duplicate copy of the regular model code. See SensitivityModels.R
# SensitivityModels.R differs from the original model script in 3 ways: (1) either CellSize, Erasize, or IntervalSize is altered each time. (2) species must have >=500 records to be analyzed. (3) MCMC sampling is reduced (only 30k iterations, to reduce computational load)
# This section handles plotting of the results only.

# Function to plot number of species in decline and growth at various resolutions:
plot.Sensitivity.TotalDeclinesGrowths <- function(type){
  # Obtain results tables at different spatial resolutions
  Results.table.high <- read.csv(ifelse(type=="Temporal (era)","SensitivityAnalysis/6eras/Results.table.csv",ifelse(type=="Temporal (interval)","SensitivityAnalysis/10intervals/Results.table.csv","SensitivityAnalysis/50x50km/Results.table.csv"))) 
  Results.table.medium <- read.csv("SensitivityAnalysis/BaseModel/Results.table.csv")
  Results.table.low <- read.csv(ifelse(type=="Temporal (era)","SensitivityAnalysis/2eras/Results.table.csv",ifelse(type=="Temporal (interval)","SensitivityAnalysis/3intervals/Results.table.csv","SensitivityAnalysis/200x200km/Results.table.csv")))
  Declines.high.recent <- setdiff(Results.table.high[-1,]$RecentDecline[Results.table.high[-1,]$RecentDecline!=""], Results.table.high[-1,]$HistoricDecline[Results.table.high[-1,]$HistoricDecline!=""]) #Remove first row which is just a description of the column.
  Declines.high.historic <- setdiff(Results.table.high[-1,]$HistoricDecline[Results.table.high[-1,]$HistoricDecline!=""], Results.table.high[-1,]$RecentDecline[Results.table.high[-1,]$RecentDecline!=""])
  Declines.high.both <- intersect(Results.table.high[-1,]$RecentDecline[Results.table.high[-1,]$RecentDecline!=""], Results.table.high[-1,]$HistoricDecline[Results.table.high[-1,]$HistoricDecline!=""])
  Declines.medium.recent <- setdiff(Results.table.medium[-1,]$RecentDecline[Results.table.medium[-1,]$RecentDecline!=""], Results.table.medium[-1,]$HistoricDecline[Results.table.medium[-1,]$HistoricDecline!=""])
  Declines.medium.historic <- setdiff(Results.table.medium[-1,]$HistoricDecline[Results.table.medium[-1,]$HistoricDecline!=""], Results.table.medium[-1,]$RecentDecline[Results.table.medium[-1,]$RecentDecline!=""])
  Declines.medium.both <- intersect(Results.table.medium[-1,]$RecentDecline[Results.table.medium[-1,]$RecentDecline!=""], Results.table.medium[-1,]$HistoricDecline[Results.table.medium[-1,]$HistoricDecline!=""])
  Declines.low.recent <- setdiff(Results.table.low[-1,]$RecentDecline[Results.table.low[-1,]$RecentDecline!=""], Results.table.low[-1,]$HistoricDecline[Results.table.low[-1,]$HistoricDecline!=""])
  Declines.low.historic <- setdiff(Results.table.low[-1,]$HistoricDecline[Results.table.low[-1,]$HistoricDecline!=""], Results.table.low[-1,]$RecentDecline[Results.table.low[-1,]$RecentDecline!=""])
  Declines.low.both <- intersect(Results.table.low[-1,]$RecentDecline[Results.table.low[-1,]$RecentDecline!=""], Results.table.low[-1,]$HistoricDecline[Results.table.low[-1,]$HistoricDecline!=""])
  Growths.high.recent <- setdiff(Results.table.high[-1,]$RecentGrowth[Results.table.high[-1,]$RecentGrowth!=""], Results.table.high[-1,]$HistoricGrowth[Results.table.high[-1,]$HistoricGrowth!=""])
  Growths.high.historic <- setdiff(Results.table.high[-1,]$HistoricGrowth[Results.table.high[-1,]$HistoricGrowth!=""], Results.table.high[-1,]$RecentGrowth[Results.table.high[-1,]$RecentGrowth!=""])
  Growths.high.both <- intersect(Results.table.high[-1,]$RecentGrowth[Results.table.high[-1,]$RecentGrowth!=""], Results.table.high[-1,]$HistoricGrowth[Results.table.high[-1,]$HistoricGrowth!=""])
  Growths.medium.recent <- setdiff(Results.table.medium[-1,]$RecentGrowth[Results.table.medium[-1,]$RecentGrowth!=""], Results.table.medium[-1,]$HistoricGrowth[Results.table.medium[-1,]$HistoricGrowth!=""])
  Growths.medium.historic <- setdiff(Results.table.medium[-1,]$HistoricGrowth[Results.table.medium[-1,]$HistoricGrowth!=""], Results.table.medium[-1,]$RecentGrowth[Results.table.medium[-1,]$RecentGrowth!=""])
  Growths.medium.both <- intersect(Results.table.medium[-1,]$RecentGrowth[Results.table.medium[-1,]$RecentGrowth!=""], Results.table.medium[-1,]$HistoricGrowth[Results.table.medium[-1,]$HistoricGrowth!=""])
  Growths.low.recent <- setdiff(Results.table.low[-1,]$RecentGrowth[Results.table.low[-1,]$RecentGrowth!=""], Results.table.low[-1,]$HistoricGrowth[Results.table.low[-1,]$HistoricGrowth!=""])
  Growths.low.historic <- setdiff(Results.table.low[-1,]$HistoricGrowth[Results.table.low[-1,]$HistoricGrowth!=""], Results.table.low[-1,]$RecentGrowth[Results.table.low[-1,]$RecentGrowth!=""])
  Growths.low.both <- intersect(Results.table.low[-1,]$RecentGrowth[Results.table.low[-1,]$RecentGrowth!=""], Results.table.low[-1,]$HistoricGrowth[Results.table.low[-1,]$HistoricGrowth!=""])
  # Create labels for each resolution
  if(type=="Temporal (era)"){
    ResolutionLabels <- c("20","30","60")
    Xlabel <- "Era Size (yrs)"
    Ylabel <- ""
  }else{
    if(type=="Temporal (interval)"){
      ResolutionLabels <- c("3","5","10")
      Xlabel <- "Interval Size (yrs)"
      Ylabel <- ""
    }else{
      ResolutionLabels <- c("50x50","100x100","200x200")
      Xlabel <- "Grid Size (km)"
      Ylabel <- "Number of species"}}
  # Create dataframe of information to plot
  sensitivity.data <- data.frame(
    Group = rep(ResolutionLabels, each = 6),
    Bar = rep(c("Decline", "Growth"), times = 3, each = 3),
    Stack = rep(c("Recent","Historic","Both"), times = 6),
    Value = c(length(Declines.high.recent), length(Declines.high.historic), length(Declines.high.both), length(Growths.high.recent), length(Growths.high.historic), length(Growths.high.both),
              length(Declines.medium.recent), length(Declines.medium.historic),length(Declines.medium.both), length(Growths.medium.recent),length(Growths.medium.historic), length(Growths.medium.both),
              length(Declines.low.recent),length(Declines.low.historic), length(Declines.low.both), length(Growths.low.recent),length(Growths.low.historic),length(Growths.low.both)))
  sensitivity.data$BarStack <- interaction(sensitivity.data$Bar, sensitivity.data$Stack)
  custom_order <- ResolutionLabels
  sensitivity.data$Group <- factor(sensitivity.data$Group, levels = custom_order)
  sensitivity.data$BarStack <- factor(sensitivity.data$BarStack, levels = c("Decline.Recent","Decline.Historic","Decline.Both","Growth.Recent","Growth.Historic","Growth.Both"))
  # Create stacked bar graph
  ggplot(sensitivity.data, aes(x = Bar, y = Value, fill = BarStack)) +
    geom_bar(stat = "identity", position = "stack", width = 0.95) +
    facet_wrap(~ Group, ncol = 3, switch="x") + 
    ylim(0,30) +
    labs(x = Xlabel, y = Ylabel) +
    theme_minimal() +
    scale_fill_manual(values= c(
      "Decline.Both" = "#7D3438",
      "Decline.Historic" = "#AC4F55",
      "Decline.Recent" = "#E65F67",
      "Growth.Both" = "#216C86",
      "Growth.Historic" = "#2F9BBF",
      "Growth.Recent" = "#3CBCE6"), 
      labels = c(
        "Decline.Both" = "Both Recent and Historical Decline",
        "Decline.Historic" = "Historical Decline only",
        "Decline.Recent" = "Recent Decline only",
        "Growth.Both" = "Both Recent and Historical Increase",
        "Growth.Historic" = "Historical Increase only",
        "Growth.Recent" = "Recent Increase only")) +
    if(type=="Temporal (interval)"){ #for multi panel plots. Only keep legend for the last plot.
      theme(legend.position = c(0.5,0.87), legend.text = element_text(size = 18), legend.background = element_rect(fill = "white", color = NA), strip.text = element_text(size = 22, margin = margin(t = -1.8)), title = element_blank(), axis.title.x = element_text(size=24,margin = margin(t = 25)), axis.title.y = element_text(size=22, margin = margin(r = 10)), axis.text.y=element_blank() , axis.text.x = element_blank(), panel.grid.major.x=element_blank(), panel.spacing = unit(0, "lines")) 
    }else{
      if(type=="Temporal (era)"){
        theme(legend.position = "none",legend.background = element_rect(fill = "white", color = NA), strip.text = element_text(size = 22, margin = margin(t = -1.8)), title = element_blank(), axis.title.x = element_text(size=24,margin = margin(t = 25)), axis.title.y = element_text(size=22, margin = margin(r = 10)), axis.text.y=element_blank(), axis.text.x = element_blank() , panel.grid.major.x=element_blank(), panel.spacing = unit(0, "lines")) 
      }else{
        theme(legend.position = "none",legend.background = element_rect(fill = "white", color = NA), strip.text = element_text(size = 20, margin = margin(t = -1.8)), title = element_blank(), axis.title.x = element_text(size=24,margin = margin(t = 25)), axis.title.y = element_text(size=22, margin = margin(r = 10)), axis.text.y=element_text(size=24), axis.text.x = element_blank() ,panel.grid.major.x=element_blank(), panel.spacing = unit(0, "lines")) 
      }}}
#plot.Sensitivity.TotalDeclinesGrowths("Spatial")
#plot.Sensitivity.TotalDeclinesGrowths("Temporal (era)")
#plot.Sensitivity.TotalDeclinesGrowths("Temporal (interval)")
#grid.arrange(plot.Sensitivity.TotalDeclinesGrowths("Spatial"),plot.Sensitivity.TotalDeclinesGrowths("Temporal (era)"),plot.Sensitivity.TotalDeclinesGrowths("Temporal (interval)"),ncol=3)

# Function to compare model estimates of baseline occupancy:
plot.Sensitivity.mu.psi <- function(type){ #Function that takes either "Spatial", "Temporal (era)", or "Temporal (interval)" and plots the associated sensitivity analysis
  ifelse(type=="Temporal (era)",load("SensitivityAnalysis/6eras/ModelSummary.RData"),ifelse(type=="Temporal (interval)",load("SensitivityAnalysis/10intervals/ModelSummary.RData"),load("SensitivityAnalysis/50x50km/ModelSummary.RData"))) #Load in data from the Sensitivity Analysis folder
  ModelSummary.high <- mutate(as.data.frame(ModelSummary), Resolution=ifelse(type=="Temporal (era)","20",ifelse(type=="Temporal (interval)","3","50x50")))
  load("SensitivityAnalysis/BaseModel/ModelSummary.RData")
  ModelSummary.medium <- mutate(as.data.frame(ModelSummary), Resolution=ifelse(type=="Temporal (era)","30",ifelse(type=="Temporal (interval)","5","100x100")))
  ifelse(type=="Temporal (era)",load("SensitivityAnalysis/2eras/ModelSummary.RData"),ifelse(type=="Temporal (interval)",load("SensitivityAnalysis/3intervals/ModelSummary.RData"),load("SensitivityAnalysis/200x200km/ModelSummary.RData")))
  ModelSummary.low <- mutate(as.data.frame(ModelSummary), Resolution=ifelse(type=="Temporal (era)","60",ifelse(type=="Temporal (interval)","10","200x200")))
  params.to.plot <- as.data.frame(rbind(ModelSummary.high[1,],ModelSummary.medium[1,],ModelSummary.low[1,])) #Choose parameters to plot
  params.to.plot$Resolution <- factor(params.to.plot$Resolution, levels=unique(params.to.plot$Resolution),ordered=T) #Do this to prevent ggplot from automatically sorting the x axis alphabetically
  ytitle <- ifelse(type=="Spatial","Baseline Occupancy Parameter (mu.psi)","") #Only give y axis title to first plot when making multi-panel plot
  ggplot(params.to.plot, aes(x= Resolution, y = Mean)) +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0.25, position = position_dodge(width = 0), linewidth=1.5, color="grey60") +
    geom_point(position = position_dodge(width = 0), size=3.5, color="black") +
    #labs(title = paste("Baseline Occupancy at Various ",type," Resolutions",sep=""), subtitle= "(black dot = mean, grey bars = 95% credible interval)", x = ifelse(type=="Temporal (era)","Era Size",ifelse(type=="Temporal (interval)","Interval Size","Grid Cell Size")), y = "Baseline Occupancy (mu.psi)") +
    labs(x = ifelse(type=="Temporal (era)","Era Size (yrs)",ifelse(type=="Temporal (interval)","Interval Size (yrs)","Grid Cell Size (km)")), y = ytitle) +
    ylim(0.15,2) +
    guides(color=F) +
    theme_minimal() +
    if(type=="Spatial"){
      theme(panel.grid.major.x=element_blank(), axis.title=element_text(size=22), axis.text = element_text(size=16), axis.title.y = element_text(margin = margin(r = 15)), axis.title.x = element_text(margin = margin(t = 15)))
    }else{
      theme(panel.grid.major.x=element_blank(), axis.title=element_text(size=22), axis.text = element_text(size=18), axis.text.y = element_blank(), axis.title.x = element_text(margin = margin(t = 15)))
    }}
#plot.Sensitivity.mu.psi("Spatial")
#plot.Sensitivity.mu.psi("Temporal (era)")
#plot.Sensitivity.mu.psi("Temporal (interval)")
#grid.arrange(plot.Sensitivity.mu.psi("Spatial"),plot.Sensitivity.mu.psi("Temporal (era)"),plot.Sensitivity.mu.psi("Temporal (interval)"),ncol=3)













############################
#### Trend Correlations ####
############################
# ANOVA was done to test effect of range size
# t-test was done to test effect of body size
# Chi square tests were done to test association between recent and historic trends, as well as effect of subfamily and larval guild

#### Add Recent Trend and Historic Trend columns to the Species Table ####
Species_table$RecentTrend <- NA #create empty columns to populate
Species_table$HistoricTrend <- NA
# Recent Trend status
for (i in 1:nrow(Species_table)){ 
  if(Species_table$genusspecies[i] %in% RecentDecline_sp){
    Species_table$RecentTrend[i] <- "Decline"
  }else{
    if(Species_table$genusspecies[i] %in% RecentGrowth_sp){
      Species_table$RecentTrend[i] <- "Increase"
    }else{
      Species_table$RecentTrend[i] <- "Stable"
    }}}
# Historical Trend status
for (i in 1:nrow(Species_table)){ 
  if(Species_table$genusspecies[i] %in% HistoricDecline_sp){
    Species_table$HistoricTrend[i] <- "Decline"
  }else{
    if(Species_table$genusspecies[i] %in% HistoricGrowth_sp){
      Species_table$HistoricTrend[i] <- "Increase"
    }else{
      Species_table$HistoricTrend[i] <- "Stable"
    }}}
Species_table$RecentTrend <- as.factor(Species_table$RecentTrend)
Species_table$HistoricTrend <- as.factor(Species_table$HistoricTrend)



#### Is range size correlated with Recent Decline/Growth? One way ANOVA ####
# Test assumptions:
#  -the residuals are normally distributed
#  -the error variance is the same for all groups (homoscedasticity)
#  -the residuals are independent.
aov.range.recent <- aov(range.size ~ RecentTrend, Species_table)
#plot(aov.range.recent)
#shapiro.test(residuals(aov.range.recent)) #Residuals are not normal
#leveneTest(range.size ~ RecentTrend, Species_table) #There is homoscedasticity
# There is non-normality and homoscedasticity; therefore we must use Kruskal-Wallis, which is non-parametric and does not rely on assumptions.
kruskal.test(range.size ~ RecentTrend, Species_table) #range size has significant effect on RecentTrend
# Plot Tukeys post-hoc
RangeSizeRecent <- glht(aov.range.recent, linfct = mcp(RecentTrend = "Tukey"))
RangeSizeRecent <- cld(RangeSizeRecent)
RangeSizeRecent$xname <- "Recent Trend (between 1961\u20131990 and 1991\u20132020)"
RangeSizeRecent$yname <- "Range Size (number of grid cells)"
#plot(RangeSizeRecent)
# Cannot use Tukey's post hoc after a non-parametric test. Use Dunn's Test instead. Check if results are identical to Tukey's, so we may still use Tukey's results for easy plotting.
dunnTest(range.size ~ RecentTrend, Species_table, method = "bonferroni")
# Species with small range sizes are more likely to have no recent trend. Perhaps this reflects a lack of data on species with small range sizes, which results in large confidence intervals and thus less likely to have a significant trend.

# Same analysis but with Historic Trends:
aov.range.historic <- aov(range.size ~ HistoricTrend, Species_table)
#plot(aov.range.historic)
#shapiro.test(residuals(aov.range.historic)) #Residuals are not normal
#leveneTest(range.size ~ HistoricTrend, Species_table) #There is homoscedasticity
kruskal.test(range.size ~ HistoricTrend, Species_table) #range size has significant effect on HistoricTrend
#plot(TukeyHSD(aov.range.historic))
RangeSizeHistoric <- glht(aov.range.historic, linfct = mcp(HistoricTrend = "Tukey"))
RangeSizeHistoric <- cld(RangeSizeHistoric)
RangeSizeHistoric$xname <- "Historical Trend (between 1900\u20131930 and 1991\u20132020)"
RangeSizeHistoric$yname <- "Range Size (number of grid cells)"
#plot(RangeSizeHistoric)
# Cannot use Tukey's post hoc after a non-parametric test. Use Dunn's Test instead. Check if results are identical to Tukey's, so we may still use Tukey's results for easy plotting.
dunnTest(range.size ~ HistoricTrend, Species_table, method = "bonferroni")
# Species with small range sizes are also more likely to have no historic trend.



#### Is subfamily associated with Recent Decline/Growth? Chi square test ####
# 80% of the cells must have an expected count greater than 5 in order to use a Chi square test
# Contingency table of recent trend and subfamily
X2.subfam.recent <- Species_table %>%
  group_by(RecentTrend,Subfamily) %>% #group by Recent Trend and Subfamily
  summarise(count = n()) %>% 
  arrange(.by_group=TRUE) %>%
  ungroup() #Ungroup
X2.subfam.recent$Subfamily <- factor(X2.subfam.recent$Subfamily, levels = names(sort(table(X2.subfam.recent$Subfamily), decreasing = TRUE))) #Treat subfamily as an ordered factor variable, ordered by frequency
X2.subfam.recent <- t(xtabs(count ~ RecentTrend + Subfamily, X2.subfam.recent)) #Convert to table format to account for all combos of Trend + Subfamily, and transpose
# Check expected values
#chisq.test(X2.subfam.recent)$expected #Only 67% of expected values are greater than 5.
# Low expected counts; therefore, use Fisher Exact Test
fisher.test(X2.subfam.recent)
# Mosaic diagram
#mosaic(X2.subfam.recent, shade=T, legend=F, rot_labels=c(0,90,0,0), just_labels=c("center","right"), offset_varnames = c(1,6), margins = c(4, 5, 4, 4), set_varnames=c(RecentTrend="Recent Trend (between 1961\u20131990 and 1991\u20132020)")) #Red and Blue shades represent deviations from the expected values if the categories were truly independent. Blue means that group has more observations than expected. Red means fewer than expected. Shading (and the presented p value) is based on a Chi square test, not Fisher's exact test, so be wary.
# Syrphines are significantly more likely, and Eristalines are significantly less likely, to be in recent decline than would be expected under independence.

# Same analysis but with Historic Trends:
X2.subfam.historic <- Species_table %>%
  group_by(HistoricTrend,Subfamily) %>% 
  summarise(count = n()) %>% 
  arrange(.by_group=TRUE) %>%
  ungroup() #Ungroup
X2.subfam.historic$Subfamily <- factor(X2.subfam.historic$Subfamily, levels = names(sort(table(X2.subfam.historic$Subfamily), decreasing = TRUE))) #Treat subfamily as an ordered factor variable, ordered by frequency
X2.subfam.historic <- t(xtabs(count ~ HistoricTrend + Subfamily, X2.subfam.historic))
#chisq.test(X2.subfam.historic)$expected #Only 67% of expected values are greater than 5.
# Low expected counts; therefore, use Fisher Exact Test
fisher.test(X2.subfam.historic)
#mosaic(X2.subfam.historic, shade=T, legend=F, rot_labels=c(0,90,0,0), just_labels=c("center","right"), offset_varnames = c(1,6), margins = c(4, 5, 4, 4), set_varnames=c(HistoricTrend="Historic Trend (between 1900\u20131930 and 1991\u20132020)")) #Red and Blue shades represent deviations from the expected values if the categories were truly independent. Blue means that group has more observations than expected. Red means fewer than expected. Shading (and the presented p value) is based on a Chi square test, not Fisher's exact test, so be wary.
# Subfamily and Historic Trends are not significantly associated.



#### Is larval guild associated with Recent Decline/Growth? Chi square test ####
# Contingency table of recent trend and larval guild
X2.larval.recent <- Species_table %>%
  group_by(RecentTrend,LarvalGuild) %>% 
  summarise(count = n()) %>% 
  arrange(.by_group=TRUE) %>%
  ungroup() #Ungroup
X2.larval.recent$LarvalGuild <- factor(X2.larval.recent$LarvalGuild, levels = c("Saprophagous","Predator","BroodParasite","Phytophagous","Unknown")) #Treat LarvalGuild as an ordered factor variable, ordered by frequency
X2.larval.recent <- t(xtabs(count ~ RecentTrend + LarvalGuild, X2.larval.recent))
#chisq.test(X2.larval.recent)$expected #Only 60% of expected values are greater than 5.
# Low expected counts; therefore, use Fisher Exact Test
fisher.test(X2.larval.recent)
#mosaic(X2.larval.recent, shade=T, legend=F, rot_labels=c(0,90,0,0), just_labels=c("center","right"), offset_varnames = c(1,6), margins = c(4, 5, 4, 4), set_varnames=c(RecentTrend="Recent Trend (between 1961\u20131990 and 1991\u20132020)")) #Red and Blue shades represent deviations from the expected values if the categories were truly independent. Blue means that group has more observations than expected. Red means fewer than expected. Shading (and the presented p value) is based on a Chi square test, not Fisher's exact test, so be wary.
# Predators are significantly more likely, and Saprophages less likely, to be in recent decline than would be expected. This is linked to subfamily traits.

# Same analysis but with Historic Trends:
X2.larval.historic <- Species_table %>%
  group_by(HistoricTrend,LarvalGuild) %>% 
  summarise(count = n()) %>% 
  arrange(.by_group=TRUE) %>%
  ungroup() #Ungroup
X2.larval.historic$LarvalGuild <- factor(X2.larval.historic$LarvalGuild, levels = c("Saprophagous","Predator","BroodParasite","Phytophagous","Unknown")) #Treat LarvalGuild as an ordered factor variable, ordered by frequency
X2.larval.historic <- t(xtabs(count ~ HistoricTrend + LarvalGuild, X2.larval.historic))
#chisq.test(X2.larval.historic)$expected #Only 60% of expected values are greater than 5.
# Low expected counts; therefore, use Fisher Exact Test
fisher.test(X2.larval.historic)
#mosaic(X2.larval.historic, shade=T, legend=F, rot_labels=c(0,90,0,0), just_labels=c("center","right"), offset_varnames = c(1,6), margins = c(4, 5, 4, 4), set_varnames=c(HistoricTrend="Historic Trend (between 1900\u20131930 and 1991\u20132020)")) #Red and Blue shades represent deviations from the expected values if the categories were truly independent. Blue means that group has more observations than expected. Red means fewer than expected. Shading (and the presented p value) is based on a Chi square test, not Fisher's exact test, so be wary.
# Larval Guild and Historic Trends are not significantly associated.



#### Is being a rot-hole specialist associated with recent decline/growth? Chi Square Test ####
# Add column to species table denoting whether a species is a rot hole specialist
Species_table$rothole <- NA
for (i in 1:nrow(Species_table)){
  if(Species_table$Genus[i] %in% c("Blera","Criorhina","Mallota","Meromacrus","Milesia","Myolepta","Spilomyia")){
    Species_table$rothole[i] <- "Rothole Specialist"
  }else{
    Species_table$rothole[i] <- "Other"}}
# Contingency table of recent trend and rot hole
X2.rothole.recent <- Species_table %>%
  group_by(RecentTrend,rothole) %>% 
  summarise(count = n()) %>% 
  arrange(.by_group=TRUE) %>%
  ungroup() #Ungroup
X2.rothole.recent$rothole <- factor(X2.rothole.recent$rothole, levels = c("Rothole Specialist","Other"))
X2.rothole.recent <- t(xtabs(count ~ RecentTrend + rothole, X2.rothole.recent))
#chisq.test(X2.rothole.recent)$expected
# Low expected counts; therefore, use Fisher Exact Test
fisher.test(X2.rothole.recent)
#mosaic(X2.rothole.recent, shade=T, legend=F, rot_labels=c(0,90,0,0), just_labels=c("center","right"), offset_varnames = c(1,6), margins = c(4, 5, 4, 4), set_varnames=c(RecentTrend="Recent Trend (between 1961\u20131990 and 1991\u20132020)")) #Red and Blue shades represent deviations from the expected values if the categories were truly independent. Blue means that group has more observations than expected. Red means fewer than expected. Shading (and the presented p value) is based on a Chi square test, not Fisher's exact test, so be wary.
# Rothole specialization and Recent Trend are not associated.

# Same analysis but with Historic Trends:
X2.rothole.historic <- Species_table %>%
  group_by(HistoricTrend,rothole) %>% 
  summarise(count = n()) %>% 
  arrange(.by_group=TRUE) %>%
  ungroup() #Ungroup
X2.rothole.historic$rothole <- factor(X2.rothole.historic$rothole, levels = c("Rothole Specialist","Other"))
X2.rothole.historic <- t(xtabs(count ~ HistoricTrend + rothole, X2.rothole.historic))
#chisq.test(X2.rothole.historic)$expected
# Low expected counts; therefore, use Fisher Exact Test
fisher.test(X2.rothole.historic)
#mosaic(X2.rothole.historic, shade=T, legend=F, rot_labels=c(0,90,0,0), just_labels=c("center","right"), offset_varnames = c(1,6), margins = c(4, 5, 4, 4), set_varnames=c(HistoricTrend="Historic Trend (between 1900\u20131930 and 1991\u20132020)")) #Red and Blue shades represent deviations from the expected values if the categories were truly independent. Blue means that group has more observations than expected. Red means fewer than expected. Shading (and the presented p value) is based on a Chi square test, not Fisher's exact test, so be wary.
# Rothole specialization and Historic Trend are not associated.



#### Are Recent Trends associated with Historic Trends? Chi square test ####
# Contingency table of recent and historic trend
X2.recent.historic <- Species_table %>%
  group_by(RecentTrend,HistoricTrend) %>% 
  summarise(count = n()) %>% 
  arrange(.by_group=TRUE) %>%
  ungroup() #Ungroup
X2.recent.historic <- t(xtabs(count ~ RecentTrend + HistoricTrend, X2.recent.historic))
#chisq.test(X2.recent.historic)$expected #Only 78% of expected values are greater than 5.
# Low expected counts; therefore, use Fisher Exact Test
fisher.test(X2.recent.historic)
#mosaic(X2.recent.historic, shade=T, gp=shading_hcl, gp_args = list(interpolate = c(2, 10)), legend=F, rot_labels=c(90,90,0,0), just_labels=c("left","right"), offset_varnames = c(3.5,3.5), margins = c(4, 10, 4, 4), set_varnames=c(HistoricTrend="Historic Trend (between 1900\u20131930 and 1991\u20132020)",RecentTrend="Recent Trend (between 1961\u20131990 and 1991\u20132020)")) #Red and Blue shades represent deviations from the expected values if the categories were truly independent. Blue means that group has more observations than expected. Red means fewer than expected. Shading (and the presented p value) is based on a Chi square test, not Fisher's exact test, so be wary. Custom shading cutoff is used to remove the different colour intensities for varying levels of deviation from expected value.
# Recent declines are significantly more likely to also be historic declines. Recent increases are significantly more likely to also be historic increases. Species with no recent trend are more likely to also have no historical trend.
# In summary, recent trends are significantly more likely to match historic trends than would be expected under independence. Species trends generally agree at both timescales.



#### Is Body Size associated with trends? t-test ####
# Create table of body sizes (for only species with significant trends) taken from Skevington et al. (2019) and BugGuide.com
AllDeclines <- union(RecentDecline_sp, HistoricDecline_sp) 
AllGrowth <- union(RecentGrowth_sp, HistoricGrowth_sp)
AllDeclinesAndGrowth <- union(AllDeclines, AllGrowth) #Find all species with significant trends
BodySize_table <- matrix(ncol=2, byrow=TRUE, dimnames = list(NULL,c("genusspecies","BodySize")),
       c("Anasimyia anausis", 8.85,
         "Anasimyia chrysostoma", 9.4,
         "Cheilosia orilliaensis (consentiens)", 8.8,
         "Cheilosia prima", 9.1,
         "Chrysogaster antitheus (ontario)", 6.0,
         "Chrysotoxum chinook", "NA",
         "Chrysotoxum coloradense (fasciatum)", "NA",
         "Chrysotoxum derivatum (columbianum,integre,minor,ventricosum)", 9.1,
         "Chrysotoxum flavifrons", 13.85,
         "Chrysotoxum plumeum (perplexum)", 10.3,                         
         "Chrysotoxum pubescens (luteopilosum)", 12.65,
         "Criorhina caudata",     15,
         "Criorhina verbosa", 15.5,
         "Dasysyrphus creper",    8.6,
         "Dasysyrphus intrudens", 9.35,
         "Dasysyrphus limatus", 9.2,                                    
         "Dasysyrphus pauxillus",    7.5,
         "Eristalis arbustorum", 10.15,                                 
         "Eristalis brousii", 11,
         "Eristalis fraterculus", 24.4,                                         
         "Eristalis hirta", 11.2,
         "Helophilus groenlandicus", 12.3,
         "Helophilus latifrons", 13.45,
         "Lejota cyanea", 7.6,                                                 
         "Melangyna arctica", 7.35,
         "Melangyna lasiophthalma", 8.5,                                       
         "Melangyna umbellatarum", 9.6,
         "Meligramma (Melangyna) triangulifera", 8.5,                          
         "Neoascia globosa (distincta)", 4.95,
         "Neoascia sandsi", "NA",                             
         "Neoascia tenur", 5,
         "Neocnemodon elongata", 6.8,                                         
         "Paragus angustifrons", 5.45,
         "Parasyrphus relictus", 8.4,                                          
         "Parasyrphus tarsatus", 8.95,
         "Parasyrphus vockerothi", 7.25,                                        
         "Pipiza femoralis", 8, 
         "Pipiza puella (nigrotibiata)", 9.25,                                  
         "Platycheirus albimanus", 7.9,
         "Platycheirus angustatus", 6.8,                                 
         "Platycheirus chilosia", 6.2,
         "Platycheirus clypeatus", 7.4,                                        
         "Platycheirus coerulescens", 7.4,
         "Platycheirus confusus", 7.1,                                         
         "Platycheirus hyperboreus (erraticus)", 7,
         "Platycheirus immarginatus (felix)", 8.05,                             
         "Platycheirus luteipennis", 9.4,
         "Platycheirus modestus", 6.75,                                         
         "Platycheirus naso (holarcticus)", 8.6,
         "Platycheirus nodosus", 7,                                          
         "Platycheirus pictipes (concinnus,rufimaculatus)", 9.1,
         "Platycheirus scutatus", 7.75,                                         
         "Platycheirus stegnus",      7.95,
         "Platycheirus striatus", 9.15,                                         
         "Pyrophaena (Platycheirus) granditarsis", 9.1,
         "Sphaerophoria abbreviata", 8.25,                                      
         "Sphaerophoria asymmetrica", 7.9,
         "Sphaerophoria contigua (cylindricus)", 7.65,                          
         "Sphaerophoria philanthus (dubia,nigritarsi,robusta)", 8.6,
         "Sphaerophoria pyrrhina (guttulata)", 6.25,                            
         "Sphecomyia brevicornis",    13.5,
         "Sphecomyia dyari",       12.8,                                       
         "Sphegina flavimana", 5.65,
         "Syrphus torvus", 10.95,                                                
         "Temnostoma balyras", 10,
         "Volucella facialis (lateralis,rufomaculata)", 12,                  
         "Xylota flavitibia", 12.75,
         "Brachypalpus femoratus", 13.6,                                        
         "Scaeva affinis (pyrastri)", 13.4,
         "Volucella evecta",  15.25,                                            
         "Allograpta obliqua", 7.2,
         "Epistrophe grossulariae", 12.7,                                       
         "Eristalis flavipes", 13.55,
         "Eristalis tenax", 13.75,                                              
         "Eristalis transversa", 11.15,
         "Eupeodes americanus", 9,                                          
         "Eupeodes fumipennis",     10,
         "Helophilus fasciatus", 13,                                         
         "Hypocritanus (Ocyptamus) fascipennis", 11.15,
         "Polydontomyia (Lejops) curvipes", 14.75,                               
         "Sericomyia militaris", 14.4,
         "Spilomyia fusca", 16.35,                                              
         "Spilomyia sayi", 13.5,
         "Toxomerus geminatus", 6.85,                                          
         "Toxomerus marginatus", 5.3,
         "Baccha elongata",       8.7,                                        
         "Chalcosyrphus (Xylota) libo", 9.5,
         "Copestylum vesicularium", 8.9,                                      
         "Eristalinus aeneus", 10.05,
         "Mallota (Imatisma) posticata (bardus,separata)",  15.25,              
         "Melanostoma mellinum (angustatum,fallax,melanderi,pallitarsis)", 7.4,
         "Meliscaeva cinctella",  9.4,
         "Merodon equestris", 14.75,
         "Milesia virginiensis",  20.05,                                        
         "Ocyptamus fuscipennis", 9.05,
         "Orthonevra pulchella",  5.95,                                    
         "Paragus haemorrhous", 5.1,
         "Parhelophilus rex", 9.7,                                      
         "Platycheirus obscurus", 8.1,
         "Pyrophaena (Platycheirus) rosarum", 8.2,                            
         "Sericomyia chrysotoxoides", 12.45,
         "Sericomyia lata", 13.4,                                   
         "Sericomyia transversa", 12.4,
         "Syritta pipiens", 8
       )) %>% data.frame()
BodySize_table$BodySize <- as.numeric(BodySize_table$BodySize)
Species_table <- right_join(BodySize_table, Species_table, by="genusspecies") %>% relocate("BodySize",.after="rothole") #add bodysize column to species table
RecentDeclineSizes <- filter(Species_table, BodySize != "NA" & RecentTrend == "Decline")$BodySize
RecentIncreaseSizes <- filter(Species_table, BodySize != "NA" & RecentTrend == "Increase")$BodySize
HistoricalDeclineSizes <- filter(Species_table, BodySize != "NA" & HistoricTrend == "Decline")$BodySize
HistoricalIncreaseSizes <- filter(Species_table, BodySize != "NA" & HistoricTrend == "Increase")$BodySize
# Check for t-test assumptions of normality and equal variance (homoscedasticity)
shapiro.test(RecentDeclineSizes)
shapiro.test(RecentIncreaseSizes)
shapiro.test(HistoricalDeclineSizes)
shapiro.test(HistoricalIncreaseSizes)
# Some groups are non-normal, so we must use a non-parametric test.
# Recent Trend
wilcox.test(RecentDeclineSizes, RecentIncreaseSizes)
boxplot(RecentDeclineSizes, RecentIncreaseSizes)
# At recent baseline, species that are increasing generally have larger body sizes!
# Historical Trend
wilcox.test(HistoricalDeclineSizes, HistoricalIncreaseSizes)
boxplot(HistoricalDeclineSizes, HistoricalIncreaseSizes)
# No relationship between body size and historical trend.



#### Are introduced species more likely to be increasing/decreasing? ####
# Create table of introduced species. "Introduced" only means there is evidence or historical knowledge of the species having been introduced to North America (typically within the past 200 years). Other species may simply be considered holarctic or cosmopolitan in distribution (ex. Melanostoma mellinum).
Introduced_table <- matrix(ncol=2, byrow=TRUE, dimnames = list(NULL,c("genusspecies","Introduced")),
                         c("Anasimyia anausis", "no",
                           "Anasimyia chrysostoma", "no",
                           "Cheilosia orilliaensis (consentiens)", "no",
                           "Cheilosia prima", "no",
                           "Chrysogaster antitheus (ontario)", "no",
                           "Chrysotoxum chinook", "NA",
                           "Chrysotoxum coloradense (fasciatum)", "NA",
                           "Chrysotoxum derivatum (columbianum,integre,minor,ventricosum)", "no",
                           "Chrysotoxum flavifrons", "no",
                           "Chrysotoxum plumeum (perplexum)", "no",                         
                           "Chrysotoxum pubescens (luteopilosum)", "no",
                           "Criorhina caudata", "no",
                           "Criorhina verbosa", "no",
                           "Dasysyrphus creper",    "no",
                           "Dasysyrphus intrudens", "no",
                           "Dasysyrphus limatus", "no",                                    
                           "Dasysyrphus pauxillus",   "no",
                           "Eristalis arbustorum", "yes",                                 
                           "Eristalis brousii", "no",
                           "Eristalis fraterculus", "no",                                         
                           "Eristalis hirta", "no",
                           "Helophilus groenlandicus", "no",
                           "Helophilus latifrons", "no",
                           "Lejota cyanea", "no",                                                 
                           "Melangyna arctica", "no",
                           "Melangyna lasiophthalma", "no",                                       
                           "Melangyna umbellatarum", "no",
                           "Meligramma (Melangyna) triangulifera", "no",                          
                           "Neoascia globosa (distincta)", "no",
                           "Neoascia sandsi", "NA",                             
                           "Neoascia tenur", "no",
                           "Neocnemodon elongata", "no",                                         
                           "Paragus angustifrons", "no",
                           "Parasyrphus relictus", "no",                                          
                           "Parasyrphus tarsatus", "no",
                           "Parasyrphus vockerothi", "no",                                        
                           "Pipiza femoralis", "no", 
                           "Pipiza puella (nigrotibiata)", "no",                                  
                           "Platycheirus albimanus", "no",
                           "Platycheirus angustatus", "no",                                 
                           "Platycheirus chilosia", "no",
                           "Platycheirus clypeatus", "no",                                        
                           "Platycheirus coerulescens", "no",
                           "Platycheirus confusus", "no",                                         
                           "Platycheirus hyperboreus (erraticus)", "no",
                           "Platycheirus immarginatus (felix)", "no",                             
                           "Platycheirus luteipennis", "no",
                           "Platycheirus modestus", "no",                                         
                           "Platycheirus naso (holarcticus)", "no",
                           "Platycheirus nodosus", "no",                                          
                           "Platycheirus pictipes (concinnus,rufimaculatus)", "no",
                           "Platycheirus scutatus", "no",                                         
                           "Platycheirus stegnus",  "no",
                           "Platycheirus striatus", "no",                                         
                           "Pyrophaena (Platycheirus) granditarsis", "no",
                           "Sphaerophoria abbreviata", "no",                                      
                           "Sphaerophoria asymmetrica", "no",
                           "Sphaerophoria contigua (cylindricus)", "no",                          
                           "Sphaerophoria philanthus (dubia,nigritarsi,robusta)", "no",
                           "Sphaerophoria pyrrhina (guttulata)", "no",                            
                           "Sphecomyia brevicornis", "no",
                           "Sphecomyia dyari", "no",                                       
                           "Sphegina flavimana", "no",
                           "Syrphus torvus", "no",                                                
                           "Temnostoma balyras", "no",
                           "Volucella facialis (lateralis,rufomaculata)", "no",                  
                           "Xylota flavitibia", "no",
                           "Brachypalpus femoratus", "no",                                        
                           "Scaeva affinis (pyrastri)", "no",
                           "Volucella evecta",  "no",                                            
                           "Allograpta obliqua", "no",
                           "Epistrophe grossulariae", "no",                                       
                           "Eristalis flavipes", "no",
                           "Eristalis tenax", "yes",                                              
                           "Eristalis transversa", "no",
                           "Eupeodes americanus", "no",                                          
                           "Eupeodes fumipennis", "no",
                           "Helophilus fasciatus", "no",                                         
                           "Hypocritanus (Ocyptamus) fascipennis", "no",
                           "Polydontomyia (Lejops) curvipes", "no",                               
                           "Sericomyia militaris", "no",
                           "Spilomyia fusca", "no",                                              
                           "Spilomyia sayi", "no",
                           "Toxomerus geminatus", "no",                                          
                           "Toxomerus marginatus", "no",
                           "Baccha elongata",   "no",                                        
                           "Chalcosyrphus (Xylota) libo", "no",
                           "Copestylum vesicularium", "no",                                      
                           "Eristalinus aeneus", "yes",
                           "Mallota (Imatisma) posticata (bardus,separata)",  "no",              
                           "Melanostoma mellinum (angustatum,fallax,melanderi,pallitarsis)", "no",
                           "Meliscaeva cinctella",  "no",
                           "Merodon equestris", "yes",
                           "Milesia virginiensis",  "no",                                        
                           "Ocyptamus fuscipennis", "no",
                           "Orthonevra pulchella", "no",                                    
                           "Paragus haemorrhous", "no",
                           "Parhelophilus rex", "no",                                      
                           "Platycheirus obscurus", "no",
                           "Pyrophaena (Platycheirus) rosarum", "no",                            
                           "Sericomyia chrysotoxoides", "no",
                           "Sericomyia lata", "no",                                   
                           "Sericomyia transversa", "no",
                           "Syritta pipiens", "yes",
                           "Eumerus funeralis", "yes",  #Add this stable but introduced species
                           "Eumerus strigatus", "yes" #Add this stable but introduced species
                         )) %>% data.frame()
Species_table <- right_join(Introduced_table, Species_table, by="genusspecies") %>% relocate("Introduced",.after="BodySize") #add Introduced column to species table
# Of the 5 introduced species with strongly supported changes in occupancy, 4 were increasing and 1 was in decline. However, at least 2 introduced species exhibited stable occupancy.



#Save the updated Species Table
#write.csv(Species_table, file="FinalModel/Species_table.csv")










#################
#### Figures ####
#################
par(mfrow=c(1,1)) #Set default plotting grid

# Pie and Bar graphs of Collection sizes
slices <- c(length(which(SyrphidDataset$Collection=="CNC")), length(which(SyrphidDataset$Collection=="LACM")), length(which(SyrphidDataset$Collection=="Guelph")), length(which(SyrphidDataset$Collection=="Spencer")), length(which(SyrphidDataset$Collection=="Alaska")), length(which(SyrphidDataset$Collection=="Strickland")), length(which(SyrphidDataset$Collection=="Texas")))
pie(slices, labels=NA, col=c("#F50E0E70","#E850FF70","#36CC0A70","#3129F970","#00F1FF70","#E6EA1670","#30B1F999")) #no labels
legend('right', legend = factor(c("CNC","LACM","UGuelph","Spencer","UAlaska","Strickland","UTexas")), col = c("#F50E0E70","#E850FF70","#36CC0A70","#3129F970","#00F1FF70","#E6EA1670","#30B1F999"), pch = 16, cex = 1.2, pt.cex = 2, inset = c(-0.4, 0.35), y.intersp = 1.1, x.intersp =0.8, bty='n')



# Map of records. Save as .png file in working directory
png(filename="RecordsMap.png",width=2500,height=2500, res=200)
plot(grid$geometry, border="lightblue"); plot(m.states$geometry, add=T); 
plot(SyrphidDataset$geometry[which(SyrphidDataset$Collection=="CNC")], add=T, cex=0.4, pch=20, col="#F50E0E30"); #CNC is red
plot(SyrphidDataset$geometry[which(SyrphidDataset$Collection=="Guelph")], add=T, cex=0.4, pch=20, col="#36CC0A30"); #Guelph is green
plot(SyrphidDataset$geometry[which(SyrphidDataset$Collection=="Strickland")], add=T, cex=0.4, pch=20, col="#E6EA1630"); #Strickland is yellow
plot(SyrphidDataset$geometry[which(SyrphidDataset$Collection=="Spencer")], add=T, cex=0.4, pch=20, col="#3129F930"); #Spencer is dark blue
plot(SyrphidDataset$geometry[which(SyrphidDataset$Collection=="Alaska")], add=T, cex=0.4, pch=20, col="#00F1FF30"); #Alaska is cyan
plot(SyrphidDataset$geometry[which(SyrphidDataset$Collection=="Texas")], add=T, cex=0.4, pch=20, col="#30B1F930"); #Texas is light blue
plot(SyrphidDataset$geometry[which(SyrphidDataset$Collection=="LACM")], add=T, cex=0.4, pch=20, col="#E850FF30"); #LACM is purple
legend('bottomleft', legend = factor(c("CNC","LACM","UGuelph","Spencer","UAlaska","Strickland","UTexas")), col = c("#F50E0E","#E850FF","#36CC0A","#3129F9","#00F1FF","#E6EA16","#30B1F9"), pch = 16, cex = 1.6, inset = c(0.1, 0.35), y.intersp = 1.2, x.intersp =0.8, bty='n')
#title(paste("Syrphidae (n = ",nrow(SyrphidDataset),")"),cex.main=3)
dev.off()


# Histogram of records per year
#ggplot(data=SyrphidDataset, aes(x=Year)) + geom_histogram(aes(y=..density..), binwidth=1, closed="left", color="black", alpha=0.3) + geom_density() + geom_rug() + scale_x_continuous(breaks = seq(1900,2020,10), labels = seq(1900,2020,10)) + theme(axis.text.x=element_text(angle=30)) + ggtitle("Density of syrphid records per year")
# Bar graph of records per year
ggplot(data=SyrphidDataset, aes(x=Year)) + 
  geom_bar(width=1, fill="grey80",color="grey20", alpha=1) + 
  theme(axis.text.x=element_text(angle=0), panel.background=element_blank(),axis.line=element_line(color="black"),panel.grid.major.x = element_blank(),axis.title = element_text(size = 20), axis.text=element_text(size=18),plot.margin = margin(10, 20, 10, 10)) + 
  scale_x_discrete(limits=seq(1900,2020,30)) + 
  xlab("Year") + 
  ylab("Number of Syrphid Records") + 
  scale_y_continuous(expand = c(0, 0))
# Specimens per era
#length(which(SyrphidDataset$Era=="1900-1930"))
#length(which(SyrphidDataset$Era=="1931-1960"))
#length(which(SyrphidDataset$Era=="1961-1990"))
#length(which(SyrphidDataset$Era=="1991-2020"))
# Histogram of grids per year
GridID.per.year <- group_by(SyrphidDataset, Year) %>% summarise(grids=n_distinct(grid.id)) %>% arrange(.by_group=TRUE) %>% ungroup() #Ungroup
GridID.per.era <- group_by(SyrphidDataset, Era) %>% summarise(grids=n_distinct(grid.id)) %>% arrange(.by_group=TRUE) %>% ungroup() #Ungroup
ggplot(data=GridID.per.year, aes(x=Year, y=grids)) + 
  geom_bar(stat="identity",width=1, fill="grey80",color="grey20", alpha=1) + 
  theme(axis.text.x=element_text(angle=0), panel.background=element_blank(),axis.line=element_line(color="black"),panel.grid.major.x = element_blank(),axis.title = element_text(size = 20), axis.text=element_text(size=18),plot.margin = margin(10, 20, 10, 10)) + 
  scale_x_discrete(limits=seq(1900,2020,30)) + 
  xlab("Year") + 
  ylab("Number of Grid Cells Sampled") + 
  scale_y_continuous(expand = c(0, 0))
# Grid Cells sampled per era
#sum(GridID.per.year$grids[1:31])
#sum(GridID.per.year$grids[32:61])
#sum(GridID.per.year$grids[62:91])
#sum(GridID.per.year$grids[92:121])


# Bar graphs of number of species in decline, increase, or none
p1 <- ggplot(Species_table, aes(x = RecentTrend, fill=RecentTrend)) +
  geom_bar() +
  labs(x = "Recent Occupancy Trend\n(between 1961\u20131990 and 1991\u20132020)", y = "Number of Species") +
  scale_fill_manual(values = c("Decline" = "#AC4F55", "Increase" = "#2F9BBF", "Stable" = "#B1B1B1")) +
  theme(legend.position="none", panel.background=element_blank(),axis.line=element_line(color="black"),panel.grid.major.x = element_blank(),axis.title = element_text(size = 18), axis.text=element_text(size=18),plot.margin = margin(10, 20, 10, 10), axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(r = 20))) +
  scale_y_continuous(breaks = seq(0, 250, by = 50),limits=c(0,250), expand = c(0, 0)) +
  scale_x_discrete(labels=c("Decline","Increase","Stable"))
p2 <- ggplot(Species_table, aes(x = HistoricTrend, fill=HistoricTrend)) +
  geom_bar() +
  labs(x = "Historical Occupancy Trend\n(between 1900\u20131930 and 1991\u20132020)", y="") +
  scale_fill_manual(values = c("Decline" = "#AC4F55", "Increase" = "#2F9BBF", "Stable" = "#B1B1B1")) +
  theme(legend.position="none", panel.background=element_blank(),axis.line=element_line(color="black"),panel.grid.major.x = element_blank(),axis.title = element_text(size = 18), axis.text=element_text(size=18),plot.margin = margin(10, 20, 10, 10), axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(r = 10))) +
  scale_y_continuous(breaks = seq(0, 250, by = 50),limits=c(0,250), expand = c(0, 0)) +
  scale_x_discrete(labels=c("Decline","Increase","Stable"))
grid.arrange(p1, p2, ncol=2)


# Recent and historic caterpillar plots
# Function to create caterpillar plot of delta occupancy of all species (between two given eras)
caterpillar_plot_delta.occ <- function(era1,era2){
  contrast.to.plot <- ifelse(era1 > era2, paste(era1,"-",era2,sep=""), paste(era2,"-",era1,sep="")) #Find the name of the contrast to plot. Make sure it is in the format [larger era]-[lesser era]
  delta.occ.to.plot <- filter(delta.occSummary, EraContrast==contrast.to.plot) %>% arrange(Mean) #Filter the delta.occupancy values to only the desired era contrast. Sort by ascending Mean
  delta.occ.to.plot$Species <- factor(delta.occ.to.plot$Species, levels=unique(delta.occ.to.plot$Species),ordered=T) #Do this to prevent ggplot from automatically sorting the species on the x axis alphabetically
  #Create caterpillar plot using ggplot2
  ggplot(delta.occ.to.plot[-which(delta.occ.to.plot$Species=="cross.sp.mean"),], aes(x = Species, y = Mean)) + #do not plot the cross.sp.mean as a vertical bar (take it out of the dataset)
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95, color=ifelse(Lower95<0 & Upper95>0,"nochange",ifelse(Upper95<=0,"decline","increase"))), width = 0, position = position_dodge(width = 0), linewidth=2) +
    geom_point(position = position_dodge(width = 0), size=0.75, color="red") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_ribbon(aes(x= 1:relevant.dat$nsp, ymin = delta.occ.to.plot[which(delta.occ.to.plot$Species=="cross.sp.mean"),]$Lower95, ymax = delta.occ.to.plot[which(delta.occ.to.plot$Species=="cross.sp.mean"),]$Upper95), fill = "red", alpha = 0.15) + #Shading for 95% CI of cross species mean
    geom_hline(yintercept = delta.occ.to.plot[which(delta.occ.to.plot$Species=="cross.sp.mean"),]$Mean, linetype = "dashed", color = "red2", linewidth = 0.5) + #cross-species mean
    scale_y_continuous(labels = percent_format()) +
    #scale_y_continuous(breaks = c(-1,-0.105,0,1,2), minor_breaks = c(-0.5,0.5,1.5), labels = c("-100%","-10.5%","0%","100%","200%")) + #Labels for eras 4-3
    scale_y_continuous(breaks = c(-1,0,0.084,1,2,3), minor_breaks = c(-0.5,0.5,1.5,2.5), labels = c("-100%","0%","8.4%","100%","200%","300%")) + #Labels for eras 4-1
    labs(x = paste("Species Identity (n=",relevant.dat$nsp,")",sep=""), y = gsub("-","\u2013",paste("Change in Occupancy between ",EraLabels[ifelse(era1>era2,era2,era1)]," and ",EraLabels[ifelse(era1>era2,era1,era2)]," (%)",sep=""))) +
    scale_color_manual(values=c("nochange" = "#B1B1B1", "decline" = "#AC4F55", "increase" = "#2F9BBF")) +
    guides(color=F) +
    theme_minimal() +
    theme(panel.grid.major.x=element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(size = 22), plot.title = element_text(face = "bold"), axis.title.y = element_text(margin = margin(r = 10), size=22), axis.title.x = element_text(size=22))
}
caterpillar_plot_delta.occ(1,4)


# Results Table showing each species and their recent and historic trends
Results.table.v2 <- Species_table[,c("genusspecies","RecentTrend","HistoricTrend")] %>%
  arrange(RecentTrend,HistoricTrend,genusspecies) %>% #sort by recent trend, then historic trend, then alphabetically
  mutate(ExceptionalTrend = ifelse(genusspecies %in% ExceptionalRecentDecline_sp | genusspecies %in% ExceptionalRecentGrowth_sp,1,0)) #add column denoting whether the trend is an exceptional one (same trend between the most recent era and all 3 prior eras)
write.csv(Results.table.v2, file="Results.table.v2.csv")



# Prior Predictive Check plot is plotted in another script. See PriorPredictiveCheck.R


# Simulation plots are plotted in another script. See Simulation.R


# Posterior Predictive Check
# Percent Similarity
par(mar = c(5, 5, 4, 2) + 0.1) #Increase left margin so y axis title doesnt get cut off
plot.PercentSimilarity.PPC(1)
par(mar = c(5, 4, 4, 2) + 0.1)  # Default margin values
# Number of Detections
par(mar = c(5, 7, 4, 2) + 0.1) #Increase left margin so y axis title doesnt get cut off
par(bty = "l")
plot.numDetections.PPC(1)
mtext("Number of Detections", cex=2, side = 2, line = 5, las = 0)
par(mar = c(5, 4, 4, 2) + 0.1)  # Default margin values


# Sensitivity Analysis
# Total Number of Species in Decline and Growth
grid.arrange(plot.Sensitivity.TotalDeclinesGrowths("Spatial"),plot.Sensitivity.TotalDeclinesGrowths("Temporal (era)"),plot.Sensitivity.TotalDeclinesGrowths("Temporal (interval)"),ncol=3)
# Baseline Occupancy
grid.arrange(plot.Sensitivity.mu.psi("Spatial"),plot.Sensitivity.mu.psi("Temporal (era)"),plot.Sensitivity.mu.psi("Temporal (interval)"),ncol=3)



# Results Tables for Association Analyses
# Plot of RangeSize for recent and historical
par(mfrow=c(1,2))
plot(RangeSizeRecent, cex.lab=1.2, cex.axis=1.25)
mtext("A", side=3, line=2, adj=-0.18, cex=2,font=2)
RangeSizeHistoric$yname <- " " #Remove y axis title on second graph
plot(RangeSizeHistoric, cex.lab=1.2, cex.axis=1.25)
mtext("B", side=3, line=2, adj=-0.18, cex=2,font=2)

# Plot of BodySize for recent and historical
par(mfrow=c(1,2))
boxplot(RecentDeclineSizes,RecentIncreaseSizes, names = c("Decline", "Increase"), xlab="Recent Trend (between 1961\u20131990 and 1991\u20132020)", ylab="Body Size (mm)", cex.lab=1.2, cex.axis=1.25)
mtext("A", side=3, line=2, adj=-0.18, cex=2,font=2)
boxplot(HistoricalDeclineSizes,HistoricalIncreaseSizes, names = c("Decline", "Increase"), xlab="Historical Trend (between 1900\u20131930 and 1991\u20132020)", cex.lab=1.2, cex.axis=1.25)
mtext("B", side=3, line=2, adj=-0.18, cex=2,font=2)

# Mosaic plots
# Subfamily
png(filename="Subfamily.png", width = 1500, height = 600)
grid.arrange(grid.grabExpr(mosaic(X2.subfam.recent, shade=T, legend=F, rot_labels=c(90,90,0,0), just_labels=c("left","right"), offset_varnames = c(3.5,6), margins = c(4, 15, 4, 4), set_varnames=c(RecentTrend="Recent Trend (between 1961\u20131990 and 1991\u20132020)"),gp_labels = gpar(fontsize = 15),gp_varnames = gpar(fontsize = 16, fontface=2))), 
             grid.grabExpr(mosaic(X2.subfam.historic, shade=T, legend=F, rot_labels=c(90,90,0,0), just_labels=c("left","right"), offset_varnames = c(3.5,6), margins = c(4, 15, 4, 4), set_varnames=c(Subfamily="",HistoricTrend="Historical Trend (between 1900\u20131930 and 1991\u20132020)"),gp_labels = gpar(fontsize = 15),gp_varnames = gpar(fontsize = 16, fontface=2) )), ncol=2, respect=T)
dev.off()

# Larval Guild
png(filename="Larval.png", width = 1500, height = 600)
grid.arrange(grid.grabExpr(mosaic(X2.larval.recent, shade=T, legend=F, rot_labels=c(90,90,0,0), just_labels=c("left","right"), offset_varnames = c(3.5,6), margins = c(4, 15, 4, 4), set_varnames=c(LarvalGuild="Larval Guild",RecentTrend="Recent Trend (between 1961\u20131990 and 1991\u20132020)"),gp_labels = gpar(fontsize = 15),gp_varnames = gpar(fontsize = 16, fontface=2))), 
             grid.grabExpr(mosaic(X2.larval.historic, shade=T, legend=F, rot_labels=c(90,90,0,0), just_labels=c("left","right"), offset_varnames = c(3.5,6), margins = c(4, 15, 4, 4), set_varnames=c(LarvalGuild="",HistoricTrend="Historical Trend (between 1900\u20131930 and 1991\u20132020)"),gp_labels = gpar(fontsize = 15),gp_varnames = gpar(fontsize = 16, fontface=2) )), ncol=2, respect=T)
dev.off()

# Rothole Specialists
png(filename="Rothole.png", width = 1500, height = 600)
grid.arrange(grid.grabExpr(mosaic(X2.rothole.recent, shade=T, legend=F, rot_labels=c(90,90,0,0), just_labels=c("left","right"), offset_varnames = c(3.5,6), margins = c(4, 15, 4, 4), set_varnames=c(rothole="",RecentTrend="Recent Trend (between 1961\u20131990 and 1991\u20132020)"),gp_labels = gpar(fontsize = 15),gp_varnames = gpar(fontsize = 16, fontface=2))), 
             grid.grabExpr(mosaic(X2.rothole.historic, shade=T, legend=F, rot_labels=c(90,90,0,0), just_labels=c("left","right"), offset_varnames = c(3.5,6), margins = c(4, 15, 4, 4), set_varnames=c(rothole="",HistoricTrend="Historical Trend (between 1900\u20131930 and 1991\u20132020)"),gp_labels = gpar(fontsize = 15),gp_varnames = gpar(fontsize = 16, fontface=2) )), ncol=2, respect=T)
dev.off()

# Recent vs. Historic
png(filename="RecentHistoric.png", width = 800, height = 600)
mosaic(X2.recent.historic, shade=T, gp=shading_hcl, gp_args = list(interpolate = c(2, 10)), legend=F, rot_labels=c(90,90,0,0), just_labels=c("left","right"), offset_varnames = c(3.5,3.5), margins = c(4, 15, 4, 4), set_varnames=c(HistoricTrend="Historical Trend (between 1900\u20131930 and 1991\u20132020)",RecentTrend="Recent Trend (between 1961\u20131990 and 1991\u20132020)"),gp_labels = gpar(fontsize = 15),gp_varnames = gpar(fontsize = 16, fontface=2)) 
dev.off()









##########################
#### APPENDIX FIGURES ####
##########################

# Appendix C Figure: Strongly Supported Declines and Increases

# Caterpillar plots of just the species in decline or just the species in growth
# Function to bold the species that show the same trend between the current era and all 3 previous eras
xlab_func <- function(sp){
  sapply(sp, function(label) {
    if (label %in% ExceptionalRecentDecline_sp | label %in% ExceptionalRecentGrowth_sp) {
      bquote(bolditalic(.(label)))  # Bold labels "B" and "D"
    } else {
      bquote(italic(.(label)))  # Regular labels for others
    }})}
# Recent Declines
caterpillar_plot_delta.occ_RecentDecline <- function(){
  delta.occ.to.plot <- filter(delta.occSummary, EraContrast==paste(relevant.dat$nera,"-",relevant.dat$nera-1,sep="") & Species %in% RecentDecline_sp) %>% arrange(Mean) #Filter the delta.occupancy values to only the desired era contrast and species in recent decline. Sort by ascending Mean
  delta.occ.to.plot$Species <- factor(delta.occ.to.plot$Species, levels=unique(delta.occ.to.plot$Species),ordered=T) #Do this to prevent ggplot from automatically sorting the species on the x axis alphabetically
  # Compute average cross-species occupancy change for this subset of species (must be done manually here, unlike for all 318 species which was explicitly estimated by the model)
  sum <- 0
  for(sp in delta.occ.to.plot$Species){ #Add all species' delta occupancy posteriors together
    sum <- matrix(as.matrix(delta.occ_list[[sp]][["4-3"]]), nrow=model.out$jags.out$sample, ncol=length(MCMC)) + sum}
  delta.occ <- sum/length(delta.occ.to.plot$Species) #take average
  delta.occ <- as.mcmc.list(lapply(seq_len(ncol(delta.occ)), function(i) as.mcmc(delta.occ[,i]))) #convert back to mcmc object
  cross.sp.Lower95 <- HPDinterval(as.mcmc(unlist(delta.occ)), prob=0.95)[[1]]
  cross.sp.Mean <- summary(delta.occ)$statistics[[1]]
  cross.sp.Upper95 <- HPDinterval(as.mcmc(unlist(delta.occ)), prob=0.95)[[2]]
  # Create caterpillar plot using ggplot2
  ggplot(delta.occ.to.plot, aes(x = Species, y = Mean)) +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95, color=ifelse(Lower95<0 & Upper95>0,"nochange",ifelse(Upper95<=0,"decline","increase"))), width = 0, position = position_dodge(width = 0), linewidth=5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_ribbon(aes(x= 1:length(Species), ymin = cross.sp.Lower95, ymax = cross.sp.Upper95), fill = "red", alpha = 0.2) + #Shading for 95% CI of cross species mean
    geom_hline(yintercept = cross.sp.Mean, linetype = "dashed", color = "red2", linewidth = 0.5) + #cross-species mean
    geom_point(position = position_dodge(width = 0), size=1.5, color="red") +
    coord_cartesian(ylim=c(-1,0)) + #Set y axis
    scale_y_continuous(labels = label_percent()) +
    scale_x_discrete( labels= xlab_func) +
    labs(y = gsub("-","\u2013",paste("Occupancy change between ",EraLabels[relevant.dat$nera-1]," and ",EraLabels[relevant.dat$nera]," (%)",sep=""))) +
    scale_color_manual(values=c("nochange" = "#B1B1B180", "decline" = "#AC4F5580", "increase" = "#2F9BBF80")) +
    guides(color=F) +
    theme_minimal() +
    theme(panel.grid.major.x=element_blank(), axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5, size=12), axis.text.y = element_text(size = 18), axis.title.y = element_text(margin = margin(r = 10), size=20, hjust=1), axis.title.x = element_text(size=20))}
caterpillar_plot_delta.occ_RecentDecline()
# Historic Declines
caterpillar_plot_delta.occ_HistoricDecline <- function(){
  delta.occ.to.plot <- filter(delta.occSummary, EraContrast==paste(relevant.dat$nera,"-",1,sep="") & Species %in% HistoricDecline_sp) %>% arrange(Mean) #Filter the delta.occupancy values to only the desired era contrast and species in recent decline. Sort by ascending Mean
  delta.occ.to.plot$Species <- factor(delta.occ.to.plot$Species, levels=unique(delta.occ.to.plot$Species),ordered=T) #Do this to prevent ggplot from automatically sorting the species on the x axis alphabetically
  sum <- 0
  for(sp in delta.occ.to.plot$Species){ #Add all species' delta occupancy posteriors together
    sum <- matrix(as.matrix(delta.occ_list[[sp]][["4-1"]]), nrow=model.out$jags.out$sample, ncol=length(MCMC)) + sum}
  delta.occ <- sum/length(delta.occ.to.plot$Species) #take average
  delta.occ <- as.mcmc.list(lapply(seq_len(ncol(delta.occ)), function(i) as.mcmc(delta.occ[,i]))) #convert back to mcmc object
  cross.sp.Lower95 <- HPDinterval(as.mcmc(unlist(delta.occ)), prob=0.95)[[1]]
  cross.sp.Mean <- summary(delta.occ)$statistics[[1]]
  cross.sp.Upper95 <- HPDinterval(as.mcmc(unlist(delta.occ)), prob=0.95)[[2]]
  ggplot(delta.occ.to.plot, aes(x = Species, y = Mean)) +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95, color=ifelse(Lower95<0 & Upper95>0,"nochange",ifelse(Upper95<=0,"decline","increase"))), width = 0, position = position_dodge(width = 0), linewidth=5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_ribbon(aes(x= 1:length(Species), ymin = cross.sp.Lower95, ymax = cross.sp.Upper95), fill = "red", alpha = 0.2) + #Shading for 95% CI of cross species mean
    geom_hline(yintercept = cross.sp.Mean, linetype = "dashed", color = "red2", linewidth = 0.5) + #cross-species mean
    geom_point(position = position_dodge(width = 0), size=1.5, color="red") +
    coord_cartesian(ylim=c(-1,0)) + #Set y axis
    scale_y_continuous(labels = label_percent()) +
    scale_x_discrete( labels= xlab_func) +
    labs(y = gsub("-","\u2013",paste("Occupancy change between ",EraLabels[1]," and ",EraLabels[relevant.dat$nera]," (%)",sep=""))) +
    scale_color_manual(values=c("nochange" = "#B1B1B180", "decline" = "#AC4F5580", "increase" = "#2F9BBF80")) +
    guides(color=F) +
    theme_minimal() +
    theme(panel.grid.major.x=element_blank(), axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5, size=12), axis.text.y = element_text(size = 18), axis.title.y = element_text(margin = margin(r = 10), size=20, hjust=1), axis.title.x = element_text(size=20))}
caterpillar_plot_delta.occ_HistoricDecline()
# Recent Growth
caterpillar_plot_delta.occ_RecentGrowth <- function(){
  delta.occ.to.plot <- filter(delta.occSummary, EraContrast==paste(relevant.dat$nera,"-",relevant.dat$nera-1,sep="") & Species %in% RecentGrowth_sp) %>% arrange(Mean) #Filter the delta.occupancy values to only the desired era contrast and species in recent decline. Sort by ascending Mean
  delta.occ.to.plot$Species <- factor(delta.occ.to.plot$Species, levels=unique(delta.occ.to.plot$Species),ordered=T) #Do this to prevent ggplot from automatically sorting the species on the x axis alphabetically
  sum <- 0
  for(sp in delta.occ.to.plot$Species){ #Add all species' delta occupancy posteriors together
    sum <- matrix(as.matrix(delta.occ_list[[sp]][["4-3"]]), nrow=model.out$jags.out$sample, ncol=length(MCMC)) + sum}
  delta.occ <- sum/length(delta.occ.to.plot$Species) #take average
  delta.occ <- as.mcmc.list(lapply(seq_len(ncol(delta.occ)), function(i) as.mcmc(delta.occ[,i]))) #convert back to mcmc object
  cross.sp.Lower95 <- HPDinterval(as.mcmc(unlist(delta.occ)), prob=0.95)[[1]]
  cross.sp.Mean <- summary(delta.occ)$statistics[[1]]
  cross.sp.Upper95 <- HPDinterval(as.mcmc(unlist(delta.occ)), prob=0.95)[[2]]
  ggplot(delta.occ.to.plot, aes(x = Species, y = Mean)) +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95, color=ifelse(Lower95<0 & Upper95>0,"nochange",ifelse(Upper95<=0,"decline","increase"))), width = 0, position = position_dodge(width = 0), linewidth=5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_ribbon(aes(x= 1:length(Species), ymin = cross.sp.Lower95, ymax = cross.sp.Upper95), fill = "red", alpha = 0.2) + #Shading for 95% CI of cross species mean
    geom_hline(yintercept = cross.sp.Mean, linetype = "dashed", color = "red2", linewidth = 0.5) + #cross-species mean
    geom_point(position = position_dodge(width = 0), size=1.5, color="red") +
    coord_cartesian(ylim=c(0,2)) + #Set y axis
    scale_y_continuous(labels = label_percent(), breaks = c(0,1,2)) +
    scale_x_discrete( labels= xlab_func) +
    labs(y = gsub("-","\u2013",paste("Occupancy change between ",EraLabels[relevant.dat$nera-1]," and ",EraLabels[relevant.dat$nera]," (%)",sep=""))) +
    scale_color_manual(values=c("nochange" = "#B1B1B180", "decline" = "#AC4F5580", "increase" = "#2F9BBF80")) +
    guides(color=F) +
    theme_minimal() +
    theme(panel.grid.major.x=element_blank(), axis.text.x = element_text(angle=90, hjust=1,vjust = 0.5, size=12), axis.text.y = element_text(size = 18), axis.title.y = element_text(margin = margin(r = 10), size=16, hjust=1), axis.title.x = element_text(size=16))}
caterpillar_plot_delta.occ_RecentGrowth()
# Historic Growth
caterpillar_plot_delta.occ_HistoricGrowth <- function(){
  delta.occ.to.plot <- filter(delta.occSummary, EraContrast==paste(relevant.dat$nera,"-",1,sep="") & Species %in% HistoricGrowth_sp) %>% arrange(Mean) #Filter the delta.occupancy values to only the desired era contrast and species in recent decline. Sort by ascending Mean
  delta.occ.to.plot$Species <- factor(delta.occ.to.plot$Species, levels=unique(delta.occ.to.plot$Species),ordered=T) #Do this to prevent ggplot from automatically sorting the species on the x axis alphabetically
  sum <- 0
  for(sp in delta.occ.to.plot$Species){ #Add all species' delta occupancy posteriors together
    sum <- matrix(as.matrix(delta.occ_list[[sp]][["4-1"]]), nrow=model.out$jags.out$sample, ncol=length(MCMC)) + sum}
  delta.occ <- sum/length(delta.occ.to.plot$Species) #take average
  delta.occ <- as.mcmc.list(lapply(seq_len(ncol(delta.occ)), function(i) as.mcmc(delta.occ[,i]))) #convert back to mcmc object
  cross.sp.Lower95 <- HPDinterval(as.mcmc(unlist(delta.occ)), prob=0.95)[[1]]
  cross.sp.Mean <- summary(delta.occ)$statistics[[1]]
  cross.sp.Upper95 <- HPDinterval(as.mcmc(unlist(delta.occ)), prob=0.95)[[2]]
  ggplot(delta.occ.to.plot, aes(x = Species, y = Mean)) +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95, color=ifelse(Lower95<0 & Upper95>0,"nochange",ifelse(Upper95<=0,"decline","increase"))), width = 0, position = position_dodge(width = 0), linewidth=5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_ribbon(aes(x= 1:length(Species), ymin = cross.sp.Lower95, ymax = cross.sp.Upper95), fill = "red", alpha = 0.2) + #Shading for 95% CI of cross species mean
    geom_hline(yintercept = cross.sp.Mean, linetype = "dashed", color = "red2", linewidth = 0.5) + #cross-species mean
    geom_point(position = position_dodge(width = 0), size=1.5, color="red") +
    coord_cartesian(ylim=c(0,3)) + #Set y axis
    scale_y_continuous(labels = label_percent()) +
    scale_x_discrete( labels= xlab_func) +
    labs(y = gsub("-","\u2013",paste("Occupancy change between ",EraLabels[1]," and ",EraLabels[relevant.dat$nera]," (%)",sep=""))) +
    scale_color_manual(values=c("nochange" = "#B1B1B180", "decline" = "#AC4F5580", "increase" = "#2F9BBF80")) +
    guides(color=F) +
    theme_minimal() +
    theme(panel.grid.major.x=element_blank(), axis.text.x = element_text(angle=90, hjust=1,vjust = 0.5, size=12), axis.text.y = element_text(size = 18), axis.title.y = element_text(margin = margin(r = 10), size=20, hjust=1), axis.title.x = element_text(size=20))}
caterpillar_plot_delta.occ_HistoricGrowth()



# Appendix E Figure: All pairwise contrasts between eras.

caterpillar_plot_delta.occ_all <- function(era1,era2){
  contrast.to.plot <- ifelse(era1 > era2, paste(era1,"-",era2,sep=""), paste(era2,"-",era1,sep="")) #Find the name of the contrast to plot. Make sure it is in the format [larger era]-[lesser era]
  delta.occ.to.plot <- filter(delta.occSummary, EraContrast==contrast.to.plot) %>% arrange(Mean) #Filter the delta.occupancy values to only the desired era contrast. Sort by ascending Mean
  delta.occ.to.plot$Species <- factor(delta.occ.to.plot$Species, levels=unique(delta.occ.to.plot$Species),ordered=T) #Do this to prevent ggplot from automatically sorting the species on the x axis alphabetically
  #Create caterpillar plot using ggplot2
  ggplot(delta.occ.to.plot[-which(delta.occ.to.plot$Species=="cross.sp.mean"),], aes(x = Species, y = Mean)) + #do not plot the cross.sp.mean as a vertical bar (take it out of the dataset)
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95, color=ifelse(Lower95<0 & Upper95>0,"nochange",ifelse(Upper95<=0,"decline","increase"))), width = 0, position = position_dodge(width = 0), linewidth=0.8) +
    geom_point(position = position_dodge(width = 0), size=0.2, color="red") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_ribbon(aes(x= 1:relevant.dat$nsp, ymin = delta.occ.to.plot[which(delta.occ.to.plot$Species=="cross.sp.mean"),]$Lower95, ymax = delta.occ.to.plot[which(delta.occ.to.plot$Species=="cross.sp.mean"),]$Upper95), fill = "red", alpha = 0.2) + #Shading for 95% CI of cross species mean
    geom_hline(yintercept = delta.occ.to.plot[which(delta.occ.to.plot$Species=="cross.sp.mean"),]$Mean, linetype = "dashed", color = "red2", linewidth = 0.5) + #cross-species mean
    scale_y_continuous(limits=c(-1,4), labels = percent_format()) +
    labs(y = ifelse(contrast.to.plot=="4-1", "Change in Occupancy (%)","")) +
    scale_color_manual(values=c("nochange" = "#B1B1B1", "decline" = "#AC4F55", "increase" = "#2F9BBF")) +
    guides(color=F) +
    theme_minimal() +
    if (contrast.to.plot=="4-3" | contrast.to.plot=="4-1" | contrast.to.plot=="3-1"){ #Remove y axis labels for certain plots
      theme(panel.grid.major.x=element_blank(), axis.title.y.left =element_text(size=25), axis.text.x=element_blank(), axis.text.y = element_text(size = 20), axis.title.y = element_text(margin = margin(r = 10), size=20), axis.title.x = element_blank())
    }else{
      theme(panel.grid.major.x=element_blank(), axis.title.y.left =element_blank(), axis.text.x=element_blank(), axis.text.y = element_blank(), axis.title.y = element_text(margin = margin(r = 10), size=20), axis.title.x = element_blank())
    }
}

spacer <- nullGrob()
plots <- grid.arrange(caterpillar_plot_delta.occ_all(3,4),caterpillar_plot_delta.occ_all(2,4),spacer,spacer,caterpillar_plot_delta.occ_all(1,4),caterpillar_plot_delta.occ_all(2,3),spacer,spacer,caterpillar_plot_delta.occ_all(1,3),caterpillar_plot_delta.occ_all(1,2),ncol=2, heights=c(1,0.08,1,0.08,1), widths=c(1,0.88))
x_axis_label <- textGrob("Species Identity (n=318)", rot = 0, vjust = 1, hjust = 0.5, gp = gpar(fontsize = 18))
grid.arrange(plots, x_axis_label, ncol = 1, heights = c(5, 0.2))




# Appendix F: Repeated Simulations is plotted in a different script. See Simulation.R




# Appendix G Figure: Repeated Post Predictive Checks

load("PosteriorPredictiveChecks.RData") #load Posterior Predictive Check results

# Stacked Barplot of Percent Similarity of observations for each species, for a given simulated dataset
plot.PercentSimilarity.PPC <- function(PPCnum){
  PostPredictiveCheck.plot <- arrange(PostPredictiveChecks[[PPCnum]][-1,],PercentSimilarity) #sort the PostPredCheck by Percent Similarity, and remove the first row that displays totals.
  StackedBars <- rbind(100*PostPredictiveCheck.plot$RightZeroes/PostPredictiveCheck.plot$NumObservations,
                       100*PostPredictiveCheck.plot$RightOnes/PostPredictiveCheck.plot$NumObservations,
                       100*PostPredictiveCheck.plot$WrongZeroes/PostPredictiveCheck.plot$NumObservations,
                       100*PostPredictiveCheck.plot$WrongOnes/PostPredictiveCheck.plot$NumObservations)
  barplot(StackedBars, ylim=c(0,100), col=c("#348FA4","#74C0D2","#FFE479","#F6FEAA"), yaxt = ifelse(PPCnum==1,"s","n"), ylab=ifelse(PPCnum==1,"Percentage of the species' dataset (%)",""), xlab= ifelse(PPCnum==2,paste("Species Identity (n=",length(splist),")",sep=""),""), space=0, border="transparent", cex.lab=1.6, las=1) #Create stacked bar plot
  #axis(side = 2, at = seq(0, 100, by = 10)) #y-axis ticks go up by 10's
  abline(a=mean(PostPredictiveChecks[[PPCnum]][-1,]$PercentSimilarity), b=0, lty=1, col="red") #average percent similarity
  #mtext("Percent Similarity between Simulated and Real Datasets", cex=1.5, font=2, line=2.5)
  #mtext(paste("(Simulated Dataset #",PPCnum,"; red line denotes average)",sep=""), cex=1, line=1.5)
  if(PPCnum==3){
  legend(c(60,310),c(50,1),c("Simulated 1, Real 0", "Simulated 0, Real 1","Matching 1's","Matching 0's"), fill = c("#F6FEAA","#FFE479","#74C0D2","#348FA4"), cex = 1.6, inset = c(0.1, 0.35), y.intersp = 1, x.intersp =0.2)}
}

plot.numDetections.PPC <- function(PPCnum){
  PostPredictiveCheck.plot <- arrange(PostPredictiveChecks[[PPCnum]][-1,],PositiveObs.real) #sort the PPC by number of real detections, and remove the first row that displays totals.
  plot(PostPredictiveCheck.plot$PositiveObs.real, pch=19, col="blue", cex=0.5, ylim=c(0,1000), ylab= ifelse(PPCnum==1,"Number of Detections",""), xlab=ifelse(PPCnum==2,paste("Species Identity (n=",length(splist),")",sep=""),""), xaxt="n", cex.lab=1.6, las=1, yaxt="n")
  points(PostPredictiveCheck.plot$PositiveObs.Simulated, pch=19, col="red", cex=0.5)
  if(PPCnum==1){
    axis(2, at = c(0,200,400,600,800,1000), las=1)}
  #mtext("Number of Detections between Simulated and Real Datasets", cex=1.5, font=2, line=2.5)
  #mtext(paste("(Simulated Dataset #",PPCnum,"; blue = real value, red = simulated value)",sep=""), cex=1, line=1.5)
}

par(mar = c(5, 5, 2, 1) + 0.1, mfrow = c(2,3)) #Increase left margin so y axis title doesnt get cut off
plot.PercentSimilarity.PPC(1)
plot.PercentSimilarity.PPC(2)
plot.PercentSimilarity.PPC(3)
plot.numDetections.PPC(1)
plot.numDetections.PPC(2)
plot.numDetections.PPC(3)
#par(mar = c(5, 4, 4, 2) + 0.1)  # Default margin values






### get the dataset
u0 <- paste0(
    'http://inla.r-inla-download.org/',
    'r-inla.org/case-studies/Cameletti2012/')
coofl <- 'coordinates.csv'
datafl <- 'Piemonte_data_byday.csv'
bordersfl <- 'Piemonte_borders.csv'

if(!file.exists(bordersfl)) 
    download.file(paste0(u0, bordersfl), bordersfl)

dim(pborders <- read.csv(bordersfl))

if(!file.exists(coofl)) 
    download.file(paste0(u0, coofl), coofl)
dim(locs <- read.csv(coofl))

if(!file.exists(datafl)) 
    download.file(paste0(u0, datafl), datafl)
dim(pdata <- read.csv(datafl))

head(pdata)

### prepare and select time 
range(pdata$Date <- as.Date(pdata$Date, '%d/%m/%y'))
pdata$time <- as.integer(difftime(
    pdata$Date, min(pdata$Date), units='days'))+1


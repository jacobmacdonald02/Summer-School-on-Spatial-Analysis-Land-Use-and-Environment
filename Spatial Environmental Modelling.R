###-------------------------------------------------------------------------------------------------
##  Applications of Spatial Analysis with Spatial Econometrics:
##    Spatial Urban-Environmental Modelling In R
##  
##  Summer School - January, 2020
##    Universidad del Pacifico 
##    Lima, Peru
## 
##  Jacob L. Macdonald
##    Jacob.Macdonald@liverpool.ac.uk
## 
## Version: January 27, 2020
###-------------------------------------------------------------------------------------------------


### TABLE OF CONTENTS - COURSE OUTLINE
###-------------------------------------------------------------------------------------------------
### SECTION 1 - WORKING AND DATABASE MANAGEMENT IN R
##            LESSON 00: SETTING THE WORKING DIRECTORY AND READING PACKAGES
##            LESSON 01: READ IN HOUSING DATASET
##            LESSON 02: GENERATING BASIC DATABASE VARIABLES
##            LESSON 03: SUMMARY STATISTICS, CORRELATIONS, SCATTER PLOTS
## 
### SECTION 2 - SPATIAL INTERPOLATION OF POLLUTION LEVELS
##            LESSON 04: LOADING AND GENERATING SPATIAL DATA (POINTS, POLYGONS, RASTERS)
##            LESSON 05: INVERSE DISTANCE AND INVERSE SQUARED DISTANCE INTERPOLATION
##            LESSON 06: NEIGHBOURHOOD LEVEL AGGREGATION
## 
### SECTION 3 - SPATIAL MODELING DIAGNOSTICS AND TECHNIQUES
##            LESSON 07: FIT BASELINE OLS MODEL
##            LESSON 08: IMPROVEMENT OF BASELINE MODEL WITH GIS VARIABLES (OLS)
##            LESSON 09: CONTROLLING FOR LOCATION PRICE DIFFERENTIALS
##            LESSON 10: SPATIAL DEPENDANCE - GENERATING AND VISUALIZING SPATIAL WEIGHTS
##            LESSON 11: SPATIAL DIAGNOSTICS AND MODEL SELECTION
##            LESSON 12: MODELS OF SPATIAL DEPENDENCE (SPATIAL LAG, ERROR)
###-------------------------------------------------------------------------------------------------



### SECTION 1 - WORKING AND DATABASE MANAGEMENT IN R
####################################################################################################

# LESSON 00: SETTING THE WORKING DIRECTORY AND READING PACKAGES
###-------------------------------------------------------------------------------------------------
# Setting the path on a Windows OS
wd <- "C:/Desktop/Spatial Environmental Modelling"
# Setting the path on a Mac OS
wd <- paste0("/Users/", Sys.info()[[7]], "/Desktop/Spatial Environmental Modelling")
# Setting the path on a Linux OS
wd <- paste0("/Users/", Sys.info()[[7]], "/Dropbox/Teaching/Spatial Environmental Modelling/Spatial Environmental Modelling")

setwd(wd)

# install.packages(c("sp", "spdep", "rgeos", "rgdal", "ggplot2", "geosphere",  
#   "gstat", "raster", "car", "lmtest", "tseries", "sandwich", "Matrix", "stargazer"))
# install.packages("zeallot")
# install.packages("backports")

library(sp)
library(spdep)
library(rgeos)
library(rgdal)
library(ggplot2)
library(geosphere)
library(gstat)
library(raster)
library(car)
library(lmtest)
library(tseries)
library(sandwich)
library(Matrix)
library(stargazer)
###-------------------------------------------------------------------------------------------------

# LESSON 01: READ IN HOUSING DATASET
###-------------------------------------------------------------------------------------------------
housing.sample <- read.csv(paste0(wd, "/Data/housing_sample.csv"), fileEncoding="latin1", stringsAsFactors=T)

names(housing.sample)
dim(housing.sample)
str(housing.sample)
summary(housing.sample)
###-------------------------------------------------------------------------------------------------

# LESSON 02: GENERATING BASIC DATABASE VARIABLES
###-------------------------------------------------------------------------------------------------
housing.sample$price.m2 <- housing.sample$price/housing.sample$area
summary(housing.sample$price.m2)

housing.sample$new.garden <- housing.sample$new*housing.sample$garden
housing.sample$area.2 <- housing.sample$area*housing.sample$area
summary(housing.sample$new.garden)
summary(housing.sample$area.2)

housing.sample$ln.area <- log(housing.sample$area)
housing.sample$ln.price <- log(housing.sample$price)
summary(housing.sample$ln.area)
summary(housing.sample$ln.price)

housing.sample$ln.near.park <- log(housing.sample$near.park)
housing.sample$ln.avgD.park <- log(housing.sample$avgD.park)
housing.sample$ln.near.shopping <- log(housing.sample$near.shopping)
housing.sample$ln.avgD.shopping <- log(housing.sample$avgD.shopping)

housing.sample$ln.dist.baixa <- log(housing.sample$dist.baixa)
housing.sample$ln.dist.tagus <- log(housing.sample$dist.tagus)
###-------------------------------------------------------------------------------------------------

# LESSON 03: SUMMARY STATISTICS, CORRELATIONS, SCATTER PLOTS
###-------------------------------------------------------------------------------------------------
summary(housing.sample)

summary(housing.sample$price)
sd(housing.sample$price)
hist(housing.sample$price)
hist(housing.sample$ln.price)

cor(housing.sample[, unlist(lapply(housing.sample, is.numeric))])
cor(housing.sample$price, housing.sample$area)
cor(housing.sample$ln.price, housing.sample$ln.area)

ggplot(data = housing.sample, mapping = aes(x = area, y = price)) +
  geom_smooth(col = "black", size = 0.5) +
  geom_point(size = 0.5) +
  labs(x = "Area m^2", y = "Price")

ggplot(data = housing.sample, mapping = aes(x = ln.area, y = ln.price)) +
  geom_smooth(col = "black", size = 0.5) +
  geom_point(size = 0.5) +
  labs(x = "ln(Area m^2)", y = "ln(Price)")

ggplot(data = housing.sample, mapping = aes(x = ln.dist.baixa, y = ln.dist.tagus)) +
  geom_smooth(col = "black", size = 0.5) +
  geom_point(size = 0.5) +
  labs(x = "ln(Distance to Baixa)", y = "ln(Distance to River)")
cor(housing.sample$ln.dist.baixa, housing.sample$ln.dist.tagus)
###-------------------------------------------------------------------------------------------------


### SECTION 2 - SPATIAL INTERPOLATION OF POLLUTION LEVELS
####################################################################################################

# LESSON 04: LOADING AND GENERATING SPATIAL DATA (POINTS, POLYGONS, RASTERS)
###-------------------------------------------------------------------------------------------------
housing.sample <- SpatialPointsDataFrame(cbind(housing.sample$longitude, housing.sample$latitude), housing.sample)
proj4string(housing.sample) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

housing.sample
summary(housing.sample)
summary(housing.sample@data)
summary(housing.sample@data$price)

plot(housing.sample)

ggplot() + geom_point(data=housing.sample@data, aes(x=longitude, y=latitude, col = price)) +
  scale_colour_gradient(low = "blue", high = "red")


pollution <- read.csv(paste0(wd, "/Data/pollution_2007.csv"), stringsAsFactors=T)
pollution$lnPM10 <- log(pollution$PM10.2007)
pollution <- SpatialPointsDataFrame(cbind(pollution$LONGITUDE, pollution$LATITUDE), pollution)
proj4string(pollution) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

int.unit <- extent(bbox(spTransform(pollution, CRS("+init=epsg:3763"))))
int.unit <- raster(ncol=abs(int.unit[1]-int.unit[2])/100, nrow=abs(int.unit[3]-int.unit[4])/100, xmn=int.unit[1], xmx=int.unit[2], ymn=int.unit[3], ymx=int.unit[4])
proj4string(int.unit) <- CRS("+init=epsg:3763")
int.unit <- projectRaster(int.unit, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

agg.unit <- spTransform(readOGR(dsn = paste0(wd, "/Data"), layer = "freguesias"), CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
agg.unit <- agg.unit[,c("FREGUESIA")]
###-------------------------------------------------------------------------------------------------

# LESSON 05: INVERSE DISTANCE AND INVERSE SQUARED DISTANCE INTERPOLATION
###-------------------------------------------------------------------------------------------------
sp.interpolation <- mask(crop(interpolate(int.unit, gstat(formula=lnPM10 ~ 1, locations=pollution, set=list(idp=1))), agg.unit), agg.unit)
plot(sp.interpolation)
plot(exp(sp.interpolation))

RMSE <- list()
for (k in 1:nrow(pollution)){
  init.val <- c(3, 2)
  RMSE[[k]] <- sqrt(as.numeric((predict(gstat(formula=as.formula("lnPM10 ~ 1"), locations=pollution[-k,], set=list(idp=1)), newdata=pollution[k,], debug.level=0)$var1.pred -
      pollution[k,][,c("lnPM10")]@data)^2))
}
mean(unlist(RMSE))

sp.interpolation <- mask(crop(interpolate(int.unit, gstat(formula=lnPM10 ~ 1, locations=pollution, set=list(idp=2))), agg.unit), agg.unit)
plot(sp.interpolation)
plot(exp(sp.interpolation))

RMSE <- list()
for (k in 1:nrow(pollution)){
  init.val <- c(3, 2)
  RMSE[[k]] <- sqrt(as.numeric((predict(gstat(formula=as.formula("lnPM10 ~ 1"), locations=pollution[-k,], set=list(idp=2)), newdata=pollution[k,], debug.level=0)$var1.pred -
      pollution[k,][,c("lnPM10")]@data)^2))
}
mean(unlist(RMSE))
###-------------------------------------------------------------------------------------------------

# LESSON 06: NEIGHBOURHOOD LEVEL AGGREGATION
###-------------------------------------------------------------------------------------------------
sp.interpolation <- mask(crop(interpolate(int.unit, gstat(formula=lnPM10 ~ 1, locations=pollution, set=list(idp=2))), agg.unit), agg.unit)
sp.aggregation <- do.call(rbind, lapply(extract(sp.interpolation, agg.unit), FUN=mean))
agg.unit <- data.frame(agg.unit, lnPM10 = sp.aggregation)

housing.sample <- merge(housing.sample, agg.unit, by.x="freguesia", by.y="FREGUESIA", all.x=T, sort=F)
###-------------------------------------------------------------------------------------------------


### SECTION 3 - SPATIAL MODELING DIAGNOSTICS AND TECHNIQUES
####################################################################################################

# LESSON 07: FIT BASELINE OLS MODEL AND RUN DIAGNOSTICS
###-------------------------------------------------------------------------------------------------
OLS.baseline <- list()
OLS.baseline[[1]] <- lm(ln.price ~ area + area.2 + new + pool + parking + fireplace + dwindows + garden + aircond + view, na.action=na.omit, data=housing.sample@data)
summary(OLS.baseline[[1]])

OLS.baseline[[2]] <- lm(ln.price ~ area + area.2 + new + parking + fireplace + aircond + view, na.action=na.omit, data=housing.sample@data)
summary(OLS.baseline[[2]])

anova(OLS.baseline[[1]], OLS.baseline[[2]])

lapply(OLS.baseline, vif)
lapply(OLS.baseline, AIC)

lapply(OLS.baseline, function(x) tail(anova(x)[,2], 1))

lapply(OLS.baseline, function(x) jarque.bera.test(x$residuals))
hist(OLS.baseline[[1]]$residuals)
hist(OLS.baseline[[2]]$residuals)

lapply(OLS.baseline, dwtest)
lapply(OLS.baseline, bgtest)
lapply(OLS.baseline, bptest)

lapply(OLS.baseline, function(y){ coeftest(y, function(x) vcovHC(x, type="HC0")) })


stargazer(OLS.baseline, type="html", title="OLS Baseline Model (Robust S.E.)", style="aer", 
  star.cutoffs=c(0.10, 0.05, 0.01), se = lapply(OLS.baseline, function(y){ coeftest(y, function(x) vcovHC(x, type="HC0"))[,2] }),
  add.lines = list(c("AIC", unlist(lapply(OLS.baseline, function(x){ tryCatch(as.character(signif(AIC(x), digits=5)), error=function(err) NA) }))),
    c("SSE", unlist(lapply(OLS.baseline, function(x){ tryCatch(as.character(signif(tail(anova(x)[,2], 1), digits=5)), error=function(err) NA) }))),
    c("Jarque-Bera", unlist(lapply(OLS.baseline, function(x){ tryCatch(ifelse(jarque.bera.test(x$residuals)$p.value > 0.01, paste0(as.character(signif(jarque.bera.test(x$residuals)$statistic, digits=5)), ""), 
      ifelse(jarque.bera.test(x$residuals)$p.value > 0.05 & jarque.bera.test(x$residuals)$p.value < 0.01, paste0(as.character(signif(jarque.bera.test(x$residuals)$statistic, digits=5)), "*"), 
	    ifelse(jarque.bera.test(x$residuals)$p.value > 0.01 & jarque.bera.test(x$residuals)$p.value < 0.05, paste0(as.character(signif(jarque.bera.test(x$residuals)$statistic, digits=5)), "**"), 
	    ifelse(jarque.bera.test(x$residuals)$p.value < 0.01, paste0(as.character(signif(jarque.bera.test(x$residuals)$statistic, digits=5)), "***"), NA)))), error=function(err) NA) }))),
    c("Breusch-Pagan", unlist(lapply(OLS.baseline, function(x){ tryCatch(ifelse(bptest(x)$p.value > 0.01, paste0(as.character(signif(bptest(x)$statistic, digits=5)), ""), 
	    ifelse(bptest(x)$p.value > 0.05 & bptest(x)$p.value < 0.01, paste0(as.character(signif(bptest(x)$statistic, digits=5)), "*"), 
	    ifelse(bptest(x)$p.value > 0.01 & bptest(x)$p.value < 0.05, paste0(as.character(signif(bptest(x)$statistic, digits=5)), "**"), 
	    ifelse(bptest(x)$p.value < 0.01, paste0(as.character(signif(bptest(x)$statistic, digits=5)), "***"), NA)))), error=function(err) NA) }))),
    c("Breusch-Godfrey", unlist(lapply(OLS.baseline, function(x){ tryCatch(ifelse(bgtest(x)$p.value > 0.01, paste0(as.character(signif(bgtest(x)$statistic, digits=5)), ""), 
	    ifelse(bgtest(x)$p.value > 0.05 & bgtest(x)$p.value < 0.01, paste0(as.character(signif(bgtest(x)$statistic, digits=5)), "*"), 
	    ifelse(bgtest(x)$p.value > 0.01 & bgtest(x)$p.value < 0.05, paste0(as.character(signif(bgtest(x)$statistic, digits=5)), "**"), 
	    ifelse(bgtest(x)$p.value < 0.01, paste0(as.character(signif(bgtest(x)$statistic, digits=5)), "***"), NA)))), error=function(err) NA) }))),
	  c("Durbin-Watson", unlist(lapply(OLS.baseline, function(x){ tryCatch(ifelse(dwtest(x)$p.value > 0.01, paste0(as.character(signif(dwtest(x)$statistic, digits=5)), ""), 
	    ifelse(dwtest(x)$p.value > 0.05 & dwtest(x)$p.value < 0.01, paste0(as.character(signif(dwtest(x)$statistic, digits=5)), "*"), 
	    ifelse(dwtest(x)$p.value > 0.01 & dwtest(x)$p.value < 0.05, paste0(as.character(signif(dwtest(x)$statistic, digits=5)), "**"), 
	    ifelse(dwtest(x)$p.value < 0.01, paste0(as.character(signif(dwtest(x)$statistic, digits=5)), "***"), NA)))), error=function(err) NA) }))),
	  c("Mean adj(VIF) of Int.", unlist(lapply(OLS.baseline, function(x){ tryCatch(as.character(signif(mean(vif(x)), digits=5)), error=function(err) NA) }))),
	  c("Median adj(VIF) of Int.", unlist(lapply(OLS.baseline, function(x){ tryCatch(as.character(signif(median(vif(x)), digits=5)), error=function(err) NA) }))),
	  c("Max adj(VIF) of Int.", unlist(lapply(OLS.baseline, function(x){ tryCatch(as.character(signif(max(vif(x)), digits=5)), error=function(err) NA) })))),
  digits = 5, out=paste0(wd, "/OLS_Baseline.htm"))
###-------------------------------------------------------------------------------------------------

# LESSON 08: IMPROVEMENT OF BASELINE MODEL WITH GIS VARIABLES
###-------------------------------------------------------------------------------------------------
OLS.gis <- list()
OLS.gis[[1]] <- lm(ln.price ~ area + area.2 + new + parking + fireplace + aircond + view + Ppop.unemp, na.action=na.omit, data=housing.sample@data)
summary(OLS.gis[[1]])

OLS.gis[[2]] <- lm(ln.price ~ area + area.2 + new + parking + fireplace + aircond + view + Ppop.unemp + ln.dist.baixa, na.action=na.omit, data=housing.sample@data)
summary(OLS.gis[[2]])
OLS.gis[[3]] <- lm(ln.price ~ area + area.2 + new + parking + fireplace + aircond + view + Ppop.unemp + ln.dist.tagus, na.action=na.omit, data=housing.sample@data)
summary(OLS.gis[[3]])

OLS.gis[[4]] <- lm(ln.price ~ area + area.2 + new + parking + fireplace + aircond + view + Ppop.unemp + ln.dist.baixa + lnPM10 + ln.avgD.park , na.action=na.omit, data=housing.sample@data)
summary(OLS.gis[[4]])
OLS.gis[[5]] <- lm(ln.price ~ area + area.2 + new + parking + fireplace + aircond + view + Ppop.unemp + ln.dist.tagus + lnPM10 + ln.avgD.park , na.action=na.omit, data=housing.sample@data)
summary(OLS.gis[[5]])


stargazer(OLS.gis, type="html", title="OLS Models with GIS Variables (Robust S.E.)", style="aer", 
  star.cutoffs=c(0.10, 0.05, 0.01), se = lapply(OLS.gis, function(y){ coeftest(y, function(x) vcovHC(x, type="HC0"))[,2] }),
  add.lines = list(c("AIC", unlist(lapply(OLS.gis, function(x){ tryCatch(as.character(signif(AIC(x), digits=5)), error=function(err) NA) }))),
    c("SSE", unlist(lapply(OLS.gis, function(x){ tryCatch(as.character(signif(tail(anova(x)[,2], 1), digits=5)), error=function(err) NA) }))),
    c("Jarque-Bera", unlist(lapply(OLS.gis, function(x){ tryCatch(ifelse(jarque.bera.test(x$residuals)$p.value > 0.01, paste0(as.character(signif(jarque.bera.test(x$residuals)$statistic, digits=5)), ""), 
      ifelse(jarque.bera.test(x$residuals)$p.value > 0.05 & jarque.bera.test(x$residuals)$p.value < 0.01, paste0(as.character(signif(jarque.bera.test(x$residuals)$statistic, digits=5)), "*"), 
	    ifelse(jarque.bera.test(x$residuals)$p.value > 0.01 & jarque.bera.test(x$residuals)$p.value < 0.05, paste0(as.character(signif(jarque.bera.test(x$residuals)$statistic, digits=5)), "**"), 
	    ifelse(jarque.bera.test(x$residuals)$p.value < 0.01, paste0(as.character(signif(jarque.bera.test(x$residuals)$statistic, digits=5)), "***"), NA)))), error=function(err) NA) }))),
    c("Breusch-Pagan", unlist(lapply(OLS.gis, function(x){ tryCatch(ifelse(bptest(x)$p.value > 0.01, paste0(as.character(signif(bptest(x)$statistic, digits=5)), ""), 
	    ifelse(bptest(x)$p.value > 0.05 & bptest(x)$p.value < 0.01, paste0(as.character(signif(bptest(x)$statistic, digits=5)), "*"), 
	    ifelse(bptest(x)$p.value > 0.01 & bptest(x)$p.value < 0.05, paste0(as.character(signif(bptest(x)$statistic, digits=5)), "**"), 
	    ifelse(bptest(x)$p.value < 0.01, paste0(as.character(signif(bptest(x)$statistic, digits=5)), "***"), NA)))), error=function(err) NA) }))),
    c("Breusch-Godfrey", unlist(lapply(OLS.gis, function(x){ tryCatch(ifelse(bgtest(x)$p.value > 0.01, paste0(as.character(signif(bgtest(x)$statistic, digits=5)), ""), 
	    ifelse(bgtest(x)$p.value > 0.05 & bgtest(x)$p.value < 0.01, paste0(as.character(signif(bgtest(x)$statistic, digits=5)), "*"), 
	    ifelse(bgtest(x)$p.value > 0.01 & bgtest(x)$p.value < 0.05, paste0(as.character(signif(bgtest(x)$statistic, digits=5)), "**"), 
	    ifelse(bgtest(x)$p.value < 0.01, paste0(as.character(signif(bgtest(x)$statistic, digits=5)), "***"), NA)))), error=function(err) NA) }))),
	  c("Durbin-Watson", unlist(lapply(OLS.gis, function(x){ tryCatch(ifelse(dwtest(x)$p.value > 0.01, paste0(as.character(signif(dwtest(x)$statistic, digits=5)), ""), 
	    ifelse(dwtest(x)$p.value > 0.05 & dwtest(x)$p.value < 0.01, paste0(as.character(signif(dwtest(x)$statistic, digits=5)), "*"), 
	    ifelse(dwtest(x)$p.value > 0.01 & dwtest(x)$p.value < 0.05, paste0(as.character(signif(dwtest(x)$statistic, digits=5)), "**"), 
	    ifelse(dwtest(x)$p.value < 0.01, paste0(as.character(signif(dwtest(x)$statistic, digits=5)), "***"), NA)))), error=function(err) NA) }))),
	  c("Mean adj(VIF) of Int.", unlist(lapply(OLS.gis, function(x){ tryCatch(as.character(signif(mean(vif(x)), digits=5)), error=function(err) NA) }))),
	  c("Median adj(VIF) of Int.", unlist(lapply(OLS.gis, function(x){ tryCatch(as.character(signif(median(vif(x)), digits=5)), error=function(err) NA) }))),
	  c("Max adj(VIF) of Int.", unlist(lapply(OLS.gis, function(x){ tryCatch(as.character(signif(max(vif(x)), digits=5)), error=function(err) NA) })))),
  digits = 5, out=paste0(wd, "/OLS_GIS.htm"))
###-------------------------------------------------------------------------------------------------

# LESSON 09: CONTROLLING FOR LOCATION PRICE DIFFERENTIALS
###-------------------------------------------------------------------------------------------------
housing.sample@data$BairroZone <- NA
housing.sample@data$BairroZone <- ifelse(as.character(housing.sample@data$freguesia)=="Sao Paulo" |
  as.character(housing.sample@data$freguesia)=="Santa Catarina" |
  as.character(housing.sample@data$freguesia)=="Merces" |
  as.character(housing.sample@data$freguesia)=="Sao Mamede" |
  as.character(housing.sample@data$freguesia)=="Coracao De Jesus" |
  as.character(housing.sample@data$freguesia)=="Sao Jose" |
  as.character(housing.sample@data$freguesia)=="Pena" |
  as.character(housing.sample@data$freguesia)=="Anjos" |
  as.character(housing.sample@data$freguesia)=="Graca" |
  as.character(housing.sample@data$freguesia)=="Sao Vicente De Fora" |
  as.character(housing.sample@data$freguesia)=="Santa Engracia" |
  as.character(housing.sample@data$freguesia)=="Castelo" |
  as.character(housing.sample@data$freguesia)=="Encarnacao" |
  as.character(housing.sample@data$freguesia)=="Madalena" |
  as.character(housing.sample@data$freguesia)=="Martires" |
  as.character(housing.sample@data$freguesia)=="Sao Cristovao e Sao Lourenco" |
  as.character(housing.sample@data$freguesia)=="Sao Miguel" |
  as.character(housing.sample@data$freguesia)=="Sao Nicolau" |
  as.character(housing.sample@data$freguesia)=="Sacramento" |
  as.character(housing.sample@data$freguesia)=="Santa Justa" |
  as.character(housing.sample@data$freguesia)=="Santiago" |
  as.character(housing.sample@data$freguesia)=="Santo Estavao" |
  as.character(housing.sample@data$freguesia)=="Se" |
  as.character(housing.sample@data$freguesia)=="Socorro", "1oBairro", housing.sample@data$BairroZone)
housing.sample@data$BairroZone <- ifelse(as.character(housing.sample@data$freguesia)=="Santa Isabel" |
  as.character(housing.sample@data$freguesia)=="Lapa" |
  as.character(housing.sample@data$freguesia)=="Santos-o-Velho" |
  as.character(housing.sample@data$freguesia)=="Prazeres" |
  as.character(housing.sample@data$freguesia)=="Santo Condestavel" |
  as.character(housing.sample@data$freguesia)=="Alcantara" |
  as.character(housing.sample@data$freguesia)=="Ajuda" |
  as.character(housing.sample@data$freguesia)=="Sao Francisco Xavier" |
  as.character(housing.sample@data$freguesia)=="Santa Maria De Belem", "2oBairro", housing.sample@data$BairroZone)
housing.sample@data$BairroZone <- ifelse(as.character(housing.sample@data$freguesia)=="Benfica" |
  as.character(housing.sample@data$freguesia)=="Campolide" |
  as.character(housing.sample@data$freguesia)=="Sao Sebastiao Da Pedreira" |
  as.character(housing.sample@data$freguesia)=="Nossa Senhora De Fatima" |
  as.character(housing.sample@data$freguesia)=="Alvalade" |
  as.character(housing.sample@data$freguesia)=="Sao Joao De Brito" |
  as.character(housing.sample@data$freguesia)=="Campo Grande" |
  as.character(housing.sample@data$freguesia)=="Sao Domingos De Benfica" |
  as.character(housing.sample@data$freguesia)=="Carnide" |
  as.character(housing.sample@data$freguesia)=="Lumiar" |
  as.character(housing.sample@data$freguesia)=="Charneca" |
  as.character(housing.sample@data$freguesia)=="Ameixoeira", "3oBairro", housing.sample@data$BairroZone)
housing.sample@data$BairroZone <- ifelse(as.character(housing.sample@data$freguesia)=="Santa Maria Dos Olivais" |
  as.character(housing.sample@data$freguesia)=="Marvila" |
  as.character(housing.sample@data$freguesia)=="Beato" |
  as.character(housing.sample@data$freguesia)=="Sao Joao" |
  as.character(housing.sample@data$freguesia)=="Penha De Franca" |
  as.character(housing.sample@data$freguesia)=="Sao Jorge De Arroios" |
  as.character(housing.sample@data$freguesia)=="Sao Joao De Deus" |
  as.character(housing.sample@data$freguesia)=="Alto Do Pina", "4oBairro", housing.sample@data$BairroZone)
housing.sample@data$BairroZone <- as.factor(housing.sample@data$BairroZone)
housing.sample@data$BairroZone <- relevel(housing.sample@data$BairroZone, ref="1oBairro")

OLS.spfe <- list()
OLS.spfe[[1]] <- lm(ln.price ~ area + area.2 + new + parking + fireplace + aircond + view + Ppop.unemp + ln.dist.baixa + BairroZone, na.action=na.omit, data=housing.sample@data)
OLS.spfe[[2]] <- lm(ln.price ~ area + area.2 + new + parking + fireplace + aircond + view + Ppop.unemp + ln.dist.tagus + BairroZone, na.action=na.omit, data=housing.sample@data)

OLS.spfe[[3]] <- lm(ln.price ~ area + area.2 + new + parking + fireplace + aircond + view + Ppop.unemp + ln.dist.baixa + ln.avgD.park + BairroZone, na.action=na.omit, data=housing.sample@data)
OLS.spfe[[4]] <- lm(ln.price ~ area + area.2 + new + parking + fireplace + aircond + view + Ppop.unemp + ln.dist.tagus + ln.avgD.park + BairroZone, na.action=na.omit, data=housing.sample@data)

OLS.spfe[[5]] <- lm(ln.price ~ area + area.2 + new + parking + fireplace + aircond + view + Ppop.unemp + ln.dist.baixa + lnPM10 + ln.avgD.park + BairroZone, na.action=na.omit, data=housing.sample@data)
OLS.spfe[[6]] <- lm(ln.price ~ area + area.2 + new + parking + fireplace + aircond + view + Ppop.unemp + ln.dist.tagus + lnPM10 + ln.avgD.park + BairroZone, na.action=na.omit, data=housing.sample@data)

stargazer(OLS.spfe, type="html", title="OLS Models with Location FE (Robust S.E.)", style="aer", 
  star.cutoffs=c(0.10, 0.05, 0.01), se = lapply(OLS.spfe, function(y){ coeftest(y, function(x) vcovHC(x, type="HC0"))[,2] }),
  add.lines = list(c("AIC", unlist(lapply(OLS.spfe, function(x){ tryCatch(as.character(signif(AIC(x), digits=5)), error=function(err) NA) }))),
    c("SSE", unlist(lapply(OLS.spfe, function(x){ tryCatch(as.character(signif(tail(anova(x)[,2], 1), digits=5)), error=function(err) NA) }))),
    c("Jarque-Bera", unlist(lapply(OLS.spfe, function(x){ tryCatch(ifelse(jarque.bera.test(x$residuals)$p.value > 0.01, paste0(as.character(signif(jarque.bera.test(x$residuals)$statistic, digits=5)), ""), 
      ifelse(jarque.bera.test(x$residuals)$p.value > 0.05 & jarque.bera.test(x$residuals)$p.value < 0.01, paste0(as.character(signif(jarque.bera.test(x$residuals)$statistic, digits=5)), "*"), 
	    ifelse(jarque.bera.test(x$residuals)$p.value > 0.01 & jarque.bera.test(x$residuals)$p.value < 0.05, paste0(as.character(signif(jarque.bera.test(x$residuals)$statistic, digits=5)), "**"), 
	    ifelse(jarque.bera.test(x$residuals)$p.value < 0.01, paste0(as.character(signif(jarque.bera.test(x$residuals)$statistic, digits=5)), "***"), NA)))), error=function(err) NA) }))),
    c("Breusch-Pagan", unlist(lapply(OLS.spfe, function(x){ tryCatch(ifelse(bptest(x)$p.value > 0.01, paste0(as.character(signif(bptest(x)$statistic, digits=5)), ""), 
	    ifelse(bptest(x)$p.value > 0.05 & bptest(x)$p.value < 0.01, paste0(as.character(signif(bptest(x)$statistic, digits=5)), "*"), 
	    ifelse(bptest(x)$p.value > 0.01 & bptest(x)$p.value < 0.05, paste0(as.character(signif(bptest(x)$statistic, digits=5)), "**"), 
	    ifelse(bptest(x)$p.value < 0.01, paste0(as.character(signif(bptest(x)$statistic, digits=5)), "***"), NA)))), error=function(err) NA) }))),
    c("Breusch-Godfrey", unlist(lapply(OLS.spfe, function(x){ tryCatch(ifelse(bgtest(x)$p.value > 0.01, paste0(as.character(signif(bgtest(x)$statistic, digits=5)), ""), 
	    ifelse(bgtest(x)$p.value > 0.05 & bgtest(x)$p.value < 0.01, paste0(as.character(signif(bgtest(x)$statistic, digits=5)), "*"), 
	    ifelse(bgtest(x)$p.value > 0.01 & bgtest(x)$p.value < 0.05, paste0(as.character(signif(bgtest(x)$statistic, digits=5)), "**"), 
	    ifelse(bgtest(x)$p.value < 0.01, paste0(as.character(signif(bgtest(x)$statistic, digits=5)), "***"), NA)))), error=function(err) NA) }))),
	  c("Durbin-Watson", unlist(lapply(OLS.spfe, function(x){ tryCatch(ifelse(dwtest(x)$p.value > 0.01, paste0(as.character(signif(dwtest(x)$statistic, digits=5)), ""), 
	    ifelse(dwtest(x)$p.value > 0.05 & dwtest(x)$p.value < 0.01, paste0(as.character(signif(dwtest(x)$statistic, digits=5)), "*"), 
	    ifelse(dwtest(x)$p.value > 0.01 & dwtest(x)$p.value < 0.05, paste0(as.character(signif(dwtest(x)$statistic, digits=5)), "**"), 
	    ifelse(dwtest(x)$p.value < 0.01, paste0(as.character(signif(dwtest(x)$statistic, digits=5)), "***"), NA)))), error=function(err) NA) }))),
	  c("Mean adj(VIF) of Int.", unlist(lapply(OLS.spfe, function(x){ tryCatch(as.character(signif(mean(vif(x)), digits=5)), error=function(err) NA) }))),
	  c("Median adj(VIF) of Int.", unlist(lapply(OLS.spfe, function(x){ tryCatch(as.character(signif(median(vif(x)), digits=5)), error=function(err) NA) }))),
	  c("Max adj(VIF) of Int.", unlist(lapply(OLS.spfe, function(x){ tryCatch(as.character(signif(max(vif(x)), digits=5)), error=function(err) NA) })))),
  digits = 5, out=paste0(wd, "/OLS_SpFE.htm"))
###-------------------------------------------------------------------------------------------------

# LESSON 10: SPATIAL DEPENDANCE - GENERATING AND VISUALIZING SPATIAL WEIGHTS
###-------------------------------------------------------------------------------------------------
dist0.5 <- dnearneigh(coordinates(housing.sample), d1=0, d2=0.5, longlat=TRUE)
neigh.dist <- nbdists(dist0.5, coordinates(housing.sample), longlat=TRUE)

SW1 <- nb2listw(dist0.5, style="W", zero.policy=TRUE)
SW2 <- nb2listw(dist0.5, glist = lapply(neigh.dist, function(x) (1/x)), style="W", zero.policy=TRUE)
SW3 <- nb2listw(knn2nb(knearneigh(coordinates(housing.sample), k=20)), zero.policy=TRUE)

summary(SW1, zero.policy=TRUE)
summary(SW2, zero.policy=TRUE)
summary(SW3, zero.policy=TRUE)

dist0.9 <- dnearneigh(coordinates(housing.sample), d1=0, d2=0.9, longlat=TRUE)
neigh.dist <- nbdists(dist0.9, coordinates(housing.sample), longlat=TRUE)

SW4 <- nb2listw(dist0.9, style="W", zero.policy=TRUE)
summary(SW4, zero.policy=TRUE)

freguesias <- spTransform(readOGR(dsn = paste0(wd, "/Data"), layer = "freguesias"), CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
plot(freguesias)
plot(SW1, coordinates(housing.sample), add=T, col="grey")
plot(freguesias)
plot(SW3, coordinates(housing.sample), add=T, col="grey")


SW1.vector <- listw2lines(SW1, coordinates(housing.sample))
writeOGR(SW1.vector, dsn = paste0(wd, "/Data"), layer = "SW1_lines", driver="ESRI Shapefile")


###-------------------------------------------------------------------------------------------------

# LESSON 11: SPATIAL DIAGNOSTICS AND MODEL SELECTION
###-------------------------------------------------------------------------------------------------
global.model <- OLS.spfe[[2]]

moran.test(housing.sample@data$ln.price, listw=SW1, alternative="greater", zero.policy=TRUE)
lm.morantest(global.model, listw=SW1, alternative="greater", zero.policy=TRUE)

moran.test(housing.sample@data$ln.price, listw=SW2, alternative="greater", zero.policy=TRUE)
lm.morantest(global.model, listw=SW2, alternative="greater", zero.policy=TRUE)

moran.test(housing.sample@data$ln.price, listw=SW3, alternative="greater", zero.policy=TRUE)
lm.morantest(global.model, listw=SW3, alternative="greater", zero.policy=TRUE)

moran.test(housing.sample@data$ln.price, listw=SW4, alternative="greater", zero.policy=TRUE)
lm.morantest(global.model, listw=SW4, alternative="greater", zero.policy=TRUE)

moran.plot(global.model$residuals, SW3, zero.policy=TRUE)
moran.plot(housing.sample@data$ln.price, SW3, zero.policy=TRUE)

lm.LMtests(global.model, listw=SW1, test="all", zero.policy=TRUE)
lm.LMtests(global.model, listw=SW2, test="all", zero.policy=TRUE)

lm.LMtests(global.model, listw=SW3, test="all", zero.policy=TRUE)
lm.LMtests(global.model, listw=SW4, test="all", zero.policy=TRUE)

###-------------------------------------------------------------------------------------------------

# LESSON 12: MODELS OF SPATIAL DEPENDENCE (SPATIAL LAG, ERROR)
###-------------------------------------------------------------------------------------------------
SEM.SW3 <- errorsarlm(global.model, data=housing.sample, listw=SW3, method="MC", zero.policy=TRUE, na.action=na.omit)
bptest.sarlm(SEM.SW3)
SEM.SW3_SE <- coeftest(lm(SEM.SW3$tary ~ SEM.SW3$tarX - 1), vcov=vcovHC(lm(SEM.SW3$tary ~ SEM.SW3$tarX - 1), type="HC0"), df=Inf)
rownames(SEM.SW3_SE) <- gsub("SEM.SW3\\$tarXI\\(x - lambda \\* WX\\)", "", rownames(SEM.SW3_SE))
summary(SEM.SW3)

SEM.SW3$SSE
SEM.SW3$s2
SEM.SW3$LL
SEM.SW3$logLik_lm.model
AIC(SEM.SW3)
SEM.SW3$AIC_lm.model


SAR.SW3 <- lagsarlm(global.model, data=housing.sample, listw=SW3, method="MC", zero.policy=TRUE, na.action=na.omit)
bptest.sarlm(SAR.SW3)
SAR.SW3_SE <- coeftest(lm(SAR.SW3$tary ~ SAR.SW3$tarX - 1), vcov=vcovHC(lm(SAR.SW3$tary ~ SAR.SW3$tarX - 1), type="HC0"), df=Inf)
rownames(SAR.SW3_SE) <- gsub("SAR.SW3\\$tarXx", "", rownames(SAR.SW3_SE))
summary(SAR.SW3)
impacts(SAR.SW3, tr=trW(forceSymmetric(as(SW3, "CsparseMatrix")), m=50, p=100, type="MC"), zstats=T)

SAR.SW3$SSE
SAR.SW3$s2
SAR.SW3$LL
SAR.SW3$logLik_lm.model
AIC(SAR.SW3)
SAR.SW3$AIC_lm.model


SARAR.SW3 <- sacsarlm(global.model, data=housing.sample, listw=SW3, method="MC", zero.policy=TRUE, na.action=na.omit) 
bptest.sarlm(SARAR.SW3)
SARAR.SW3_SE <- coeftest(lm(SARAR.SW3$tary ~ SARAR.SW3$tarX - 1), vcov=vcovHC(lm(SARAR.SW3$tary ~ SARAR.SW3$tarX - 1), type="HC0"), df=Inf)
summary(SARAR.SW3)
impacts(SARAR.SW3, tr=trW(forceSymmetric(as(SW3, "CsparseMatrix")), m=50, p=100, type="MC"), zstats=T)

SARAR.SW3$SSE
SARAR.SW3$s2
SARAR.SW3$LL
SARAR.SW3$logLik_lm.model
AIC(SARAR.SW3)
SARAR.SW3$AIC_lm.model


stargazer(list(global.model, SEM.SW3, SAR.SW3), model.numbers=FALSE, type="html", title="Spatial Results (Robust S.E.)", 
  style="aer", dep.var.labels.include=FALSE, df = FALSE, star.cutoffs=c(0.10, 0.05, 0.01),
  se = list(coeftest(global.model, function(x) vcovHC(x, type="HC0"))[,2], SEM.SW3_SE[,2], SAR.SW3_SE[,2]),  
  digits = 5, out=paste0(wd, "/SP_Models.htm"))


residuals <- housing.sample[, c("latitude", "longitude")]
residuals$residuals.olsFE <- global.model$residuals
residuals$residuals.ols <- OLS.gis[[5]]$residuals
residuals$residuals.semSW3 <- SEM.SW3$residuals
residuals$residuals.sarSW3 <- SAR.SW3$residuals
residuals$residuals.sararSW3 <- SARAR.SW3$residuals

ggplot() + geom_point(data=residuals@data, aes(x=longitude, y=latitude, col = residuals.ols)) +
  scale_colour_gradient2(low = "red", mid = "white", high = "green", midpoint = mean(residuals@data$residuals.ols))
ggplot() + geom_point(data=residuals@data, aes(x=longitude, y=latitude, col = residuals.olsFE)) +
  scale_colour_gradient2(low = "red", mid = "white", high = "green", midpoint = mean(residuals@data$residuals.olsFE))
ggplot() + geom_point(data=residuals@data, aes(x=longitude, y=latitude, col = residuals.sarSW3)) +
  scale_colour_gradient2(low = "red", mid = "white", high = "green", midpoint = mean(residuals@data$residuals.sarSW3))
###-------------------------------------------------------------------------------------------------



# APPENDIX: CODE FOR GENERATING GIS VARIABLES IN R 
#    - Requires downloading additional open-source data
###-------------------------------------------------------------------------------------------------
# census.freguesias <- spTransform(readOGR(dsn = paste0(wd, "/Data"), layer = "freguesias"), CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# 
# census.freguesias@data$BAIXA <- ifelse(as.character(census.freguesias@data$FREGUESIA)=="Santa Justa" |
#   as.character(census.freguesias@data$FREGUESIA)=="Sao Nicolau" |
#   as.character(census.freguesias@data$FREGUESIA)=="Se" |
#   as.character(census.freguesias@data$FREGUESIA)=="Sacramento" |
#   as.character(census.freguesias@data$FREGUESIA)=="Sao Cristovao e Sao Lourenco" |
#   as.character(census.freguesias@data$FREGUESIA)=="Madalena" |
#   as.character(census.freguesias@data$FREGUESIA)=="Encarnacao" |
#   as.character(census.freguesias@data$FREGUESIA)=="Castelo" |
#   as.character(census.freguesias@data$FREGUESIA)=="Santo Estavao" |
#   as.character(census.freguesias@data$FREGUESIA)=="Martires" |
#   as.character(census.freguesias@data$FREGUESIA)=="Santiago" |
#   as.character(census.freguesias@data$FREGUESIA)=="Graca" |
#   as.character(census.freguesias@data$FREGUESIA)=="Sao Miguel" |
#   as.character(census.freguesias@data$FREGUESIA)=="Sao Vicente De Fora" |
#   as.character(census.freguesias@data$FREGUESIA)=="Socorro", 1, 0)
# 
# plot(census.freguesias, col = "light grey")
# plot(subset(census.freguesias, census.freguesias@data$BAIXA==1), add = T, col = "dark grey")
# plot(housing.sample, add = T, pch = 20, size = 0.3)
# 
# census.data <- read.csv(paste0(wd, "/Data/Census Freguesia.csv"), fileEncoding="latin1", header = T, sep = ";")
# summary(census.data)
# head(census.data)
# 
# ROWS <- grepl("'1106", census.data$GEO_COD, fixed=TRUE) & census.data$NIVEL_DSG=="Freguesia"
# COLUMNS <- colnames(census.data) %in% c("ANO", "GEO_COD", "N_EDIFICIOS_CLASSICOS",
#   "N_EDIFICIOS_EXCLUSIV_RESID", "N_EDIFICIOS_PRINCIP_NAO_RESID", "N_INDIVIDUOS_RESIDENT",
#   "N_INDIVIDUOS_RESIDENT_65", "N_IND_RESID_DESEMP_PROC_1EMPRG", "N_IND_RESID_DESEMP_PROC_EMPRG",
#   "N_IND_RESID_EMPREGADOS", "N_IND_RESIDENT_ENSINCOMP_SEC", "N_IND_RESIDENT_ENSINCOMP_POSEC",
#   "N_IND_RESIDENT_ENSINCOMP_SUP")
# 
# census.data <- census.data[ROWS , COLUMNS]
# rm(ROWS, COLUMNS)
# 
# colnames(census.data)[colnames(census.data)=="ANO"] <- "YEAR"
# colnames(census.data)[colnames(census.data)=="GEO_COD"] <- "ID"
# colnames(census.data)[colnames(census.data)=="N_EDIFICIOS_CLASSICOS"] <- "Nbldg"
# colnames(census.data)[colnames(census.data)=="N_EDIFICIOS_EXCLUSIV_RESID"] <- "Nbldg.res"
# colnames(census.data)[colnames(census.data)=="N_EDIFICIOS_PRINCIP_NAO_RESID"] <- "Nbldg.non"
# colnames(census.data)[colnames(census.data)=="N_INDIVIDUOS_RESIDENT"] <- "Npop"
# colnames(census.data)[colnames(census.data)=="N_INDIVIDUOS_RESIDENT_65"] <- "Npop.65"
# colnames(census.data)[colnames(census.data)=="N_IND_RESIDENT_ENSINCOMP_SEC"] <- "Npop.secondary"
# colnames(census.data)[colnames(census.data)=="N_IND_RESIDENT_ENSINCOMP_POSEC"] <- "Npop.Psecondary"
# colnames(census.data)[colnames(census.data)=="N_IND_RESIDENT_ENSINCOMP_SUP"] <- "Npop.superior"
# colnames(census.data)[colnames(census.data)=="N_IND_RESID_DESEMP_PROC_1EMPRG"] <- "Npop.unemp1"
# colnames(census.data)[colnames(census.data)=="N_IND_RESID_DESEMP_PROC_EMPRG"] <- "Npop.unemp"
# colnames(census.data)[colnames(census.data)=="N_IND_RESID_EMPREGADOS"] <- "Npop.empl"
# 
# census.data$Pbldg.res <- census.data$Nbldg.res/census.data$Nbldg
# census.data$Pbldg.non <- census.data$Nbldg.non/census.data$Nbldg
# 
# census.data$Ppop.65 <- census.data$Npop.65/census.data$Npop
# 
# census.data$Ppop.secondary <- census.data$Npop.secondary/census.data$Npop
# census.data$Ppop.higher <- (census.data$Npop.Psecondary + census.data$Npop.superior)/census.data$Npop
# 
# census.data$Ppop.empl <- census.data$Npop.empl/census.data$Npop
# census.data$Ppop.unemp <- (census.data$Npop.unemp + census.data$Npop.unemp1)/census.data$Npop
# 
# census.data$ID <- as.factor(as.character(gsub("\\'", "", census.data$ID)))
# 
# census.freguesias <- merge(census.freguesias, census.data, by.x="DICOFRE", by.y="ID", all.x=T, sort=F)
# rm(census.data)
# 
# census.freguesias@data$pop.ha <- census.freguesias@data$Npop/census.freguesias@data$AREA_T_HA
# census.freguesias@data$bldg.ha <- census.freguesias@data$Nbldg/census.freguesias@data$AREA_T_HA
# 
# proj4string(housing.sample)
# proj4string(census.freguesias)
# 
# housing.sample@data$baixa <- over(housing.sample, census.freguesias)$BAIXA
# housing.sample@data$freguesia <- over(housing.sample, census.freguesias)$FREGUESIA
# housing.sample@data$Pbldg.res <- over(housing.sample, census.freguesias)$Pbldg.res
# housing.sample@data$Pbldg.non <- over(housing.sample, census.freguesias)$Pbldg.non
# housing.sample@data$Ppop.65 <- over(housing.sample, census.freguesias)$Ppop.65
# housing.sample@data$Ppop.secondary <- over(housing.sample, census.freguesias)$Ppop.secondary
# housing.sample@data$Ppop.higher <- over(housing.sample, census.freguesias)$Ppop.higher
# housing.sample@data$Ppop.empl <- over(housing.sample, census.freguesias)$Ppop.empl
# housing.sample@data$Ppop.unemp <- over(housing.sample, census.freguesias)$Ppop.unemp
# housing.sample@data$pop.ha <- over(housing.sample, census.freguesias)$pop.ha
# housing.sample@data$bldg.ha <- over(housing.sample, census.freguesias)$bldg.ha
# 
# baixa <- gUnaryUnion(subset(census.freguesias, census.freguesias@data$BAIXA==1))
# plot(baixa, add=T, col="blue")
# 
# proj4string(housing.sample)
# proj4string(baixa)
# 
# housing.sample@data$dist.baixa <- spDistsN1(coordinates(housing.sample), coordinates(gCentroid(baixa)), longlat=TRUE)
# summary(housing.sample@data$dist.baixa)
# 
# parks <- spTransform(readOGR(dsn = paste0(wd, "/Data"), layer = "Jardins_Parques_Urbanos"), CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# shopping.center <- spTransform(readOGR(dsn = paste0(wd, "/Data"), layer = "Centros_Comerciais"), CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# 
# plot(census.freguesias, col = "light grey")
# plot(baixa, add = T, col = "dark grey")
# plot(housing.sample, add = T, pch = 23, size = 0.3)
# plot(parks, add = T, pch = 20, size = 0.2, col="springgreen3")
# plot(shopping.center, add = T, pch = 20, size = 0.2, col="purple")
# 
# proj4string(housing.sample)
# proj4string(shopping.center)
# proj4string(parks)
# 
# parks.distances <- spDists(coordinates(housing.sample), coordinates(parks), longlat=TRUE)
# shopping.distances <- spDists(coordinates(housing.sample), coordinates(shopping.center), longlat=TRUE)
# housing.sample@data$near.park <- apply(parks.distances, 1, min)
# housing.sample@data$avgD.park <- apply(parks.distances, 1, mean)
# housing.sample@data$near.shopping <- apply(shopping.distances, 1, min)
# housing.sample@data$avgD.shopping <- apply(shopping.distances, 1, mean)
# 
# housing.sample@data$park.100 <- apply(parks.distances, 1, function(x) sum(ifelse(x <= 0.1, 1, 0)))
# housing.sample@data$park.500 <- apply(parks.distances, 1, function(x) sum(ifelse(x <= 0.5, 1, 0)))
# 
# housing.sample@data$shopping.100 <- apply(shopping.distances, 1, function(x) sum(ifelse(x <= 0.1, 1, 0)))
# housing.sample@data$shopping.500 <- apply(shopping.distances, 1, function(x) sum(ifelse(x <= 0.5, 1, 0)))
# 
# summary(housing.sample@data$near.park)
# summary(housing.sample@data$avgD.park)
# summary(housing.sample@data$near.shopping)
# summary(housing.sample@data$avgD.shopping)
# 
# summary(housing.sample@data$park.100)
# summary(housing.sample@data$park.500)
# summary(housing.sample@data$shopping.100)
# summary(housing.sample@data$shopping.500)
# 
# tagus <- spTransform(readOGR(dsn = paste0(wd, "/Data"), layer = "Tagus_Riverfront"), CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# 
# tagus.distances <- as.data.frame(dist2Line(housing.sample, tagus, distfun=distHaversine))
# housing.sample@data$dist.tagus1 <- tagus.distances$distance
# summary(housing.sample@data$dist.tagus1)
# housing.sample@data$dist.tagus1 <- housing.sample@data$dist.tagus1/1000
# 
# tagus <- gLineMerge(tagus)
# tagus <- spsample(tagus, n=round(SpatialLinesLengths(tagus), digits=0), type="regular")
# 
# tagus.distances <- spDists(coordinates(housing.sample), coordinates(tagus), longlat=TRUE)
# housing.sample@data$dist.tagus2 <- apply(tagus.distances, 1, min)
# summary(housing.sample@data$dist.tagus2)
# 
# rm(tagus, tagus.distances)
###-------------------------------------------------------------------------------------------------














OLS.spfe <- list()
OLS.spfe[[1]] <- lm(ln.price ~ area + area.2 + new + parking + fireplace + aircond + view + Ppop.unemp + ln.dist.baixa + BairroZone, na.action=na.omit, data=housing.sample@data)
OLS.spfe[[2]] <- lm(ln.price ~ ln.area + new + parking + fireplace + aircond + view + Ppop.unemp + ln.dist.baixa + BairroZone, na.action=na.omit, data=housing.sample@data)
OLS.spfe[[3]] <- lm(ln.price ~ area + area.2 + new + parking + fireplace + aircond + view + Ppop.unemp + ln.dist.baixa, na.action=na.omit, data=housing.sample@data)
OLS.spfe[[4]] <- lm(ln.price ~ ln.area + new + parking + fireplace + aircond + view + Ppop.unemp + ln.dist.baixa, na.action=na.omit, data=housing.sample@data)
OLS.spfe[[5]] <- lm(ln.price ~ area + area.2 + new + parking + fireplace + aircond + view + Ppop.unemp + baixa + BairroZone, na.action=na.omit, data=housing.sample@data)
OLS.spfe[[6]] <- lm(ln.price ~ ln.area + new + parking + fireplace + aircond + view + Ppop.unemp + baixa + BairroZone, na.action=na.omit, data=housing.sample@data)
OLS.spfe[[7]] <- lm(ln.price ~ area + area.2 + new + parking + fireplace + aircond + view + Ppop.unemp + baixa, na.action=na.omit, data=housing.sample@data)
OLS.spfe[[8]] <- lm(ln.price ~ ln.area + new + parking + fireplace + aircond + view + Ppop.unemp + baixa, na.action=na.omit, data=housing.sample@data)


stargazer(OLS.spfe, type="html", title="OLS Models (Robust S.E.)", style="aer", 
  star.cutoffs=c(0.10, 0.05, 0.01), se = lapply(OLS.spfe, function(y){ coeftest(y, function(x) vcovHC(x, type="HC0"))[,2] }),
  add.lines = list(c("AIC", unlist(lapply(OLS.spfe, function(x){ tryCatch(as.character(signif(AIC(x), digits=5)), error=function(err) NA) }))),
    c("SSE", unlist(lapply(OLS.spfe, function(x){ tryCatch(as.character(signif(tail(anova(x)[,2], 1), digits=5)), error=function(err) NA) }))),
    c("Jarque-Bera", unlist(lapply(OLS.spfe, function(x){ tryCatch(ifelse(jarque.bera.test(x$residuals)$p.value > 0.01, paste0(as.character(signif(jarque.bera.test(x$residuals)$statistic, digits=5)), ""), 
      ifelse(jarque.bera.test(x$residuals)$p.value > 0.05 & jarque.bera.test(x$residuals)$p.value < 0.01, paste0(as.character(signif(jarque.bera.test(x$residuals)$statistic, digits=5)), "*"), 
	    ifelse(jarque.bera.test(x$residuals)$p.value > 0.01 & jarque.bera.test(x$residuals)$p.value < 0.05, paste0(as.character(signif(jarque.bera.test(x$residuals)$statistic, digits=5)), "**"), 
	    ifelse(jarque.bera.test(x$residuals)$p.value < 0.01, paste0(as.character(signif(jarque.bera.test(x$residuals)$statistic, digits=5)), "***"), NA)))), error=function(err) NA) }))),
    c("Breusch-Pagan", unlist(lapply(OLS.spfe, function(x){ tryCatch(ifelse(bptest(x)$p.value > 0.01, paste0(as.character(signif(bptest(x)$statistic, digits=5)), ""), 
	    ifelse(bptest(x)$p.value > 0.05 & bptest(x)$p.value < 0.01, paste0(as.character(signif(bptest(x)$statistic, digits=5)), "*"), 
	    ifelse(bptest(x)$p.value > 0.01 & bptest(x)$p.value < 0.05, paste0(as.character(signif(bptest(x)$statistic, digits=5)), "**"), 
	    ifelse(bptest(x)$p.value < 0.01, paste0(as.character(signif(bptest(x)$statistic, digits=5)), "***"), NA)))), error=function(err) NA) }))),
    c("Breusch-Godfrey", unlist(lapply(OLS.spfe, function(x){ tryCatch(ifelse(bgtest(x)$p.value > 0.01, paste0(as.character(signif(bgtest(x)$statistic, digits=5)), ""), 
	    ifelse(bgtest(x)$p.value > 0.05 & bgtest(x)$p.value < 0.01, paste0(as.character(signif(bgtest(x)$statistic, digits=5)), "*"), 
	    ifelse(bgtest(x)$p.value > 0.01 & bgtest(x)$p.value < 0.05, paste0(as.character(signif(bgtest(x)$statistic, digits=5)), "**"), 
	    ifelse(bgtest(x)$p.value < 0.01, paste0(as.character(signif(bgtest(x)$statistic, digits=5)), "***"), NA)))), error=function(err) NA) }))),
	  c("Durbin-Watson", unlist(lapply(OLS.spfe, function(x){ tryCatch(ifelse(dwtest(x)$p.value > 0.01, paste0(as.character(signif(dwtest(x)$statistic, digits=5)), ""), 
	    ifelse(dwtest(x)$p.value > 0.05 & dwtest(x)$p.value < 0.01, paste0(as.character(signif(dwtest(x)$statistic, digits=5)), "*"), 
	    ifelse(dwtest(x)$p.value > 0.01 & dwtest(x)$p.value < 0.05, paste0(as.character(signif(dwtest(x)$statistic, digits=5)), "**"), 
	    ifelse(dwtest(x)$p.value < 0.01, paste0(as.character(signif(dwtest(x)$statistic, digits=5)), "***"), NA)))), error=function(err) NA) }))),
	  c("Mean adj(VIF)", unlist(lapply(OLS.spfe, function(x){ tryCatch(as.character(signif(mean(vif(x)), digits=5)), error=function(err) NA) }))),
	  c("Median adj(VIF)", unlist(lapply(OLS.spfe, function(x){ tryCatch(as.character(signif(median(vif(x)), digits=5)), error=function(err) NA) }))),
	  c("Max adj(VIF)", unlist(lapply(OLS.spfe, function(x){ tryCatch(as.character(signif(max(vif(x)), digits=5)), error=function(err) NA) })))),
  digits = 5, out=paste0(wd, "/OLS_Final.htm"))
###-------------------------------------------------------------------------------------------------


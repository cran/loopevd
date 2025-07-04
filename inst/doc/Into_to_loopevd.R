## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(loopevd)
library(ncdf4)      
library(terra)
library(raster)
library(sp)
library(ozmaps)
library(evd)
library(zoo)


## -----------------------------------------------------------------------------
cp = colorRampPalette(c("light yellow","yellow","red","brown4"))
rwb = colorRampPalette(c("red","white","blue"))

coastsp = methods::as(ozmaps::ozmap_data("country", quiet =TRUE),"Spatial")
coastp = vect(ozmaps::ozmap_data("country", quiet =TRUE))
coast = as.lines(coastp)
coast.layer = list("sp.lines",methods::as(coast,"Spatial"),col="black")

## -----------------------------------------------------------------------------
r = terra::rast(system.file("extdata/50km_AnnMax_agcd_v1_tmax_mean_r005_daily_1980-2019.nc"
                            ,package = "loopevd"))
#r = terra::rast("../inst/extdata/50km_AnnMax_agcd_v1_tmax_mean_r005_daily_1980-2019.nc")

## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
secondHighest_r = app(r,function(x) sort(x, decreasing =TRUE)[2])
at=c(-Inf,seq(27.5,50,2.5),Inf) #the Inf bookends add the hats/trinagles to the colour scale plot 
sp::spplot(secondHighest_r,at=at,col.regions =c("#FFFFFF",cp(9),"#000000"),main="2nd highest value in 40 year dataset")

## ----include=FALSE------------------------------------------------------------

system.time(gumbel_r <- raster_fevd(r,"fgumbel",silent=TRUE,seed=1))

system.time(gev_r <- raster_fevd(r,"fgev",silent=TRUE,seed=1))

system.time(try(gumbelx_r <- raster_fevd(r,"fgumbelx",silent=TRUE,ntries = 3,seed=1),silent=TRUE))
gev_r

## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
gev_rl20 = raster_qevd(x = gev_r,p = 1-.05, evd_mod_str = "fgev") # terra::app(x = gev_r,fun = qevd_vector,p = 1-.05, evd_mod_str = "fgev")
at=c(-Inf,seq(27.5,50,2.5),Inf)
cp = colorRampPalette(c("light yellow","yellow","red","brown4"))
#pdf("../plots/Empirical_20_year_ARI_daily_max_temp.pdf")
sp::spplot(gev_rl20,at=at,col.regions =c("#FFFFFF",cp(9),"#000000"),main="GEV 20 year ARI")

## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
r_years = as.numeric(format(time(r),"%Y"))

#https://psl.noaa.gov/data/climateindices/list/
#dat = read.table("https://psl.noaa.gov/data/correlation/gmsst.data",skip=1,fill = TRUE,header = FALSE)
#dat = read.table("https://psl.noaa.gov/data/correlation/soi.data",skip=1,fill = TRUE,header = FALSE)
#n = dim(dat)[1]
#dat = cbind(dat[,1][-n],dat[,8:13][-n,],dat[,2:7][-1,]) 
#dat[dat == -99.99] = NA
#CI = suppressWarnings(data.frame(year = dat[,1],yearMean = apply(dat[,-1],MARGIN = 1,FUN = function(x) mean(as.numeric(x),na.rm=TRUE))))
#https://www.ipcc.ch/sr15/faq/faq-chapter-1/ 
#CI = CI[CI$year >= min(r_years) & CI$year <= max(r_years) & !is.na(CI$year),]
#CI$yearMeanRoll = zoo::rollmean(CI$yearMean,20,fill=NA)
#nsloc = data.frame(yearMean = CI$yearMean) #climate index

nt = length(time(r))
nsloc = data.frame(year = seq(-1, 1, length.out = nt)) # linear trend

nsloc  = centredAndScaled(data.frame(year = r_years))
nsloctab = data.frame(year = r_years, year_cs = nsloc[,1])

plot(r_years,nsloc[,1],xlab = "Year",ylab = "Nonstationary Location covariate")

zero_year = approx(nsloctab[,2],r_years,0)$y
abline(v = zero_year)
grid()

#print("gridded non stationary gumbel fit:")
system.time(gumbelns_r <- suppressWarnings(raster_fevd(r,"fgumbel",nsloc = nsloc,seed=1)))

#print("gridded non stationary gev fit:")
system.time(gevns_r <- suppressWarnings(raster_fevd(r,"fgev",nsloc = nsloc,ntries = 3,seed=1)))


## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
AICs_r = c(gumbel_r$AIC,gev_r$AIC,gumbelns_r$AIC,gevns_r$AIC,gumbelx_r$AIC) 
best = app(AICs_r,which.min,na.rm=TRUE)
plot(best,col = c("black","red","grey","yellow","green"),plg=list(legend=c("gumbel", "gev", "gumbel_ns","gev_ns","gumbelx")),main = "Which EVD is best (lowest AIC)")
lines(coastp)
#writeRaster(best,paste0(datadir,"best_EVD.tif"))

## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
#find where gev shape is significant, defined as (1 - p-value) × 100

mle_sig = raster_se_sig(c(gevns_r$shape,gevns_r$cov_16))

at = c(-Inf,seq(50,100,5))

#where mle(par) +-1.95 se(par) is entirely on one side of zero"
sp::spplot(mle_sig,at=at,col.regions =rev(rwb(20)),main = "Confidence Value For a GEV shape parameter in a non-stationary EVD",
           sp.layout = coast.layer)


## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
#find where gev non-stationary location is significant, defined as (1 - p-value) × 100

mle_sig = raster_se_sig(c(gevns_r$locyear,gevns_r$cov_6))

at = c(-Inf,seq(50,100,5))

#where mle(par) +-1.95 se(par) is entirely on one side of zero"
sp::spplot(mle_sig,at=at,col.regions =rev(rwb(14)),main = "Confidence Value For a non-stationary covariate location (locyear) parameter",
           sp.layout = coast.layer)


## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
AEP_year = 1981
gevns_RL5pc_1981 = raster_qevd(x = gev_r,p = 1-.05, evd_mod_str = "fgev") + nsloctab$year_cs[nsloctab$year==AEP_year]*gevns_r$locyear
AEP_year = 2019
gevns_RL5pc_2019 = raster_qevd(x = gev_r,p = 1-.05, evd_mod_str = "fgev") + nsloctab$year_cs[nsloctab$year==AEP_year]*gevns_r$locyear

plotit = c(gevns_RL5pc_1981,gevns_RL5pc_2019)
names(plotit) = c("gevns_RL5pc_1981","gevns_RL5pc_2019")
at=c(-Inf,seq(27.5,50,2.5),Inf)
cp = colorRampPalette(c("light yellow","yellow","red","brown4"))
#pdf("../plots/Empirical_20_year_ARI_daily_max_temp.pdf")
sp::spplot(plotit,at=at,col.regions =c("#FFFFFF",cp(9),"#000000"),main=paste0("GEV_ns 5% AEP for the years 1981 and 2019"))

## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
at = c(-Inf,seq(-4,4,.5),Inf)
sp::spplot(gevns_RL5pc_2019-gevns_RL5pc_1981,at=at,col.regions =c("#FFFFFF",rev(rwb(17)),"#000000"),main=paste0("Non-stationary change in 5% AEP GEV extremes 1981 to 2019"))

## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
poi = vect(rbind(c(142.91964371685395,-35.873803703002594)) ,"points")
#extract a time series for the point of interest
crs(poi) = crs(gevns_r)
crs(r) = crs(poi)
valz =  as.numeric(extract(r,poi,ID = FALSE))
poi_ts = data.frame(date = time(r),annMax = valz)

gumbel_v = extract(gumbel_r,poi,xy=FALSE,ID=FALSE)
gumbel_v = as.data.frame(t(sapply(gumbel_v, function(x) x[[1]])))

gev_v = extract(gev_r,poi,xy=FALSE,ID=FALSE)
gev_v = as.data.frame(t(sapply(gev_v, function(x) x[[1]])))

gumbelx_v = extract(gumbelx_r,poi,xy=FALSE,ID=FALSE)
gumbelx_v = as.data.frame(t(sapply(gumbelx_v, function(x) x[[1]])))

gumbelns_v = extract(gumbelns_r,poi,xy=FALSE,ID=FALSE)
gumbelns_v = as.data.frame(t(sapply(gumbelns_v, function(x) x[[1]])))

gevns_v = extract(gevns_r,poi,xy=FALSE,ID=FALSE)
gevns_v = as.data.frame(t(sapply(gevns_v, function(x) x[[1]])))

gevns_63AEP = df_qevd(x = gev_v,p = 1-0.63,evd_mod_str = "fgev") + nsloctab$year_cs*gevns_v$locyear
gevns_05AEP = df_qevd(x = gev_v,p = 1-0.05,evd_mod_str = "fgev") + nsloctab$year_cs*gevns_v$locyear
gevns_01AEP = df_qevd(x = gev_v,p = 1-0.01,evd_mod_str = "fgev") + nsloctab$year_cs*gevns_v$locyear

YLIM = c(min(poi_ts$annMax),max(gevns_01AEP))
plot(r_years,valz,type = "l",main = "Non stationary EVD",ylim = YLIM);grid()

poi_ts$annMaxns = poi_ts$annMax-nsloctab$year_cs*gevns_v$locyear
lines(r_years,poi_ts$annMaxns,col="grey")

lines(r_years,gevns_63AEP,col = 4)
lines(r_years,gevns_05AEP,col = 2)
lines(r_years,gevns_01AEP,col = 6)

legend("bottomright",c("1% AEP","5% AEP","63% AEP","AnnMax","AnnMax_stat"),col = c(6,2,4,1,"grey"),lty=1,ncol = 2)




## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
ndat = length(poi_ts$annMax)

empi_AEP = -1/log((1:ndat)/(ndat + 1))
AEP = 1-exp(-1/empi_AEP)

gumbel_RL = qevd_vector(x = gumbel_v,p = 1-AEP,evd_mod_str = "fgumbel")
gev_RL = qevd_vector(x = gev_v,p = 1-AEP,evd_mod_str = "fgev")
gumbelns_RL = qevd_vector(x = gumbelns_v,p = 1-AEP,evd_mod_str = "fgumbel")
gevns_RL = qevd_vector(x = gevns_v,p = 1-AEP,evd_mod_str = "fgev")
gumbelx_RL = qevd_vector(x = gumbelx_v,p = 1-AEP,evd_mod_str = "fgumbelx",interval = c(0,20))

plot_empirical(x=poi_ts$annMax,xns = poi_ts$annMaxns)
lines(empi_AEP,gumbel_RL)
lines(empi_AEP,gev_RL,col=2)
#lines(empi_AEP,gumbelx_RL,col="green")
lines(empi_AEP,gumbelns_RL,col="grey")
lines(empi_AEP,gevns_RL,col="yellow")

#legend('topleft',legend = paste0(c("fgumbelx","fgev","fgumbel","fgumbelns"),rep(" AIC = ",4),round(as.numeric(c(rev(gumbelx)[1],rev(gev)[1],rev(gumbel)[1],rev(gumbelns)[1])),1)),col = c(4,2,1,"grey"),lty = 1)
 legend('bottomright',legend = paste0(c("fgev","fgumbel","fgumbelns","fgevns"),rep(" AIC = ",4),round(as.numeric(c(rev(gev_v)[1],rev(gumbel_v)[1],rev(gumbelns_v)[1],rev(gevns_v)[1])),1)),col = c(2,1,"grey","yellow"),lty = 1)


## -----------------------------------------------------------------------------
# = terra::rast(system.file("extdata/50km_AnnMax_agcd_v1_precip_calib_r005_daily_1980-2019.nc",package = "loopevd"))




#' Fit EVD Parameters to a List of Annual Maxima
#'
#' @description
#' Fits an extreme value distribution to each element of a list of annual maxima series, optionally using non-stationary covariates, and returns a table of fitted parameters.
#'
#' @param lst A \code{list} of \code{data.frame}s, each as returned by \code{annual_max()}, containing at least:
#'   \itemize{
#'     \item \code{annMax}: annual maximum values
#'     \item \code{annMean}: annual mean values
#'     \item \code{datestr}: timestamp strings for the maxima
#'     \item \code{date}: POSIX timestamps for the maxima
#'     \item \code{pc_complete}: completeness fraction per year
#'     \item \code{zero}: the zero level
#'     \item \code{nsloc}: optional covariate matrix used in non-stationary models.
#'   }
#'   If non-stationary fitting is required, each element may also include an \code{nsloc} matrix of covariates.
#' @param evd_mod_str A \code{character} string specifying which fitting function from \pkg{evd} to use: \code{"fgumbel"}, \code{"fgumbelx"} or \code{"fgev"}.
#' @param nsloc Optional \code{matrix} of covariates for non-stationary location modelling. Must have the same number of rows as years retained after filtering.
#' @param outfile Optional \code{character} path to a NetCDF file in which to write the results (not currently implemented).
#' @param pc_complete Numeric scalar (0-1). Minimum completeness fraction for a year to be included. Defaults to \code{0.8}.
#' @param minyear Numeric. Minimum calendar year to include. Defaults to \code{1800}.
#' @param maxyear Numeric. Maximum calendar year to include. Defaults to \code{2100}.
#' @param mslAtt \code{character}. Name of the attribute to be removed from annMax in each data.frame (e.g.\ \code{"annMean"} or \code{"zero"}). Defaults to \code{"annMean"}.
#'
#' @return A \code{data.frame} with one row per list element, containing the parameters returned by \code{evd_params()} for each annual-max series.
#'
#' @examples
#' dates = seq.Date(as.Date("1990-01-01"),as.Date("2019-12-31"), "day")
#' lst = lapply(1:10,function(x) loopevd::annual_max(data.frame(date = dates,
#'                   sea_level = stats::rnorm(length(dates),mean=x/10,sd = x),
#'                   zero = rep(0,length(dates)))))
#' loopevd::list_fevd(lst,"fgumbel",pc_complete=0)
#' @export
list_fevd = function(lst,evd_mod_str,nsloc=NULL,outfile=NULL,pc_complete=0.8,minyear=1800,maxyear=2100,mslAtt = "annMean"){
  if(!is.null(lst[[1]]$nsloc[1])) {
    nsloc = "not empty"
  }

  if( is.null(nsloc[1])) lst_evd_params = lapply(X = lst,FUN = function(x) {s = x$pc_complete > pc_complete & as.numeric(substr(x$datestr,1,4)) <= maxyear & as.numeric(substr(x$datestr,1,4)) >= minyear;
  evd_params((x$annMax-x[[mslAtt]])[s], evd_mod_str = evd_mod_str,nsloc = NULL)})

  if(!is.null(nsloc[1])) lst_evd_params = lapply(X = lst,FUN = function(x) {s = x$pc_complete > pc_complete & as.numeric(substr(x$datestr,1,4)) <= maxyear & as.numeric(substr(x$datestr,1,4)) >= minyear & !is.na(x$nsloc[,1]);
  evd_params((x$annMax-x[[mslAtt]])[s], evd_mod_str = evd_mod_str,nsloc = x$nsloc[s,, drop = FALSE])})
  lst_evd_params
  #tab_evd_params = sapply(lst_evd_params)
  #lngs = lengths(lst_evd_params)
  #lst_evd_params[[which(lngs == 1)]] = rep(NA,max(lngs))
  out = t(sapply(lst_evd_params,function(x) x))
  as.data.frame(out)
}


#' Centre and Scale Numeric Data
#'
#' \code{centredAndScaled} centres (subtracts the mean) and scales (divides by
#' the standard deviation) each column of a numeric vector or matrix.
#'
#' @param nsloc data.frame.  If \code{NULL} or of length zero,
#'   the function returns \code{nsloc} unchanged.
#'
#' @return A numeric vector or matrix of the same dimensions as the input,
#'   with each column centred to mean zero and scaled to unit variance.  If
#'   \code{nsloc} is \code{NULL}, returns \code{NULL}.
#'
#' @details
#' If \code{nsloc} has only one column, the function computes the mean and
#' standard deviation of the entire vector.  If \code{nsloc} has multiple
#' columns, each column is centred and scaled independently.
#'
#' @examples
#' # Centre and scale a simple vector
#' centredAndScaled(data.frame(1:10))
#'
#' # Centre and scale each column of a matrix
#' mat <- as.data.frame(matrix(stats::rnorm(30), nrow = 10, ncol = 3))
#' centredAndScaled(mat)
#'
#' @export
centredAndScaled = function(nsloc = NULL){
  if(!is.null(nsloc[1])){
    #centre and scale the data
    if(dim(nsloc)[2] == 1){
      cnt = mean(nsloc[,])
      scl = sd(nsloc[,])
      nsloc = (nsloc-cnt)/scl

    }
    if(dim(nsloc)[2] > 1){
      cnt = apply(nsloc,2,mean)
      scl = apply(nsloc,2,sd)
      nsloc = (nsloc - t(array(cnt,dim = rev(dim(nsloc)))))/t(array(scl,dim = rev(dim(nsloc))))
    }
  }
  return(nsloc)
}

#' Turn Raster of Annual Maximums into Extreme Value Distributions parameters for Netcdf Output
#'
#' @param r SpatRaster
#' @param evd_mod_str either a string "fgumbel", "fgev" or "fgumbelx" from the extreme value distribution (evd) in the evd package
#' @param nsloc A data frame with the same number of rows as the length of x, for linear modelling of the location parameter. The data frame is treated as a covariate matrix (excluding the intercept). A numeric vector can be given as an alternative to a single column data frame.
#' @param outfile filename to write to netcdf
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used. You can also supply a cluster object. Ignored for functions that are implemented by terra in C++ (see under fun)
#' @param ntries integer number of attempts at fitting fgumbelx
#' @param silent logical: should the report of error messages be suppressed?
#' @param seed set the seed for fitting.
#'
#' @return the parameters of the evd in a SpatRasterDataset
#' @export
#' @seealso [evd::fgev()], [evd::fgumbelx()]
#'
#' @examples
#' require(terra)
#' r = rast(system.file("extdata/50km_AnnMax_agcd_v1_tmax_mean_r005_daily_1980-2019.nc"
#' ,package = "loopevd"))
#' r2 = aggregate(r,4) #lower the resolution for a fast example
#' gumbel_r = raster_fevd(r2,"fgumbel",seed = 1)
#' plot(gumbel_r$loc,main = "location")
raster_fevd = function(r,evd_mod_str,nsloc=NULL,outfile=NULL,cores = 1,ntries=1,silent = FALSE,seed=NULL){

  if(!is.null(seed)){
    set.seed(seed)
  }

  uvdata <- evd::rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)

  if( is.null(nsloc)) empty_evd_params = evd_params(uvdata,evd_mod_str)*NA
  if(!is.null(nsloc)) {
    nns = dim(nsloc)[2]
    uvnsloc <- as.data.frame(array(stats::rnorm(1000),dim = c(100,nns)))
    names(uvnsloc) <- names(nsloc)
    empty_evd_params = evd_params(x = uvdata,evd_mod_str = evd_mod_str,nsloc = uvnsloc,returncs=FALSE,silent = FALSE)*NA
  }

  cl = cores
  if(cores > 1){  #need to export function variables to all cpus
    cl = parallel::makeCluster(cores)
    parallel::clusterExport(cl = cl,varlist = c("evd_mod_str","nsloc"),envir = environment())
  }
  if( is.null(nsloc)) r_evd_params = terra::app(x=r,fun = evd_params, evd_mod_str = evd_mod_str,nsloc = NULL,empty_evd_params=empty_evd_params,cores=cl,ntries=ntries,returncs=FALSE,silent = FALSE)
  if(!is.null(nsloc)) r_evd_params = terra::app(x=r,fun = evd_params, evd_mod_str = evd_mod_str,nsloc = nsloc,empty_evd_params=empty_evd_params,cores=cl,ntries=ntries,returncs=FALSE,silent = FALSE)
  if(cores > 1) parallel::stopCluster(cl)

  terra::crs(r_evd_params) = "+init=epsg:4326 +proj=longlat"
  #need to convert to a spatial dataset to add netcdf variable params
  r_evd_params_sds =  terra::sds(lapply(r_evd_params,function(x) x))
  nc_names = names(r_evd_params)
  terra::varnames(r_evd_params_sds) = nc_names
  #terra::longnames(r_evd_params_sds) = nc_longnames
  #terra::units(r_evd_params_sds) = nc_units
  if(is.null(outfile)) return(r_evd_params_sds)
  if(!is.null(outfile)) {
    terra::writeCDF(r_evd_params_sds,filename = outfile, overwrite=TRUE, zname="Parameter")
  }

}


#' Function to return the EVD (Gumbel or GEV) parameters as a vector.
#'
#' @param x  numeric vector of data to be fitted.
#' @param evd_mod_str either a string "fgumbel", "fgumbelx" or "fgev" from the extreme value distribution (evd) in the evd package
#' @param nsloc A data frame with the same number of rows as the length of x, for linear modelling of the location parameter. The data frame is treated as a covariate matrix (excluding the intercept). A numeric vector can be given as an alternative to a single column data frame.
#' @param empty_evd_params A preallocated vector or array used to store the return value when fitting fails
#' @param ntries number of tries to fit the evd
#' @param silent logical: should the report of error messages be suppressed?
#' @param returncs logical: should the centered and scaled values be returned
#' @importFrom stats sd AIC loess pnorm rnorm uniroot
#' @return a vector of estimate, var.cov, AIC, centered and scaled values
#' @export
#'
#' @examples
#' # Ten records of 20 random data generated from the fgumbel EVD
#' am = lapply(1:10, function(x) evd::rgumbel(20))
#' tab = as.data.frame(t(sapply(am,function(x) evd_params(x,"fgumbel"))))
evd_params = function(x,evd_mod_str,nsloc = NULL,empty_evd_params,ntries = 3,silent = FALSE,returncs=TRUE){
  #https://scientistcafe.com/ids/centering-and-scaling.html
  x = as.numeric(x)
  if(evd_mod_str == "fgumbel") {
    evd_mod = evd::fgumbel}
  if(evd_mod_str == "fgev") {
    evd_mod = evd::fgev}
  if (evd_mod_str == "fgumbelx") {
    evd_mod = evd::fgumbelx}

  m <- NULL
  m1 <- NULL
  if(!is.na(x[1])){
    if(!is.null(nsloc[1])){
      #centre and scale the data
      if(dim(nsloc)[2] == 1){
        cnt = mean(nsloc[,])
        scl = sd(nsloc[,])
        nsloc = (nsloc-cnt)/scl

      }
      if(dim(nsloc)[2] > 1){
        cnt = apply(nsloc,2,mean)
        scl = apply(nsloc,2,sd)
        nsloc = (nsloc - t(array(cnt,dim = rev(dim(nsloc)))))/t(array(scl,dim = rev(dim(nsloc))))
      }
    }
    if(evd_mod_str != "fgumbelx"){
      try(m1 <- evd::fgumbel(x),silent = silent)
      if(!is.null(m1[1])){
        if(!is.null(nsloc[1])) try(m <- evd_mod(x = x,nsloc = nsloc),silent = silent)
        if(is.null(nsloc[1])) {

            if(is.null(m[1])) try(m <- evd_mod(x = x, method = "BFGS"),silent = silent)
            if(is.null(m[1])) try(m <- evd_mod(x = x, method = "Nelder-Mead"),silent = silent)
            if(is.null(m[1])) try(m <- evd_mod(x = x, method = "CG"),silent = silent)
          for(ri in 1:ntries){
            if(is.null(m[1])) try(m <- evd_mod(x = x, method = "SANN"),silent = silent)
          }
        }
        #if(!is.null(m[1])) out = c(m$estimate,m$var.cov,m$deviance/2)
        if(!is.null(m[1])) {
          out = c(m$estimate,m$var.cov,AIC(m))
          names(out) = c(names(m$estimate),paste0("cov_",1:length(m$var.cov)),"AIC")
          if(!is.null(nsloc[1])){
            if(returncs){
              out = c(m$estimate,m$var.cov,AIC(m),cnt,scl)
              names(out) = c(names(m$estimate),paste0("cov_",1:length(m$var.cov)),"AIC",
                             paste0("Centred_",names(nsloc)),paste0("Scaled_",names(nsloc)))
            }
            if(!returncs){
              out = c(m$estimate,m$var.cov,AIC(m))
              names(out) = c(names(m$estimate),paste0("cov_",1:length(m$var.cov)),"AIC")
            }
          }
        }
      }
      if( is.null(m[1])) out = empty_evd_params
    }

    if(evd_mod_str == "fgumbelx"){
      try(m1 <- evd::fgumbel(x),silent = silent)
      if(!is.null(m1[1])){
        #try(m <- fgumbelx(x,loc1 = m1$estimate[1],scale1 = m1$estimate[2]))
        #try(m <- fgumbelx(x))
        #if(!is.null(m[1])) out = c(m$estimate,m$var.cov,m$deviance/2)

        if(is.null(m[1])) try(m <- evd::fgumbelx(x, method = "BFGS",warn.inf=!silent),silent = silent)
        if(is.null(m[1])) try(m <- evd::fgumbelx(x, method = "Nelder-Mead",warn.inf=!silent),silent = silent)
        if(is.null(m[1])) try(m <- evd::fgumbelx(x, method = "CG",warn.inf=!silent),silent = silent)
        for(ri in 1:ntries){
          if(is.null(m[1])) try(m <- evd::fgumbelx(x, method = "SANN",warn.inf=!silent),silent = silent)
        }
        if(is.null(m[1])) message("fgumbelx fit not found")
        if(!is.null(m[1])){
          out = c(m$estimate, m$var.cov, AIC(m))
          names(out) = c(names(m$estimate),paste0("cov_",1:length(m$var.cov)),"AIC")
        }
      }
      if( is.null(m[1])) out = empty_evd_params
    }
  }
  if(is.na(x[1])) out = empty_evd_params
  return(out)
}




#' Calculate One-Sided Confidence Level (%)
#'
#' @description
#' Computes the one-sided confidence level, defined as (1 - p-value) x 100, for testing whether each mean (`mu`) differs from zero under a normal approximation.
#'
#' @param muvari   Numeric array, of mean (location) values, variances corresponding to each `mu` to test against zero.
#'
#' @return A numeric vector of confidence levels (0-100%), each rounded to one decimal place.
#'
#' @details
#' For each element:
#' 1. Calculate the standard error:
#'    se = sqrt(vari).
#' 2. Compute the absolute z-score:
#'    z = abs(mu / se).
#' 3. The one-sided p-value is 1 - phi(z), where phi is the CDF of the standard normal.
#' 4. The confidence level is (1 - p-value) x 100 = phi(z) x 100.
#'
#' @examples
#' # Single value
#' se_sig(muvari = cbind(2,1))
#'
#' # Vector of values
#' se_sig(muvari = cbind(c(-1, 0, 1),c(1, 2, 3)))
#'
#' @importFrom stats pnorm
#' @export
se_sig = function(muvari) {
  mu = muvari[1]
  vari = muvari[2]
  se = sqrt(vari)                # Standard error
  z = abs(mu / se)               # Absolute z-score
  p_value = 1 - pnorm(z)         # One-tailed p-value
  confidence_percent = (1 - p_value) * 100  # Confidence level as percentage
  return(round(confidence_percent, 2))      # Rounded for clarity
}


#' Calculate One-Sided Confidence Level (%)
#'
#' @description
#' Computes the one-sided confidence level, defined as (1 - p-value) x 100, for testing whether each mean (`mu`) differs from zero under a normal approximation.
#'
#' @param muvari SpatRaster, of mean (location) values, variances corresponding to each `mu` to test against zero.
#'
#' @return A SpatRaster of confidence levels (0-100%), each rounded to one decimal place.
#'
#' @details
#' For each element:
#' 1. Calculate the standard error:
#'    se = sqrt(vari).
#' 2. Compute the absolute z-score:
#'    z = abs(mu / se).
#' 3. The one-sided p-value is 1 - phi(z), where phi is the CDF of the standard normal.
#' 4. The confidence level is (1 - p-value) x 100 = phi(z) x 100.
#'
#' @examples
#' require(terra)
#' r = rast(system.file("extdata/50km_AnnMax_agcd_v1_tmax_mean_r005_daily_1980-2019.nc",
#'                      package = "loopevd"))
#' r2 = aggregate(r,4) #lower the resolution for a fast example
#' gev_r = raster_fevd(r2,"fgev")
#' raster_se_sig(c(gev_r$shape,gev_r$cov_9))
#' @importFrom stats pnorm
#' @export
raster_se_sig = function(muvari) terra::app(muvari,se_sig)


#' Plot the Empirical Return Level Data
#'
#' @param x A numeric vector, which may contain missing values.
#' @param xns A numeric vector, corrected for the non-stationary change in location, which may contain missing values.
#' @param unitz y-label
#' @param ... parameters sent to base::plot
#'
#' @return r
#' @export
#'
#' @examples
#' ns = seq(-1,1,,50)
#' x = evd::rgev(50,loc=3)+ns
#' xns = x-ns
#' plot_empirical(x,xns)
plot_empirical = function(x,xns=NULL,unitz = "-",...){
  ndat = length(x)
  empi_AEP = -1/log((1:ndat)/(ndat + 1))
  #rng = range(x,xns,na.rm=TRUE)
  plot(empi_AEP,sort(x),col = 1,pch = 1,log='x',xlab = "Average Recurrence Interval [years]",ylab = unitz,...)
  if(!is.null(xns)) graphics::points(empi_AEP,sort(xns),col="grey")
  graphics::grid()
}

#' Compute Annual Maximum and Mean of On-the-Hour Records
#'
#' @description
#' annual_max() takes a data frame of daily (or sub-daily) observations and returns a summary of the annual maximum and mean values, the date/time of each annual maximum, and the fraction of "on-the-hour" samples (data completeness) for each calendar year.
#'
#' @param DF A data.frame containing at least:
#'   * date: a Date or POSIXlt column of observation timestamps
#'   * a numeric column of values
#' @param record_attribute A character string giving the name of the column in DF containing the values.  Defaults to "sea_level".
#'
#' @return a data.frame containing a date column and a numeric column (specified in record_attribute) for years where at least one nonNA value is present, containing:
#'   * annMax - the annual maximum
#'   * annMean - the annual mean (calendar year)
#'   * datestr - the date/time of the annual maximum, formatted "YYYYmmddHH"
#'   * date - the POSIXlt timestamp of the annual maximum
#'   * pc_complete - the fraction (0 to 1) of hourly-timestamped samples available in that year
#'
#' @details
#' For each year, only observations exactly on the hour (minute == 0 & second == 0) are counted toward completeness.  If no valid data exist for a year, that year is dropped from the output.
#'
#' @examples
#' # generate example daily data
#' dates <- seq.Date(as.Date("1990-01-01"), as.Date("1995-12-31"), by = "day")
#' DF <- data.frame(
#'   date      = dates,
#'   sea_level = rnorm(length(dates), mean = 0, sd = 1)
#' )
#' # compute annual summary
#' am <- annual_max(DF)
#' head(am)
#'
#' @export
annual_max = function(DF,record_attribute = "sea_level"){
  tm = DF$date
  dat = DF[[record_attribute]]
  YR = as.numeric(format(tm,"%Y"))
  YRS = unique(YR)
  MINS = as.numeric(format(tm,"%M"))
  ny = length(YRS)
  nt = length(YR)

  annMax = array(-Inf,dim = ny)	#annual maximum
  annMean = annMax				#annual mean (calander year)
  annMaxTm = tm[1:ny]			#date of the annual max (first)
  pc_complete = annMax
  for(i in 1:ny) {
    ixy = YR==YRS[i]
    annMax[i] = max(dat[ixy],na.rm=T)
    annMaxTm[i] = tm[which(dat == annMax[i] & YR==YRS[i])[1]]

    mint = format(tm[ixy],"%M%S")
    l1 = length(dat[ixy][mint == "0000" & !is.na(dat[ixy])]) #only use on the hour
    t1 = strptime(paste0(YRS[i],"0101"),"%Y%m%d",tz = "UTC");t2 = strptime(paste0(YRS[i]+1,"0101"),"%Y%m%d",tz = "UTC")
    l2 = (as.numeric(t2)-as.numeric(t1))/3600 #number of continuos samples

    pc_complete[i] = round(l1/l2,2)

    annMean[i] = round(mean(dat[YR==YRS[i]],na.rm=T),2)
  }
  annMax[as.numeric(annMax) == -Inf] = NaN
  notNAN = !is.nan(annMax)
  zero = annMax*0
  out = data.frame(annMax = annMax[notNAN],annMean = annMean[notNAN],datestr = format(annMaxTm[notNAN],"%Y%m%d%H"),
                   date = annMaxTm[notNAN],pc_complete=pc_complete[notNAN],zero = zero[notNAN])#,runMaxAtTm = runMaxAtTm[notNAN],runAnnMaxAtTm = runAnnMaxAtTm[notNAN])
  return(out)
}

#' Return a Vector of EVD Quantiles
#'
#' @param x vector of EVD parameters
#' @param p vector of probabilities.
#' @param evd_mod_str either a string "fgumbel", "fgev" or "fgumbelx" from the extreme value distribution (evd) in the evd package
#' @param interval A length two vector containing the end-points of the interval to be searched for the quantiles, passed to the uniroot function.
#' @param lower.tail Logical; if TRUE (default), P (\eqn{x \le y}), otherwise P (X > x).
#' @param nams names of the values of x (optional)
#' @return gives the quantile function corresponding to p
#' @export
#' @seealso [evd::qgev()], [evd::qgumbelx()]
#'
#' @examples
#' qevd_vector(c(1,0.5),1-0.05,"fgumbel",nams = c("loc","scale"))
#' df = data.frame(loc = 1,scale = 0.5)
#' qevd_vector(df,1-0.05,"fgumbel")
#'
qevd_vector = function(x,p,evd_mod_str,interval = NULL,lower.tail = TRUE,nams = NULL){
  if(is.null(nams[1])) nams = names(x)
  loci = which(nams == "loc")
  scalei = which(nams == "scale")
  shapei = which(nams == "shape")
  loc1i = which(nams == "loc1")
  scale1i = which(nams == "scale1")
  loc2i = which(nams == "loc2")
  scale2i = which(nams == "scale2")

  x = as.numeric(x)
  #if(is.null(dim(x)[1])) stop("x must have evd parameters for one site only")
  out = NA
  if(evd_mod_str == "fgumbel"  ) try(out <- evd::qgev(p, loc = x[loci], scale = x[scalei], shape = 0, lower.tail = lower.tail),silent = TRUE)
  if(evd_mod_str == "fgev"     ) try(out <- evd::qgev(p, loc = x[loci], scale = x[scalei], shape = x[shapei], lower.tail = lower.tail),silent = TRUE)
  if(evd_mod_str == "fgumbelx" ) {
    if(is.null(interval[1])) {
      try(gumq <- evd::qgev(p, loc = x[loci], scale = x[scalei], shape = 0, lower.tail = TRUE),silent = lower.tail)
      try(interval <- gumq+(evd::qgev(p, loc = 0, scale = x[scalei], shape = 0, lower.tail = TRUE))*c(-2,2),silent = TRUE)
    }
    if(!is.null(interval[1])) try(out <- evd::qgumbelx(p = p,interval = interval, loc1 = x[loc1i], scale1 = x[scale1i], loc2 = x[loc2i],scale2=x[scale2i], lower.tail = lower.tail),silent = TRUE)
  }
  return(out)
}


#' Return a raster of EVD Quantiles
#'
#' @param x SpatRasterDataset of EVD parameters, e.g. loc, scale, shape
#' @param p probability value.
#' @param evd_mod_str either a string "fgumbel", "fgev" or "fgumbelx" from the extreme value distribution (evd) in the evd package
#' @param interval A length two vector containing the end-points of the interval to be searched for the quantiles, passed to the uniroot function.
#' @param lower.tail Logical; if TRUE (default), probabilities are P \eqn{x \le y}), otherwise P (X > x).
#' @return gives the quantile function corresponding to p
#' @export
#'
#' @seealso [evd::qgev()], [evd::qgumbelx()]
#'
#' @examples
#' require(terra)
#' r = rast(system.file("extdata/50km_AnnMax_agcd_v1_tmax_mean_r005_daily_1980-2019.nc",
#'                      package = "loopevd"))
#' r2 = aggregate(r,4) #lower the resolution for a fast example
#' gumbel_r = raster_fevd(r2,"fgumbel")
#' AEP_10pc = raster_qevd(gumbel_r,1-0.1,"fgumbel") # 10% Annual Exceedance Probability.
raster_qevd = function(x,p,evd_mod_str,interval = NULL,lower.tail = TRUE){
  if(class(x) != class(sds(rast()))) stop("x must be of the class SpatRasterDataset")
  if(length(p) > 1) stop("provide one value p")
  terra::app(x = x,fun = loopevd::qevd_vector, p = p, evd_mod_str = evd_mod_str,interval=interval,lower.tail=lower.tail,nams = names(x))
}


#' Return a vector of EVD Quantiles
#'
#' @param x A data.frame of EVD parameters
#' @param p Probability value (e.g. 0.99)
#' @param evd_mod_str One of "fgumbel", "fgev", or "fgumbelx"
#' @param interval Optional; interval passed to uniroot when needed
#' @param lower.tail Logical; if TRUE (default), computes P(\eqn{x \le y})
#' @param cores Number of cores to use for parallel processing (default is 1)
#'
#' @return A numeric vector of quantiles
#' @importFrom parallel makeCluster clusterExport parLapply clusterEvalQ stopCluster
#' @export
#' @example
#' df <- data.frame(loc = 1,scale = 0.5)
#' qevd_vector(df,1-0.05,"fgumbel")
df_qevd <- function(x, p, evd_mod_str, interval = NULL, lower.tail = TRUE, cores = 1) {
  if (!is.data.frame(x)) stop("x must be a data.frame")
  if (length(p) != 1) stop("provide one value for p")

  if (cores > 1) {
    cl <- makeCluster(cores)
    clusterExport(cl, varlist = c("x", "p", "evd_mod_str", "interval", "lower.tail", "qevd_vector"), envir = environment())
    clusterEvalQ(cl, library(loopevd))

    out <- parLapply(cl, 1:nrow(x), function(i) {
      tryCatch(
        qevd_vector(x[i, ], p = p, evd_mod_str = evd_mod_str, interval = interval, lower.tail = lower.tail),
        error = function(e) {
          message(sprintf("Error in row %d: %s", i, e$message))
          return(NA)
        }
      )
    })

    stopCluster(cl)
    return(as.numeric(out))
  } else {
    out <- apply(X = x, MARGIN = 1, FUN = loopevd::qevd_vector,
                 p = p, evd_mod_str = evd_mod_str, interval = interval, lower.tail = lower.tail)
    return(out)
  }
}

#' Internal function to compute EVD quantile confidence intervals
#'
#' @param xrow A single-row data.frame of EVD parameters and covariance matrix
#' @param p Probability value (e.g. 0.99)
#' @param ci Confidence level (default 0.95)
#'
#' @return Numeric vector of confidence interval width(s)
#' @keywords internal
#'
#' @example
#' xrow <- data.frame(loc = 1, scale = 0.5,
#'       cov_1 = 0.01, cov_2 = 0.001, cov_3 = 0.002,
#'       cov_4 = 0.001)
#' cievd_vector(xrow,1-0.1,0.95)
#' xrow <- data.frame(loc = 1,cale = 0.5,shape = .01,
#'     cov_1 = 0.01, cov_2 = 0.001, cov_3 = 0.002,
#'     cov_4 = 0.001, cov_5 = 0.01,  cov_6 = 0.001,
#'     cov_7 = 0.002, cov_8 = 0.001, cov_9 = 0.02
#'     )
#' cievd_vector(xrow,1-0.1,0.95)
cievd_vector <- function(xrow, p, ci = 0.95) {
  if (!is.data.frame(xrow)) xrow <- as.data.frame(xrow)
  if (is.na(xrow[[1]])) return(NA)

  if(!any(names(xrow) == "shape")) a <- matrix(as.numeric(xrow[c("loc", "scale")]))
  if( any(names(xrow) == "shape")) a <- as.numeric(xrow[, c("loc", "scale", "shape")])
  cov_idx <- which(substr(names(xrow), 1, 3) == "cov")

  if (length(cov_idx) == 0 || (sqrt(length(cov_idx)) %% 1 != 0)) {
    stop("Covariance matrix elements not found or malformed")
  }

  sqncov <- sqrt(length(cov_idx))
  mat <- matrix(as.numeric(xrow[, cov_idx]), nrow = sqncov, ncol = sqncov)

  if(!any(names(xrow) == "shape")){
    a[3] = 0
    eps = 1e-06
    a1 <- a
    a2 <- a
    a1[1] <- a[1] + eps
    a2[2] <- a[2] + eps
    q <- ismev::gevq(a, 1 - p)
    d1 <- (ismev::gevq(a1, 1 - p) - q)/eps
    d2 <- (ismev::gevq(a2, 1 - p) - q)/eps
    d <- cbind(d1, d2)
  }
  if( any(names(xrow) == "shape")) d <- t(ismev::gev.rl.gradient(a = matrix(a), p = 1 - p))
  v <- apply(d, 1, ismev::q.form, m = mat)
  z <- stats::qnorm(1 - (1 - ci) / 2)
  cv <- z * sqrt(v)

  return(cv)
}

#' Return a numeric vector of confidence intervals for EVD quantiles
#'
#' @param x A data.frame of EVD parameters and associated covariance matrix. Must include 'loc', 'scale', 'shape', and 'cov*' column names.
#' @param p A single probability value for the quantile (e.g. 0.99).
#' @param ci Confidence level for the interval (default is 0.95).
#' @param cores Number of parallel cores to use. If >1, parallel processing is used via \code{parallel::parLapply()}.
#'
#' @return A numeric vector giving the confidence interval widths for each row in \code{x}.
#'
#' @details This function calculates confidence intervals for extreme value quantiles using the delta method. The required input is a row-wise data frame with EVD parameters and their variance-covariance elements.
#' Internally uses \code{ismev::gev.rl.gradient()} and \code{ismev::q.form()}.
#'
#' @seealso [ismev::gev.rl.gradient()], [ismev::q.form()]
#' @export
#' @importFrom parallel makeCluster clusterExport parLapply clusterEvalQ stopCluster
#' @examples
#' \dontrun{
#'   df <- data.frame(loc = 1, scale = 0.5, shape = 0.1,
#'                    cov_1 = 0.01, cov_2 = 0.001, cov_3 = 0.002,
#'                    cov_4 = 0.001, cov_5 = 0.01, cov_6 = 0.001,
#'                    cov_7 = 0.002, cov_8 = 0.001, cov_9 = 0.02)
#'   df_cievd(df, p = 0.99, ci = 0.95, cores = 2)
#' }
#'
df_cievd <- function(x, p, ci = 0.95, cores = 8) {
  if (!is.data.frame(x)) stop("x must be a data.frame")
  if (length(p) != 1) stop("p must be a single value")

  if (cores > 1) {
    cl <- makeCluster(cores)
    clusterExport(cl, varlist = c("x", "p", "ci", "cievd_vector"), envir = environment())
    clusterEvalQ(cl, library(ismev))

    result <- parLapply(cl, 1:nrow(x), function(i) {
      tryCatch(
        cievd_vector(x[i, ], p, ci),
        error = function(e) {
          message(sprintf("Error in row %d: %s", i, e$message))
          return(NA)
        }
      )
    })

    stopCluster(cl)
  } else {
    result <- lapply(1:nrow(x), function(i) {
      tryCatch(
        cievd_vector(x[i, ], p, ci),
        error = function(e) {
          message(sprintf("Error in row %d: %s", i, e$message))
          return(NA)
        }
      )
    })
  }

  return(as.numeric(result))
}

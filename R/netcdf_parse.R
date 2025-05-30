#' Add Global Attribute Metadata to Netcdf File
#'
#' Improve FAIR metadata record with consideration of CF conventions https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html
#'
#' @param outfile character the file to be edited
#' @param r SpatRaster the dataset in raster format e.g. r = terra::rast(outfile)
#' @param creator_name character, optional. The name of the person (or other creator type specified by the creator_type attribute) principally responsible for creating this data.
#' @param creator_email character, optional. The email address of the person (or other creator type specified by the creator_type attribute) principally responsible for creating this data.
#' @param references character, optional. Published or web-based references that describe the data or methods used to produce it. Recommend URIs (such as a URL or DOI) for papers or other references. This attribute is defined in the CF conventions.
#' @param title character, optional. A short phrase or sentence describing the dataset. In many discovery systems, the title will be displayed in the results list from a search, and therefore should be human readable and reasonable to display in a list of such names. This attribute is also recommended by the NetCDF Users Guide and the CF conventions.
#' @param summary character, optional. A paragraph describing the dataset, analogous to an abstract for a paper.
#' @param keywords  character, optional. A comma-separated list of key words and/or phrases. Keywords may be common words or phrases, terms from a controlled vocabulary (GCMD is often used), or URIs for terms from a controlled vocabulary (see also "keywords_vocabulary" attribute).
#' @param history  character, optional. Provides an audit trail for modifications to the original data. This attribute is also in the NetCDF Users Guide: 'This is a character array with a line for each invocation of a program that has modified the dataset. Well-behaved generic netCDF applications should append a line containing: date, time of day, user name, program name and command arguments.' To include a more complete description you can append a reference to an ISO Lineage entity; see NOAA EDM ISO Lineage guidance.
#' @param licence  character, optional.
#' @param Disclaimer character, optional.
#' @return see https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3, http://cfconventions.org/cf-conventions/cf-conventions.html
#' @export
#' @examples
#' ## Not run:
#' tf = tempfile("test.nc")
#' tf
#' file.copy(system.file("extdata/50km_AnnMax_agcd_v1_tmax_mean_r005_daily_1980-2019.nc",
#' package = "loopevd"),tf)
#' r = terra::rast(tf)
#' add_nc_atts(tf,r,creator_name="add_nc_atts examples")
add_nc_atts = function(outfile,r,creator_name="",creator_email="",references="",title="",summary="",keywords="",history="",licence="",Disclaimer=""){
  nc = ncdf4::nc_open(outfile,write = TRUE)

  ncdf4::ncatt_put(nc,0,"Conventions","CF-1.7")
  if(creator_name!="")  ncdf4::ncatt_put(nc,0,"creator_name",creator_name)
  if(creator_email!="") ncdf4::ncatt_put(nc,0,"creator_email",creator_email)
  if(references!="")    ncdf4::ncatt_put(nc,0,"references",references)
  if(title!="")         ncdf4::ncatt_put(nc,0,"title",title)
  if(summary!="")       ncdf4::ncatt_put(nc,0,"summary",summary)
  if(keywords!="")      ncdf4::ncatt_put(nc,0,"keywords",keywords)
  if(licence!="")        ncdf4::ncatt_put(nc,0,"licence",licence)
  if(history!="")        ncdf4::ncatt_put(nc,0,"history",history)
  if(Disclaimer!="")     ncdf4::ncatt_put(nc,0,"Disclaimer",Disclaimer)

  ncdf4::ncatt_put(nc,0,"grid_mapping_name","latitude_longitude")
  ncdf4::ncatt_put(nc,0,"geospatial_bounds_crs","EPSG:4326")
  ee = as.vector(terra::ext(r))
  ncdf4::ncatt_put(nc,0,"geospatial_bounds",paste0(ee,collapse = ","))
  ncdf4::ncatt_put(nc,0,"geospatial_lat_min",ee[3])
  ncdf4::ncatt_put(nc,0,"geospatial_lat_max",ee[4])
  ncdf4::ncatt_put(nc,0,"geospatial_lon_min",ee[1])
  ncdf4::ncatt_put(nc,0,"geospatial_lon_max",ee[2])
  ncdf4::nc_close(nc)
}


#' Convert Data Frame to NetCDF
#'
#' This function converts a data frame to a NetCDF file for a list of points (rows are different stations).
#' It attempts to identify CF-compliant coordinate variables, such as latitude and longitude, using default or specified column names.
#'
#' @param df A data frame containing the data to be converted.
#' @param output_file The path to the output NetCDF file.
#' @param lat_var The name of the latitude variable in the data frame. Default is "lat".
#' @param lon_var The name of the longitude variable in the data frame. Default is "lon".
#' @param global_atts A list of global attributes to add to the NetCDF file. Default is an empty list.
#' @param units A list of units for each variable in the data frame. Default is an empty list.
#' @param longnames A list of long names for each variable in the data frame. Default is an empty list.
#' @param cf_standard_names A list of CF standard names for each variable in the data frame. Default is an empty list.
#'
#' @return None. The function writes the data to a NetCDF file.
#'
#' @import ncdf4
#' @import terra
#' @importFrom utils modifyList
#'
#' @examples
#' # Example data frame
#' example_df <- data.frame(
#'   lat = c(-35.0, -34.5, -34.0),
#'   lon = c(150.0, 150.5, 151.0),
#'   variable1 = c(0.1, 0.2, 0.3),
#'   variable2 = c(0.4, 0.5, 0.6)
#' )
#'
#' # Define units, longnames, and CF standard names
#' units_list <- list(variable1 = "m", variable2 = "cm")
#' longnames_list <- list(variable1 = "Variable 1 Longname",
#'                        variable2 = "Variable 2 Longname")
#' cf_standard_names_list <- list(variable1 = "sea_surface_height",
#'                                 variable2 = "sea_water_temperature")
#' tf <- tempfile("test.nc")
#' # Convert example data frame to NetCDF
#' df_to_netcdf(example_df, tf, global_atts = list(
#'   title = "Example NetCDF",
#'   summary = "This is a test NetCDF file created from an example data frame.",
#'   source = "Example data",
#'   references = "N/A"
#' ), units = units_list, longnames = longnames_list,
#'   cf_standard_names = cf_standard_names_list)
#'
#' @export
df_to_netcdf <- function(df, output_file, lat_var = "lat", lon_var = "lon",
                         global_atts = list(), units = list(), longnames = list(), cf_standard_names = list()) {
  # Check if latitude and longitude columns exist
  if (!(lat_var %in% names(df))) {
    stop("latitude variable not found in the data frame.")
  }
  if (!(lon_var %in% names(df))) {
    stop("longitude variable not found in the data frame.")
  }

  # Extract latitude and longitude
  latitudes <- df[[lat_var]]
  longitudes <- df[[lon_var]]

  # Remove latitude and longitude from the data frame for variable definitions
  df <- df[, !(names(df) %in% c(lat_var, lon_var))]

  # Define dimensions
  dim_station <- ncdf4::ncdim_def("station", "station index", seq_len(nrow(df)))

  # Define variables
  var_list <- list(
    ncdf4::ncvar_def(lat_var, "degrees_north", dim_station, longname = "latitude"),
    ncdf4::ncvar_def(lon_var, "degrees_east", dim_station, longname = "longitude")
  )

  for (var_name in names(df)) {
    unit <- ifelse(var_name %in% names(units), units[[var_name]], "unitless")
    longname <- ifelse(var_name %in% names(longnames), longnames[[var_name]], var_name)
    var <- ncdf4::ncvar_def(var_name, unit, dim_station, longname = longname)
    var_list <- c(var_list, list(var))  # Use list() to add a new element to var_list
  }

  # Create the NetCDF file
  nc_out <- ncdf4::nc_create(output_file, var_list)

  # Write the data to the NetCDF file
  ncdf4::ncvar_put(nc_out, var_list[[1]], latitudes)
  ncdf4::ncvar_put(nc_out, var_list[[2]], longitudes)

  for (i in 1:ncol(df)) {
    ncdf4::ncvar_put(nc_out, var_list[[i + 2]], df[[i]])
    var_name <- names(df)[i]
    if (var_name %in% names(cf_standard_names)) {
      ncdf4::ncatt_put(nc_out, var_list[[i + 2]]$name, "standard_name", cf_standard_names[[var_name]])
    }
  }

  # Add global attributes (CF-compliant metadata)
  default_global_atts <- list(
    title = "Converted Data Frame to NetCDF",
    Conventions = "CF-1.8",
    institution = "Your Institution",
    source = "Data frame",
    history = paste(Sys.Date(), "Created"),
    references = "N/A",
    comment = "Converted using df_to_netcdf function in R"
  )

  global_atts <- modifyList(default_global_atts, global_atts)

  for (att_name in names(global_atts)) {
    ncdf4::ncatt_put(nc_out, 0, att_name, global_atts[[att_name]])
  }

  # Close the NetCDF file
  ncdf4::nc_close(nc_out)
}

#' Read Station-based evd NetCDF File into a Data Frame
#'
#' This function reads a NetCDF file with EVD parameters (e.g. loc and scale) with a station dimension and converts it
#' back into a data frame. The function assumes that the NetCDF file has a
#' 'station' dimension with associated variables such as latitude, longitude,
#' and other station-specific data.
#'
#' @param netcdf_file The path to the NetCDF file.
#' @param exclude_cov Logical, if TRUE, variables starting with 'cov' will be excluded. Default is FALSE.
#' @return A data frame containing the data from the NetCDF file.
#' @import ncdf4
#'
#' @examples
#' tf = tempfile("test.nc")
#' # Ten records of 20 random data generated from the fgumbel EVD
#' am = lapply(1:10, function(x) evd::rgumbel(20))
#' tab = as.data.frame(t(sapply(am,function(x) evd_params(x,"fgumbel"))))
#' tab$lon = rnorm(10,sd=10) #station latitude
#' tab$lat = rnorm(10,sd=20) #station longitude
#' loopevd::df_to_netcdf(df = tab,output_file = tf)
#' tab2 = loopevd::netcdf_to_df(tf)
#'
#' @export
netcdf_to_df <- function(netcdf_file, exclude_cov = FALSE) {

  # Open the NetCDF file
  nc <- nc_open(netcdf_file)

  # Get the number of stations
  num_stations <- nc$dim$station$len

  # Initialize a list to store the data
  data_list <- list()

  # Extract the data for each variable associated with the station dimension
  for (var_name in names(nc$var)) {
    if (!(exclude_cov && startsWith(var_name, "cov"))) {
      data <- ncvar_get(nc, var_name)
      data_list[[var_name]] <- data
    }
  }

  # Convert the list to a data frame
  df <- as.data.frame(data_list)

  # Close the NetCDF file
  nc_close(nc)

  return(df)
}

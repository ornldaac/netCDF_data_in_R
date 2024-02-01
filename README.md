# How to Open and Work with NetCDF Data in R  

### *Presented by the ORNL DAAC*  <https://daac.ornl.gov>  

### *March 26, 2018*  

### *Keywords: R, netCDF*  

***
**UPDATE:** The repo was updated on 2021-12-03 and again on 2024-01-02. The recent update includes a revision of the tutorial to replace deprecated packages (e.g., rgdal and raster).

# 1. Overview

In this tutorial we will open some geospatial data that is stored in a netCDF file. Then we will select the variable and time range of interest and export the data as a GeoTIFF file so that we can continue the analysis in R or other geospatial software.

![Arctic Growing Season NDVI, 1982](GIMMS3g_1982_NDVI.png)

# 2. Dataset

**Long-Term Arctic Growing Season NDVI Trends from GIMMS 3g, 1982-2012**

The dataset used in this tutorial provides normalized difference vegetation index (NDVI) data for the Arctic growing season derived primarily with data from Advanced Very High Resolution Radiometer (AVHRR) sensors on board several NOAA satellites over the years 1982 through 2012. The NDVI data, which show vegetation activity, were averaged annually for the arctic growing season (June, July, and August) in each year. The data are circumpolar in coverage at 8-kilometer resolution and limited to greater than 20 degrees North.

Guay, K.C., P.S.A. Beck, and S.J. Goetz. 2015. Long-Term Arctic Growing Season NDVI Trends from GIMMS 3g, 1982-2012. ORNL DAAC, Oak Ridge, Tennessee, USA. <a href="<<<https://doi.org/10.3334/ORNLDAAC/1275>>>">https://doi.org/10.3334/ORNLDAAC/1275</a>

Specifically, we will use the file “gimms3g_ndvi_1982-2012.nc4”. With a single click, download the data here <https://daac.ornl.gov/daacdata/global_vegetation/GIMMS3g_NDVI_Trends/data/gimms3g_ndvi_1982-2012.nc4> before beginning the tutorial.

# 3. Prerequisites

Participants should have a basic understanding of R.  

## 3.1 R

1. [Download R](https://cran.r-project.org/)  
2. [Download RStudio](https://www.rstudio.com/products/rstudio/download/#download)  *Recommended*  
3. [Review R Manuals](https://cran.r-project.org/manuals.html)  *Recommended*  

# 4. Procedure

## 4.1 Tutorial  

1. [R Markdown](netCDF_in_r_ornldaac_tutorial.Rmd)  
2. [Markdown](netCDF_in_r_ornldaac_tutorial.md)  

# 5. Credits

* [R](https://www.r-project.org/) - 4.3.2 (2023-10-31) -- "Eye Holes"  
* [RStudio](https://www.rstudio.com/products/rstudio/) - IDE and notebook construction  
* [ggplot2](https://CRAN.R-project.org/package=ggplot2) - plot figures (R package)
* [ncdf4](https://CRAN.R-project.org/package=ncdf4) - netcdf manipulation (R package)  
* [terra](https://CRAN.R-project.org/package=terra) - geospatial data manipulations (R package)

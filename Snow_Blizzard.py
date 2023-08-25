#| Risk workflow for snow and blizzard hazards 

#| Blizzard  

#| A blizzard is a severe storm condition defined by low temperature, sustained wind or frequent wind gust and considerable precipitating or blowing snow. For blizzard conditions we propose the use of following impact indicator: Tmean  ≤ 0 oC, Rs (snow amount) ≥ 10 cm and Wg (wind gust) ≥ 17 m/s ( Vajda et al., 2014). This impact indicator was defined taking into account the exposure of critical infrastructure, i.e., roads, rails, power lines, telecommunication to the hazard and is based on an extensive literature review, media reports, surveys conducted with European CI operators and case studies. 
#| Heavy snowfall may cause many disruptions and impacts in various sectors; however, the impacts and consequences of this hazard depend on the affected sector, infrastructure and also preparedness of society that varies over Europe.  For example, already a few centimeters of snow can disrupt road traffic, but doesn’t normally cause any harm to energy infrastructure. Many weather services have three warning levels based on the severity of expected impacts, which are typically different for different sectors of infrastructure. There is a large variation in the national warning criteria or thresholds.

#| Heavy Snow 
#| Similarly to blizzard, the impact indicators for heavy snowfall were defined taking into account the exposure of critical infrastructure, i.e., roads, rails, power lines, telecommunication to the hazard and is based on an extensive literature review, media reports, surveys conducted with European CI operators and case studies. The qualitative description of the two-level thresholds are:
#| 1st threshold: Some adverse impacts are expected, their severity depends on the resilience of the system, transportation is mainly affected.
#| 2nd threshold: The weather phenomena are so severe that is likely that adverse impact will occur, CI system is seriously impacted.

#| This code calculates the Annual probability (%) of a blizzard and heavy snowfall during the specified period and a region of interest.
#| Please note that we used ERA5-Land data here, which is locally stored. This daily dataset emanated from the hourly datasets available on the Climate Data Store. 

#| Load libraries
import os
from glob import glob
import numpy as np
import rasterio
import xarray as xr
import rioxarray as rxr
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xclim as xc



#| Define thresholds 
lim_tas=0.       # Deg C
lim_pr=10.       # mm
lim_gust=17.     # m/s
lim_snow = 10.   # cm
lim_snow6 = 6.   # cm
lim_snow25 = 25. # cm


input_dat = "/Users/suraj_fmi/D_Drive/Prog_FMI/dmin_CLIMAAX/"
output_dar = "/Users/suraj_fmi/D_Drive/Prog_FMI/dmin_CLIMAAX/"

#| Read a Datasets 
#|
#| Read the U and V components of Wind: Unit (m/s)

u10 = xr.open_mfdataset(os.path.join(input_dat, "10m_u_component_of_wind_DMIN_era5Land_202*.nc"),chunks='auto').u10
v10 = xr.open_mfdataset(os.path.join(input_dat, "10m_v_component_of_wind_DMIN_era5Land_202*.nc"),chunks='auto').v10

#| Calculate wind speed from U and V component (m/s)
wspd = np.sqrt((u10 * u10) + (v10 * v10))
wspd = wspd.assign_attrs(units="m/s", description="Wind Speed")
del u10, v10

#| Read Temperature: Unit (K)

tas = xr.open_mfdataset(os.path.join(input_dat, "2m_temperature_DMIN_era5Land_202*.nc"),chunks='auto').t2m
tas = tas.assign_attrs(units="K", description="Temperature data")
tas = xc.units.convert_units_to(tas, "degC")


#| Read snow data:  Unit (M)
snow = xr.open_mfdataset(os.path.join(input_dat, "snow_depth_DMIN_era5Land_202*.nc"),chunks='auto').sde
snow = xc.units.convert_units_to(snow, "cm")

#! Annual sum of blizzard days 

BdayCount_annual = ((tas < lim_tas) * (snow > lim_snow) * (wspd > lim_gust) * 1).groupby('time.year').sum(dim='time')
BdayCount_annual = BdayCount_annual.where(BdayCount_annual != 0.)
BdayCount_annual = BdayCount_annual.assign_attrs(units="number", long_name="Annual number of blizzard days")
BdayCount_annual = BdayCount_annual.to_dataset(name='blizzard_days')
BdayCount_annual_mean = BdayCount_annual.mean('year')
BdayCount_annual_mean.to_netcdf(path=os.path.join(input_dat, "BdayCount_annual_mean.nc"))
del BdayCount_annual

#| Probability of occurrence of blizzard days (%)
BdayCount_anaProb = ((tas < lim_tas) * (snow > lim_snow) * (wspd > lim_gust)).groupby('time.year').mean('time')
del tas, wspd, pr
BdayCount_anaProb = BdayCount_anaProb.where(BdayCount_anaProb != 0.)
BdayCount_anaProb = BdayCount_anaProb.assign_attrs(units="%", long_name="Annual probability of blizzard days")
BdayCount_anaProb = BdayCount_anaProb.to_dataset(name='blizzard_days')
BdayCount_anaProb_mean = BdayCount_anaProb.mean('year')
BdayCount_anaProb_mean = BdayCount_anaProb_mean * 100.
BdayCount_anaProb_mean.to_netcdf(path=os.path.join(input_dat, "BdayCount_AnaProb_mean.nc"))
del BdayCount_anaProb

#| Annual sum of Heavy snowfall days for snowfall > 6cm

snow6Count_annual = ((snow > lim_snow6) * 1).groupby('time.year').sum('time')
snow6Count_annual = snow6Count_annual.where(snow6Count_annual != 0.)
snow6Count_annual = snow6Count_annual.assign_attrs(units="number", long_name="Annual number of snow days")
snow6Count_annual = snow6Count_annual.to_dataset(name='snow_days')
snow6Count_annual_mean = snow6Count_annual.mean('year')
snow6Count_annual_mean.to_netcdf(path=os.path.join(input_dat, "snow6Count_annual_mean.nc"))
del snow6Count_annual


#| Annual sum of Heavy snowfall days for snowfall > 25cm

snow25Count_annual = ((snow > lim_snow25) * 1).groupby('time.year').sum('time')
snow25Count_annual = snow25Count_annual.where(snow25Count_annual != 0.)
snow25Count_annual = snow25Count_annual.assign_attrs(units="number", long_name="Annual number of snow days")
snow25Count_annual = snow25Count_annual.to_dataset(name='snow_days')
snow25Count_annual_mean = snow25Count_annual.mean('year')
snow25Count_annual_mean.to_netcdf(path=os.path.join(input_dat, "snow25Count_annual_mean.nc"))
del snow25Count_annual

#| Probability of occurrence in % for snowfall > 6cm

snow6Prob_annual = (snow > 6.).groupby('time.year').mean('time')
snow6Prob_annual = snow6Prob_annual.where(snow6Prob_annual != 0.)
snow6Prob_annual = snow6Prob_annual.assign_attrs(units="%", long_name="Annual probability of snow days")
snow6Prob_annual = snow6Prob_annual.to_dataset(name='snow_days')
snow6Prob_annual = snow6Prob_annual * 100.
snow6Prob_annual_mean = snow6Prob_annual.mean('year')
snow6Prob_annual_mean.to_netcdf(path=os.path.join(input_dat, "snow6Prob_annual_mean.nc"))
del snow6Prob_annual

#| Probability of occurrence in % for snowfall > 25cm

snow25Prob_annual = (snow > 25.).groupby('time.year').mean('time')
snow25Prob_annual = snow25Prob_annual.where(snow25Prob_annual != 0.)
snow25Prob_annual = snow25Prob_annual.assign_attrs(units="%", long_name="Annual probability of snow days")
snow25Prob_annual = snow25Prob_annual.to_dataset(name='snow_days')
snow25Prob_annual = snow25Prob_annual * 100.
snow25Prob_annual_mean = snow25Prob_annual.mean('year')
snow25Prob_annual_mean.to_netcdf(path=os.path.join(input_dat, "snow25Prob_annual_mean.nc"))
del snow25Prob_annual

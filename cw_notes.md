### Notes for Chris's benefit.  
I estimated the mean annual rainfall and temperature for the LSC and Dobsons Sites using the following:

in QGIS for rainfall (quite a chunky raster but better than any of the available BOM temperature gauges)
/CWmacbook/data/spatial/Australia/BOM rainfall grids/meanAnnualRainfall.tif

and for mean temperature from the mwstr database for the two most downstream reaches 
`SELECT reach, meant_y30_2022 FROM subc_env WHERE site IN (61638, 71416);`
(see [the stream manual](https://tools.thewerg.unimelb.edu.au/mwstr_manual/c4-env_tabs.html#air-temperature) for the derivation of meant_y30_2022: I figured the mean for the last 30 years is most appropriate, given the 1.5 increase in mean temperature over the last century.


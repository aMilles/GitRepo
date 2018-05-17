CmbTransects_segmentize.R		combines all transects into one shapefile and creates equally spaced polygons parallel to the transect lines
					output: segments.shp and transects.shp

combine_FSO.R 				tries to harmonize FSOD, however they are currently not usable
					output: none

createafrica_shape.R			creates a continent and country level shape of africa
					output: africa_countries.shp and africa_continent.shp

GLW_creation.R				uses Gridded Livestock of the World data to calculate Tropical Livestock Units (TLU) as an Indicator of Livestock Density (LD)
					output: tlu.tif

harmonize_RSOD.R			handles issues in single GEC datasets, harmonizes observation codes, times (converts to utc) and position (EPSG 4326), the final product is a combined set of spatial points with information of RSOD
					output: GEC_points.shp

PA.R					CURRENTLY NOT APPLICABLE: creates a .shp of African protected areas (PA) by using WDPA data
					output: PA_africa.shp

SgmtExtract_PA-GLW.R			uses segments.shp to extract information about LD and PA as well as observation counts and observation types by including nearby points in GEC_points.shp
					output: segments.shp with added information
# Japan_Bat_Supplemental
### Supplementary data of "Scale-dependent influences of environmental, historical, and spatial processes on beta diversity of Japanese bat assemblages"
Takahiro Maki, Nozomi Sannomiya, Toshihide Hirao and Dai Fukui*

This data package contains the original data and code of statistical analysis in "Scale-dependent influences of environmental, historical, and spatial processes on beta diversity of Japanese bat assemblages"

# CONTENTS
### coo_Entire.shp and coo_Main.shp
##### Coordinates of centroids in land area within each mesh or island at Entire archipelago scale and 4 main islands scale.

### env_Entire.csv and env_Main.csv
##### Dataset of environmental factors at Entire archipelago scale and 4 main islands scale.
"Rplantf" means the ratio of plantation in forest area.
"Ranthrof" means the ratio of artificial environment including agricultural land and human residential area in terrestrial area.
"temp" means the annual mean temperature.
"pre" means the annual precipitation. 
"mele" means the mean altitude.
"mesl" means the mean slope.
"area 1" means land area of belonging island in common logarithm.

### his_Entire.csv and his_Main.csv
##### Dataset of historical factors at Entire archipelago scale and 4 main islands scale.
"BA_YDSt" means difference of mean annual temperature between Bølling-Allerød interstadial and Younger Dryas stadial. 
"BA_YDSp" means difference of annual precipitation between Bølling-Allerød interstadial and Younger Dryas stadial. 
"YDS_Nt" means difference of mean annual temperature between Younger Dryas stadial and current time. 
"YDS_Np" means difference of annual precipitation between Younger Dryas stadial and current time. 
"tokara" was assigned a value of 1 to islands south of the Tokara Strait, while the rest were assigned a value of 0 
"tsugaru" was assigned a value of 1 to grids and islands north of the Tsugaru Strait, while the rest were assigned a value of 0

### bat_dist_Entire and bat_dist_Main
##### Distribution data of Japanese bats in presence/absence data at Entire archipleggo scale and 4 main islands scale.
Rows mean the assemblages of grids or islands. 
columns mean the species of bats.

### VarPar_Entire and Varpar_Main
##### R code of the statistical analysis

# Marine fish climate responses

Marine fish abundance responses to climate change paper

Data and code for Mangnall et al. (unpublished)

***

**assets**

***

The following original datasets were used for analyses: See Table S1 for details:

Alcatrazes_monitoring_program_660.csv.gz  
CRED_Fish_Surveys_in_the_Pacific_Ocean.csv.gz  
CSIRO_Marine_Data_Warehouse.csv.gz  
ECNASAP_East_Coast_North_America_Strategic_Assessment.csv.gz  
FishGlob_.csv.gz  
IMOS_Cryptobenthic_Abundance.csv.gz  
IMOS_NZ_Fish.csv.gz  
IMOS_Reef_Abundance.csv.gz  
Marine_Fish_Underwater_Surveys_in_The_Israeli_Mediterranean.csv.gz  
MARMAP_Chevron_Trap_Survey_1990_2009.csv.gz  
Pacific_Shrimp_Trawl_Survey(OBIS_Canada).csv.gz  
Pelagic_Fish_Observations_1968-1999.csv.gz  
Previous_fisheries_REVIZEE_Program_135.csv.gz  
Previous_fisheries_REVIZEE_Program_284.csv.gz  
Previous_fisheries_REVIZEE_Program_285.csv.gz  
St_Croix_USVI_Fish_Assessment_and_Monitoring_Data.csv.gz

The following processed data were used for analyses"

Data.csv.gz Bound survey data, with corrected names and filted by depth (0-200m)  
HADISST.csv.gz Processed relevant temperature data sourced from the Hadley Centre https://www.metoffice.gov.uk/hadobs/hadisst/  
gbif_data_northern.csv.gz Downloaded records for surveyed species from GBIF  
gbif_data_southern.csv.gz Downloaded survey for surveyed species from GBIF abundances.csv.gz  
Intermediary processing file, abundance of each species, per survey, grid cell and year, with temperature data  
final_abundance_changes.csv.gz Final analysed abundance time series data, filtered per survey, grid cell and year, with temperature data.

***

**scripts**

***

1_Get_Data.R Code for compiling abundance data  
2_Get_Climate.R Code for linking temperature data to abundance data  
3_Get_GBIF.R Code for obtaining the distributional records and range limits of species  
4_Get_Abundances.R Code for obtaining abundance trends of species  
5_Get_AnalysisA.R Code for plotting time series data, and modelling overall temperature changes  
6_Get_AnalysisB.R Code for plotting and modelling abundance x temperature responses  
7_Get_AnalysisC.R Code for plotting and analyses of the influence of life history traits  
8_Get_AnalysisD.R Code for exploring patterns of of change in local communities

***

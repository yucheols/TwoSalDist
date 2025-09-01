# Ecological niche modeling of two distantly related lungless salamanders  
Ecological niche modeling and associated analyses for the two species of Korean endemic salamanders.
*Some of the custom R scripts used in this project have been updated and included in the R package ENMwrap (https://github.com/yucheols/ENMwrap)*

## Software and package dependencies
#### Core analyses
- R (version 4.2.2)
- ENMeval (version 2.0.4)
- dismo (version 1.3.14)
- ntbox (version 0.6.6.6)
- humboldt (version 1.0.0.420121)
- ENMTools (version 1.0.5)
- megaSDM (version 2.0.0)
- blockCV (version 3.0.0)
- plyr (version 1.8.8)
- dplyr (version 1.1.0)
- readr (version 2.1.4)
- raster (version 3.6.14)
- terra (version 1.7.65)
- rgdal (version 1.6.4)
- sf (version 1.0.14)
- tmap (version 3.3.3)
- MASS (version 7.3.58.2)
- rasterVis (version 0.51.5)
- ggplot2 (version 3.4.1)
- ggpubr (version 0.5.0)
- patchwork (version 1.1.2)

#### Additional analyses
- R (version 4.4.2)
- ENMeval (version 2.0.5)
- ENMwrap (version 1.0.1)
- dplyr (version 1.1.4)
- terra (version 1.8.21)
- rgdal (version 1.6.7)
- ggplot2 (version 3.5.1)



## Study background
![Figure_1_compressed](https://github.com/user-attachments/assets/74ca0052-7ef8-49df-a91f-62df9a5ab8a8)

- The Korean Clawed Salamander (*Onychodactylus koreanus*) and Korean Crevice Salamander (*Karsenia koreana*) are endemic to the Korean Peninsula.
- *O. koreanus* is in the family Hynobiidae and *K. koreana* is in the family Plethodontidae (it is the only Asian representative of this family).
- Although distantly related, these two species are ecologically similar and thus share similar habitat requirements. Both species are lungless and are strict habitat specialists found exclusively in heavily forested areas close to clean mountain streams.
- They are also found in broad sympatry across their known ranges.
- As these two species are both lungless, their ranges are closely tied to very speecific envieonmental conditions, and their current ranges likely have been shaped by the Pleistocene climatic fluctuations.

## Occurrence dataset
- *O. koreanus*: Occurrence points retrieved from a previous study (Shin et al., 2021. Ecology and Evolution 11: 14669-14688), which were compiled from survey records, GBIF, museum records, literature etc.
- *K. koreana*: Occurrence points compiled from a previous study (Jeon et al., 2021. Scientific Reports 11: 9193), survey records, GBIF, etc.
- Occurrence points for both species were spatially thinned to have only one point in every 1km spatial resolution raster pixel.

## Environmental dataset
- WorldClim bioclimatic variables
- CHELSA bioclimatic variables
- EarthEnv global consensus forest cover layers, with further processing in ArcGIS Pro v2.6.0
- PaleoClim data for Plio-Pleistocene climates

## Background dataset
- Two levels of spatial bias correction
- Three different sample sizes (*n*=5000, *n*=10000, *n*=15000)

## Model tuning experiments and model evaluation
- Model parameter tuning with ENMeval
- Model cross-validation using the hierarchical spatial checkerboard method
- Model evaluation using statistical metrics (e.g. AUC, CBI, omission rate) and testing against null niche models

## Current habitat suitability
![current](https://github.com/yucheols/TwoSalDist/assets/85914125/edf19032-9a3d-46d1-b9cc-84d67219a6e2)

## Hindcasting ENMs
![mapping](F:/Projects_completed/TwoSalDist/MS/Submission_BMC_Ecology_and_Evolution/Proof/Fig_3_revised_compressed.png)
- ENMs hindcasted to five different time periods of Plio-Pleistocene

![MESS](<img width="12598" height="4724" alt="Figure_S2_revised" src="https://github.com/user-attachments/assets/7edf0651-89f2-476b-971c-1bf2bf0a438c" />)
- Extrapolation risk assessed with Multivariate Environmental Similarity Surface (MESS)

## Niche analyses
![Niche analyses](https://github.com/user-attachments/assets/475537ef-6ea4-4ce9-b734-f16a641b1c48)

## Citation
```
Y Shin, A Borzée, D Park. 2025. Climatic data sources and limitations of ecological niche models impact the estimations of historical ranges and niche overlaps in distantly related Korean salamanders. BMC Ecology and Evolution. In press.
```


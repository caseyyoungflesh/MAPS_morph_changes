# MAPS_morph_changes 


[![Zenodo DOI](https://zenodo.org/badge/393847511.svg)](https://zenodo.org/badge/latestdoi/393847511)


Code for Youngflesh et al. 2022 _Nature Ecology and Evolution_

This repository contains code to characterize avian morphological change over time, latitude, elevation, and in response to temperature.

**Associated publications:**

Youngflesh, C, JF Saracco, RB Siegel, MW Tingley. 2022. [Abiotic conditions shape spatial and temporal morphological variation in North American birds.](https://www.nature.com/articles/s41559-022-01893-x) **_Nature Ecology and Evolution_** 6:1860-1870.


**Repository structure:**

* `Scripts/`
  - `1-process-data.R`
  - `2-process-elev.R`
  - `3-process-daymet/`
      - `3-proces-daymet.R`
      - `daymet_query.txt`
      - `daymet_spt_batch.sh`
  - `4-si-tle.R`
  - `5-wi-tle.R`
  - `6-si-temp-space.R`
  - `7-si-temp-lag.R`
  - `8-temp-cov.R`
  - `9-backtransform.R`
  - `10-analyze-results.R`
  - `11-haldanes.R`
  - `model_files`
    - `morph-temp-space.stan`
    - `morph-temp-ss2.stan`
    - `morph-tle.stan`
    - `temp_sp_cov.stan`
* `Data/`
  - `bird_phylo/`
    - `phylo_names_key-2021-05-19.csv`
    - `tree-pruner-XXXX/` [ignored - available [here](http://www.birdtree.org)]
  - `GMTED2010_mean/` [ignored - available [here](https://www.usgs.gov/core-science-systems/eros/coastal-changes-and-impacts/gmted2010?qt-science_support_page_related_con=0#qt-science_support_page_related_con)]
  - `MAPS_data.csv` [ignored - available [here](DRYAD LINK)]
  - `Gen_time_Bird_et_al_2020_table_4.csv` [ignored - available [here](https://conbio.onlinelibrary.wiley.com/doi/abs/10.1111/cobi.13486)]
  - `Hendry_et_al_2008_Mol_Eco.csv` [ignored - available [here](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-294X.2007.03428.x?casa_token=xSDmUC85ERkAAAAA%3Ap5kaUUhK9w2wV7MiVqOvU52Eo-2sjFA1_h977u6L0YiXTcEuhPIfPwU8G1H8Q3CgoCSdPtJBAwOHFEdp)]
* `Results/` [ignored]
* `Figures/`
  - `PPO/` - prior posterior overlap

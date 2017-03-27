# cFOS OXT Image Analysis
This GitHub repository provides the scripts to processing and analysis scripts for OXT/cFOS/DAPI-stained immunofluorescent image stacks acquired from Zeis LSM700 confocal microscope. Data from this analysis is published by **Maynard *et al.*, in prep**, *insert link to publication*. 

## Processing 
* Project phenotype data is generated using the **make_pheno.R** script by parsing the precisely named file names.

* Primary processing and data extraction is done using a custom MATLAB script **cellSegCFOS.m**. This custom function was ran with MATLAB 2016b and requires the *CellSegM* Matlab toolbox (https://github.com/ehodneland/cellsegm). 

* Batch execution and summarization of extracted data into the R environment requires running the R script **preprocess_cFOS_v2.R** on an SGE compute cluster. Processed data files can be requested from the authors.

## Analysis
Main statistical analysis for this dataset occurs in the Rscript **analyze_cFos_v2.R**. 

### Project design
* Nested image dataset (3 wild-type vs. 3 mutant animals, 6-9 image stacks per animal, N=48 total)
* Image are tile-stitched in x & y, retaining z information
* PVN area are image-dependent and must be computed to adjust for # of OXT (+) cells in view

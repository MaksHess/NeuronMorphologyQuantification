# NeuronMorphologyQuantification
Code accompanying the manuscript "In-vivo quantitative image analysis of age-related morphological changes of C. elegans neurons reveals a correlation between neurite bending and novel neurite outgrowths" by Hess et al.

# Data
You can download the full dataset from [Zenodo](https://zenodo.org/record/2350066#.XBd7YWhKiUk)

# Instructions
1. Clone the repository and download the full data from [Zenodo](https://zenodo.org/record/2350066#.XBd7YWhKiUk)
2. In batchProcessing.py set the parameter root to the directory containing the data
3. Run the script batchProcessing.py

# Output
## Measurement Files:
- root/yyyy-mm-dd_hh-mm-ss_IndividualMeasurements.csv:\
  Contains individual length measurements of structures (i. e. somaoutgrowths)
- root/yyyy-mm-dd_hh-mm-ss_SummaryMeasurements.csv:\
  Contains mesaruements corresponding to individual neurons (i. e. sharp-bend counts)
- root/yyyy-mm-dd_hh-mm-ss_Parameters.txt:\
  Contains a copy of the parameters used
## Visualisations
- root/(ALM|PLM)/classifiedtrees/\
  This folder contains .swc files with classified nodes.
- root/(ALM|PLM)/wavytrees/\
  This folder contains .swc files with sharp bends and beads visualised.

# DeepRead
The Rscript and jupyter Notebooks associated with the paper "Convolutional neural networks assisted peak classification in targeted UHPLC-HRMS/MS for equine doping control screening analyses".

The Rscript shows the reprocess of raw data, using the package Rawrr, to extract chromatograms and produce inputs for the CNN and LDA models. It generates chromatograms as JPEG images for the two most intense transitions of each targeted drug, a CSV file with peak information (number of points, area, relative retention time deviation) and a PDF report for each sample containing all chromatograms.

One jupyter notebook illustrates the train and test process for the CNN and LDA models, and the other one presents an example of prediction based on a set of anonymized inputs.

![workflow](https://github.com/user-attachments/assets/d64378de-50ce-4cab-ab1d-abd53377a5c8)

*Figure: Worflow of the AI-based tool from raw data to CNN scoring and LDA classification.*

## R Dependencies
* tidyverse
* Rawrr (https://github.com/fgcz/rawrr)
* ggplot2
* cardidates
* cowplot
* doParallel
* foreach

## Python Dependencies
* sklearn
* numpy/pandas
* keras
* tensorflow
* cuda
* matplotlib

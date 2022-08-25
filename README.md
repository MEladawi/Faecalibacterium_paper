# Differential analysis code

This repository contains the differential analysis R code used in the paper titled: "Human commensal Faecalibacterium prausnitzii ameliorates chronic kidney disease in mice through the gut microbiota-butyrate-renal GPR43 axis"


## Installing the required packages
To install the required packages, use the following code:
```
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
```

## To reproduce the results
* Replace CodeDirectory in the code line 23 with the directory of the Diff_analysis.R file on your machine
* Prepare two csv files for each dataset named as "[ds_name]_data.csv" & "[ds_name]_meta.csv" for count and metadata respectively. For example, for the dataset "CKD-AmericanHuman," the required files are CKD-AmericanHuman_data.csv & CKD-AmericanHuman_meta.csv
* Place the data files in the same directory of the Diff_analysis.R file.

## Output results:
A folder for each dataset (with the same dataset name) will be created, and the corresponding output will be generated inside these folders.
For example, for the dataset "CKD-AmericanHuman," a folder with the name "CKD-AmericanHuman" will be created. The output files will be as follows:
* normalized counts.csv: for normalized counts using DSeq2 normalization factors.
* for each contrast, a file will be generated named [contrast]_de_ALL.csv that contains the differential analysis data. For instance, for a contrast between "CKD" & "Control," the file will be named "CKD_vs_Control_de.ALL.csv"



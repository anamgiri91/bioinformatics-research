# OutlierMeth


The OutlierMeth R package was designed to standardize an approach to flag aberrant outliers using Illumina Human 450K array DNA methylation data.
---

## Changelog
[20220831] v0.01.0:OutlierMeth.
[20240814] v0.01.1:OutlierMeth.
[20250330] v1.02.1:OutlierMeth.

---

## Installation

The preferred method for downloading this script file is through github.


#### Step 1: Download the OutlierMeth R package.

Please download the R Package [OutlierMeth_1.2.1.tar.gz] 

https://github.com/bdowns4/BinaryMeth/blob/master/OutlierMeth_1.2.1.tar.gz


#### Step 2: Install the R the package.

In R, install the package with the function:

install.packages("~/Path/OutlierMeth_1.2.1.tar.gz", repos = NULL, type = "source")


#### Step 3:Load the R the package.

library(OutlierMeth)


## Usage

The four functions in the OutlierMeth package are **flagMeth**, **deltMeth**, **relMeth** and **referenceMeth**.

-The function **flagMeth** binarize flagged outliers using normal tissue reference thresholds from Illumina Human 450K array DNA methylation data.

-The function **deltMeth** calculates the ∆Meth of flagged outliers using beta values from normal tissue reference thresholds.

-The function **relMeth** calculates the relative ∆Meth of flagged outliers using beta values from normal tissue reference thresholds.

-The function **referenceMeth** makes new reference threshold datasets.

When flagging the outliers the DNA methylation data, the user can select the significance levels **0.01**, **0.001**, **0.0001**, and **0.00001**, for four collections of reference samples: **tcga**, **geo**, **all** (tcga-geo) combined, and **blood**.

##Note: The references geo and all no longer include dataset GSE127857

## Citation
Downs BM, Thursby SJ, Cope L. Detecting aberrant DNA methylation in Illumina DNA methylation arrays: a toolbox and recommendations for its use. Epigenetics. 2023 Dec;18(1):2213874. doi: 10.1080/15592294.2023.2213874. PubMed PMID: 37218167.

-https://doi.org/10.1080/15592294.2023.2213874

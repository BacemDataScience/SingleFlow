SingleFlow
==============
[![SingleFlow license](https://img.shields.io/github/license/hernet/SingleFlow.svg?colorB=26af64&style=popout)](https://github.com/hernet/SingleFlow/blob/master/LICENSE)

SingleFlow is a pipeline for analyzing scRNA-seq data built using the workflow
management system Nextflow.

## Installation
### Dependencies
* Compiler Java 8
* Runtime Java 8 or later
* [Nextflow](https://github.com/nextflow-io/nextflow)

#### Install Nextflow
[Nextflow](https://github.com/nextflow-io/nextflow) can be downloaded through the following command

```sh
curl -fsSL get.nextflow.io | bash
```

and subsequently moved to a location in $PATH.

Nextflow can also be installed from Conda

```sh
conda install -c bioconda nextflow
```

### Python and R packages
The Python and R packages required for the specific analysis being performed have to be
installed. The environment used when developing *SingleFlow* can be found in the
[requirements file](https://github.com/hernet/SingleFlow/blob/master/requirements.txt).
To include glm-PCA in the anlaysis see https://github.com/willtownes/scrna2019.

## Usage

```sh
$ nextflow run SingleFlow.nf --input <folder(s)>

```

### Parameters

#### Mandatory arguments
* ```--input```: Folder or comma-separated list of folders containing the scRNA-seq data for analysis

#### Optional arguments
* ```--outdir```: Folder for the results of the anlalysis to be put in, default: 'results'.
* ```--colors```: Comma-separated colors for use in plotting the results of the analysis.
* ```--custom_colors```: Comma-separated colors for use when defining custom cell clusters.
* ```--ccc```: Specify whether to perform cell cycle correction of the data set or not.
* ```--custom_clusters```: Specify whether to define custom cell clusters or not
* ```--deg```: Specify whether to perform differentially expressed gene analysis or not.
* ```--rna_velocity```: Specify whether to embed RNA velocity or not, if yes input a .loom file similar to the output of velocyto.
* ```--phenograph_clusters```: Specify whether to calculate and visualize Phenograph clusters or not.
* ```--gene_trends```: Specify whether to calculate the gene trends.
* ```--gene_trends_late```: If set, the gene trends will be calculated and clustered from the given pseudotime.
* ```--gene_expression```: Specify a file with genes whose expression levels should be visualized.
* ```--glmpca```: Use the glm-PCA method or the conventional normalization and subsequent PCA.
* ```--imputation```: Specify which imputation method to use, default: 'MAGIC'.
* ```--exclude```: Specify a file containing cells that should be excluded form the analysis.
* ```--degWin```: Show the result of DEG analysis in the specified file.
* ```--factorAnalysis```: File containing the factors used as input to the factor analysis.
* ```--help```: Show helper menu.


## Examples
If we have three samples we want to include in our analysis and they are placed in
three separate folders 'Sample 14', 'Sample 15' and 'Sample 16', within a folder
called 'Donor 2', we can include all these samples in our analysis by specifying
the input as in the following example

```sh
$ nextflow run SingleFlow.nf --input "Donor 2/Sample 14,Donor 2/Sample 15,Donor 2/Sample 16"

```

*SingleFlow* will in the subsequent analysis recognize that these comma-separated folders
contain data from different samples.

For a specific example of the application of *SingleFlow* to obtain new biological insight,
see our preprint on [NK cell differentiation](https://www.biorxiv.org/content/10.1101/630657v1).
Some of the results of the analysis is available and can be explored at http://herknet.io/nk.

## License
*SingleFlow* is released under the MIT license.

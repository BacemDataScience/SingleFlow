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

### Python packages

### R packages


## Examples
If we have three samples we want to include in our analysis and they are placed in
three separate folders 'Sample 14', 'Sample 15' and 'Sample 16', within a folder
called 'Donor 2', we can include all these samples in our analysis by specifying
the input as in the following example

```sh
$ nextflow run SingleFlow.nf --input "Donor 2/Sample 14,Donor 2/Sample 15,Donor 2/Sample 16"

```

SingleFlow will in the subsequent analysis recognize that these comma-separated folders
contain data from different samples.

## Parameters

## License
*SingleFlow* is released under the MIT license.

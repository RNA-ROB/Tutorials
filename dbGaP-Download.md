# Downloading files from dbGaP

dbGaP (database of Genotypes and Phenotypes) is a database created by NCBI to archive and distribute the results of studies that have investigated the interactions of genotype and phenotype. The data archived on dbGaP either involves (1) open access or (2) controlled access. To access dbGaP controlled-access data, researchers need to comply with the NIH's security best practices. 

## **Obtaining a dbGaP repository key:**

Once you have been granted formal access to a particular dataset on dbGaP (the PI needs to fill out a request form describing how the requested data will be used), you will first need to get a dbGaP repository key (this is a file with a <tt>.ngc</tt> extension).

Login to dbGaP and select the **'My Projects'** tab. Select your project of interest, and then select **'Get dbGaP repository key'**. The dbGaP repository key file will be subsequently downloaded.

## **Installing SRA-Toolkit:**

The SRA-Toolkit is required when downloading/unpacking any files from dbGaP. To install and compile the SRA-Toolkit, please refer to the instructions found on their GitHub Wiki, linked [here](https://github.com/ncbi/sra-tools/wiki).

## **Working with the Run Selector**

To download raw sequencing data associated with a project on dbGaP, select the **'Run Selector'** tab on the project webpage. This will take you to a new window containing a data matrix of all available sequencing data. Options for dataset filters (e.g., WES/WGS/RNA-seq) can be found on the left side of the webpage. Once you have selected your sequencing datasets of interest, you will be able to download an encrypted cart file (file with extension <tt>.krt</tt>). This cart file can then be used for retrieving the actual raw sequencing datasets with the help of the *prefetch* module in the SRA-Toolkit as follows:

```
[/path/to/prefetch/binary] --ngc [/path/to/dbGaP/repository/key] [/path/to/cart/file]
```

Once we finish executing the *prefetch* command using our cart file, we should have <tt>.sra</tt> files downloaded in the current directory. We can next use the *fastq-dump* module to convert the <tt>.sra</tt> files into compressed <tt>FASTQ</tt> files as follows:

```
# --split-files will separate forward and reverse reads for paired-end data
# -F will preserve the original sequence names in the FASTQ defline

[/path/to/fastq-dump/binary] --ngc [/path/to/dbGaP/repository/key] --split-files \
    -F --gzip -O [/path/to/output/directory] [accession ID]
```

## **Downloading the SRA Run Table**

It may be useful to download the SRA run table (a data matrix characterizing every sequencing dataset) displayed for a dbGaP project. Unfortunately, this cannot be done in the current version of the Run Selector and will require you to revert to the *Old Run Selector* (this should be a directly accessible option when you are on the dbGaP Run Selector webpage). On the Old Run Selector, you should be able to download the run table as a data matrix file by selecting **'RunInfo Table'** under **'Download'**.

## **Working with the File Selector**

To download non-sequencing but protected data from a project on dbGaP (e.g., metadata containing information on patient covariates), select the **'File Selector'** tab on the project webpage. This will take you to a new window containing another data matrix of all protected non-sequencing data files. Once you have selected your files of interest, you will be able to download an encrypted cart file (file with extension <tt>.krt</tt>). We will again use the *prefetch* command in the SRA-Toolkit (refer to the previous section) on this cart file. Once our desired files have been prefetched, we will have several encrypted files (file extension <tt>*.enc*</tt>) in the current directory. To decrypt these files, we will need to use the *vdb-decrypt* command in the SRA-Toolkit as follows:

```
[/path/to/vdb-decrypt/binary] --ngc [/path/to/dbGaP/repository/key] [/path/to/directory/with/encrypted/files]
```
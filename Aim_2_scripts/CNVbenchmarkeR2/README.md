# CNVbenchmarkeR2 #

CNVbenchmarkeR2 is a framework to benchmark germline copy number variant (CNV) calling tools on multiple NGS datasets. Current version supports DECoN, CoNVaDING, panelcn.MOPS, ExomeDepth, CODEX2, ClinCNV, clearCNV, GATK-gCNV, Atlas-CNV, Cobalt, CNVkit and VisCap tools.

Previous version, CNVbenchmarkeR, is available [here](https://github.com/TranslationalBioinformaticsIGTP/CNVbenchmarkeR).

### Citation ###

Please cite [last publication](https://academic.oup.com/bib/article/26/1/bbae645/7922578) when using CNVbenchmarker2:

Elisabet Munté, Carla Roca, Jesús Del Valle, Lidia Feliubadaló, Marta Pineda, Bernat Gel, Elisabeth Castellanos, Barbara Rivera, David Cordero, Víctor Moreno, Conxi Lázaro, José Marcos Moreno-Cabrera, *Detection of germline CNVs from gene panel data: benchmarking the state of the art*, Briefings in Bioinformatics, Volume 26, Issue 1, January 2025, bbae645, https://doi.org/10.1093/bib/bbae645


### Prerequisites ###

Tools should be properly installed. Links for tools installation:

- [Panelcn.mops](https://github.com/bioinf-jku/panelcn.mops)
- [CoNVaDING](https://molgenis.gitbooks.io/convading)
- [DECoN](https://github.com/RahmanTeam/DECoN)
- [CODEX2](https://github.com/yuchaojiang/CODEX2)
- [ExomeDepth](https://cran.r-project.org/web/packages/ExomeDepth/index.html)
- [ClinCNV](https://github.com/imgag/ClinCNV)
- [clearCNV](https://github.com/bihealth/clear-cnv)
- [Atlas-CNV](https://github.com/theodorc/Atlas-CNV)
- [Cobalt](https://github.com/ARUP-NGS/cobalt)
- [CNVkit](https://github.com/etal/cnvkit)
- [VisCap](https://github.com/pughlab/VisCap)
- [GATK-gCNV (GermlineCNVCaller)](https://hub.docker.com/r/broadinstitute/gatk) (At the moment, it only works using singularity).

Also, R/Bioconductor should be installed including these packages: GenomicRanges, biomaRt, regioneR, vcfR, optparse.

### How to use
1. Get Code
```
git clone https://github.com/jpuntomarcos/CNVbenchmarkeR2 
```

2. **Configure tools.yaml** to set which tools will be benchmarked.
  
3. **Set tool parameter values** located at tools/[name_of_the_tool]/[name_of_the_tool]Params.yaml. Please, note that **tool paths must be set** for tools that are not R packages (DECoN, Convading, ClinCNV, clearCNV, GATK-gCNV, Atlas-CNV, Cobalt, CNVkit and Viscap).

4. **Configure datasets.yaml** to define on which datasets the tools will be executed. Within this file, it is important to provide files with the exact expected format (**special attention** to `validated_results_file` and `bed_file` that are **tab-delimited** files). To do so, please **check the [examples](https://github.com/jpuntomarcos/CNVbenchmarkeR2/tree/master/examples) folder**.


5. Launch CNVbenchmarkeR2
```
cd CNVbenchmarkerR2
Rscript runBenchmark.R [-t tools_yaml] [-d datasets_yaml] [-f include_temp_files]
```


### Output ###

A summary file and a .csv results file will be generated at output/summary folder. Stats include sensitivity, specificity, no-call rate, precision (PPV), NPV, F1, MCC and kappa coefficient.

Statistics are calculated per ROI, per gene and at whole strategy level.

Logs files will be generated in the logs folder. Output for each tool and dataset will be generated at output folder.


### Troubleshooting  ###

Two important checks to ensure that metrics are computed correctly:

- The **sample names in the `validated_results_file` should match the file names of your bam files** (excluding the .bam extension). For example, if the `validated_results_file` contains sample names like mySample2312, your bam files should have file names like mySample2312.bam .
- Provide and use chromosomes names with the same format, that is, do not use "chr5" and "5" in you bed and `validated_results_file` files, for example.


## Extra feature: evaluate parameters ##
A parameter evaluator is also included in the framework. It executes each tool parameter over a broad range of values to assess its impact on tool performance. Up to 15 values are evaluated for each numerical param and all the available options for categorical ones. 
The parameter evaluator facilitates understanding how individual parameters influence the overall performance of CNV calling tools.


### Prerequisites ###

An SGE cluster system has to be available.

### Run evaluate parameters
1. Configure all the steps needed to run the benckmark (expained in section how to use)

2. Before running any evaluations, you need to create subfolders  that will contain the YAML files. These files define the settings values for each execution. Execute the setUpFolders script from the CNVbenchmarkeR2 directory.
``` 
Rscript evaluate_parameters/setUpFolders.R [-t tools_yaml] [-d datasets_yaml]
```
3. Modify [jobs.sh](https://github.com/jpuntomarcos/CNVbenchmarkeR2/blob/master/evaluate_parameters/job.sh) to set your SGE parameters and paths.
   
4. Execute the runEvaluate script. This process involves calculating the values for each modified parameter within the selected tools and datasets. For parameters not explicitly modified, the script will apply default values. Execute the setUpFolders script from the CNVbenchmarkeR2 directory.

```
Rscript evaluate_parameters/runEvaluate.R [-t tools_file] [-d datasets_file] [-f keepTempFiles]

```
For space optimization, it is recommended to set the -f parameter to false, which deletes all intermediate files.


5. After the evaluation completes, generate summary CSV files for a comprehensive overview. Again, ensure you are in the CNVbenchmarkeR2 directory and run the summaryEvaluate script:
```
Rscript evaluate_parameters/summaryEvaluate.R [-t tools_file] [-d datasets_file]

```
Each CSV file will contain a summary for every parameter in each dataset, stored in the following path: evaluate_parameters/tool/dataset/parameter/results-dataset_tool_param.csv.


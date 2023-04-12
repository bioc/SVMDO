#' @title SVMDO
#' @return Providing package-level manual page
#' 
#' @description Package Description:
#' It is an easy-to-use GUI using disease information for detecting tumor/normal sample 
#' discriminating gene sets from differentially expressed genes. Our approach is based on an 
#' iterative algorithm filtering genes with disease ontology enrichment analysis and wilk’s lambda 
#' criterion connected to SVM classification model construction. Along with gene set extraction, 
#' SVMDO also provides individual prognostic marker detection among the discriminative genes. The algorithm is designed 
#' for FPKM and RPKM normalized RNA-Seq transcriptome datasets. To provide experience about the GUI usage, a test section 
#' involving dummy example using SummarizedExperiment objects of transcriptome (small form) and clinical datasets is also included.  

#' @section Package Sections:
#' 1. **Analysis:** Acquiring discriminative gene sets and further detecting the gene subset with prognostic characteristics 
#' 2. **Result:** Visualization and Download of discriminative gene sets and survival plot list of prognostic genes

#' @section Steps of Analysis Screen:

#'  1. To search your transcriptome dataset, use the file detection in Choose Your Expression Dataset section. The file will be automatically uploaded into the GUI.
#'  2. To prevent clashing with test datasets, "None" option has to be selected from the radio button section.
#'  3. By clicking on DEG Analysis button you further apply differential expression analysis. Labels of tissue_type column in dataset must contain “Nor” and “Tum” for determining normal/tumour (or tumor) samples. A message window saying Process Completed will appear if there is not any problem.
#'  4. When the differential expression process is completed, a user-defined input size (n) is selected to filter the initial gene list (i.e., n number of upregulated and downregulated genes) by entering a number in Input Size section. It is predetermined as 50 in GUI which can be changed based on the user. If there is problem with the value of input size, you will get a warning about inappropriate input size selection. If the input size remains, algorithm selects all of the differentiallly expressed genes to be used in the next process.
#'  5. To apply disease ontology-based gene filtration, click on DO Analysis button. A message window saying process completed will appear if there is not any problem.
#'  6. To further apply the following feature selection and classification processes, click on the Classification button. A message window saying process completed will appear if there is not any problem.
#'  7. Acquired discriminative gene set can be further used for survival analysis to detect individual prognostic genes. To apply this process, use the file detection in Choose Clinical Data section for searching clinical data about patient survival followed by clicking on Survival Analysis button.

#' @section Steps of Result Screen:

#'  1. To visualize discriminative gene sets inside GUI screen, click on Show Gene Results button. When you click this button, a table of gene set will appear. If there is a problem in the analysis, an error message will appear.
#'  2. To visualize survival plots of individual genes, two steps have to be applied. First of all, click on Prepare Plot Lists button to feed plot information to the visualization system. After that, click on Show Plots button to visualize survival plots.
#'  3. Before downloading files, you can adjust the output directory with Choose Directory button. It can be used for separating files by selecting a destination before clicking download buttons. If it is desired, files can be downloaded to the same folder by selecting an output directory just one time before the download steps.If you do not select any output directory, files will be downloaded to your working directory.
#'  4. To download the resulting discriminative gene set, it is obligatory to define a filename in the Enter Final Gene Set Filename section. After that, you can click on Download Gene List button to complete the process.
#'  5. To download survival plots, you have to click on Download Plot List button. Names of plot files are automatically done by assigning gene names.

#' @section Application of Test Datasets:

#' SVMDO includes test datasets providing dummy examples for gaining experience on the GUI usage.
#' Test datasets consist of simplified forms of TCGA-COAD (COAD) and TCGA-LUSC (LUSC) with 400 genes along with clinical datasets loaded into summarized experiment objects. 
#' When test datasets are used, predetermined expression and clinical datasets are automatically uploaded into the GUI.
#' A test-based analysis is done with predefined input size (n=50). Therefore, users have to continue with DO Analysis after DEG Analysis.

#' @section Workspace Clearance:

#'When the user task is completed, click on the Clear Environment button to remove the global variables created during the algorithm sections. 
#'To prevent error in the next usages of GUI, it is a necessary process. 
#'It can be applied at any moment without the necessity of completing all of the steps of algorithm.


#' @docType package
#' @name SVMDO

NULL
#> NULL
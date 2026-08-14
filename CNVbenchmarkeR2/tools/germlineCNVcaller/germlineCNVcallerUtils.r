# Utils to call GermlineCNVcaller steps

# Creates Interval list if not exists
createIntervalList <- function(picard, bedFile, destPath, fastaDict){
  
  # Create if not exists
  if(!file.exists(interval_list)){
    cmd_interval <- paste("java -jar", picard,
                          "BedToIntervalList",
                          "-I", bedFile,
                          "-O", destPath,
                          "-SD", fastaDict)
    print(cmd_interval);system(cmd_interval);
  }  
}


# Calls GATK PreprocessIntervals
preprocessIntervals <- function(gatk, fastaFile, outputFolder, interval_list, destPath){
  
  # Create if not exists
  if(!file.exists(preprocessInterval_list)){
    cmd_preprocessIntervals <- paste("singularity exec -B /fs03/vh83 -B /fs04/vh83 -B /fs04/scratch2/vh83 -B", paste0(dirname(fastaFile), ",", outputFolder),
                                     gatk, 'gatk PreprocessIntervals',
                                     "-R", fastaFile,
                                     "-L", interval_list,
                                     "-imr OVERLAPPING_ONLY", 
                                     "-O", destPath)
    print(cmd_preprocessIntervals);system(cmd_preprocessIntervals);
  }
}


# Calls GATK AnnotateIntervals
annotateIntervals <- function(gatk, fastaFile, outputFolder, preprocessInterval_list, destPath){
  
  # Create if not exists
  if(!file.exists(annotatedIntervals_list)){
    cmd_annotateIntervals <- paste("singularity exec -B /fs03/vh83 -B /fs04/vh83 -B /fs04/scratch2/vh83 -B", paste0(dirname(fastaFile), ",", outputFolder),
                                   gatk, 'gatk AnnotateIntervals',
                                   "-R", fastaFile,
                                   "-L", preprocessInterval_list,
                                   "-imr OVERLAPPING_ONLY", 
                                   "-O", destPath)
    print(cmd_annotateIntervals);system(cmd_annotateIntervals);
  }
}


# Calls GATK CollectReadCounts for all bam files
collectReadCounts <- function(gatk, bamFiles, outputFolder, preprocessInterval_list, readCountsFolder){
  
  dir.create(readCountsFolder, showWarnings = FALSE)
  print(paste("Starting at", Sys.time(), "collect Read counts"))
  
  
  for (bam in bamFiles){
    print(basename(bam))
    bamsDir <- dirname(bam)
    sample_name <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(bam))
    cmd_read_counts <- paste("singularity exec -B /fs03/vh83 -B /fs04/vh83 -B /fs04/scratch2/vh83 -B", paste0(bamsDir, ",", outputFolder),
                             gatk, 'gatk CollectReadCounts',
                             "-I", paste0(bamsDir,"/",sample_name, ".bam"),
                             "-L", preprocessInterval_list,
                             "--interval-merging-rule OVERLAPPING_ONLY",
                             "-O", paste0(readCountsFolder,"/", sample_name, ".counts.hdf5"))
    print(cmd_read_counts);system(cmd_read_counts);
  }
}


# Calls GATK FilterIntervals
filterIntervals <- function(gatk, outputFolder, preprocessInterval_list, hdf5sFiles, annotatedIntervals_list, filteredIntervals_list){
  if(!file.exists(filteredIntervals_list)){
    cmd_filterIntervals_list <- paste("singularity exec -B /fs03/vh83 -B /fs04/vh83 -B /fs04/scratch2/vh83 -B", outputFolder,
                                      gatk, 'gatk FilterIntervals',
                                      " -L ",preprocessInterval_list,
                                      " -I ", hdf5sFiles,
                                      " --annotated-intervals ", annotatedIntervals_list,
                                      " --interval-merging-rule OVERLAPPING_ONLY",
                                      " -O ", filteredIntervals_list )
    print(cmd_filterIntervals_list);system(cmd_filterIntervals_list);
  }  
}



# Calls GATK DetermineGermlineContigPloidy
determineContigPloidy <- function(gatk, contigPloidyDir, bamsDir, outputFolder, contigFile, all_hdf5_names, filteredIntervals_list){
  if(!dir.exists(contigPloidyDir)){
    dir.create(contigPloidyDir)
    print(paste("Starting at", Sys.time(), "contigPloidy"))
    cmd_det_contig_ploidy <- paste("singularity exec -B /fs03/vh83 -B /fs04/vh83 -B /fs04/scratch2/vh83 -B", paste0(bamsDir,",",outputFolder, ",", dirname(contigFile)),
                                   gatk, 'gatk DetermineGermlineContigPloidy',
                                   "-I", all_hdf5_names,
                                   "-L", filteredIntervals_list,
                                   "--interval-merging-rule OVERLAPPING_ONLY",
                                   "--contig-ploidy-priors", contigFile,
                                   "-O", contigPloidyDir,
                                   "--output-prefix contig")
    print(cmd_det_contig_ploidy);system(cmd_det_contig_ploidy);
  }
}
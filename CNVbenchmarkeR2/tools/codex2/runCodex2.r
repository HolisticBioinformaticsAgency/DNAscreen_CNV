# Runs codex2 over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runCodex2.r [codex2_params_file] [datasets_params_file] [keepTempFiles]
# keepTempFiles: if true, temp files will not be removed (Default: true)
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(CODEX2))
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions

#Functions----
# translates del/dup to deletion/duplication
auxCNname <- function(x) {
  if (x == "del") return("deletion")
  else if (x == "dup") return("duplication")
}

# Functions that need to be adapted
#available at https://github.com/yuchaojiang/CODEX2/blob/master/targeted_sequencing/segment_targeted.R
segment_targeted <- function(yi, yhati, sampname_qc, refi, genei, lmax, mode) {
  chri <- as.matrix(seqnames(refi))[1]
  finalcall <- matrix(ncol = 10)
  lmax <- max(1, lmax - 1)
  if(is.vector(yi)){
    yi=t(yi)
    yhati=t(yhati)
  }
  for (sampno in 1:ncol(yi)) {
    #message("Segmenting sample ", sampno, ": ", sampname_qc[sampno], ".")
    y <- yi[, sampno]
    yhat <- yhati[, sampno]
    num <- length(y)
    y <- c(y, rep(0, lmax))
    yhat <- c(yhat, rep(0, lmax))
    i <- rep(1:num, rep(lmax, num))
    j <- rep(1:lmax, num) + i
    yact <- rep(0, length(i))
    lambda <- rep(0, length(i))
    for (k in 1:num) {
      yact[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(y[k:(k +
                                                                lmax)])[-1]
      lambda[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(yhat[k:(k +
                                                                     lmax)])[-1]
    }
    i <- i[j <= num]
    yact <- yact[j <= num]
    lambda <- lambda[j <= num]
    j <- j[j <= num]
    yact[lambda<20]=20
    lambda[lambda<20]=20
    if (mode == "integer") {
      chat <- round(2 * (yact/lambda))
    } else if (mode == "fraction") {
      chat <- 2 * (yact/lambda)
    }
    lratio <- (1 - chat/2) * lambda + log((chat + 1e-04)/2.0001) * yact
    # chat[chat > 5] <- 5
    if (sum(lratio > 0) > 0) {
      if (sum(lratio > 0) >= 2) {
        finalmat <- (cbind(i, j, yact, lambda, chat, lratio))[lratio >
                                                                0, ]
        finalmat <- finalmat[order(-finalmat[, 6]), ]
        s <- 1
        while (s <= (nrow(finalmat))) {
          rowstart <- finalmat[s, 1]
          rowend <- finalmat[s, 2]
          rowsel <- (finalmat[, 1] <= rowend & finalmat[, 2] >=
                       rowstart)
          rowsel[s] <- FALSE
          finalmat <- finalmat[!rowsel, ]
          if (is.vector(finalmat)) {
            finalmat <- t(as.matrix(finalmat))
          }
          s <- s + 1
        }
      }
      if (sum(lratio > 0) == 1) {
        finalmat <- (cbind(i, j, yact, lambda, chat, lratio))[lratio >
                                                                0, ]
        finalmat <- t(as.matrix(finalmat))
      }
      finalmat <- round(finalmat, digits = 3)
      loglikeij <- cumsum(finalmat[, 6])
      mBIC <- rep(NA, length(loglikeij))
      for (s in 1:nrow(finalmat)) {
        tau <- sort(unique(c(as.vector(finalmat[1:s, 1:2]), 1, num)))
        P <- length(tau) - 2
        mbic <- loglikeij[s]
        mbic <- mbic - 0.5 * sum(log(tau[2:length(tau)] -
                                       tau[1:(length(tau) - 1)]))
        mbic <- mbic + (0.5 - P) * log(num)
        mBIC[s] <- mbic
      }
      mBIC <- round(mBIC, digits = 3)
      if (mBIC[1] > 0) {
        finalmat <- cbind(rep(sampname_qc[sampno], nrow(finalmat)),
                          rep(chri, nrow(finalmat)),rep(genei, nrow(finalmat)), finalmat)
        finalmat <- (cbind(finalmat, mBIC)[1:which.max(mBIC), ])
        finalcall <- rbind(finalcall, finalmat)
      }
    }
  }
  finalcall <- finalcall[-1, ]
  if (is.vector(finalcall)) {
    finalcall <- t(as.matrix(finalcall))
  }
  st <- start(refi)[as.numeric(finalcall[, 4])]
  ed <- end(refi)[as.numeric(finalcall[, 5])]
  cnvtype <- rep(NA, length(st))
  cnvtype[as.numeric(finalcall[, 8]) < 2] <- "del"
  cnvtype[as.numeric(finalcall[, 8]) > 2] <- "dup"
  if (nrow(finalcall) == 1) {
    finalcall <- t(as.matrix(c(finalcall[, 1:3], cnvtype, st, ed, (ed -
                                                                     st + 1)/1000, finalcall[, 4:10])))
  } else {
    finalcall <- cbind(finalcall[, 1:3], cnvtype, st, ed, (ed - st +
                                                             1)/1000, finalcall[, 4:10])
  }
  colnames(finalcall) <- c("sample_name", "chr", "gene" ,"cnv", "st_bp", "ed_bp",
                           "length_kb", "st_exon", "ed_exon", "raw_cov",
                           "norm_cov", "copy_no", "lratio", "mBIC")
  rownames(finalcall) <- rep("", nrow(finalcall))
  finalcall
}


getgc =function (ref, genome = NULL) {
  if(is.null(genome)){genome=BSgenome.Hsapiens.UCSC.hg19}
  gc=rep(NA,length(ref))
  for(chr in unique(seqnames(ref))){
    message("Getting GC content for chr ", chr, sep = "")
    chr.index=which(as.matrix(seqnames(ref))==chr)
    ref.chr=IRanges(start= start(ref)[chr.index] , end = end(ref)[chr.index])
    if (chr == "X" | chr == "x" | chr == "chrX" | chr == "chrx") {
      chrtemp <- 'X'
    } else if (chr == "Y" | chr == "y" | chr == "chrY" | chr ==
               "chry") {
      chrtemp <- 'Y'
    } else {
      chrtemp <- as.numeric(mapSeqlevels(as.character(chr), "NCBI")[1])
    }
    if (length(chrtemp) == 0) message("Chromosome cannot be found in NCBI Homo sapiens database!")
    chrm <- unmasked(genome[[paste('chr',chrtemp,sep='')]])
    seqs <- Views(chrm, ref.chr)
    af <- alphabetFrequency(seqs, baseOnly = TRUE, as.prob = TRUE)
    gc[chr.index] <- round((af[, "G"] + af[, "C"]) * 100, 2)
  }
  gc
}


#Get parameters----
## Read args----
args <- commandArgs(TRUE)
print(args)
if(length(args)>0) {
  codex2ParamsFile <- args[1]
  datasetsParamsFile <- args[2]
  keepTempFiles <- args[3]
} else {
  codex2ParamsFile <- "codex2Params.yaml"
  datasetsParamsFile <- "../../datasets.yaml"
  keepTempFiles <- "true"
}

## Load the parameters file----
params <- yaml.load_file(codex2ParamsFile)
datasets <- yaml.load_file(datasetsParamsFile)
print(paste("Params for this execution:", list(params)))
print(paste("Datasets for this execution:", list(datasets)))

# Dataset iteration ----
# go over datasets and run codex2 for those which are active
for (name in names(datasets)) {
  dataset <- datasets[[name]]
  if (dataset$include){
    print(paste("Starting codex2 for", name, "dataset", sep=" "))

    # extract fields
    bamsDir <- file.path(dataset$bams_dir)
    bedFile <- file.path(dataset$bed_file)
    fastaFile <- file.path(dataset$fasta_file)
    print(getwd())
    # Create output folder
    if (!is.null(params$outputFolder)) {
      outputFolder <- params$outputFolder
    } else{
      outputFolder <- file.path(getwd(), "output", paste0("codex2-", name)) }
    if (is.null(params$execution) || params$execution != "skipPrecalcPhase") {
      unlink(outputFolder, recursive = TRUE);
      dir.create(outputFolder)
    }

    files <- list.files(bamsDir, pattern = '*.bam$', full.names = TRUE)

    # Do pre-calc part of the algorithm
    if (is.null(params$execution) || params$execution != "skipPrecalcPhase") {


      # get bam directories, read in bed file, get sample names
      sampname <- as.matrix(unlist(strsplit(files,"\\.bam")))
      bambedObj <- CODEX2::getbambed(bamdir = files,
                                     bedFile = bedFile,
                                     sampname = sampname,
                                     projectname = "projectname"
      )
      ref <- bambedObj$ref

      ## Get GC content and mappability----
      gc <- getgc(ref)
      mapp <- getmapp(ref)

      #Getting gene names, needed for targeted sequencing
      gene <- read.csv2(bedFile, sep="\t", header = F)$V4
      values(ref) <- cbind(values(ref), DataFrame(gc, mapp, gene))

      ## Get depth of coverage----
      coverageObj <- getcoverage(bambedObj, mapqthres = 20)
      Y <- coverageObj$Y

      ## Quality control----

      qcObj <- qc(Y, sampname, ref,
                  cov_thresh = c(params$cov_thresh_down, params$cov_thresh_up),
                  length_thresh =c(params$length_thresh_down, params$length_thresh_up),
                  mapp_thresh = params$mapp_thresh,
                  gc_thresh = c(params$gc_thresh_down, params$gc_thresh_up))

      Y_qc <- qcObj$Y_qc; sampname_qc <- qcObj$sampname_qc
      ref_qc <- qcObj$ref_qc; qcmat <- qcObj$qcmat; gc_qc <- ref_qc$gc


      ## Estimating library size factor for each sample----
      Y.nonzero <- Y_qc[apply(Y_qc, 1, function(x){!any(x==0)}),]
      pseudo.sample <- apply(Y.nonzero,1,function(x){exp(1/length(x)*sum(log(x)))})
      N <- apply(apply(Y.nonzero, 2, function(x){x/pseudo.sample}), 2, median)
      #plot(N, apply(Y,2,sum), xlab='Estimated library size factor', ylab='Total sum of reads')

      ## Genome-wide normalization using normalize_null----
      # If there are negative control samples, use normalize_codex2_ns()
      # If there are negative control regions, use normalize_codex2_nr()
      normObj.null <- normalize_null(Y_qc = Y_qc,
                                     gc_qc = gc_qc,
                                     K = 1:4, N = N)
      Yhat <- normObj.null$Yhat
      AIC <- normObj.null$AIC; BIC <- normObj.null$BIC
      RSS <- normObj.null$RSS


      ## CBS segmentation per gene: optinmal for targeted seq----
      #source('segment_targeted.R')
      # Available at: https://github.com/yuchaojiang/CODEX2/blob/master/targeted_sequencing/segment_targeted.R
      optK = which.max(BIC)
      finalcall = matrix(ncol = 14, nrow = 0)
      colnames(finalcall)=c('sample_name','chr','gene','cnv',
                            'st_bp','ed_bp','length_kb',
                            'st_exon','ed_exon','raw_cov',
                            'norm_cov','copy_no','lratio',
                            'mBIC')

      for(genei in unique(ref_qc$gene)){
        cat('Segmenting gene',genei,'\n')
        geneindex = which(ref_qc$gene == genei)
        yi = Y_qc[geneindex,, drop = FALSE]
        yhati = Yhat[[optK]][geneindex,, drop=FALSE]
        refi = ref_qc[geneindex]
        finalcalli = segment_targeted(yi, yhati, sampname_qc, refi, genei, lmax=length(geneindex), mode='fraction')
        finalcall = rbind(finalcall,finalcalli)
      }

      cn <- (as.numeric(as.matrix(finalcall[,'copy_no'])))
      cn.filter <- (cn <= params$cn_del_threshold) | (cn >= params$cn_dup_threshold)
      finalcall <- finalcall[cn.filter,]

      ## Set right sample name and cnv naming----
      for (i in 1:nrow(finalcall)){
        parts <- strsplit(finalcall[i, "sample_name"], "/")[[1]]
        finalcall[i, "sample_name"] <- parts[length(parts)]
        finalcall[i, "cnv"] <- auxCNname(finalcall[i, "cnv"])
      }

      # Save results----
      ##CSV file----
      write.table(finalcall, file = file.path(outputFolder, "cnvFounds.csv"), sep='\t', quote=F, row.names=F)

      ##GRanges format----
      message("Saving CNV GenomicRanges")
      saveResultsFileToGR(outputFolder, "cnvFounds.csv",
                          geneColumn = "gene",
                          sampleColumn = "sample_name",
                          chrColumn = "chr",
                          startColumn = "st_bp",
                          endColumn = "ed_bp",
                          cnvTypeColumn = "cnv")

      print(paste("CODEX2 for", name, "dataset finished", sep=" "))
      cat("\n\n\n")
      ##Temporary files----
      #Delete temporary files if specified
      if(keepTempFiles == "false"){
        filesAll <- list.files(outputFolder, full.names = TRUE)
        filesToKeep <- c("failedROIs.csv", "grPositives.rds", "cnvs_summary.tsv", "cnvFounds.csv", "cnvFounds.txt", "all_cnv_calls.txt", "calls_all.txt", "failures_Failures.txt", "cnv_calls.tsv")
        filesToRemove <- list(filesAll[!(filesAll %in% grep(paste(filesToKeep, collapse = "|"), filesAll, value = TRUE))])
        do.call(unlink, filesToRemove)
      }
    }
  }
}

print(paste("Finishing at", endTime <- Sys.time()))
cat("\nElapsed time:")
print(endTime - startTime)

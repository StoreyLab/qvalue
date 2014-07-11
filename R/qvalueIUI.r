qvalueIUI <- function(outname="qvalue.obj", return.object=TRUE, ...) {
  # This is a function that interacts with the user to produce q values.
  #
  # Args:
  #   outname: The name of the q value object created.
  #   return.object: If TRUE, there is a q object returned. Default is FALSE.
  #
  # Returns:
  #   An interactive feature that helps the user produce q values from the
  #   data with the option to produce various plots and returns q value object:
  #   call: Gives the function call.
  #   pi0: An estimate of the proportion of null p-values.
  #   qvalues: A vector of the estimated q-values.
  #   pvalues: A vector of the original p-values.
  #   lfdr: A vector of the local FDR values.
  #   significant: If fdr.level is TRUE, an indicator of whether the q-values
  #                fell below fdr.level (taking all such q-values to be 
  #                significant controls FDR at level fdr.level).
  #   lambda: A vector of lambda value(s) choosen.
  #   pi0.lambda: A vector of the proportion of null values at lambda.
  VarPval <- function(evar, varname) {
  # Finds p value from environment
  #
  # Args:
  #   evar: text for user to enter pvalue name
  #
  # Returns:
  #   p-value variable from user input
    if (varname == 9) {
      qvalueIUI()
      return(invisible(NULL))
    } else{
        if (grepl('\\$',varname)) {
	   pLis <- strsplit(varname, '\\$')
           if (!(pLis[[1]][1] %in% ls(envir=.GlobalEnv))) {
             cat("\n \t Please select an appropriate value, 9 to restart or \n\t press Ctrl+C (Cmd+C) to abort\n")
             varname <- readline(evar)
             pval <- VarPval(evar, varname)
           } else if (!(pLis[[1]][2] %in% ls(get(pLis[[1]][1])))) {
             cat("\n \t Please select an appropriate value, 9 to restart or \n\t press Ctrl+C (Cmd+C) to abort\n")
             varname <- readline(evar)
             pval <- VarPval(evar, varname) 
           } else {  
             pval <- unlist(as.vector(get(pLis[[1]][1])[pLis[[1]][2]]))
           } 
        } else {
          pval <- CheckValues(varname, ls(envir=.GlobalEnv), evar)
	  cat(typeof(pval))
          if (typeof(pval)=="character") {
            pval <- get(pval)
          } 
        }
      }
    while (!is.vector(pval)) {
      ErrMsg("qvalueIUI", "The selected variable is not a numerical vector.")
      VarPval(evar)
    }
  return(pval)
  }
  CheckValues <- function(x, values, msg) {
    # Determines if user input in correct
    #
    # Args:
    #   x: User input.
    #   values: String of accepted input values.
    #   msg: Output string asking for an input.
    #
    # Returns:
    #   An accepted user input value
    if (x %in% values) {
      invisible(x)
      } else {
         if (grepl('\\$', x)) {
           x <- VarPval(msg, x)
         } else if (!(x %in% values)) {
	    cat("\n \t Please select an appropriate value, 9 to restart or \n\t press Ctrl+C (Cmd+C) to abort\n")
            x <- readline(msg)
            x <- CheckValues(x, values, msg)
         }    
     }
  return(x) 
  }
  ErrMsg <- function(err.func = "qvalueIUI", msg){
    # Error functions
    #
    # Args:
    #   err.func: A string that identifies the source of error.
    #   msg: A string of the options for input.
    # 
    # Returns:
    #   An error message along with options to input again.
    cat('\n')
    cat('\t')
    cat('ERROR in the', err.func, 'function:')
    cat('\n')
    cat('\t')
    cat(msg,'\n\n')
  }
  # Initializations
  err.func <- "qvalueIUI"
  enum <- "\n Enter number: "
  efil <- "\n Enter file name: "
  evar <- "\n Enter variable name: "
  fvar <- "\n Enter new object name: "
  cat("\n------\n")
  cat("
          Welcome to the QVALUE interactive user interface,
          part of the QVALUE software package. This tool
          will ask you a series of questions to set up 
          your QVALUE analysis. Before using this tool, 
          you must either have stored p-values in a variable
          or have a file with the p-values in the correct
          format. See the QVALUE manual for directions
          on how to formate your p-value file. Press Ctrl+C (Cmd+C)
          at any time to abort this interactive user interface. 
		  \n")
  cat("\n------\n")
  # Choice of p values stored in session or selecting a file
  msg <- "\nAre the p-values stored as a variable in this R session?\n 1 -- Yes \n 2 -- No, I would like to load a p-value file now\n"
  cat(msg)
  infile <- readline(enum)
  infile <- CheckValues(infile, c("1","2"), enum)
  cat("\n------\n")
  if (infile=="1") {
    msg2 <- "\nWhat is the variable name where the p-values are stored? Press 9 to restart.\n"
    cat(msg2)
    varname <- readline(evar)
    pval <- VarPval(evar, varname)
    cat(paste("\n \t **Read ", length(pval), "p-values \n"))
  } else {
    msg2 <- "\nWhat is the file name where the p-values are stored? Press 9 to restart.\n\n"
    cat(msg2)
    pvalfile <- file.choose()
    if (pvalfile == 9) {
      qvalue()
      return(invisible(NULL))
    } else{
      pval <- try(scan(pvalfile, what=double(0), quiet=TRUE), silent=TRUE)
    }
    while (!is.vector(pval)) {
      cat("The selected file is not in the QVALUE format.\nEach p-value should be numeric and separated by a new line in the file.\nFor example:\n      0.23\n      0.012\n      0.92\n      ...\n Please select another file.\n")
      pvalfile <- file.choose()
      pval <- try(scan(pvalfile, what=double(), quiet=TRUE), silent=TRUE)
    }
    cat(paste("\n \t **Read ", length(pval), "p-values \n"))
  }
  # Calculates q values
  cat("\nCalculating q-values...")
  qobj <- qvalue(p=pval, ...)
  cat(" Done\n")
  cat("\nEstimated proportion of true nulls, pi0:", round(qobj$pi0, 3), "\n")
  cat("\n")
  cat("Cumulative number of significant calls:\n")
  cat("\n")
  cuts <- c(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 1)
  counts <- sapply(cuts, function(x) c("p-value"=sum(qobj$pvalues < x), 
                   "q-value"=sum(qobj$qvalues < x), 
                   "local FDR"=sum(qobj$lfdr < x)))
  colnames(counts) <- paste("<", cuts, sep="")
  print(counts)
  cat("\n")
  readline("\nPress return to continue... ")
  cat("\n------\n")
  # Plotting option
  msg3 <- "\nWould you like to view or save q-value plots?\n 1 -- Yes, view plots in R \n 2 -- Yes, save plots as a PDF file \n 3 -- Yes to both\n 4 -- No\n"
  cat(msg3)
  makeplots <- readline(enum)
  makeplots <- CheckValues(makeplots, c("1", "2", "3", "4"), enum)
  if (makeplots==1 || makeplots==3) {
    cat("\n------\n")
    cat("\nSee the first plot in the other window\n")
    plot(qobj, ...)
    readline("\nPress return to see the next plot... ")
    hist(qobj)
    cat("\nSee the second plot in the other window\n")
    readline("\nPress return to continue... ")
  }
  if (makeplots==2 || makeplots==3) {
    plotfile <- file.choose(TRUE)
    pdf(file=plotfile)
    plot(qobj, ...)
    hist(qobj)
    dev.off()
  }
  cat("\n------\n")
  # Save options
  msg4 <- "\nWould you like to save your results to a tab-delimited text file?\n 1 -- Yes \n 2 -- No \n"
  cat(msg4)
  makeresults <- readline(enum)
  makeresults <- CheckValues(makeresults, c("1", "2"), enum)
  cat("\n------\n")
  if (makeresults==1) {
    resultsfile <- file.choose(TRUE)
    qwrite(qobj, filename=resultsfile, ...)
  }
  if (return.object) {
    assign(outname, qobj, envir=globalenv())
    cat(paste("\n
          Thank you for using the QVALUE interactive user interface.\n
          A qvalue object called: '",outname,"' containing all of the
          above results has been defined in your global R environment. \n\n", sep=""))
    cat("\n------\n\n\n\n")
    return(invisible(NULL))
  } else {
    cat("\n
          Thank you for using the QVALUE interactive user interface.\n\n")
    cat("\n------\n\n\n\n")
    return(qobj)
  }
}

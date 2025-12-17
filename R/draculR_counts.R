#' The geometric mean of vector
#'
#' This function has been extracted from psych v2.4.6.26
#' in order to reduce the amount of libraries required.
#' It calculates the geometric mean of a given vector.
#'
#' @param x vector with values
#' @param na.rm Boolean to remove NAs
#' @return A numeric value 
#' @export
#' @examples
#' geometric.mean(c(2, 3, 4, NA))
#' geometric.mean(c(2, 3, 4, NA), na.rm=TRUE)
geometric.mean <- function(x,na.rm=TRUE){ 
  ## extracted from psych v2.4.6.26
  if (is.null(nrow(x))) {
    exp(mean(log(x),na.rm=na.rm)) 
  } else {
    exp(apply(log(x),2,mean,na.rm=na.rm))
  } 
}

#' Parse miRNA counts to calculate draculR haemolysis score
#'
#' This is the minimun working code to calculate the haemolysis score for each sample.
#'
#' @param counts_df Dataframe containing samples as columns and miRNAs as rows with raw counts
#' @param drop_miRs Vector of miRNAs to remove from the original 20 miRNA used in the calculation
#' @param verbose Boolean to show more/less information
#' @param filterNum Minimun number of counts for each sample
#' @return A dataframe with the classification of each sample as Clear or Caution followin the original calculation by draculR
#' @export
draculR_parse_counts <- function(counts_df, drop_miRs=c(''), verbose=FALSE, filterNum = 1, 
                                 listOfMiRs = c('')) {
  
  library(plyr)
  library(edgeR)
  library(magrittr)
  library(dplyr)
  
  listOfMiRs = unique(c(listOfMiRs, 
                        "hsa-miR-106b-3p", "hsa-miR-140-3p", "hsa-miR-142-5p",
                        "hsa-miR-532-5p",  "hsa-miR-17-5p",  "hsa-miR-19b-3p",
                        "hsa-miR-30c-5p",  "hsa-miR-324-5p", "hsa-miR-192-5p",
                        "hsa-miR-660-5p",  "hsa-miR-186-5p", "hsa-miR-425-5p",
                        "hsa-miR-25-3p",   "hsa-miR-363-3p", "hsa-miR-183-5p",
                        "hsa-miR-451a",    "hsa-miR-182-5p", "hsa-miR-191-5p",
                        "hsa-miR-194-5p",  "hsa-miR-20b-5p"))
  
  # global objects for imported data calculations
  classifier_miRs <- data.frame(
    SYMBOL = listOfMiRs
  )
  
  print(head(classifier_miRs))
  print(dim(classifier_miRs))
  
  
  # negate %in%
  `%notin%` <- Negate(`%in%`)
  
  #-------------------
  ## Rank counts
  #-------------------
  # rank the samples by read counts and by unique miRs
  # this table will be joined downstream with the distribution difference table
  rank <- base::as.data.frame(base::colSums(counts_df)) %>%
    magrittr::set_colnames(., "readCounts") %>% 
    dplyr::arrange(., -(readCounts)) %>% 
    tibble::rownames_to_column("samplename") %>% 
    dplyr::mutate(., rank_readCounts = 1:nrow(.)) %>% 
    dplyr::full_join(.,
                     as.data.frame(t(numcolwise(nonzero)(as.data.frame(counts_df)))) %>%
                       tibble::rownames_to_column() %>%
                       magrittr::set_colnames(., c("samplename", "unique_miRs")) %>%
                       arrange(., desc(unique_miRs)) %>%
                       mutate(., rank_unique = 1:nrow(.)),
                     by = "samplename")
  
  if (verbose) { 
    print("rank")
    print(rank) 
  }
  #-------------------
  
  #-------------------
  # create the (super) minimal metadata table
  #-------------------
  meta <- data.frame(samplename = rank$samplename,
                     readCounts = rank$readCounts) 
  
  if (verbose) { 
    print("meta")
    print(meta) 
  }
  #-------------------
  
  #-------------------
  # Create a DGEList object
  #-------------------
  print("Calculate normalisation factors and apply to the DGEList object")
  DGEList_public <- edgeR::DGEList(counts = counts_df,
                                   samples = rank)
  
  # calculate normalisation factors and apply to the DGEList object
  DGEList_public <- edgeR::calcNormFactors(DGEList_public, method = "TMM")
  #-------------------
  
  if (verbose) { 
    print("DGEList_public")
    print(DGEList_public) 
  }
  #-------------------
  
  #-------------------
  # calculate the CPMs and reduce
  #-------------------
  print("+ Calculate the CPMs")
  rawCPM <- edgeR::cpm(DGEList_public, log = FALSE)
  
  if (verbose) { 
    print("rawCPM")
    print(rawCPM) 
  }
  
  # remove low expressed genes
  print("+ Remove low expressed genes")
  keep.exprs <- rowSums(rawCPM > 40) >= filterNum
  DGEList_public <- DGEList_public[keep.exprs,, keep.lib.sizes = FALSE]
  
  if (verbose) { 
    
    print("keep.exprs")
    print(keep.exprs)
    print("DGEList_public")
    print(DGEList_public) 
  }
  #-------------------
  
  #-------------------
  ## Distributions differences
  #-------------------
  print("+ Calculate the difference between the geometric mean of the distributions:  1 = classifier, 0 = other, 2 = dropped")
  
  ## Calculate the difference between the geometric mean of the distributions
  ### here we calculate the geometric mean of the classifier distribution and the
  ### "other" ensuring those taken from the classifier list are not included in
  ### other.
  
  # create a vector of sample names for use in the lapply
  varc <- dplyr::select(DGEList_public$samples, samplename) %>%
    tibble::remove_rownames() %>% 
    dplyr::pull(., samplename)
  #-------------------
  
  
  #---------------
  # drop miRNAs if necessary
  #---------------
  
  if (length(drop_miRs)>0) {
    print("+ Drop miRNAs if necessary")  
  }
  
  # define the dropped classifiers as input from the groupCheckboxInput
  dropped <- subset(classifier_miRs, SYMBOL %in% drop_miRs)
  
  # define the final set of classifiers
  final_classifiers <- subset(classifier_miRs, SYMBOL %notin% drop_miRs)
  #---------------
  
  #-------------------
  # calculate the geometric mean of the two distributions (1 = classifier, 0 = other, 2 = dropped)
  #-------------------
  distributionDifference <- lapply(varc, function(x) {
    dtmp <- dplyr::select(as.data.frame(edgeR::cpm(DGEList_public$counts, log = TRUE)), x) %>%
      tibble::rownames_to_column("mirna") %>% 
      mutate(., classifier = as.factor(ifelse(mirna %in% final_classifiers$SYMBOL, 1,
                                              ifelse(mirna %in% dropped$SYMBOL, 2,
                                                     ifelse(mirna %notin% classifier_miRs$SYMBOL, 0, NA)))))
    
    cdat_tmp <- with(dtmp,tapply(get(x), classifier, geometric.mean,na.rm=T))
    
    cdat <- data.frame("classifier"=rownames(cdat_tmp),"geometric.mean"=cdat_tmp)
    
    # calculate the difference between the two geometric means (classifier-other)  
    cdat_out <- dplyr::filter(cdat, classifier == 1)$geometric.mean - dplyr::filter(cdat, classifier == 0)$geometric.mean
    return(cdat_out)
  })
  
  ## set names
  names(distributionDifference) <- varc
  
  unlist_distributionDifference <- do.call(cbind.data.frame, distributionDifference) %>% 
    t() %>%
    magrittr::set_colnames("distributionDifference") %>%
    base::as.data.frame() %>% 
    tibble::rownames_to_column("samplename") %>% 
    dplyr::mutate(., haemoResult = ifelse(distributionDifference < 1.9, "Clear",
                                          ifelse(distributionDifference >= 1.9, "Caution", NA)))
  
  unlist_distributionDifference$haemoResult <- as.factor(unlist_distributionDifference$haemoResult)
  
  #-------------------
  ## join and rename
  #-------------------
  distributionDifference.full <- dplyr::full_join(rank, unlist_distributionDifference, by = "samplename")
  colnames(distributionDifference.full)[6] <- "draculR.score"
  
  return(distributionDifference.full)  
}

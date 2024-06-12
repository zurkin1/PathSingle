#' @author Yaru Zhang
#'
#' @name pathway activity score
#' @details This code is motified from R package GSVA, more sufficient for big data set
#' GSVA url: https://bioconductor.org/packages/release/bioc/html/GSVA.html
#'
#' description: Calculates pathway activity score matrix
#'
#'
#' import foreach bigstatsr
#' importFrom bigstatsr as_FBM, big_parallelize
#' importFrom Rcpp sourceCpp
#' importFrom parallel makeCluster
#' importFrom doParallel registerDoParallel
#'
#'
#'
#' @param expr: gene expression profile
#' The expr should be provided as FBM format(see `bigstatsr` package)
#' @param gset.idx.list:
#' @param genes_name: gene names of gene expression profile.
#' Length of gene names are equal to number of columns of gene expression profile
#' @param cells_name: cell names of gene_expression profile.
#' Length of cell names are equal to number of row of gene expression profile
#' @param method:
#' @param min.sz:
#' @param max.sz:
#' @param parallel.sz:
#' @param mx.diff:
#' @param tau:
#' @param ssgsea.norm:
#' @param src_folder:
#' @param verbose:
#'
#'
#'
#' @return
#'
#'
#' @example
#'
library(AUCell)

print('source PAS success')
setGeneric("PAS", function(expr, gset.idx.list, ...) standardGeneric("PAS"))


setMethod("PAS", signature(expr="FBM",
                           gset.idx.list = "list"),
          function(expr,
                   gset.idx.list,
                   genes_name,
                   cells_name,
                   src_folder,
                   method=c("gsva", "aucell", "zscore", "plage"),
                   min.sz=3,
                   max.sz=Inf,
                   parallel.sz=1,
                   mx.diff=TRUE,
                   abs.ranking=FALSE,
                   tau=switch(method, gsva=1, aucell=0.25, NA),
                   ssgsea.norm=TRUE,
                   verbose=TRUE){

            rnaseq = FALSE
            kernel = TRUE
            kcdf = "Gaussian"

            if(method == "plage"){
              min.sz = max(2, min.sz)
            }
            if (nrow(expr) < 2)
              stop("Less than two genes in the input expression data matrix\n")



            ### filter genes according to sd
            sdGenes = big_apply(expr,
                                a.FUN=function(expr, ind)
                                  apply(expr[ind,],1,sd),
                                ind=rows_along(expr),
                                a.combine='c')
            if (any(sdGenes == 0) || any(is.na(sdGenes))) {
              warning(sum(sdGenes == 0 | is.na(sdGenes)),
                      " genes with constant expression values throuhgout the samples.")
              if (method != "ssgsea") {
                warning("Since argument method!=\"ssgsea\", genes with constant expression values are discarded.")
                expr = expr[sdGenes > 0 & !is.na(sdGenes), ]
                genes_name = genes_name[sdGenes > 0 & !is.na(sdGenes)]
                print(dim(expr))
                print(length(genes_name))
              }
            }



            ###################################################
            ######        mapping pathway sets           ######
            ###################################################
            mapped.gset.idx.list = lapply(gset.idx.list,
                                          function(x ,y) na.omit(match(x, y)),
                                          genes_name)
            gc()
            #print("min cells ....................")
            #print(min.sz)

            ###################################################
            ######        filtering pathway sets         ######
            ###################################################
            filterGeneSets = function(gSets,
                                      min.sz=1,
                                      max.sz=Inf) {
              gSetsLen = sapply(gSets,length)
              return (gSets[gSetsLen >= min.sz & gSetsLen <= max.sz])
            }
            mapped.gset.idx.list = filterGeneSets(mapped.gset.idx.list,
                                                  min.sz=max(3, min.sz),
                                                  max.sz=max.sz)
            rm(filterGeneSets,
               min.sz,
               max.sz,
               gset.idx.list)
            gc()

            if (length(unlist(mapped.gset.idx.list, use.names=FALSE)) == 0)
              stop("No identifiers in the gene sets could be matched to the identifiers in the expression data.")


            if(length(mapped.gset.idx.list) == 0){
              stop("The gene set list is empty!  Filter may be too stringent.")
            }


            ###################################################
            ######        calculate PAS using method     ######
            ###################################################
            pas = switch(
              EXPR = method,
              'aucell' = aucell_noWeight(
                expr = expr,
                gset.idx.list = mapped.gset.idx.list,
                cells_name = cells_name,
                src_folder = src_folder,
                tau = tau,
                parallel.sz = parallel.sz,
                normalization = ssgsea.norm,
                verbose = verbose
              ),
              'gsva' = gsva_noWeight(
                expr = expr,
                gset.idx.list = mapped.gset.idx.list,
                cells_name = cells_name,
                src_folder = src_folder,
                parallel.sz = parallel.sz,
                abs.ranking = abs.ranking,
                mx.diff = mx.diff,
                tau = tau,
                kernel = kernel,
                rnaseq = rnaseq,
                verbose = verbose
              ),
              'zscore' = zscore_noWeight(
                expr = expr,
                gset.idx.list = mapped.gset.idx.list,
                cells_name = cells_name,
                parallel.sz = parallel.sz,
                verbose = verbose
              ),
              'plage' = plage_noWeight(
                expr = expr,
                gset.idx.list = mapped.gset.idx.list,
                cells_name = cells_name,
                parallel.sz = parallel.sz,
                verbose = verbose
              ),
              stop("Unknown PAS method: ", method)
            )

            return(pas)

          })


###################################################
######             ssGSEA                    ######
###################################################
## step1: rank
## step2: order
## step3: cumulate

aucell_noWeight = function(expr,
                           gset.idx.list,
                           cells_name,
                           tau,
                           parallel.sz,
                           normalization,
                           src_folder,
                           verbose) {

  if (verbose)
    cat("Computing pathway activity scores using AUCell\n")

  n_cores = min(parallel::detectCores(),parallel.sz)
  print(n_cores)
  # Calculate enrichment scores
  pas_= AUCell_run(expr, gset.idx.list) #, aucMaxRank=nrow(cells_rankings)*0.05)
  # pas = big_parallelize(expr,
  #                       p.FUN = function(expr,
  #                                        ind,
  #                                        gset.idx.list,
  #                                        tau){
  #                         Rcpp::sourceCpp(file.path(src_folder,'ks_ssgsea.cpp'))
  #                         #Rcpp::sourceCpp('ks_test_2.cpp')
  #                         foreach(xc = expr[,ind],.combine = cbind) %do% {
  #                           R = rank(xc)
  #                           O = order(R,decreasing = T)
  #                           ks_gset_ssgsea(geneset_idxs_lis = gset.idx.list,
  #                                          expr = R,
  #                                          sort_idxs = O,
  #                                          n_gsets = length(gset.idx.list),
  #                                          tau = tau)
  #                         }
  #                       },
  #                       p.combine = 'cbind',
  #                       ncores = n_cores,
  #                       gset.idx=gset.idx.list,
  #                       tau=tau)

  # if(nrow(pas) == length(gset.idx.list) && ncol(pas) == expr$ncol){
  #   rownames(pas) = names(gset.idx.list)
  #   colnames(pas) = cells_name
  # }else{
  #   stop("dimensions of pas matrix are not correct, may be the memory cannot allocate\n")
  # }

  if(normalization){
    pas = pas / (range(pas)[2] - range(pas)[1])
  }

  return(pas)
}


###################################################
######             GSVA                      ######
###################################################
## step1: calculate gene density
## step2: order gene density
## step3: custom-rank
## step4: cumulate
gsva_noWeight = function(expr,
                         gset.idx.list,
                         cells_name,
                         src_folder,
                         parallel.sz,
                         abs.ranking,
                         mx.diff,
                         tau,
                         kernel,
                         rnaseq,
                         verbose){

  if (verbose)
    cat("Computing pathway activity scores using GSVA\n")

  n_cores = min(parallel::detectCores(),parallel.sz)
  n_cores2 = min(5, n_cores)

  compute.gene.density = function(expr,
                                  parallel.sz,
                                  src_folder,
                                  rnaseq=FALSE,
                                  kernel=TRUE){

    gene.density = NA
    if (kernel) {
      gene.density = big_parallelize(expr,
                                     p.FUN = function(expr,
                                                      ind,
                                                      src_folder,
                                                      rnaseq){
                                       Rcpp::sourceCpp(file.path(src_folder,'ker.cpp'))
                                       matrix_d(expr[ind,], as.integer(rnaseq))
                                     },
                                     ind = 1:nrow(expr),
                                     p.combine = 'rbind',
                                     src_folder = src_folder,
                                     ncores = n_cores2,
                                     rnaseq = rnaseq)
      print("gensity success")
      rm(expr)
      gc()
      if(nrow(gene.density) == nrow(expr) && ncol(gene.density) == ncol(expr)){
        gene.density = as_FBM(gene.density)
      }else{
        stop("dimensions of gene.density are not correct.\n")
      }

    } else {
      gene.density = t(apply(expr, 1, function(x) {
        f = ecdf(x[sample.idxs])
        f(x)
      }, sample.idxs))
      gene.density = log(gene.density / (1-gene.density))
    }
    #rownames(gene.density) = rownames(expr)
    #colnames(gene.density) = colnames(expr)

    if(is.na(gene.density)){
      stop("Error in calculating gene density!\n")
    }

    return(gene.density)    ## return a list
  }


  # compute_rank_score = function(sort_idx_vec, num_genes){
  #   tmp = rep(0, num_genes)
  #   tmp[sort_idx_vec] = abs(seq(from=num_genes,to=1) - num_genes/2)
  #   return (tmp)
  # }

  gene.density = compute.gene.density(expr,
                                      parallel.sz,
                                      src_folder,
                                      rnaseq, kernel)
  gc()

  if(verbose){
    cat("Gene density success, dimensions of gene.density: \n")
    cat(dim(gene.density))
    cat("\n")
  }

  pas = big_parallelize(gene.density,
                        p.FUN = function(gene.density,
                                         ind,
                                         gset.idx.list,
                                         src_folder,
                                         num_genes,
                                         tau,
                                         mx.diff,
                                         abs.ranking){
                          Rcpp::sourceCpp(file.path(src_folder,'ks_gsva.cpp'))
                          foreach(xc = gene.density[,ind],.combine = cbind) %do% {
                            O = order(xc,decreasing = T)
                            R = rep(0, num_genes)
                            R[O] = abs(seq(from=num_genes,to=1) - num_genes/2)
                            ks_gset_gsva(geneset_idxs_lis = gset.idx.list,
                                         expr = R,
                                         sort_idxs = O,
                                         n_gsets = length(gset.idx.list),
                                         tau = tau,
                                         mx_diff = as.integer(mx.diff),
                                         abs_rank = as.integer(abs.ranking))
                          }
                        },
                        p.combine = 'cbind',
                        ncores = n_cores,
                        gset.idx.list = gset.idx.list,
                        num_genes = nrow(gene.density),
                        src_folder = src_folder,
                        tau = tau,
                        mx.diff = mx.diff,
                        abs.ranking = abs.ranking)


  if(nrow(pas) == length(gset.idx.list) && ncol(pas) == expr$ncol){
    rownames(pas) = names(gset.idx.list)
    colnames(pas) = cells_name
  }else{
    cat("Error in cumulation of GSVA, dimensions of PAS matrix:\n")
    cat(dim(pas))
    stop("\nDimensions of pas matrix are not correct, may be the memory cannot allocate\n")
  }

  return(pas)
}


###################################################
######            zscore                     ######
###################################################
## step1: zscore-convert
## step4: sum
zscore_noWeight = function(expr,
                           gset.idx.list,
                           cells_name,
                           parallel.sz,
                           verbose){

  if (verbose)
    cat("Computing pathway activity scores using zscore\n")

  library(foreach)

  n_cores = min(parallel::detectCores(),parallel.sz)

  Z = big_parallelize(expr,
                      p.FUN = function(expr,
                                       ind,
                                       gset.idx.list){
                        t(scale(t(expr[ind,])))
                      },
                      ind = 1:nrow(expr),
                      p.combine = 'rbind',
                      ncores = n_cores)

  #Z = as_FBM(Z)

  cl=NA
  tryCatch({
    cl = parallel::makeCluster(n_cores,type="FORK")
    print(cl)},
    error = function(e){
             warning("type of FORK is only legal for unix-based platform")
             cl = parallel::makeCluster(n_cores)
           })
  parallel::clusterExport(cl, c("Z"))
  #print(1)
  doParallel::registerDoParallel(cl)
  #print(2)

  pas = foreach(gset = gset.idx.list,
                .combine = rbind) %dopar% {
                  sqrt_len_gset = sqrt(length(gset))
                  if(length(gset) == 1){
                    Z[gset,]
                  }else{
                    colSums(Z[gset,],na.rm = T) / sqrt_len_gset
                  }
                }
  doParallel::stopImplicitCluster()

  if(nrow(pas) == length(gset.idx.list) && ncol(pas) == expr$ncol){
    rownames(pas) = names(gset.idx.list)
    colnames(pas) = cells_name
  }else{
    stop("dimensions of pas matrix are not correct, may be the memory cannot allocate")
  }


  return(pas)

}


###################################################
######               plage                   ######
###################################################
## step1: zscore-convert
## step4: calculate right singular vector
plage_noWeight = function(expr,
                          gset.idx.list,
                          cells_name,
                          parallel.sz,
                          verbose) {

  library(foreach)
  if (verbose)
    cat("Computing pathway activity scores using plage\n")

  n_cores = 4 #min(parallel::detectCores(),parallel.sz)

  Z = big_parallelize(expr,
                      p.FUN = function(expr, ind, gset.idx.list){
                        t(scale(t(expr[ind,])))
                      },
                      ind = 1:nrow(expr),
                      p.combine = 'rbind',
                      ncores = n_cores)

  # Check if Z is defined and accessible
  if (exists("Z")) {
    print("Z is defined in the global environment")
    print(Z)  # Print Z to ensure it's correctly defined
    str(Z)    # Inspect the structure of Z
  } else {
    stop("Z is not defined in the global environment")
  }
  #Z = as_FBM(Z)

  cl=NA
  cl <- tryCatch({
  parallel::makeCluster(n_cores, type = "PSOCK")
  }, error = function(e) {
    warning("Failed to create cluster with PSOCK: ", e$message)
    NULL
  })

  # Ensure the cluster was created successfully
  if (!is.null(cl)) {
    print(0)
    print(ls(.GlobalEnv))  # Print all objects in the global environment to check the presence of Z
    parallel::clusterExport(cl = cl, varlist = c("Z"), envir = .GlobalEnv)
  } else {
    warning("Cluster was not created successfully.")
  }

  parallel::clusterExport(cl = cl, c("Z"), envir = .GlobalEnv)
  print(1)
  doParallel::registerDoParallel(cl)
  print(2)

  pas = foreach(gset = gset.idx.list,
                .combine = rbind) %dopar% {
                  sqrt_len_gset = sqrt(length(gset))
                  if(length(gset) == 1){
                    Z[gset,]
                  }else{
                    svd(Z[gset,])$v[,1]
                  }
                }
  doParallel::stopImplicitCluster()

  if(nrow(pas) == length(gset.idx.list) && ncol(pas) == expr$ncol){
    rownames(pas) = names(gset.idx.list)
    colnames(pas) = cells_name
  }else{
    stop("dimensions of pas matrix are not correct, may be the memory cannot allocate")
  }


  return(pas)
}
#'
#' parallelize version of FindAllMarkers
#' import parallel
#'
#'

FindAllMarkers_me <- function(
  object,
  assay = NULL,
  features = NULL,
  logfc.threshold = 0.25,
  test.use = 'wilcox',
  slot = 'data',
  min.pct = 0.1,
  min.diff.pct = -Inf,
  node = NULL,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  return.thresh = 1e-1,
  para_size,
  ...
) {
  MapVals <- function(vec, from, to) {
    vec2 <- setNames(object = to, nm = from)[as.character(x = vec)]
    vec2[is.na(x = vec2)] <- vec[is.na(x = vec2)]
    return(unname(obj = vec2))
  }
  if ((test.use == "roc") && (return.thresh == 1e-2)) {
    return.thresh <- 0.7
  }
  if (is.null(x = node)) {
    idents.all <- sort(x = unique(x = Idents(object = object)))
  } else {
    tree <- Tool(object = object, slot = 'BuildClusterTree')
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' before finding markers on nodes")
    }
    descendants <- DFT(tree = tree, node = node, include.children = TRUE)
    all.children <- sort(x = tree$edge[, 2][!tree$edge[, 2] %in% tree$edge[, 1]])
    descendants <- MapVals(
      vec = descendants,
      from = all.children,
      to = tree$tip.label
    )
    drop.children <- setdiff(x = tree$tip.label, y = descendants)
    keep.children <- setdiff(x = tree$tip.label, y = drop.children)
    orig.nodes <- c(
      node,
      as.numeric(x = setdiff(x = descendants, y = keep.children))
    )
    tree <- drop.tip(phy = tree, tip = drop.children)
    new.nodes <- unique(x = tree$edge[, 1, drop = TRUE])
    idents.all <- (tree$Nnode + 2):max(tree$edge)
  }

  cal_markers = function(ident){
    tryCatch(
    FindMarkers(
      object = object,
      assay = assay,
      ident.1 = ident,
      ident.2 = NULL,
      features = features,
      logfc.threshold = logfc.threshold,
      test.use = test.use,
      slot = slot,
      min.pct = min.pct,
      min.diff.pct = min.diff.pct,
      verbose = F,
      only.pos = only.pos,
      max.cells.per.ident = max.cells.per.ident,
      random.seed = random.seed,
      latent.vars = latent.vars,
      min.cells.feature = min.cells.feature,
      min.cells.group = min.cells.group,
      pseudocount.use = pseudocount.use,
      ...
    ),
    error = function(e) {
      #print("error in cal markers")
      NULL
    }
    )
  }

  #print("aaaaaaaaaaaaaaaaaaaaaaaaaa")

  mclapp = get('mclapply', envir=getNamespace('parallel'))
  options(mc.cores = para_size)

   genes.de = mclapp(idents.all, cal_markers)


  if(length(genes.de) == 0){
    print("length of genes.de is 0")
    return(data.frame())
  }

  print("all diff genes")
  print(length(genes.de))
  gde.all <- data.frame()
  for (i in 1:length(x = idents.all)) {
    if (is.null(x = unlist(x = genes.de[i]))) {
      next
    }
    print(i)
    gde <- genes.de[[i]]
    #print(dim(gde))
    #print(gde)
    #print(is.null(gde))
    #print(!is.null(gde) && nrow(gde) > 0)
    if (!is.null(gde) && nrow(gde) > 0) {
      if (test.use == "roc") {
        gde <- subset(
          x = gde,
          subset = (myAUC > return.thresh | myAUC < (1 - return.thresh))
        )
      } else if (is.null(x = node) || test.use %in% c('bimod', 't')) {
        gde <- gde[order(gde$p_val_adj,-gde[, 2]), ]
        gde <- subset(x = gde, subset = p_val_adj < return.thresh)
      }
      if (nrow(x = gde) > 0) {
        gde$cluster <- idents.all[i]
        gde$gene <- rownames(x = gde)
      }
      if (nrow(x = gde) > 0) {
        gde.all <- rbind(gde.all, gde)
      }
    }
  }
  if ((only.pos) && nrow(x = gde.all) > 0) {
    return(subset(x = gde.all, subset = gde.all[, 2] > 0))
  }
  rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
  if (nrow(x = gde.all) == 0) {
    warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
  }

  return(gde.all)

}

#' import ggplot2
#' import cowplot
#' import scles
#' import pheatmap



DoHeatmap_me <- function(
  object,
  features = NULL,
  cells = NULL,
  group.by = 'ident',
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = 'scale.data',
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
) {
  cells <- colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- rev(x = unique(x = features))
  disp.max <- ifelse(
    test = slot == 'scale.data',
    yes = 2.5,
    no = 6
  )
  # make sure features are present
  possible.features <- rownames(x = GetAssayData(object = object, slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if(length(x = features) == 0) {
      stop("No requested features found in the ", slot, " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", slot,
            " slot for the ", assay, " assay: ", paste(bad.features, collapse = ", "))
  }
  #print(cells[1:3])
  #print(features[1:3])
  #print(dim(as.matrix(GetAssayData(
  #  object = object,
  #  slot = slot)[features, cells, drop = FALSE])))

  data <- as.data.frame(x = as.matrix(x = t(x = as.matrix(GetAssayData(
    object = object,
    slot = slot)[features, cells, drop = FALSE]))))

  object <- suppressMessages(expr = StashIdent(object = object, save.name = 'ident'))
  group.by <- 'ident'
  groups.use <- object[[group.by]][cells, , drop = FALSE]
  # group.use <- switch(
  #   EXPR = group.by,
  #   'ident' = Idents(object = object),
  #   object[[group.by, drop = TRUE]]
  # )
  # group.use <- factor(x = group.use[cells])
  plots <- vector(mode = 'list', length = ncol(x = groups.use))
  for (i in 1:ncol(x = groups.use)) {
    #print(dim(data))
    #print(i)
    data.group <- data
    group.use <- groups.use[, i, drop = TRUE]
    if (!is.factor(x = group.use)) {
      group.use <- factor(x = group.use)
    }
    names(x = group.use) <- cells
    if (draw.lines) {
      # create fake cells to serve as the white lines, fill with NAs
      lines.width <- ceiling(x = nrow(x = data.group) * 0.0025)
      placeholder.cells <- sapply(
        X = 1:(length(x = levels(x = group.use)) * lines.width),
        FUN = function(x) {
          return(RandomName(length = 20))
        }
      )
      placeholder.groups <- rep(x = levels(x = group.use), times = lines.width)
      group.levels <- levels(x = group.use)
      names(x = placeholder.groups) <- placeholder.cells
      group.use <- as.vector(x = group.use)
      names(x = group.use) <- cells
      group.use <- factor(x = c(group.use, placeholder.groups), levels = group.levels)
      na.data.group <- matrix(
        data = NA,
        nrow = length(x = placeholder.cells),
        ncol = ncol(x = data.group),
        dimnames = list(placeholder.cells, colnames(x = data.group))
      )
      data.group <- rbind(data.group, na.data.group)
    }
    lgroup <- length(levels(group.use))
    #print(data.group[1:2,])
    #print(class(data.group))
    #print(length(group.use))
    #print(names(x = sort(x = group.use))[1:3])
    plot <- SingleRasterMap_me(
      data = data.group,
      raster = raster,
      disp.min = disp.min,
      disp.max = disp.max,
      feature.order = features,
      cell.order = names(x = sort(x = group.use)),
      group.by = group.use
    )
    if (group.bar) {
      # TODO: Change group.bar to annotation.bar
      default.colors <- c(hue_pal()(length(x = levels(x = group.use))))
      cols <- default.colors
      if (any(is.na(x = cols))) {
        cols[is.na(x = cols)] <- default.colors[is.na(x = cols)]
        cols <- Col2Hex(cols)
        col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(
          x = cols,
          start = 1,
          stop = 7
        )))))
        through <- length(x = default.colors)
        while (length(x = col.dups) > 0) {
          pal.max <- length(x = col.dups) + through
          cols.extra <- hue_pal()(pal.max)[(through + 1):pal.max]
          cols[col.dups] <- cols.extra
          col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(
            x = cols,
            start = 1,
            stop = 7
          )))))
        }
      }
      #print(group.use[1:3])
      group.use2 <- sort(x = group.use)
      #print(group.use2[1:3])
      if (draw.lines) {
        na.group <- RandomName(length = 20)
        levels(x = group.use2) <- c(levels(x = group.use2), na.group)
        group.use2[placeholder.cells] <- na.group
        cols <- c(cols, "#FFFFFF")
      }
      pbuild <- ggplot_build(plot = plot)
      names(x = cols) <- levels(x = group.use2)
      # scale the height of the bar
      y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
      y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
      y.max <- y.pos + group.bar.height * y.range
      plot <- plot +
        annotation_raster(
          raster = t(x = cols[group.use2]),
          xmin = -Inf,
          xmax = Inf,
          ymin = y.pos,
          ymax = y.max
        ) +
        coord_cartesian(ylim = c(0, y.max), clip = 'off') +
        scale_color_manual(values = cols)
      if (label) {
        x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
        x.divs <- pbuild$layout$panel_params[[1]]$x.major
        x <- data.frame(group = sort(x = group.use), x = x.divs)
        label.x.pos <- tapply(X = x$x, INDEX = x$group, FUN = median) * x.max
        label.x.pos <- data.frame(group = names(x = label.x.pos), label.x.pos)
        plot <- plot + geom_text(
          stat = "identity",
          data = label.x.pos,
          aes_string(label = 'group', x = 'label.x.pos'),
          y = y.max + y.max * 0.03 * 0.2,
          angle = angle,
          hjust = hjust,
          size = size
        )
        plot <- suppressMessages(plot + coord_cartesian(
          ylim = c(0, y.max + y.max * 0.003 * max(nchar(x = levels(x = group.use))) * size),
          clip = 'off')
        )
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  #print(length(plots))
  plots <- plot_grid(plotlist = plots)
  return(plots)
}

SingleRasterMap_me <- function(
  data,
  raster = TRUE,
  cell.order = NULL,
  feature.order = NULL,
  colors = colorRampPalette(c("forestgreen", "gray90", "orange"))(300),
  disp.min = -2.5,
  disp.max = 2.5,
  limits = NULL,
  group.by = NULL
) {
  #print(dim(data))
  #print(disp.max)
  #print(disp.min)
  data <- MinMax(data = data, min = disp.min, max = disp.max)
  #print(dim(data))
  data <- Melt(x = t(x = data))
  #print(dim(data))
  colnames(x = data) <- c('Feature', 'Cell', 'Expression')
  if (!is.null(x = feature.order)) {
    data$Feature <- factor(x = data$Feature, levels = unique(x = feature.order))
  }
  if (!is.null(x = cell.order)) {
    data$Cell <- factor(x = data$Cell, levels = unique(x = cell.order))
  }
  if (!is.null(x = group.by)) {
    data$cell_type <- group.by[data$Cell]
  }
  limits <- c(min(data$Expression), max(data$Expression))
  if (length(x = limits) != 2 || !is.numeric(x = limits)) {
    stop("limits' must be a two-length numeric vector")
  }
  my_geom <- ifelse(test = raster, yes = geom_raster, no = geom_tile)
  plot <- ggplot(data = data) +
    my_geom(mapping = aes_string(x = 'Cell', y = 'Feature', fill = 'Expression')) +
    scale_fill_gradientn('value',
                         limits = limits,
                         colors = colors,
                         na.value = "white")  +
    labs(x = NULL, y = NULL, fill = 'Expression') +
    WhiteBackground() + NoAxes(keep.text = TRUE)
  if (!is.null(x = group.by)) {
    plot <- plot + geom_point(
      mapping = aes_string(x = 'Cell', y = 'Feature', color = 'cell_type'),
      size=10,
      alpha = 0
    ) +
      guides(color = guide_legend(override.aes = list(alpha = 1)))+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 15),
            legend.text = element_text(size = 16, face='bold'),
            legend.title = element_text(size = 20,face='bold'))
  }
  return(plot)
}

RandomName <- function(length = 5L, ...) {
  CheckDots(..., fxns = 'sample')
  return(paste(sample(x = letters, size = length, ...), collapse = ''))
}

CheckDots <- function(..., fxns = NULL) {
  args.names <- names(x = list(...))
  if (length(x = list(...)) == 0) {
    return(invisible(x = NULL))
  }
  if (is.null(x = args.names)) {
    stop("No named arguments passed")
  }
  if (length(x = fxns) == 1) {
    fxns <- list(fxns)
  }
  for (f in fxns) {
    if (!(is.character(x = f) || is.function(x = f))) {
      stop("CheckDots only works on characters or functions, not ", class(x = f))
    }
  }
  fxn.args <- suppressWarnings(expr = sapply(
    X = fxns,
    FUN = function(x) {
      x <- tryCatch(
        expr = if (isS3stdGeneric(f = x)) {
          as.character(x = methods(generic.function = x))
        } else {
          x
        },
        error = function(...) {
          return(x)
        }
      )
      x <- if (is.character(x = x)) {
        sapply(X = x, FUN = argsAnywhere, simplify = FALSE, USE.NAMES = TRUE)
      } else if (length(x = x) <= 1) {
        list(x)
      }
      return(sapply(
        X = x,
        FUN = function(f) {
          return(names(x = formals(fun = f)))
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  ))
  fxn.args <- unlist(x = fxn.args, recursive = FALSE)
  fxn.null <- vapply(X = fxn.args, FUN = is.null, FUN.VALUE = logical(length = 1L))
  if (all(fxn.null) && !is.null(x = fxns)) {
    stop("None of the functions passed could be found")
  } else if (any(fxn.null)) {
    warning(
      "The following functions passed could not be found: ",
      paste(names(x = which(x = fxn.null)), collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
    fxn.args <- Filter(f = Negate(f = is.null), x = fxn.args)
  }
  dfxns <- vector(mode = 'logical', length = length(x = fxn.args))
  names(x = dfxns) <- names(x = fxn.args)
  for (i in 1:length(x = fxn.args)) {
    dfxns[i] <- any(grepl(pattern = '...', x = fxn.args[[i]], fixed = TRUE))
  }
  if (any(dfxns)) {
    dfxns <- names(x = which(x = dfxns))
    if (any(nchar(x = dfxns) > 0)) {
      fx <- vapply(
        X = Filter(f = nchar, x = dfxns),
        FUN = function(x) {
          if (isS3method(method = x)) {
            x <- unlist(x = strsplit(x = x, split = '\\.'))
            x <- x[length(x = x) - 1L]
          }
          return(x)
        },
        FUN.VALUE = character(length = 1L)
      )
      message(
        "The following functions and any applicable methods accept the dots: ",
        paste(unique(x = fx), collapse = ', ')
      )
      if (any(nchar(x = dfxns) < 1)) {
        message(
          "In addition, there is/are ",
          length(x = Filter(f = Negate(f = nchar), x = dfxns)),
          " other function(s) that accept(s) the dots"
        )
      }
    } else {
      message("There is/are ", length(x = dfxns), 'function(s) that accept(s) the dots')
    }
  } else {
    unused <- Filter(
      f = function(x) {
        return(!x %in% unlist(x = fxn.args))
      },
      x = args.names
    )
    if (length(x = unused) > 0) {
      msg <- paste0(
        "The following arguments are not used: ",
        paste(unused, collapse = ', ')
      )
      switch(
        EXPR = getOption(x = "Seurat.checkdots"),
        "warn" = warning(msg, call. = FALSE, immediate. = TRUE),
        "stop" = stop(msg),
        "silent" = NULL,
        stop("Invalid Seurat.checkdots option. Please choose one of warn, stop, silent")
      )
      unused.hints <- sapply(X = unused, FUN = OldParamHints)
      names(x = unused.hints) <- unused
      unused.hints <- na.omit(object = unused.hints)
      if (length(x = unused.hints) > 0) {
        message(
          "Suggested parameter: ",
          paste(unused.hints, "instead of", names(x = unused.hints), collapse = '; '),
          "\n"
        )
      }
    }
  }
}

Melt <- function(x) {
  if (!is.data.frame(x = x)) {
    x <- as.data.frame(x = x)
  }
  return(data.frame(
    rows = rep.int(x = rownames(x = x), times = ncol(x = x)),
    cols = unlist(x = lapply(X = colnames(x = x), FUN = rep.int, times = nrow(x = x))),
    vals = unlist(x = x, use.names = FALSE)
  ))
}

DimPlot_me <- function(
  data,
  dims = c(1, 2),
  cells = NULL,
  cols = NULL,
  pt.size = NULL,
  reduction = NULL,
  group.by = NULL,
  split.by = NULL,
  shape.by = NULL,
  order = NULL,
  label = FALSE,
  label.size = 4,
  repel = FALSE,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50',
  combine = TRUE,
  ncol = NULL,
  ...
) {

  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  reduction <- reduction %||% DefaultDimReduc(object = object)
  cells <- cells %||% colnames(x = object)
  data <- Embeddings(object = object[[reduction]])[cells, dims]
  data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)
  object[['ident']] <- Idents(object = object)
  group.by <- group.by %||% 'ident'
  data[, group.by] <- object[[group.by]][cells, , drop = FALSE]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    data[, split.by] <- object[[split.by, drop = TRUE]]
  }
  plots <- lapply(
    X = group.by,
    FUN = function(x) {
      plot <- SingleDimPlot(
        data = data[, c(dims, x, split.by, shape.by)],
        dims = dims,
        col.by = x,
        cols = cols,
        pt.size = pt.size,
        shape.by = shape.by,
        order = order,
        label = FALSE,
        cells.highlight = cells.highlight,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        na.value = na.value
      )
      if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = x,
          repel = repel,
          size = label.size,
          split.by = split.by
        )
      }
      if (!is.null(x = split.by)) {
        plot <- plot + FacetTheme() +
          facet_wrap(
            facets = vars(!!sym(x = split.by)),
            ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
              length(x = unique(x = data[, split.by]))
            } else {
              ncol
            }
          )
      }
      return(plot)
    }
  )
  if (combine) {
    plots <- CombinePlots(
      plots = plots,
      ncol = if (!is.null(x = split.by) && length(x = group.by) > 1) {
        1
      } else {
        ncol
      },
      ...
    )
  }
  return(plots)
}

SingleDimPlot_me <- function(
  data,
  dims,
  col.by = NULL,
  cols = NULL,
  pt.size = NULL,
  shape.by = NULL,
  order = NULL,
  label = FALSE,
  repel = FALSE,
  label.size = 4,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50'
) {
  pt.size <- pt.size %||% AutoPointSize(data = data)
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
    stop("Cannot find dimensions to plot in data")
  } else if (is.numeric(x = dims)) {
    dims <- colnames(x = data)[dims]
  }
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = sizes.highlight %||% pt.size,
      cols.highlight = cols.highlight,
      col.base = cols[1] %||% '#C3C3C3',
      pt.size = pt.size
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    if (typeof(x = order) == "logical") {
      if (order) {
        data <- data[order(data[, col.by]), ]
      }
    } else {
      order <- rev(x = c(
        order,
        setdiff(x = unique(x = data[, col.by]), y = order)
      ))
      data[, col.by] <- factor(x = data[, col.by], levels = order)
      new.order <- order(x = data[, col.by])
      data <- data[new.order, ]
      if (length(x = pt.size) == length(x = new.order)) {
        pt.size <- pt.size[new.order]
      }
    }
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  } else {
    # col.index <- grep(pattern = col.by, x = colnames(x = data), fixed = TRUE)
    col.index <- match(x = col.by, table = colnames(x = data))
    if (grepl(pattern = '^\\d', x = col.by)) {
      # Do something for numbers
      col.by <- paste0('x', col.by)
    } else if (grepl(pattern = '-', x = col.by)) {
      # Do something for dashes
      col.by <- gsub(pattern = '-', replacement = '.', x = col.by)
    }
    colnames(x = data)[col.index] <- col.by
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  plot <- ggplot(data = data) +
    geom_point(
      mapping = aes_string(
        x = dims[1],
        y = dims[2],
        color = paste0("`", col.by, "`"),
        shape = shape.by
      ),
      size = pt.size
    ) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(color = NULL)
  if (label && !is.null(x = col.by)) {
    plot <- LabelClusters(
      plot = plot,
      id = col.by,
      repel = repel,
      size = label.size
    )
  }
  if (!is.null(x = cols)) {
    plot <- plot + if (length(x = cols) == 1 && (is.numeric(x = cols) || cols %in% rownames(x = brewer.pal.info))) {
      scale_color_brewer(palette = cols, na.value = na.value)
    } else {
      scale_color_manual(values = cols, na.value = na.value)
    }
  }
  plot <- plot + theme_cowplot()
  return(plot)
}

LabelClusters <- function(
  plot,
  id,
  clusters = NULL,
  labels = NULL,
  split.by = NULL,
  repel = TRUE,
  ...
) {
  xynames <- unlist(x = GetXYAesthetics(plot = plot), use.names = TRUE)
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[, id])))
  groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[, id])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", paste(groups[!groups %in% possible.clusters], collapse = ","))
  }
  labels.loc <- lapply(
    X = groups,
    FUN = function(group) {
      data.use <- data[data[, id] == group, , drop = FALSE]
      data.medians <- if (!is.null(x = split.by)) {
        do.call(
          what = 'rbind',
          args = lapply(
            X = unique(x = data.use[, split.by]),
            FUN = function(split) {
              medians <- apply(
                X = data.use[data.use[, split.by] == split, xynames, drop = FALSE],
                MARGIN = 2,
                FUN = median,
                na.rm = TRUE
              )
              medians <- as.data.frame(x = t(x = medians))
              medians[, split.by] <- split
              return(medians)
            }
          )
        )
      } else {
        as.data.frame(x = t(x = apply(
          X = data.use[, xynames, drop = FALSE],
          MARGIN = 2,
          FUN = median,
          na.rm = TRUE
        )))
      }
      data.medians[, id] <- group
      return(data.medians)
    }
  )
  labels.loc <- do.call(what = 'rbind', args = labels.loc)
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop("Length of labels (", length(x = labels),  ") must be equal to the number of clusters being labeled (", length(x = labels.loc), ").")
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  geom.use <- ifelse(test = repel, yes = geom_text_repel, no = geom_text)
  plot <- plot + geom.use(
    data = labels.loc,
    mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
    ...
  )
  return(plot)
}


bubble_plot = function(obSeurat,
                       all_markers,
                       file_path) {
  #n_top = 5
 # tops = all_markers %>% group_by(cluster) %>% top_n(n = n_top, wt = avg_logFC)

  #while( nrow(tops) > 30 ){
  #  n_top = n_top - 1
  #  tops = all_markers %>% group_by(cluster) %>% top_n(n = n_top, wt = avg_logFC)
  #}
  shown_markers = unique(all_markers$gene)

  data.features <- FetchData(object = obSeurat, vars = shown_markers)
  data.features$id = Idents(object = obSeurat)
  data.features$id <- factor(x = data.features$id)
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  max_p = max(-log10(all_markers$p_val_adj[which(all_markers$p_val_adj > 0)]))
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = x))
        }
      )
      pct.exp <- sapply(colnames(data.use), function(p){
        pval = -log10(all_markers[intersect(which(all_markers$gene == p), which(all_markers$cluster == ident)), 'p_val_adj'])
        pval[which(pval == Inf)] = max_p * 1.1
        pval
      })
      return(list(avg.exp = avg.exp, pct.exp = as.numeric(pct.exp)))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = rev(x = shown_markers)
  )

  plot <- ggplot(data = data.plot, mapping = aes_string(x = 'id', y = 'features.plot')) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp')) +
    scale_radius(range = c(0, 8)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = '-log10 pval')) +
    labs(x = 'cell types',y = 'pathways') +
    theme_cowplot()+
    scale_color_gradient(low = "lightgrey", high = 'red')+
    theme(axis.text.x = element_text(angle = 90, size=16, face = "bold"),
          axis.text.y = element_text(angle = 0, size=15),
          axis.title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 16,face='bold'),
          legend.title = element_text(size = 20,face='bold'))+
    guides(color = guide_colorbar(title = 'Average PAS'))
  return(plot)
}

feature_vio = function(obSeurat,
                       feature,
                       path_dir,
                       feature_name,
                       x_angle,
                       pic_type){

  #file_path = paste0(feature_name,'.vil.png')
  x_angle=0
  p = VlnPlot(obSeurat, features = feature, pt.size=0)+
    theme(axis.text.x = element_blank(),
          #axis.text.x = element_text(angle = x_angle, size=15, hjust=0.5),
          axis.text.y = element_text(angle = 0, size=15),
          axis.title.x =element_text(size=18, face = "bold",margin = margin(0, 0, 21, 0, "pt"),vjust=-3.5),
          axis.title.y =element_text(size=16, face = "bold"))+
    labs(x="cell types",y="pathway activity score",title=NULL)

  if(pic_type == 'png'){
    png(file.path(path_dir,paste0(feature_name,'.vil.png')))
    print(p)
    dev.off()
  }else if(pic_type == 'pdf'){
    pdf(paste0(file.path(path_dir,feature_name,'.vil.pdf')))
    print(p)
    dev.off()
  }else{
    stop("--pic_type must be png or pdf")
  }

}

feature_highlight = function(obSeurat,
                             feature,
                             path_dir,
                             feature_name,
                             pic_type){

  p = FeaturePlot(obSeurat, features = feature,cols = c("lightgrey", "red"))+
    theme(axis.text.x = element_text( size=15),
          axis.text.y = element_text( size=15),
          axis.title.x =element_text(size=16, face = "bold"),
          axis.title.y =element_text(size=16, face = "bold"))+
    labs(title=NULL)

  if(pic_type == 'png'){
    png(file.path(path_dir,paste0(feature_name,'.fet.png')))
    print(p)
    dev.off()
  }else if(pic_type == 'pdf'){
    pdf(file.path(path_dir,paste0(feature_name,'.fet.pdf')))
    print(p)
    dev.off()
  }else{
    stop("--pic_type must be png or pdf")
  }

}

genes_heatmap = function(expr,
                         genes,
                         cell_type,
                         path_dir,
                         feature_name,
                         pic_type){
  #max_genes = 20
  expr_plot = expr[intersect(genes,rownames(expr)), ]
  #expr_plot = na.omit(expr_plot)
  #print(dim(expr_plot))
  #print(dim(cell_type))
  #print(head(cell_type))
  #print(expr_plot[1:3,1:3])
  #if(dim(expr_plot)[1] > max_genes){
  #  gene_sd = apply(expr_plot,1, sd)
  #  expr_plot = expr_plot[names(sort(gene_sd,decreasing = T)[0:max_genes]),]
  #}
  #print(dim(expr_plot))
  #write.csv(expr_plot, mat_path, quote=F)
  p = pheatmap(expr_plot,
           color = colorRampPalette(c("cyan2", "gray90", "purple"))(100),
           scale = 'row',
           show_rownames = T,
           show_colnames = F,
           annotation_col = cell_type,
           clustering_method = "ward.D2")

  if(pic_type == 'png'){
    png(file.path(path_dir,paste0(feature_name,'.het.png')))
    print(p)
    dev.off()
  }else if(pic_type == 'pdf'){
    pdf(file.path(path_dir,paste0(feature_name,'.het.pdf')))
    print(p)
    dev.off()
  }else{
    stop("--pic_type must be png or pdf")
  }

}


ReadData <- function(samples, ST_paths) {
  cl <- makeCluster(11)
  registerDoParallel(cl)
  seurat_st <-
    foreach(
      i = 1:length(samples),
      .packages = c("Seurat", "tidyverse", "magrittr")
    ) %dopar% {
      sample = samples[i]
      tmp <- Load10X_Spatial(ST_paths[i],
                             slice = sample)
      tmp <-
        RenameCells(tmp, new.names =  str_c(sample, 
                                            rownames(tmp@meta.data), sep ="@")) # globalKey
      tmp@meta.data$sample = sample
      tmp <-
      PercentageFeatureSet(tmp, pattern = "^MT-", col.name = "percent.mt")
      tmp <- tmp[!grepl("^MT-", rownames(tmp)), ]
      tmp <- SCTransform(tmp,
                         assay = "Spatial",
                         verbose = FALSE)

      # RunPCA(tmp, assay = "SCT", verbose = FALSE) %>%
      #   RunUMAP(reduction = "pca", dims = 1:30) %>%
      #   FindNeighbors(reduction = "pca",dims = 1:50)
    }
  stopCluster(cl)
  seurat_st
}


integratedData <- function(seurat_st,nfeatures = 2000,gene_list=NULL) {
  # library(Seurat)
  features <- SelectIntegrationFeatures(seurat_st,
                                        nfeatures = nfeatures)
  
  features <- c(features,gene_list[!gene_list %in% features],"DRD1")
  
  seurat_st <- PrepSCTIntegration(object.list = seurat_st,
                                  anchor.features = features)
  
  st.anchors <- FindIntegrationAnchors(object.list = seurat_st, 
                                       # assay = "SCT",
                                       normalization.method = "SCT",
                                       anchor.features = features)
  
  # This is much faster...
  st.combined.sct <-
    IntegrateData(anchorset = st.anchors, normalization.method = "SCT")
}


sign2star <- function(x){
  cut(x,breaks=c(0,0.0001,0.001,0.01,0.05,1),
      labels=c("****", "***", "**", "*", ""))
}

addAnno <- function(st.combined.sct, Anno_paths, samples) {
  for (i in 1:length(samples)) {
    ManualAnno <- read.csv(Anno_paths[i])
    names(ManualAnno) <- c("Barcode", "Layer")
    ManualAnno$globalKey <-
      str_c(samples[i], ManualAnno$Barcode, sep = "@")
    
    st.combined.sct@meta.data[ManualAnno$globalKey, "Layer"] <-
      ManualAnno$Layer
  }
  st.combined.sct
}

alphaImg <- function(seuratobject, sample = NULL) {
  image_size <-
    seuratobject@images[[sample]]@image %>% {
      list(width = ncol(.), height = nrow(.))
    }
  function(alpha) {
    annotation_custom(
      rasterGrob(
        GetImage(seuratobject, "raw", image = sample)  %>% {
          matrix(rgb(.[, , 1],
                     .[, , 2],
                     .[, , 3], alpha),
                 nrow = nrow(.))
        },
        width = unit(1, "npc"),
        height = unit(1, "npc")
      ),
      xmin = 0,
      xmax = image_size$width,
      ymin = image_size$height,
      ymax = 0
    )
  }
}

LoadOutline <-
  function(seuratobject, sample = NULL, line_paths = NULL) {
    image_size <-
      seuratobject@images[[sample]]@image %>% {
        list(width = ncol(.), height = nrow(.))
      }
    function() {
      annotation_custom(
        rasterGrob(
          read.bitmap(line_paths),
          width = unit(1, "npc"),
          height = unit(1, "npc")
        ),
        xmin = 0,
        xmax = image_size$width,
        ymin = image_size$height,
        ymax = 0
      )
    }
  }

SpatialPlot2 <- function (object,
                          group.by = NULL,
                          features = NULL,
                          images = NULL,
                          cols = NULL,
                          image.alpha = 1,
                          crop = TRUE,
                          slot = "data",
                          min.cutoff = NA,
                          max.cutoff = NA,
                          cells.highlight = NULL,
                          cols.highlight = c("#DE2D26", "grey50"),
                          facet.highlight = FALSE,
                          label = FALSE,
                          label.size = 5,
                          label.color = "white",
                          label.box = TRUE,
                          repel = FALSE,
                          ncol = NULL,
                          combine = TRUE,
                          pt.size.factor = 1.6,
                          alpha = c(1, 1),
                          stroke = 0.25,
                          interactive = FALSE,
                          do.identify = FALSE,
                          identify.ident = NULL,
                          do.hover = FALSE,
                          information = NULL,
                          legend = T,
                          byrow = T)
{
  if (isTRUE(x = do.hover) || isTRUE(x = do.identify)) {
    warning(
      "'do.hover' and 'do.identify' are deprecated as we are removing plotly-based interactive graphics, use 'interactive' instead for Shiny-based interactivity",
      call. = FALSE,
      immediate. = TRUE
    )
    interactive <- TRUE
  }
  if (!is.null(x = group.by) & !is.null(x = features)) {
    stop("Please specific either group.by or features, not both.")
  }
  images <-
    images %||% Images(object = object, assay = DefaultAssay(object = object))
  if (length(x = images) == 0) {
    images <- Images(object = object)
  }
  if (length(x = images) < 1) {
    stop("Could not find any spatial image information")
  }
  if (is.null(x = features)) {
    if (interactive) {
      return(
        ISpatialDimPlot(
          object = object,
          image = images[1],
          group.by = group.by,
          alpha = alpha
        )
      )
    }
    group.by <- group.by %||% "ident"
    object[["ident"]] <- Idents(object = object)
    data <- object[[group.by]]
    for (group in group.by) {
      if (!is.factor(x = data[, group])) {
        data[, group] <- factor(x = data[, group])
      }
    }
  }
  else {
    if (interactive) {
      return(
        ISpatialFeaturePlot(
          object = object,
          feature = features[1],
          image = images[1],
          slot = slot,
          alpha = alpha
        )
      )
    }
    data <- FetchData(object = object,
                      vars = features,
                      slot = slot)
    features <- colnames(x = data)
    min.cutoff <- mapply(
      FUN = function(cutoff, feature) {
        return(ifelse(
          test = is.na(x = cutoff),
          yes = min(data[,
                         feature]),
          no = cutoff
        ))
      },
      cutoff = min.cutoff,
      feature = features
    )
    max.cutoff <- mapply(
      FUN = function(cutoff, feature) {
        return(ifelse(
          test = is.na(x = cutoff),
          yes = max(data[,
                         feature]),
          no = cutoff
        ))
      },
      cutoff = max.cutoff,
      feature = features
    )
    check.lengths <- unique(x = vapply(
      X = list(features,
               min.cutoff, max.cutoff),
      FUN = length,
      FUN.VALUE = numeric(length = 1)
    ))
    if (length(x = check.lengths) != 1) {
      stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    data <- sapply(
      X = 1:ncol(x = data),
      FUN = function(index) {
        data.feature <- as.vector(x = data[, index])
        min.use <- SetQuantile(cutoff = min.cutoff[index],
                               data.feature)
        max.use <- SetQuantile(cutoff = max.cutoff[index],
                               data.feature)
        data.feature[data.feature < min.use] <- min.use
        data.feature[data.feature > max.use] <- max.use
        return(data.feature)
      }
    )
    colnames(x = data) <- features
    rownames(x = data) <- Cells(x = object)
  }
  features <- colnames(x = data)
  colnames(x = data) <- features
  rownames(x = data) <- colnames(x = object)
  facet.highlight <-
    facet.highlight && (!is.null(x = cells.highlight) &&
                          is.list(x = cells.highlight))
  if (do.hover) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning(
        "'do.hover' requires only one image, using image ",
        images,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (length(x = features) > 1) {
      features <- features[1]
      type <- ifelse(test = is.null(x = group.by),
                     yes = "feature",
                     no = "grouping")
      warning(
        "'do.hover' requires only one ",
        type,
        ", using ",
        features,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (facet.highlight) {
      warning(
        "'do.hover' requires no faceting highlighted cells",
        call. = FALSE,
        immediate. = TRUE
      )
      facet.highlight <- FALSE
    }
  }
  if (facet.highlight) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning(
        "Faceting the highlight only works with a single image, using image ",
        images,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    ncols <- length(x = cells.highlight)
  }
  else {
    ncols <- length(x = images)
  }
  plots <- vector(mode = "list", length = length(x = features) *
                    ncols)
  for (i in 1:ncols) {
    plot.idx <- i
    image.idx <- ifelse(test = facet.highlight, yes = 1,
                        no = i)
    image.use <- object[[images[[image.idx]]]]
    coordinates <- GetTissueCoordinates(object = image.use)
    highlight.use <- if (facet.highlight) {
      cells.highlight[i]
    }
    else {
      cells.highlight
    }
    for (j in 1:length(x = features)) {
      cols.unset <- is.factor(x = data[, features[j]]) &&
        is.null(x = cols)
      if (cols.unset) {
        cols <- hue_pal()(n = length(x = levels(x = data[,
                                                         features[j]])))
        names(x = cols) <- levels(x = data[, features[j]])
      }
      plot <- SingleSpatialPlot(
        data = cbind(coordinates,
                     data[rownames(x = coordinates), features[j],
                          drop = FALSE]),
        image = image.use,
        image.alpha = image.alpha,
        col.by = features[j],
        cols = cols,
        alpha.by = if (is.null(x = group.by)) {
          features[j]
        }
        else {
          NULL
        },
        pt.alpha = if (!is.null(x = group.by)) {
          alpha[j]
        }
        else {
          NULL
        },
        geom = if (inherits(x = image.use, what = "STARmap")) {
          "poly"
        }
        else {
          "spatial"
        },
        cells.highlight = highlight.use,
        cols.highlight = cols.highlight,
        pt.size.factor = pt.size.factor,
        stroke = stroke,
        crop = crop
      )
      if (is.null(x = group.by)) {
        plot <- plot + scale_fill_gradientn(name = features[j],
                                            colours = SpatialColors(n = 100)) +
          theme(legend.position = "none",
                legend.title = element_blank()) +
          scale_alpha(range = alpha) + guides(alpha = FALSE)
      }
      else if (label) {
        plot <-
          LabelClusters(
            plot = plot,
            id = ifelse(
              test = is.null(x = cells.highlight),
              yes = features[j],
              no = "highlight"
            ),
            geom = if (inherits(x = image.use,
                                what = "STARmap")) {
              "GeomPolygon"
            }
            else {
              "GeomSpatial"
            },
            repel = repel,
            size = label.size,
            color = label.color,
            box = label.box,
            position = "nearest"
          )
      }
      plot$info <-
        list(image = images[[image.idx]], feature = features[j])
      if (i == 1 && length(x = images) > 1 && !facet.highlight) {
        # plot <- plot + ggtitle(label = features[j]) +
        #   theme(plot.title = element_text(hjust = 0.5))
      }
      if (facet.highlight) {
        plot <- plot + ggtitle(label = names(x = cells.highlight)[i]) +
          theme(plot.title = element_text(hjust = 0.5)) +
          NoLegend()
      }
      plots[[plot.idx]] <- plot
      plot.idx <- plot.idx + ncols
      if (cols.unset) {
        cols <- NULL
      }
    }
  }
  if (length(x = images) > 1 && combine) {
    plots <-
      wrap_plots(plots = plots,
                 ncol = ncol,
                 byrow = byrow)
  }
  else if (length(x = images) == 1 && combine) {
    plots <- wrap_plots(plots = plots, ncol = ncol)
  }
  return(plots)
}

environment(SpatialPlot2) <- environment(SpatialPlot)

Dimplot2 <-
  function (object,
            dims = c(1, 2),
            cells = NULL,
            cols = NULL,
            pt.size = NULL,
            reduction = NULL,
            group.by = NULL,
            split.by = NULL,
            shape.by = NULL,
            order = NULL,
            shuffle = FALSE,
            seed = 1,
            label = FALSE,
            label.size = 4,
            label.color = "black",
            label.box = FALSE,
            repel = FALSE,
            cells.highlight = NULL,
            cols.highlight = "#DE2D26",
            sizes.highlight = 1,
            na.value = "grey50",
            ncol = NULL,
            combine = TRUE,
            raster = NULL)
  {
    if (length(x = dims) != 2) {
      stop("'dims' must be a two-length vector")
    }
    reduction <- reduction %||% DefaultDimReduc(object = object)
    cells <- cells %||% colnames(x = object)
    data <- Embeddings(object = object[[reduction]])[cells, dims]
    data <- as.data.frame(x = data)
    dims <- paste0(Key(object = object[[reduction]]), dims)
    object[["ident"]] <- Idents(object = object)
    orig.groups <- group.by
    group.by <- group.by %||% "ident"
    data <- cbind(data, object[[group.by]][cells, , drop = FALSE])
    group.by <- colnames(x = data)[3:ncol(x = data)]
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
    if (isTRUE(x = shuffle)) {
      set.seed(seed = seed)
      data <- data[sample(x = 1:nrow(x = data)),]
    }
    plots <- lapply(
      X = group.by,
      FUN = function(x) {
        plot <- SingleDimPlot(
          data = data[, c(dims, x, split.by,
                          shape.by)],
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
          na.value = na.value,
          raster = raster
        )
        if (label) {
          plot <- LabelClusters(
            plot = plot,
            id = x,
            repel = repel,
            size = label.size,
            split.by = split.by,
            box = label.box,
            color = label.color
          )
        }
        if (!is.null(x = split.by)) {
          plot <-
            plot + FacetTheme() + facet_wrap(facets = vars(!!sym(x = split.by)),
                                             ncol = if (length(x = group.by) > 1 ||
                                                        is.null(x = ncol)) {
                                               length(x = unique(x = data[, split.by]))
                                             }
                                             else {
                                               ncol
                                             })
        }
        plot <- if (is.null(x = orig.groups)) {
          plot + labs(title = NULL)
        }
        else {
          plot + CenterTitle()
        }
      }
    )
    if (!is.null(x = split.by)) {
      ncol <- 1
    }
    
    if (combine) {
      plots <- wrap_plots(plots, ncol = orig.groups %iff% ncol)
    }
    return(plots)
  }

environment(Dimplot2) <- environment(DimPlot)

FeaturePlot2 <- function(object,
                         features,
                         dims = c(1, 2),
                         cells = NULL,
                         cols = if (blend) {
                           c('lightgrey', '#ff0000', '#00ff00')
                         } else {
                           c('lightgrey', 'blue')
                         },
                         pt.size = NULL,
                         order = FALSE,
                         min.cutoff = NA,
                         max.cutoff = NA,
                         reduction = NULL,
                         split.by = NULL,
                         keep.scale = "feature",
                         shape.by = NULL,
                         slot = 'data',
                         blend = FALSE,
                         blend.threshold = 0.5,
                         label = FALSE,
                         label.size = 4,
                         label.color = "black",
                         repel = FALSE,
                         ncol = NULL,
                         coord.fixed = FALSE,
                         by.col = TRUE,
                         sort.cell = NULL,
                         interactive = FALSE,
                         combine = TRUE,
                         raster = NULL,
                         raster.dpi = c(512, 512)) {
  # TODO: deprecate fully on 3.2.0
  if (!is.null(x = sort.cell)) {
    warning(
      "The sort.cell parameter is being deprecated. Please use the order ",
      "parameter instead for equivalent functionality.",
      call. = FALSE,
      immediate. = TRUE
    )
    if (isTRUE(x = sort.cell)) {
      order <- sort.cell
    }
  }
  if (interactive) {
    return(
      IFeaturePlot(
        object = object,
        feature = features[1],
        dims = dims,
        reduction = reduction,
        slot = slot
      )
    )
  }
  # Check keep.scale param for valid entries
  if (!(is.null(x = keep.scale)) &&
      !(keep.scale %in% c("feature", "all"))) {
    stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
  }
  # Set a theme to remove right-hand Y axis lines
  # Also sets right-hand Y axis text label formatting
  no.right <- theme(
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.y.right = element_text(
      face = "bold",
      size = 14,
      margin = margin(r = 7)
    )
  )
  # Get the DimReduc to use
  reduction <- reduction %||% DefaultDimReduc(object = object)
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  # Figure out blending stuff
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  # Set color scheme for blended FeaturePlots
  if (blend) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(
      EXPR = as.character(x = length(x = cols)),
      '0' = {
        warning(
          "No colors provided, using default colors",
          call. = FALSE,
          immediate. = TRUE
        )
        default.colors
      },
      '1' = {
        warning(
          "Only one color provided, assuming specified is double-negative and augmenting with default colors",
          call. = FALSE,
          immediate. = TRUE
        )
        c(cols, default.colors[2:3])
      },
      '2' = {
        warning(
          "Only two colors provided, assuming specified are for features and agumenting with '",
          default.colors[1],
          "' for double-negatives",
          call. = FALSE,
          immediate. = TRUE
        )
        c(default.colors[1], cols)
      },
      '3' = cols,
      {
        warning(
          "More than three colors provided, using only first three",
          call. = FALSE,
          immediate. = TRUE
        )
        cols[1:3]
      }
    )
  }
  if (blend && length(x = cols) != 3) {
    stop("Blending feature plots only works with three colors; first one for negative cells")
  }
  # Name the reductions
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)
  # Get plotting data
  data <- FetchData(
    object = object,
    vars = c(dims, 'ident', features),
    cells = cells,
    slot = slot
  )
  # Check presence of features/dimensions
  if (ncol(x = data) < 4) {
    stop(
      "None of the requested features were found: ",
      paste(features, collapse = ', '),
      " in slot ",
      slot,
      call. = FALSE
    )
  } else if (!all(dims %in% colnames(x = data))) {
    stop("The dimensions requested were not found", call. = FALSE)
  }
  features <- colnames(x = data)[4:ncol(x = data)]
  # Determine cutoffs
  min.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = min.cutoff,
    feature = features
  )
  max.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = max(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = max.cutoff,
    feature = features
  )
  check.lengths <- unique(x = vapply(
    X = list(features, min.cutoff, max.cutoff),
    FUN = length,
    FUN.VALUE = numeric(length = 1)
  ))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(
    test = length(x = cols) == 1,
    yes = brewer.pal.info[cols,]$maxcolors,
    no = length(x = cols)
  )
  # Apply cutoffs
  data[, 4:ncol(x = data)] <- sapply(
    X = 4:ncol(x = data),
    FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <-
        SetQuantile(cutoff = min.cutoff[index - 3], data.feature)
      max.use <-
        SetQuantile(cutoff = max.cutoff[index - 3], data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      if (brewer.gran == 2) {
        return(data.feature)
      }
      data.cut <- data.feature
      return(data.cut)
    }
  )
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells
  # Figure out splits (FeatureHeatmap)
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  } else {
    switch(EXPR = split.by,
           ident = Idents(object = object)[cells, drop = TRUE],
           object[[split.by, drop = TRUE]][cells, drop = TRUE])
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  # Set shaping variable
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  # Make list of plots
  plots <- vector(mode = "list",
                  length = ifelse(
                    test = blend,
                    yes = 4,
                    no = length(x = features) * length(x = levels(x = data$split))
                  ))
  # Apply common limits
  xlims <-
    c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, dims[1]])))
  ylims <-
    c(floor(min(data[, dims[2]])), ceiling(x = max(data[, dims[2]])))
  # Set blended colors
  if (blend) {
    ncol <- 4
    color.matrix <- BlendMatrix(
      two.colors = cols[2:3],
      col.threshold = blend.threshold,
      negative.color = cols[1]
    )
    cols <- cols[2:3]
    colors <- list(color.matrix[, 1],
                   color.matrix[1,],
                   as.vector(x = color.matrix))
  }
  # Make the plots
  for (i in 1:length(x = levels(x = data$split))) {
    # Figure out which split we're working with
    ident <- levels(x = data$split)[i]
    data.plot <-
      data[as.character(x = data$split) == ident, , drop = FALSE]
    # Blend expression values
    if (blend) {
      features <- features[1:2]
      no.expression <-
        features[colMeans(x = data.plot[, features]) == 0]
      if (length(x = no.expression) != 0) {
        stop(
          "The following features have no value: ",
          paste(no.expression, collapse = ', '),
          call. = FALSE
        )
      }
      data.plot <-
        cbind(data.plot[, c(dims, 'ident')], BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    # Make per-feature plots
    for (j in 1:length(x = features)) {
      feature <- features[j]
      # Get blended colors
      if (blend) {
        cols.use <-
          as.numeric(x = as.character(x = data.plot[, feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      } else {
        cols.use <- NULL
      }
      data.single <-
        data.plot[, c(dims, 'ident', feature, shape.by)]
      # Make the plot
      plot <- SingleDimPlot(
        data = data.single,
        dims = dims,
        col.by = feature,
        order = order,
        pt.size = pt.size,
        cols = cols.use,
        shape.by = shape.by,
        label = FALSE,
        raster = raster
      ) +
        scale_x_continuous(limits = xlims) +
        scale_y_continuous(limits = ylims) +
        theme_cowplot() +
        CenterTitle()
      # theme(plot.title = element_text(hjust = 0.5))
      # Add labels
      if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = 'ident',
          repel = repel,
          size = label.size,
          color = label.color
        )
      }
      # Make FeatureHeatmaps look nice(ish)
      if (length(x = levels(x = data$split)) > 1) {
        plot <-
          plot + theme(panel.border = element_rect(fill = NA, colour = 'black'))
        # Add title
        plot <- plot + if (i == 1) {
          labs(title = feature)
        } else {
          labs(title = NULL)
        }
        # Add second axis
        if (j == length(x = features) && !blend) {
          suppressMessages(
            expr = plot <- plot +
              scale_y_continuous(
                sec.axis = dup_axis(name = ident),
                limits = ylims
              ) +
              no.right
          )
        }
        # Remove left Y axis
        if (j != 1) {
          plot <- plot + theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y.left = element_blank()
          )
        }
        # Remove bottom X axis
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank()
          )
        }
      } else {
        plot <- plot + labs(title = feature)
      }
      # Add colors scale for normal FeaturePlots
      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        } else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warning("All cells have the same value (",
                    unique.feature.exp,
                    ") of ",
                    feature,
                    ".")
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            } else{
              cols.grad <- cols
            }
          }
          plot <-
            suppressMessages(expr = plot + scale_color_gradientn(colors = cols.grad,
                                                                 guide = "colorbar"))
        }
      }
      if (!(is.null(x = keep.scale)) &&
          keep.scale == "feature" && !blend) {
        max.feature.value <- max(data[, feature])
        min.feature.value <- min(data[, feature])
        plot <-
          suppressMessages(plot &
                             scale_color_gradientn(
                               colors = cols,
                               limits = c(min.feature.value, max.feature.value)
                             ))
      }
      # Add coord_fixed
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      # I'm not sure why, but sometimes the damn thing fails without this
      # Thanks ggplot2
      plot <- plot
      # Place the plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  # Add blended color key
  if (blend) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(
        x = plots,
        values = list(
          blend.legend +
            scale_y_continuous(
              sec.axis = dup_axis(name = ifelse(
                test = length(x = levels(x = data$split)) > 1,
                yes = levels(x = data$split)[ii],
                no = ''
              )),
              expand = c(0, 0)
            ) +
            labs(
              x = features[1],
              y = features[2],
              title = if (ii == 1) {
                paste('Color threshold:', blend.threshold)
              } else {
                NULL
              }
            ) +
            no.right
        ),
        after = 4 * ii - 1
      ))
    }
  }
  # Remove NULL plots
  plots <- Filter(f = Negate(f = is.null), x = plots)
  # Combine the plots
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  ncol <- ifelse(
    test = is.null(x = split.by) || blend,
    yes = ncol,
    no = length(x = features)
  )
  legend <- if (blend) {
    'none'
  } else {
    split.by %iff% 'none'
  }
  # Transpose the FeatureHeatmap matrix (not applicable for blended FeaturePlots)
  if (combine) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(
        X = plots,
        FUN = function(x) {
          return(
            suppressMessages(
              expr = x +
                theme_cowplot() +
                ggtitle("") +
                scale_y_continuous(sec.axis = dup_axis(name = ""), limits = ylims) +
                no.right
            )
          )
        }
      )
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) + 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(expr = plots[[i]] +
                                         scale_y_continuous(
                                           sec.axis = dup_axis(name = features[[idx]]),
                                           limits = ylims
                                         ) +
                                         no.right)
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots) %% length(x = features) == 1)) {
        plots[[i]] <- plots[[i]] +
          ggtitle(levels(x = data$split)[[idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] +
            ggtitle(levels(x = data$split)[[idx]]) +
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- 1
        nrow <- nsplits
      } else {
        nrow <- split.by %iff% length(x = levels(x = data$split))
      }
      plots <- plots[c(do.call(
        what = rbind,
        args = split(
          x = 1:length(x = plots),
          f = ceiling(x = seq_along(along.with = 1:length(x = plots)) / length(x = features))
        )
      ))]
      # Set ncol to number of splits (nrow) and nrow to number of features (ncol)
      plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
      if (!is.null(x = legend) && legend == 'none') {
        plots <- plots & NoLegend()
      }
    } else {
      plots <-
        wrap_plots(plots,
                   ncol = ncol,
                   nrow = split.by %iff% length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == 'none') {
      plots <- plots & NoLegend()
    }
    if (!(is.null(x = keep.scale)) &&
        keep.scale == "all" && !blend) {
      max.feature.value <- max(data[, features])
      min.feature.value <- min(data[, features])
      plots <-
        suppressMessages(plots &
                           scale_color_gradientn(
                             colors = cols,
                             limits = c(min.feature.value, max.feature.value)
                           ))
    }
  }
  return(plots)
}

FeaturePlot_size <- function(object,
                             features,
                             dims = c(1, 2),
                             cells = NULL,
                             cols = if (blend) {
                               c('lightgrey', '#ff0000', '#00ff00')
                             } else {
                               c('lightgrey', 'blue')
                             },
                             pt.size = NULL,
                             order = FALSE,
                             min.cutoff = NA,
                             max.cutoff = NA,
                             reduction = NULL,
                             split.by = NULL,
                             keep.scale = "feature",
                             shape.by = NULL,
                             slot = 'data',
                             blend = FALSE,
                             blend.threshold = 0.5,
                             label = FALSE,
                             label.size = 4,
                             label.color = "black",
                             repel = FALSE,
                             ncol = NULL,
                             coord.fixed = FALSE,
                             by.col = TRUE,
                             sort.cell = NULL,
                             interactive = FALSE,
                             combine = TRUE,
                             raster = NULL,
                             raster.dpi = c(512, 512)) {
  # TODO: deprecate fully on 3.2.0
  if (!is.null(x = sort.cell)) {
    warning(
      "The sort.cell parameter is being deprecated. Please use the order ",
      "parameter instead for equivalent functionality.",
      call. = FALSE,
      immediate. = TRUE
    )
    if (isTRUE(x = sort.cell)) {
      order <- sort.cell
    }
  }
  if (interactive) {
    return(
      IFeaturePlot(
        object = object,
        feature = features[1],
        dims = dims,
        reduction = reduction,
        slot = slot
      )
    )
  }
  # Check keep.scale param for valid entries
  if (!(is.null(x = keep.scale)) &&
      !(keep.scale %in% c("feature", "all"))) {
    stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
  }
  # Set a theme to remove right-hand Y axis lines
  # Also sets right-hand Y axis text label formatting
  no.right <- theme(
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.y.right = element_text(
      face = "bold",
      size = 14,
      margin = margin(r = 7)
    )
  )
  # Get the DimReduc to use
  reduction <- reduction %||% DefaultDimReduc(object = object)
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  # Figure out blending stuff
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  # Set color scheme for blended FeaturePlots
  if (blend) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(
      EXPR = as.character(x = length(x = cols)),
      '0' = {
        warning(
          "No colors provided, using default colors",
          call. = FALSE,
          immediate. = TRUE
        )
        default.colors
      },
      '1' = {
        warning(
          "Only one color provided, assuming specified is double-negative and augmenting with default colors",
          call. = FALSE,
          immediate. = TRUE
        )
        c(cols, default.colors[2:3])
      },
      '2' = {
        warning(
          "Only two colors provided, assuming specified are for features and agumenting with '",
          default.colors[1],
          "' for double-negatives",
          call. = FALSE,
          immediate. = TRUE
        )
        c(default.colors[1], cols)
      },
      '3' = cols,
      {
        warning(
          "More than three colors provided, using only first three",
          call. = FALSE,
          immediate. = TRUE
        )
        cols[1:3]
      }
    )
  }
  if (blend && length(x = cols) != 3) {
    stop("Blending feature plots only works with three colors; first one for negative cells")
  }
  # Name the reductions
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)
  # Get plotting data
  data <- FetchData(
    object = object,
    vars = c(dims, 'ident', features),
    cells = cells,
    slot = slot
  )
  # Check presence of features/dimensions
  if (ncol(x = data) < 4) {
    stop(
      "None of the requested features were found: ",
      paste(features, collapse = ', '),
      " in slot ",
      slot,
      call. = FALSE
    )
  } else if (!all(dims %in% colnames(x = data))) {
    stop("The dimensions requested were not found", call. = FALSE)
  }
  features <- colnames(x = data)[4:ncol(x = data)]
  # Determine cutoffs
  min.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = min.cutoff,
    feature = features
  )
  max.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = max(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = max.cutoff,
    feature = features
  )
  check.lengths <- unique(x = vapply(
    X = list(features, min.cutoff, max.cutoff),
    FUN = length,
    FUN.VALUE = numeric(length = 1)
  ))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(
    test = length(x = cols) == 1,
    yes = brewer.pal.info[cols,]$maxcolors,
    no = length(x = cols)
  )
  # Apply cutoffs
  data[, 4:ncol(x = data)] <- sapply(
    X = 4:ncol(x = data),
    FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <-
        SetQuantile(cutoff = min.cutoff[index - 3], data.feature)
      max.use <-
        SetQuantile(cutoff = max.cutoff[index - 3], data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      if (brewer.gran == 2) {
        return(data.feature)
      }
      data.cut <- data.feature
      return(data.cut)
    }
  )
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells
  # Figure out splits (FeatureHeatmap)
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  } else {
    switch(EXPR = split.by,
           ident = Idents(object = object)[cells, drop = TRUE],
           object[[split.by, drop = TRUE]][cells, drop = TRUE])
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  # Set shaping variable
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  # Make list of plots
  plots <- vector(mode = "list",
                  length = ifelse(
                    test = blend,
                    yes = 4,
                    no = length(x = features) * length(x = levels(x = data$split))
                  ))
  # Apply common limits
  xlims <-
    c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, dims[1]])))
  ylims <-
    c(floor(min(data[, dims[2]])), ceiling(x = max(data[, dims[2]])))
  # Set blended colors
  if (blend) {
    ncol <- 4
    color.matrix <- BlendMatrix(
      two.colors = cols[2:3],
      col.threshold = blend.threshold,
      negative.color = cols[1]
    )
    cols <- cols[2:3]
    colors <- list(color.matrix[, 1],
                   color.matrix[1,],
                   as.vector(x = color.matrix))
  }
  # Make the plots
  for (i in 1:length(x = levels(x = data$split))) {
    # Figure out which split we're working with
    ident <- levels(x = data$split)[i]
    data.plot <-
      data[as.character(x = data$split) == ident, , drop = FALSE]
    # Blend expression values
    if (blend) {
      features <- features[1:2]
      no.expression <-
        features[colMeans(x = data.plot[, features]) == 0]
      if (length(x = no.expression) != 0) {
        stop(
          "The following features have no value: ",
          paste(no.expression, collapse = ', '),
          call. = FALSE
        )
      }
      data.plot <-
        cbind(data.plot[, c(dims, 'ident')], BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    # Make per-feature plots
    for (j in 1:length(x = features)) {
      feature <- features[j]
      # Get blended colors
      if (blend) {
        cols.use <-
          as.numeric(x = as.character(x = data.plot[, feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      } else {
        cols.use <- NULL
      }
      data.single <-
        data.plot[, c(dims, 'ident', feature, shape.by)]
      # Make the plot
      plot <- SingleDimPlot(
        data = data.single,
        dims = dims,
        col.by = feature,
        order = order,
        pt.size = pt.size,
        cols = cols.use,
        shape.by = shape.by,
        label = FALSE,
        raster = raster
      ) +
        scale_x_continuous(limits = xlims) +
        scale_y_continuous(limits = ylims) +
        theme_cowplot() +
        CenterTitle()
      # theme(plot.title = element_text(hjust = 0.5))
      # Add labels
      if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = 'ident',
          repel = repel,
          size = label.size,
          color = label.color
        )
      }
      # Make FeatureHeatmaps look nice(ish)
      if (length(x = levels(x = data$split)) > 1) {
        plot <-
          plot + theme(panel.border = element_rect(fill = NA, colour = 'black'))
        # Add title
        plot <- plot + if (i == 1) {
          labs(title = feature)
        } else {
          labs(title = NULL)
        }
        # Add second axis
        if (j == length(x = features) && !blend) {
          suppressMessages(
            expr = plot <- plot +
              scale_y_continuous(
                sec.axis = dup_axis(name = ident),
                limits = ylims
              ) +
              no.right
          )
        }
        # Remove left Y axis
        if (j != 1) {
          plot <- plot + theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y.left = element_blank()
          )
        }
        # Remove bottom X axis
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank()
          )
        }
      } else {
        plot <- plot + labs(title = feature)
      }
      # Add colors scale for normal FeaturePlots
      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        } else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warning("All cells have the same value (",
                    unique.feature.exp,
                    ") of ",
                    feature,
                    ".")
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            } else{
              cols.grad <- cols
            }
          }
          plot <-
            suppressMessages(expr = plot + scale_color_gradientn(colors = cols.grad,
                                                                 guide = "colorbar"))
        }
      }
      if (!(is.null(x = keep.scale)) &&
          keep.scale == "feature" && !blend) {
        max.feature.value <- max(data[, feature])
        min.feature.value <- min(data[, feature])
        plot <-
          suppressMessages(plot &
                             scale_color_gradientn(
                               colors = cols,
                               limits = c(min.feature.value, max.feature.value)
                             ))
      }
      # Add coord_fixed
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      # I'm not sure why, but sometimes the damn thing fails without this
      # Thanks ggplot2
      plot <- plot
      # Place the plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  # Add blended color key
  if (blend) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(
        x = plots,
        values = list(
          blend.legend +
            scale_y_continuous(
              sec.axis = dup_axis(name = ifelse(
                test = length(x = levels(x = data$split)) > 1,
                yes = levels(x = data$split)[ii],
                no = ''
              )),
              expand = c(0, 0)
            ) +
            labs(
              x = features[1],
              y = features[2],
              title = if (ii == 1) {
                paste('Color threshold:', blend.threshold)
              } else {
                NULL
              }
            ) +
            no.right
        ),
        after = 4 * ii - 1
      ))
    }
  }
  # Remove NULL plots
  plots <- Filter(f = Negate(f = is.null), x = plots)
  # Combine the plots
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  ncol <- ifelse(
    test = is.null(x = split.by) || blend,
    yes = ncol,
    no = length(x = features)
  )
  legend <- if (blend) {
    'none'
  } else {
    split.by %iff% 'none'
  }
  # Transpose the FeatureHeatmap matrix (not applicable for blended FeaturePlots)
  if (combine) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(
        X = plots,
        FUN = function(x) {
          return(
            suppressMessages(
              expr = x +
                theme_cowplot() +
                ggtitle("") +
                scale_y_continuous(sec.axis = dup_axis(name = ""), limits = ylims) +
                no.right
            )
          )
        }
      )
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) + 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(expr = plots[[i]] +
                                         scale_y_continuous(
                                           sec.axis = dup_axis(name = features[[idx]]),
                                           limits = ylims
                                         ) +
                                         no.right)
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots) %% length(x = features) == 1)) {
        plots[[i]] <- plots[[i]] +
          ggtitle(levels(x = data$split)[[idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] +
            ggtitle(levels(x = data$split)[[idx]]) +
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- 1
        nrow <- nsplits
      } else {
        nrow <- split.by %iff% length(x = levels(x = data$split))
      }
      plots <- plots[c(do.call(
        what = rbind,
        args = split(
          x = 1:length(x = plots),
          f = ceiling(x = seq_along(along.with = 1:length(x = plots)) / length(x = features))
        )
      ))]
      # Set ncol to number of splits (nrow) and nrow to number of features (ncol)
      plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
      if (!is.null(x = legend) && legend == 'none') {
        plots <- plots & NoLegend()
      }
    } else {
      plots <-
        wrap_plots(plots,
                   ncol = ncol,
                   nrow = split.by %iff% length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == 'none') {
      plots <- plots & NoLegend()
    }
    if (!(is.null(x = keep.scale)) &&
        keep.scale == "all" && !blend) {
      max.feature.value <- max(data[, features])
      min.feature.value <- min(data[, features])
      plots <-
        suppressMessages(plots &
                           scale_color_gradientn(
                             colors = cols,
                             limits = c(min.feature.value, max.feature.value)
                           ))
    }
  }
  return(plots)
}

environment(FeaturePlot2) <- environment(Seurat::FeaturePlot)

SingleDimPlot2 <-
  function (data,
            dims,
            col.by = NULL,
            cols = NULL,
            pt.size = NULL,
            shape.by = NULL,
            alpha.by = NULL,
            order = NULL,
            label = FALSE,
            repel = FALSE,
            label.size = 4,
            cells.highlight = NULL,
            cols.highlight = "#DE2D26",
            sizes.highlight = 1,
            na.value = "grey50",
            raster = NULL) {
    pt.size <- pt.size %||% AutoPointSize(data = data, raster = raster)
    if ((nrow(x = data) > 1e+05) & !isFALSE(raster)) {
      message(
        "Rasterizing points since number of points exceeds 100,000.",
        "\nTo disable this behavior set `raster=FALSE`"
      )
    }
    raster <- raster %||% (nrow(x = data) > 1e+05)
    if (length(x = dims) != 2) {
      stop("'dims' must be a two-length vector")
    }
    if (!is.data.frame(x = data)) {
      data <- as.data.frame(x = data)
    }
    if (is.character(x = dims) &&
        !all(dims %in% colnames(x = data))) {
      stop("Cannot find dimensions to plot in data")
    }
    else if (is.numeric(x = dims)) {
      dims <- colnames(x = data)[dims]
    }
    if (!is.null(x = cells.highlight)) {
      highlight.info <- SetHighlight(
        cells.highlight = cells.highlight,
        cells.all = rownames(x = data),
        sizes.highlight = sizes.highlight %||%
          pt.size,
        cols.highlight = cols.highlight,
        col.base = cols[1] %||%
          "#C3C3C3",
        pt.size = pt.size
      )
      order <- highlight.info$plot.order
      data$highlight <- highlight.info$highlight
      col.by <- "highlight"
      pt.size <- highlight.info$size
      cols <- highlight.info$color
    }
    if (!is.null(x = order) && !is.null(x = col.by)) {
      if (typeof(x = order) == "logical") {
        if (order) {
          data <- data[order(!is.na(x = data[, col.by]),
                             data[, col.by]),]
        }
      }
      else {
        order <- rev(x = c(order, setdiff(
          x = unique(x = data[,
                              col.by]), y = order
        )))
        data[, col.by] <- factor(x = data[, col.by], levels = order)
        new.order <- order(x = data[, col.by])
        data <- data[new.order,]
        if (length(x = pt.size) == length(x = new.order)) {
          pt.size <- pt.size[new.order]
        }
      }
    }
    if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
      warning("Cannot find ", col.by, " in plotting data, not coloring plot")
      col.by <- NULL
    }
    else {
      col.index <- match(x = col.by, table = colnames(x = data))
      if (grepl(pattern = "^\\d", x = col.by)) {
        col.by <- paste0("x", col.by)
      }
      else if (grepl(pattern = "-", x = col.by)) {
        col.by <- gsub(pattern = "-",
                       replacement = ".",
                       x = col.by)
      }
      colnames(x = data)[col.index] <- col.by
    }
    if (!is.null(x = shape.by) &&
        !shape.by %in% colnames(x = data)) {
      warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
    }
    if (!is.null(x = alpha.by) &&
        !alpha.by %in% colnames(x = data)) {
      warning(
        "Cannot find alpha variable ",
        alpha.by,
        " in data, setting to NULL",
        call. = FALSE,
        immediate. = TRUE
      )
      alpha.by <- NULL
    }
    plot <- ggplot(data = data)
    plot <- if (isTRUE(x = raster)) {
      plot + geom_scattermore(
        mapping = aes_string(
          x = dims[1],
          y = dims[2],
          color = paste0("`", col.by, "`"),
          shape = shape.by,
          alpha = alpha.by
        ),
        pointsize = pt.size
      )
    }
    else {
      plot + geom_point(
        mapping = aes_string(
          x = dims[1],
          y = dims[2],
          color = paste0("`", col.by, "`"),
          shape = shape.by,
          alpha = alpha.by,
          size = paste0("`", col.by, "`")
        ),
      )
    }
    plot <-
      plot + guides(color = guide_legend(override.aes = list(size = 3))) +
      labs(color = NULL, title = col.by) + CenterTitle()
    if (label && !is.null(x = col.by)) {
      plot <- LabelClusters(
        plot = plot,
        id = col.by,
        repel = repel,
        size = label.size
      )
    }
    if (!is.null(x = cols)) {
      if (length(x = cols) == 1 && (is.numeric(x = cols) ||
                                    cols %in% rownames(x = brewer.pal.info))) {
        scale <- scale_color_brewer(palette = cols, na.value = na.value)
      }
      else if (length(x = cols) == 1 && (cols %in% c(
        "alphabet",
        "alphabet2",
        "glasbey",
        "polychrome",
        "stepped"
      ))) {
        colors <- DiscretePalette(length(unique(data[[col.by]])),
                                  palette = cols)
        scale <-
          scale_color_manual(values = colors, na.value = na.value)
      }
      else {
        scale <- scale_color_manual(values = cols, na.value = na.value)
      }
      plot <- plot + scale
    }
    plot <- plot + theme_cowplot()
    plot <- scale_size_continuous(range = c(0, 1.5))
    return(plot)
  }
environment(SingleDimPlot2) <- environment(Seurat::SingleDimPlot)

function (plot,
          id,
          clusters = NULL,
          labels = NULL,
          split.by = NULL,
          repel = TRUE,
          box = FALSE,
          geom = "GeomPoint",
          position = "median",
          ...)
{
  xynames <- unlist(x = GetXYAesthetics(plot = plot, geom = geom),
                    use.names = TRUE)
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) &&
      !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <-
    as.character(x = na.omit(object = unique(x = data[,
                                                      id])))
  groups <-
    clusters %||% as.character(x = na.omit(object = unique(x = data[,
                                                                    id])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ",
         paste(groups[!groups %in%
                        possible.clusters], collapse = ","))
  }
  pb <- ggplot_build(plot = plot)
  if (geom == "GeomSpatial") {
    xrange.save <- layer_scales(plot = plot)$x$range$range
    yrange.save <- layer_scales(plot = plot)$y$range$range
    data[, xynames["y"]] = max(data[, xynames["y"]]) - data[,
                                                            xynames["y"]] + min(data[, xynames["y"]])
    if (!pb$plot$plot_env$crop) {
      y.transform <- c(0, nrow(x = pb$plot$plot_env$image)) -
        pb$layout$panel_params[[1]]$y.range
      data[, xynames["y"]] <-
        data[, xynames["y"]] + sum(y.transform)
    }
  }
  data <- cbind(data, color = pb$data[[1]][[1]])
  labels.loc <- lapply(
    X = groups,
    FUN = function(group) {
      data.use <- data[data[, id] == group, , drop = FALSE]
      data.medians <- if (!is.null(x = split.by)) {
        do.call(what = "rbind",
                args = lapply(
                  X = unique(x = data.use[,
                                          split.by]),
                  FUN = function(split) {
                    medians <- apply(
                      X = data.use[data.use[, split.by] ==
                                     split, xynames, drop = FALSE],
                      MARGIN = 2,
                      FUN = median,
                      na.rm = TRUE
                    )
                    medians <-
                      as.data.frame(x = t(x = medians))
                    medians[, split.by] <-
                      split
                    return(medians)
                  }
                ))
      }
      else {
        as.data.frame(x = t(x = apply(
          X = data.use[, xynames,
                       drop = FALSE],
          MARGIN = 2,
          FUN = median,
          na.rm = TRUE
        )))
      }
      data.medians[, id] <- group
      data.medians$color <- data.use$color[1]
      return(data.medians)
    }
  )
  if (position == "nearest") {
    labels.loc <- lapply(
      X = labels.loc,
      FUN = function(x) {
        group.data <- data[as.character(x = data[, id]) ==
                             as.character(x[3]), ]
        nearest.point <-
          nn2(data = group.data[, 1:2],
              query = as.matrix(x = x[c(1,
                                        2)]),
              k = 1)$nn.idx
        x[1:2] <- group.data[nearest.point, 1:2]
        return(x)
      }
    )
  }
  labels.loc <- do.call(what = "rbind", args = labels.loc)
  labels.loc[, id] <-
    factor(x = labels.loc[, id], levels = levels(data[,
                                                      id]))
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop(
      "Length of labels (",
      length(x = labels),
      ") must be equal to the number of clusters being labeled (",
      length(x = labels.loc),
      ")."
    )
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  if (box) {
    geom.use <- ifelse(test = repel, yes = geom_label_repel,
                       no = geom_label)
    plot <-
      plot + geom.use(
        data = labels.loc,
        mapping = aes_string(
          x = xynames["x"],
          y = xynames["y"],
          label = id,
          fill = id
        ),
        show.legend = FALSE,
        ...
      ) + scale_fill_manual(values = labels.loc$color[order(labels.loc[,
                                                                       id])])
  }
  else {
    geom.use <- ifelse(test = repel, yes = geom_text_repel,
                       no = geom_text)
    plot <-
      plot + geom.use(
        data = labels.loc,
        mapping = aes_string(x = xynames["x"],
                             y = xynames["y"], label = id),
        show.legend = FALSE,
        ...
      )
  }
  if (geom == "GeomSpatial") {
    plot <-
      suppressMessages(expr = plot + coord_fixed(xlim = xrange.save,
                                                 ylim = yrange.save))
  }
  return(plot)
}

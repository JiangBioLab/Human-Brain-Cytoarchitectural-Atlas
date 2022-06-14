# Fig1D Global Umap of cell depth ----------------------------------
activedir <- "1-ANALY/Fig1" %T>%  dir.create()

{
  ggplot(umap_all ,
         aes(x = UMAP_1,
             y = UMAP_2,
             colour = aver.depth)) +
    geom_point(size = 0.01) +
    scale_color_gradientn(
      colours = ColorBar2(100),
      limits = c(1, 7),
      breaks = 1:7,
      labels = Layers,
      name = ""
    ) +
    theme_nothing() +
    theme(
      legend.text = element_text(family = "ArialMT", size = 20),
      legend.key.size = unit(30, units = "pt"),
      legend.position = c(0.9, 0.2)
    ) +
    guides(color = guide_colorbar(reverse = T))
} %>%
  ggsave(
    file.path(activedir,  "Fig1D.png"),
    .,
    height = 10,
    width = 10,
    dpi = 600,
    limitsize = F
  )

# FigS3B Non-neural Umap of cell depth ----------------------------------
activedir <- "1-ANALY/FigS3" %T>%  dir.create()
{
  ggplot(
    cbind(EXC_UMAP,
          umap_all[EXC_UMAP$cells, c("hicat_cluster_merge_anno", "sample")]) %>%
      select(UMAP1, UMAP2, hicat_cluster_merge_anno) %>%
      left_join(
        cellType_stat2,
        by = c("hicat_cluster_merge_anno" = "cellType")
      ),
    aes(x = UMAP1,
        y = UMAP2,
        colour = m_)
  ) +
    geom_point(size = 0.05) +
    scale_color_gradientn(
      colours = ColorBar(11),
      limits = c(1, 7),
      breaks = 1:7,
      labels = Layers,
      name = ""
    ) +
    theme_nothing() +
    theme(legend.position = c(0.9, 0.2)) +
    guides(color = guide_colorbar(reverse = T))
} %>%
  ggsave(
    file.path(activedir,  "FigS3B.png"),
    .,
    height = 5,
    width = 5,
    limitsize = F
  )

# Fig6D Non-neural Umap of cell depth ----------------------------------
activedir <- "1-ANALY/Fig6" %T>%  dir.create()
{
  ggplot(
    cbind(NonNeuron_UMAP,
          umap_all[NonNeuron_UMAP$cells, c("hicat_cluster_merge_anno", "sample")]) %>%
      select(UMAP1, UMAP2, hicat_cluster_merge_anno) %>%
      left_join(
        cellType_stat2,
        by = c("hicat_cluster_merge_anno" = "cellType")
      ),
    aes(x = UMAP1,
        y = UMAP2,
        colour = m_)
  ) +
    geom_point(size = 0.05) +
    scale_color_gradientn(
      colours = ColorBar(11),
      limits = c(1, 7),
      breaks = 1:7,
      labels = Layers,
      name = ""
    ) +
    theme_nothing() +
    theme(legend.position = c(0.9, 0.2)) +
    guides(color = guide_colorbar(reverse = T))
} %>%
  ggsave(
    file.path(activedir,  "Fig6D.png"),
    .,
    height = 5,
    width = 5,
    limitsize = F
  )


# Fig2A&S2A Layer annotation -------------------------
activedir <- "1-ANALY/Fig2" %T>%  dir.create()

plots <-   purrr::map(
  Seurat::SpatialDimPlot(
    st.combined.sct,
    cols = color_layer,
    combine = F,
    group.by = "Layer2",
    label = T,
    repel = T,
    label.size = 20
  ),
  ~ . #+ guides(fill = guide_legend(override.aes = list(size = 5)))
)

purrr::walk2(plots,
             samples,
             ~ ggsave(
               file.path(activedir, str_c(.y, ".png")),
               # file.path(activedir, "legend.pdf"),
               .x,
               height = 10,
               width = 10,
               limitsize = F,
             ))

# Fig2B -- Vlnplot for nCount sample--------------------------------
{
  activedir <- "1-ANALY/Fig2" %T>%  dir.create()
  # st.combined.sct$Layer = as.factor(st.combined.sct$Layer)
  # st.combined.sct$Layer_2 = as.factor(st.combined.sct$Layer)
  {
    VlnPlot(
      st.combined.sct,
      # sort = T,
      cols = color_sample,
      features = c("nFeature_Spatial"),
      group.by = c("sample"),
      pt.size = 0
    ) + geom_boxplot() + theme(legend.position = "none") +
      labs(y = "Number of Genes", title = "", x = "Sample")
  } %>%
    {
      ggsave(
        file.path("1-ANALY/Fig2", "Fig2B.pdf"),
        .,
        height =  5,
        width = 5,
        limitsize = F
      )
    }
}

# Fig2C -- Vlnplot for nCount layer -------------------------------
{
  activedir <- "1-ANALY/Fig2" %T>%  dir.create()
  # st.combined.sct$Layer = as.factor(st.combined.sct$Layer)
  # st.combined.sct$Layer_2 = as.factor(st.combined.sct$Layer)
  {
    VlnPlot(
      st.combined.sct,
      # sort = T,
      cols = color_layer,
      features = c("nFeature_Spatial"),
      group.by = c("Layer2"),
      pt.size = 0
    ) + geom_boxplot(fill = "white") + theme(legend.position = "none") +
      labs(y = "Number of Genes", title = "", x = "Layer")
  } %>%
    {
      ggsave(
        file.path("1-ANALY/Fig2" , "Fig2C.pdf"),
        .,
        height =  5,
        width = 5,
        limitsize = F
      )
    }
}


# Fig2D -- Umap layer ------------------------------
{
  activedir <- "1-ANALY/Fig2" %T>%  dir.create()
  # st.combined.sct$Layer = as.factor(st.combined.sct$Layer)
  # st.combined.sct$Layer_2 = as.factor(st.combined.sct$Layer)
  {
    DimPlot(
      st.combined.sct,
      order  = T,
      group.by = "Layer2",
      cols = alpha(color_layer, alpha = 0.5),
      shuffle = T,
      # cells.highlight = subset(st.combined.sct, Layer =="L4") %>% Cells(),
      # group.by = ,
      # ncol = 4,
      label = T
    ) + theme(legend.position = "none", plot.title = element_blank())
  } %>%
    {
      ggsave(
        file.path(activedir, "Fig2D.pdf"),
        .,
        height =  3,
        width = 3,
        limitsize = F
      )
    }
  
  # st.combined.sct$Layer_2 = as.factor(st.combined.sct$Layer)
  #
  # DimPlot(
  #   st.combined.sct,
  #   order  = T,
  #   split.by = "Layer",
  #   # cols = mypal,
  #   group.by = "sample",
  #   # cells.highlight = subset(st.combined.sct, Layer =="L4") %>% Cells(),
  #   # group.by = ,
  #   ncol = 4,
  #   label = T
  # ) %>%
  #   {
  #     ggsave(
  #       file.path(activedir, "Fig2B_2.pdf"),
  #       .,
  #       height =  16,
  #       width = 16,
  #       limitsize = F
  #     )
  #   }
  #
  # DimPlot(
  #   st.combined.sct,
  #   order  = T,
  #   split.by = "sample",
  #   group.by = "Layer",
  #   # cols = mypal,
  #   # cells.highlight = subset(st.combined.sct, Layer =="L4") %>% Cells(),
  #   # group.by = ,
  #   ncol = 4,
  #   label = T
  # ) %>%
  #   {
  #     ggsave(
  #       file.path(activedir, "Fig2B_3.pdf"),
  #       .,
  #       height =  8,
  #       width = 8,
  #       limitsize = F
  #     )
  #   }
}


# Fig2D -- Umap each layer --------
{
  activedir <- "1-ANALY/Fig2" %T>%  dir.create()
  
  for (i in spotAll %>%
       group_split(Layer2) %>%
       purrr::map( ~ slice_max(., order_by = nCount_Spatial, prop = .9))) {
    ggsave(
      file.path(activedir, str_c(i[[1, "Layer2"]], ".png")),
      ggplot(spotAll,
             aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(color = "grey", size = 0.05) +
        ggnewscale::new_scale_color() +
        geom_point(
          data = i,
          mapping = aes(x = UMAP_1, y = UMAP_2),
          color = color_layer[i[[1, "Layer2"]]],
          size = 0.05
        ) +
        scale_color_manual(values = color_layer,
                           guide = guide_legend(override.aes = list(size = 3))) +
        cowplot::theme_nothing() +
        geom_label(
          data = labels.loc %>% filter(ident == i[[1, "Layer2"]]),
          mapping = aes(
            x = x,
            y = y,
            label = ident,
            fill = ident
          ),
          colour = "white",
          size = 7,
          label.size = 0.25
        ) + scale_fill_manual(values = alpha(color_layer[Layers], alpha = 0.5)) +
        theme(legend.position = "none"),
      height =  3,
      width = 3,
      dpi = 300,
      limitsize = F
    )
  }
}

# Fig2D -- Umap gene each layer --------
{
  activedir <- "1-ANALY/Fig2" %T>%  dir.create()
  
  # plot <-
  FeaturePlot2(
    st.combined.sct,
    features = gene_list,
    order = T,
    combine = F,
    cols = ColorBar(10)
  ) %>%
    purrr::map(theme_umap(.)) %>%
    # + geom_label(
    #   data = labels.loc,
    #   mapping = aes(
    #     x = x,
    #     y = y,
    #     label = ident,
    #     fill = ident
    #   ),
    #   colour = "white",
    #   size = 7,
    #   label.size = 0.25
  # ) + scale_fill_manual(values = alpha(color_layer[Layers], alpha = 0.5)) + theme(legend.position = "none"))
  purrr::walk(~ ggsave(
    file.path(activedir, str_c(.$labels$title, ".png")),
    .,
    height = 3,
    width = 3,
    limitsize = F
  ))
  
}


# Fig2G -- Spatial Dimplot Gene expression for 6 regions ---------------------------------
activedir <- "1-ANALY/Fig2" %T>%
  dir.create()
DefaultAssay(st.combined.sct) <- "SCT"

plot <-
  SpatialPlot2(
    st.combined.sct,
    images = sample_rep,
    # image.alpha = 0.2,
    # slot = "scale.data",
    alpha = c(1, 1),
    cols = ColorBar(100),
    features = gene_list,
    byrow = F,
    legend = F,
    ncol = length(gene_list)
  )

# plot <- plot + scale_fill_continuous(trans = "log2")

plot %>% {
  ggsave(
    file.path("1-ANALY/Fig2", "Fig2G.png"),
    .,
    height = length(sample_rep) * 3,
    width = length(gene_list) * 3,
    dpi = 70,
    limitsize = F
  )
}


# FigS2B -- Spatial Dimplot Gene expression for rest 5 regions ---------------------------------
activedir <- "1-ANALY/FigS2" %T>%
  dir.create()

plot <-
  SpatialPlot2(
    st.combined.sct,
    images = setdiff(samples, sample_rep),
    # image.alpha = 0.2,
    # slot = "scale.data",
    alpha = c(1, 1),
    cols = ColorBar(100),
    features = gene_list,
    byrow = F,
    legend = F
  )

# plot <- plot + scale_fill_continuous(trans = "log2")

plot %>% {
  ggsave(
    file.path(activedir, "FigS2B.png"),
    .,
    height = length(sample_rep) * 3,
    width = length(gene_list) * 3,
    dpi = 70,
    limitsize = F
  )
}
# Fig2H gene expression boxplot for each layer-----------------------
activedir <- "1-ANALY/Fig2"  %T>%  dir.create()

tmp <-
  FetchData(st.combined.sct,
            vars = c("sample", gene_list, "Layer2")) %>%
  pivot_longer(
    cols = !c(sample, Layer2),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  group_by(sample, gene, Layer2) %>%
  summarise(aver = median(expression))

plots = list()

for (.gene in gene_list) {
  tmp2 <- tmp %>% filter(gene == .gene)
  tmp2$CN <-
    tmp2$Layer2 %>% factor(levels = sort(Layers)) %>% as.numeric
  sgs = Layers #spots groups to compare
  # each group compare to the median of the rest
  signs <- map(sgs, function(sg) {
    comps <-
      tmp2 %>% mutate(split = Layer2 == sg) %>% group_by(split) %>% group_split()
    t <-
      wilcox.test(comps[[1]]$aver, comps[[2]]$aver)[["p.value"]]
  }) %>% unlist
  signs_star <- signs %>% sign2star() %>% as.character()
  # comb <- paste(format(signs, digits = 3),signs_star,  sep="\n")
  print(.gene)
  print(signs)
  signLab = data.frame(
    star_x = sgs,
    star_y = rep(max(tmp2$aver) * 1.1, length(sgs)),
    star = signs_star
  )
  
  # tmp2_med <- tmp2 %>% group_by(CN) %>% summarise(median = median(aver_prob))
  {
    plots[[.gene]] <-
      ggboxplot(tmp2,
                x = "Layer2",
                y = "aver",
                color = "Layer2")  +
      geom_smooth(
        data = tmp2,
        mapping = aes(x = CN),
        se = 0,
        color = "grey"
      ) +
      geom_point(mapping = aes(color  = Layer2)) +
      stat_compare_means(
        method = "anova",
        label.x = 2.5,
        label.y = max(tmp2$aver) * 1.2
      ) +      # Add global p-value
      geom_text(
        data = signLab,
        mapping = aes(
          x = star_x,
          y = star_y,
          label = star,
          color = star_x
        ),
        size = 5
      ) +
      # geom_text(data = signLab,
      #           mapping = aes(x = p_x,
      #                         y = p_y,
      #                         label = p,
      #                         color = star_x),
      #           size = 5) +
      scale_color_manual(values = color_layer) +
      # stat_compare_means(label = "p.signif",
      #                    method = "wilcox",
      #                    ref.group = ".all.")  +
      labs(y = "Average Probability Per Spot", x = "Cluster", title = .gene) +
      guides(color = F)
  }
}

plots %>% wrap_plots(nrow = 2) %>% {
  ggsave(
    file.path("1-ANALY/Fig2", "Fig2H.pdf"),
    .,
    height = 2 * 3,
    width = length(gene_list) / 2 * 3,
    dpi = 300,
    limitsize = F
  )
}
# Fig2H -- Box and Vlnplot Gene expression for each layer ---------------------------------
activedir <- "1-ANALY/Fig2" %T>%
  dir.create()

st.combined.sct@active.assay <- "integrated"

tmp <-
  FetchData(st.combined.sct,
            vars = c(gene_list, "nCount_SCT", "sample", "Layer2")) %>%
  group_by(Layer2, sample) %>%
  summarise(across(everything(), median))

tmp$Layer2 <-
  factor(tmp$Layer2, levels = sort(unique(tmp$Layer2)))
tmp$LayerN <- as.numeric(tmp$Layer2)

Fig2F_list <- list()

for (gene in gene_list) {
  print(gene)
  # tmp$layer %<>% factor(levels = rev(annotations))
  p <-
    summary(aov(as.formula(str_c(
      gene, "~", "Layer2"
    )), tmp))[[1]] %>% `[`(1, 5) %>% signif(digits = 3)
  
  Fig2F_list[[gene]] <-
    ggplot(tmp,
           aes_string(x = "Layer2",
                      y = gene)) +
    theme(legend.position = "none") +
    labs(
      title = gene,
      subtitle = str_c("P = ", p),
      y = "Normalized Expression",
      x = "Layer"
    ) +
    # geom_jitter()+
    geom_violin(mapping = aes(fill = Layer2), adjust = 4) +
    geom_smooth(mapping = aes(x = LayerN),
                color = "grey50",
                se = F)   +
    geom_boxplot(fill = "white") +
    scale_fill_manual(values = color_layer) +
    theme_cowplot() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 15),
      legend.position = "none"
    )
  
  #
  # ggsave(
  #   str_c("1-ANALY/Fig2/","Fig2F_", gene, ".pdf"),
  #   plot,
  #   height = 3.5,
  #   width = 3.5,
  #   limitsize = F
  # )
  #
}

Fig2F_list %>% wrap_plots(nrow = 2) %>% {
  ggsave(
    file.path("1-ANALY/Fig2", "Fig2H.png"),
    .,
    height = 2 * 3.5,
    width = length(gene_list) / 2 * 3.5,
    dpi = 300,
    limitsize = F
  )
}



# Fig* -- Umap with Cell Types Probability ----------------------------------
{
  activedir <- "1-ANALY/other" %T>%  dir.create()
  
  theme_my <-
    theme_nothing() +
    theme(
      text = element_text(family = "ArialMT"),
      plot.title = element_text(size = 15, hjust = 0.5),
      plot.background = element_rect(fill = "white",
                                     colour = NA),
      panel.background = element_rect(fill = "white",
                                      colour = NA),
      # panel.border = element_line(colour = "black"),
      legend.key = element_rect(fill = "white",
                                colour = NA),
      rect = element_rect(
        fill = "white",
        colour = "white",
        linetype = 1
      ),
      legend.position = c(0.9, 0.15),
      legend.key.width = unit(8, 'points'),
      legend.key.height =  unit(8, 'points'),
      legend.margin = margin(c(0, 0, 0, 0), unit = "mm"),
      plot.margin = margin(c(0, 0.2, 0.3, 0), unit = "in"),
      # legend.title = element_text(size = 3),
      legend.text = element_text(size = 15, margin = unit(c(0, 0, 0, 0), "mm")),
      legend.text.align = 1,
      legend.box.margin = margin(c(0, 0, 0, 0), unit = "mm")
    )
  
  {  nCut = 50
    spotAll$cutX <- cut_interval(spotAll$umap_x,n = nCut)
    centexX <- seq(min(spotAll$umap_x),max(spotAll$umap_x), length.out = nCut*2 +1 )[(1:(nCut*2 +1 )) %%2 ==0]
    spotAll$cutY <- cut_interval(spotAll$umap_y,n = nCut)
    centexY <- seq(min(spotAll$umap_y),max(spotAll$umap_y), length.out = nCut*2 +1 )[(1:(nCut*2 +1 )) %%2 ==0]}

  
  
  for (i in 1:nrow(cellType.meta)) {
    ct <- cellType.meta[i, ] %$% anno
    
   { 
     tmp = spotAll[c("cutX","cutY",ct)]
    names(tmp) <- c("cutX","cutY","prob")
    average_prob <- tmp %>% group_by(cutX, cutY) %>% summarise(z = mean(prob))
    levels(average_prob$cutX) <- centexX
    levels(average_prob$cutY) <- centexY
    average_prob$cutX %<>% as.character %>% as.numeric()
    average_prob$cutY %<>% as.character %>% as.numeric()
    
    
    density_pz = kde2d.weighted(
      average_prob$cutX,
      average_prob$cutY,
      w = average_prob$z,
      n = 100,
      h = 3,
      # lims = c(
      #   min(imagecol),
      #   max(imagecol),
      #   min(imagerow),
      #   max(imagerow)
      # ) 
    ) # row is x ,col is y  =.=#
    
    colnames(density_pz$z) <- density_pz$y
    rownames(density_pz$z) <- density_pz$x
    density_pz_ <-
      density_pz$z %>% as.data.frame %>%
      rownames_to_column("x") %>%
      pivot_longer(cols = !x,
                   names_to = "y",
                   values_to = "z") %>%
      mutate(y = as.numeric(y),
             x = as.numeric(x))
    
    selected_spots = spotAll[c("umap_x", "umap_y",ct)] %>% slice_tail( prop = 0.05)
   }
    
    {
      ggplot() +
        geom_point(
          spotAll,
          mapping = aes(x = umap_x, y = umap_y),
          size = 0.1,
          color = "grey70"
        )     +
        geom_raster(data = density_pz_,
                    mapping = aes(
                      x = x,
                      y = y,
                      fill = z,
                      alpha = z
                    ),interpolate = T) +
        scale_alpha(range = c(0,0.7))+
        scale_fill_gradientn(colors = ColorBar(20)) +
        guides(alpha  =F) +
        new_scale("alpha") +
        new_scale_color() +
        geom_point(
          data = selected_spots,
          mapping = aes_string(
            x = "umap_x",
            y = "umap_y",
            # color = str_c("`", ct, "`"),
            alpha = str_c("`", ct, "`")
            # size =  str_c("`", ct, "`")
          ),
          color = "red",
          #,
          size = 0.1
        )   +
        # scale_alpha_continuous(n.breaks = 2) +
        # scale_color_gradientn(colors = Red,
        #                       limits = c(0,1),
        #                       breaks = c(0,1)) +
        # scale_size_continuous(range = c(0.05,0.5)) +
        theme_my + guides(
          # size = guide_legend("Probabiliy"),
          fill = F,
          alpha = guide_legend(
            "Probabiliy",
            override.aes = list(size = 5),
            title.position = "left",
            title.theme = element_text(
              family = "ArialMT",
              size = 13,
              angle = 90,
              hjust = 0.5
            )
          )
          # alpha = F
        ) +
        labs(title = ct)
    } %>%
      {
        ggsave(
          file.path(activedir, "umap", str_c(
            str_replace_all(str_c(cellType.meta[i, c(4:7, 3)], collapse = "_"), "/", "_"), ".png"
          )),
          .,
          height =  3,
          width = 3,
          dpi = 300,
          limitsize = F
        )
      }
  }
}

# Fig* -- Global Scatterpie -----------------------------------------------------------
activedir <- "1-ANALY/other/scatterpie" %T>%  dir.create()
for (sample in c("M1", "V1_1")) {
  coord <-
    GetTissueCoordinates(st.combined.sct, image = sample) %>% rownames_to_column("globalKey")
  
  tmp <- spotAll_long %>%
    filter(sample.1 == .env$sample) %>%
    # filter(level2 == "EXC") %>%
    group_by(globalKey)   %>%
    slice_max(order_by = probability.norm, n = 5)   %>%
    group_by(globalKey, level4) %>%
    summarise(sumP = sum(probability.norm)) %>%
    left_join(coord, by = c("globalKey" = "globalKey"))
  
  image <- st.combined.sct@images[[sample]]@image
  
  plot <-  ggplot() +
    images[[sample]](0.5) +
    coord_cartesian(
      expand = FALSE,
      ylim = c(images_lim[sample, "up"], images_lim[sample, "down"]),
      xlim = c(images_lim[sample, "left"], images_lim[sample, "right"])
    ) +
    geom_scatterpie(
      data = tmp %>% pivot_wider(
        values_from = sumP,
        names_from = level4,
        values_fill = 0
      ),
      mapping = aes(y = imagerow,
                    x = imagecol),
      pie_scale = 0.4,
      cols = unique(tmp$level4),
      color = NA
    ) +
    outline[[sample]]() +
    scale_fill_manual(values = color_subclass[tmp$level4 %>% unique]) +
    theme_set(theme_bw(base_size = 10)) +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "mm"),
      plot.title  = element_text(hjust  = 0.5),
      axis.title = element_blank(),
      text = element_text(family = "ArialMT")
    )
  
  ggsave(
    file.path(activedir, str_c("Globalscatterpie_", sample, ".pdf")),
    plot,
    # type = "cairo",
    height = unit(nrow(image) / 150, "in"),
    width = unit(ncol(image) / 150, "in"),
    limitsize = F,
    dpi = 800
  )
}

# Fig4D&S4D -- EXC Scatterpie -----------------------------------------------------------
activedir <- "1-ANALY/other/scatterpie_exc" %T>%  dir.create()
for (sample in c("V1_1", "M1")) {
  coord <-
    GetTissueCoordinates(st.combined.sct, image = sample) %>% rownames_to_column("globalKey")
  
  tmp <- spotAll_long %>%
    filter(sample.1 == .env$sample) %>%
    mutate(target = ifelse(level2 == "EXC", yes = level4, no = "other")) %>%
    group_by(globalKey)   %>%
    slice_max(order_by = probability.norm, n = 20) %>%
    group_by(globalKey, target) %>%
    summarise(sumP = sum(probability.norm)) %>%
    left_join(coord, by = c("globalKey" = "globalKey"))
  
  image <- st.combined.sct@images[[sample]]@image
  
  plot <-  ggplot() +
    images[[sample]](0.5) +
    coord_cartesian(
      expand = FALSE,
      ylim = c(images_lim[sample, "up"], images_lim[sample, "down"]),
      xlim = c(images_lim[sample, "left"], images_lim[sample, "right"])
    ) +
    geom_scatterpie(
      data = tmp %>% pivot_wider(
        values_from = sumP,
        names_from = target,
        values_fill = 0
      ),
      mapping = aes(y = imagerow,
                    x = imagecol),
      pie_scale = 0.4,
      cols = unique(tmp$target),
      color = NA
    ) +
    outline[[sample]]() +
    scale_fill_manual(values = color_subclass[tmp$target %>% unique], na.value = NA) +
    theme_set(theme_bw(base_size = 10)) +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "mm"),
      plot.title  = element_text(hjust  = 0.5),
      axis.title = element_blank(),
      text = element_text(family = "ArialMT")
    )
  
  
  ggsave(
    file.path(activedir, str_c("EXCScatterpie_", sample, ".pdf")),
    plot,
    # type = "cairo",
    height = unit(nrow(image) / 150, "in"),
    width = unit(ncol(image) / 150, "in"),
    limitsize = F,
    dpi = 800
  )
}


# Fig* -- Non-neuron Scatterpie -----------------------------------------------------------
activedir <-
  "1-ANALY/other/scatterpie_non_neuron" %T>%  dir.create()
for (sample in samples) {
  coord <-
    GetTissueCoordinates(st.combined.sct, image = sample) %>% rownames_to_column("globalKey")
  
  tmp <- spotAll_long %>%
    filter(sample.1 == .env$sample) %>%
    mutate(target = ifelse(level1 == "Non-neuron", yes = level4, no = "other")) %>%
    group_by(globalKey)   %>%
    slice_max(order_by = probability.norm, n = 20) %>%
    group_by(globalKey, target) %>%
    summarise(sumP = sum(probability.norm)) %>%
    left_join(coord, by = c("globalKey" = "globalKey"))
  
  image <- st.combined.sct@images[[sample]]@image
  
  plot <-  ggplot() +
    images[[sample]](0.5) +
    coord_cartesian(
      expand = FALSE,
      ylim = c(images_lim[sample, "up"], images_lim[sample, "down"]),
      xlim = c(images_lim[sample, "left"], images_lim[sample, "right"])
    ) +
    geom_scatterpie(
      data = tmp %>% pivot_wider(
        values_from = sumP,
        names_from = target,
        values_fill = 0
      ),
      mapping = aes(y = imagerow,
                    x = imagecol),
      pie_scale = 0.4,
      cols = unique(tmp$target),
      color = NA
    ) +
    outline[[sample]]() +
    scale_fill_manual(values = color_subclass[tmp$target %>% unique], na.value = NA) +
    theme_set(theme_bw(base_size = 10)) +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "mm"),
      plot.title  = element_text(hjust  = 0.5),
      axis.title = element_blank(),
      text = element_text(family = "ArialMT")
    )
  
  
  ggsave(
    file.path(activedir, str_c(
      "Non_Neuron_Scatterpie_", sample, ".pdf"
    )),
    plot,
    # type = "cairo",
    height = unit(nrow(image) / 150, "in"),
    width = unit(ncol(image) / 150, "in"),
    limitsize = F,
    dpi = 800
  )
}




# Fig1A  -- Ridge plot --------------------------------------
activedir <- "1-ANALY/Fig1" %T>%  dir.create()

# cellType_layer_density <-
#   spotAll_long %>%
#   group_by(Layer2,cellType) %>%
#   mutate(Layer_n = n()) %>%
#   group_by(cellType) %>%
#   top_frac(n = 0.05,wt = probability)%>%
#   group_by(Layer2,cellType,Layer_n) %>%
#   summarise(p = n()) %>%
#   mutate(average = p/Layer_n)%>%
#   left_join(cellType.meta, by = c("cellType" = "anno"))

cellType_layer_density <-
  spotAll_long %>%
  group_by(sample, cellType) %>%
  mutate(probability.norm = probability / sum(probability)) %>%
  group_by(Layer2, cellType) %>%
  summarise(average = mean(probability.norm)) %>%
  left_join(cellType.meta, by = c("cellType" = "anno"))


cellType_layer_density$Layer2 %<>% factor(levels = sort(unique(cellType_layer_density$Layer2)))

cellType_layer_density$cellType %<>% factor(levels = rev(ct.order))

p_ridge_all <-
  ggplot(
    cellType_layer_density,
    aes(
      x = Layer2,
      y = cellType,
      height = average,
      group = cellType,
      fill = cellType,
      colour = cellType
    )
  ) + geom_ridgeline(scale = 1,
                     size = 0) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(),
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major.x = element_line(colour = "grey20",
                                      linetype = "dotted"),
    axis.text.x = element_text(angle = 70, hjust = 1)
  ) +
  scale_fill_manual(values = color_celltype) +
  scale_colour_manual(values = color_celltype)

ggsave(
  file.path(activedir, "all-ridge.pdf"),
  p_ridge_all,
  height = 35,
  width = 3
)


# Fig3B  -- Heatmap for INH CGE ----------------------------------------------

{
  activedir <- "1-ANALY/Fig3" %T>%  dir.create()
  
  cellType_layer_density <-
    spotAll_long %>%
    filter(str_detect(level3, "CGE")) %>%
    group_by(sample, cellType) %>%
    mutate(probability.norm = probability) %>%
    group_by(Layer2, cellType) %>%
    summarise(average_layer = mean(probability.norm)) %>%
    group_by(cellType) %>%
    mutate(
      m_ = mean(average_layer),
      sd_ = sd(average_layer),
      p.z = (average_layer - m_) / sd_
    ) %>%
    left_join(cellType.meta, by = c("cellType" = "anno"))
  
  
  cellType_layer_density$Layer2 %<>% factor(levels = sort(unique(cellType_layer_density$Layer2)))
  
  cellType_layer_density$cellType %<>% factor(levels = ct.order)
  
  {
    ggplot(cellType_layer_density,
           aes(
             x = cellType,
             y = Layer2,
             # height = average,
             # group = cellType,
             color = p.z,
             # fill = p.z,
             size = p.z
           )) +
      geom_point() +
      # scale_x_discrete(position = "top") +
      theme_nothing() +
      theme(
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "right",
        panel.background = element_blank(),
        axis.text.x = element_text(
          family = "ArialMT",
          angle = -45,
          # vjust = 0.5,
          hjust = 0
        ),
        axis.text.y = element_text(family = "ArialMT",
                                   # angle = 270,
                                   hjust = 0.5),
        plot.margin = margin(20, 5, 20, 5)
      ) +
      scale_colour_gradientn(colours = ColorBar(100)) +
      guides(
        color = guide_legend("Normalized \nProbability"),
        size = guide_legend("Normalized \nProbability")
      ) +
      coord_fixed(ratio = 1)
  } %>% ggsave(file.path(activedir, "Fig3C.pdf"),
               .,
               height = 4,
               width = 12)
}

# Fig3G  -- Heatmap for INH MGE -------------------------------------------------

{
  activedir <- "1-ANALY/Fig3" %T>%  dir.create()
  
  cellType_layer_density <-
    spotAll_long %>%
    filter(str_detect(level3, "MGE")) %>%
    group_by(sample, cellType) %>%
    mutate(probability.norm = probability) %>%
    group_by(Layer2, cellType) %>%
    summarise(average_layer = mean(probability.norm)) %>%
    group_by(cellType) %>%
    mutate(
      m_ = mean(average_layer),
      sd_ = sd(average_layer),
      p.z = (average_layer - m_) / sd_
    ) %>%
    left_join(cellType.meta, by = c("cellType" = "anno"))
  
  
  cellType_layer_density$Layer2 %<>% factor(levels = sort(unique(cellType_layer_density$Layer2)))
  
  cellType_layer_density$cellType %<>% factor(levels = ct.order)
  
  {
    ggplot(cellType_layer_density,
           aes(
             x = cellType,
             y = Layer2,
             # height = average,
             # group = cellType,
             color = p.z,
             # fill = p.z,
             size = p.z
           )) +
      geom_point() +
      scale_x_discrete(position = "top") +
      theme_nothing() +
      theme(
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "right",
        panel.background = element_blank(),
        axis.text.x = element_text(
          family = "ArialMT",
          angle = -45,
          # vjust = 0.5,
          hjust = 0
        ),
        axis.text.y = element_text(family = "ArialMT",
                                   # angle = 270,
                                   hjust = 0.5),
        plot.margin = margin(20, 5, 20, 5)
      ) +
      scale_colour_gradientn(colours = ColorBar(100)) +
      guides(
        color = guide_legend("Normalized \nProbability"),
        size = guide_legend("Normalized \nProbability")
      ) +
      coord_fixed(ratio = 1)
  } %>% ggsave(file.path(activedir, "Fig3G.pdf"),
               .,
               height = 4,
               width = 12)
}
# Fig4B EXC Umap of cell depth ----------------------------------
activedir <- "1-ANALY/Fig4" %T>%  dir.create()

{
  ggplot(umap_EXC ,
         aes(x = UMAP1,
             y = UMAP2,
             colour = aver.depth)) +
    geom_point(size = 0.01) +
    scale_color_gradientn(
      colours = ColorBar2(100),
      limits = c(1, 7),
      breaks = 1:7,
      labels = Layers,
      name = ""
    ) +
    theme_nothing() +
    theme(
      legend.text = element_text(family = "ArialMT", size = 20),
      legend.key.size = unit(30, units = "pt"),
      legend.position = c(0.9, 0.2)
    ) +
    guides(color = guide_colorbar(reverse = T))
} %>%
  ggsave(
    file.path(activedir,  "Fig4B.png"),
    .,
    height = 10,
    width = 10,
    dpi = 600,
    limitsize = F
  )

# Fig4C(down)  -- Heatmap for EXC -------------------------------------------------

{
  activedir <- "1-ANALY/Fig4" %T>%  dir.create()
  
  cellType_layer_density <-
    spotAll_long %>%
    filter(str_detect(level2, "EXC")) %>%
    group_by(sample.x, cellType) %>%
    mutate(probability.norm = probability / sum(probability)) %>%
    group_by(Layer2, cellType) %>%
    summarise(average_layer = mean(probability.norm)) %>%
    group_by(cellType) %>%
    mutate(
      m_ = mean(average_layer),
      sd_ = sd(average_layer),
      p.z = (average_layer - m_) / sd_
    ) %>%
    left_join(cellType.meta, by = c("cellType" = "anno"))
  
  
  cellType_layer_density$Layer2 %<>% factor(levels = sort(unique(cellType_layer_density$Layer2)))
  
  cellType_layer_density$cellType %<>% factor(levels = ct.order)
  
  {
    ggplot(cellType_layer_density,
           aes(
             x = cellType,
             y = Layer2,
             # height = average,
             # group = cellType,
             color = p.z,
             # fill = p.z,
             size = p.z
           )) +
      geom_point() +
      theme_nothing() +
      theme(
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "right",
        panel.background = element_blank(),
        axis.text.x = element_text(
          size = 10,
          family = "ArialMT",
          angle = -90,
          # vjust = 0.5,
          hjust = 0
        ),
        axis.text.y = element_text(family = "ArialMT",
                                   # angle = 270,
                                   hjust = 0.5),
        plot.margin = margin(20, 5, 20, 5)
      ) +
      scale_colour_gradientn(colours = ColorBar(100)) +
      guides(
        color = guide_legend("Probability \n(z-score)"),
        size = guide_legend("Probability \n(z-score)")
      ) +
      coord_fixed(ratio = 1)
  } %>% ggsave(file.path(activedir, "Fig4C_down.pdf"),
               .,
               height = 4,
               width = 14)
}


# Fig5C  -- Heatmap for Non-neural -------------------------------------------------

{
  activedir <- "1-ANALY/Fig5" %T>%  dir.create()
  
  cellType_layer_density <-
    spotAll_long %>%
    filter(str_detect(level1, "Non")) %>%
    group_by(sample, cellType) %>%
    mutate(probability.norm = probability / sum(probability)) %>%
    group_by(Layer2, cellType) %>%
    summarise(average_layer = mean(probability.norm)) %>%
    group_by(cellType) %>%
    mutate(
      m_ = mean(average_layer),
      sd_ = sd(average_layer),
      p.z = (average_layer - m_) / sd_
    ) %>%
    left_join(cellType.meta, by = c("cellType" = "anno"))
  
  
  cellType_layer_density$Layer2 %<>% factor(levels = sort(unique(cellType_layer_density$Layer2)))
  
  cellType_layer_density$cellType %<>% factor(levels = rev(ct.order))
  
  {
    ggplot(cellType_layer_density,
           aes(
             x = cellType,
             y = Layer2,
             # height = average,
             # group = cellType,
             color = p.z,
             # fill = p.z,
             size = p.z
           )) +
      geom_point() +
      theme_nothing() +
      theme(
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "right",
        panel.background = element_blank(),
        axis.text.x = element_text(
          size = 10,
          family = "ArialMT",
          angle = -45,
          # vjust = 0.5,
          hjust = 0
        ),
        axis.text.y = element_text(family = "ArialMT",
                                   # angle = 270,
                                   hjust = 0.5),
        plot.margin = margin(20, 5, 20, 5)
      ) +
      scale_colour_gradientn(colours = ColorBar(100)) +
      guides(
        color = guide_legend("Normalized \nProbability"),
        size = guide_legend("Normalized \nProbability")
      ) +
      coord_fixed(ratio = 1)
  }  %>%
    ggsave(file.path(activedir, "Fig5C.pdf"),
           .,
           height = 4,
           width = 8)
}


# Fig4E -- EXC Deconvolution plot for each tissue ---------------------------------------------
activedir <- "1-ANALY/Fig4"  %T>%  dir.create()

selected_cellTypes <-  rev(c(
  # "EXC FEZF2 FSHR",
  # "EXC THEMIS ITGB8",
  "EXC LINC00507 LINC01331",
  "EXC RORB ALDH1A1",
  "EXC FEZF2 CHST8",
  "EXC THEMIS LINC00343",
  "EXC FEZF2 TSPAN18",
  "EXC FEZF2 NREP"
))

plist <-
  cross2(selected_cellTypes,
         sample_rep[1:3])


plots = list()

{
  for (i in 1:length(plist)) {
    sample_ = plist[[i]][[2]]
    celltype_ = plist[[i]][[1]]
    pz <- filter(spotAll_long,
                 sample == sample_,
                 cellType == celltype_)
    density_pz = kde2d.weighted(
                        pz$imagecol,
                        pz$imagerow,
                        w = pz$probability,
                        n = 100,
                        # h = 20,
                        # lims = c(
                        #   min(imagecol),
                        #   max(imagecol),
                        #   min(imagerow),
                        #   max(imagerow)
                        # ) 
                      ) # row is x ,col is y  =.=#
    colnames(density_pz$z) <- density_pz$y
    rownames(density_pz$z) <- density_pz$x
    density_pz_ <-
      density_pz$z %>% as.data.frame %>%
      rownames_to_column("x") %>%
      pivot_longer(cols = !x,
                   names_to = "y",
                   values_to = "z") %>%
      mutate(y = as.numeric(y),
             x = as.numeric(x))
    
    plots[[i]] <-
    ggplot(pz,
           aes(x = imagecol,
               y = imagerow)) +
      images[[sample_]](1)+
      geom_raster(data = density_pz_,
                  mapping = aes(
                    x = x,
                    y = y,
                    fill = z,
                    alpha = z
                  ),interpolate = T) +
      scale_fill_gradientn(colors = ColorBar(20)) +
      scale_alpha(range = c(0, 1)) +
      new_scale("alpha") +
      new_scale_fill() +
      geom_point(#shape = 21,
                 #stroke = 0.3,
                 #color = "black",
        aes(color = probability,
            size = probability,
            alpha = probability
                     )) +
      scale_alpha(range = c(0, 1)) +
      scale_size_continuous(range = c(0, 2)) +
      scale_color_gradientn(colors = ColorBar(20)) +
      outline[[sample_]]() +
      coord_cartesian(
        expand = FALSE,
        # ylim = c(images_lim[sample_, "y_max"], images_lim[sample_, "y_min"]),
        # xlim = c(images_lim[sample_, "x_min"], images_lim[sample_, "x_max"])
        ylim = c(images_lim[sample_, "up"], images_lim[sample_, "down"]),
        xlim = c(images_lim[sample_, "left"], images_lim[sample_, "right"])
      ) +
      # shape = 21,
      # stroke = 0.25,
      # colour = "black") +
      # scale_colour_gradient2(low = "white",mid = "yellow", high = "red")+    # scale_color_gradientn(colors = ColorBar(100)) +
      
      # scale_fill_manual(values = ColorBar(20)) +
      # ggsci::scale_fill_aaas()+
      # scale_alpha_continuous() +
      # theme_set(theme_bw(base_size = 10)) +
    theme_cowplot() +
      theme(
        legend.position = "none",
        panel.grid.major = element_line(),
        panel.grid.minor = element_line(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        plot.title  = element_text(hjust  = 0.5),
        axis.title = element_blank(),
        text = element_text(family = "ArialMT")
      ) +
      guides(fill = guide_legend("Probability"),
             size = guide_legend("Probability")) +
    Seurat::SpatialTheme()
  }
  
  
  wrap_plots(plots,
             nrow  = length(sample_rep[1:3]),
             byrow = T) %>%
    ggsave(
      file.path(activedir, "Fig4E.png"),
      .,
      type = "cairo",
      height = length(sample_rep[1:3]) * 3,
      width =  length(selected_cellTypes) * 3,
      limitsize = F,
      dpi = 300
    )
  
}

# for (sample in samples) {
#   for (ct in selected_ct) {
#     # values_name = "proportion"
#
#
#
#
#   }
# }



# Fig5F -- Non-neural Deconvolution plot for each tissue ---------------------------------------------
activedir <- "1-ANALY/Fig5"  %T>%  dir.create()

selected_cellTypes <- c(#non-neural
  "ASTRO FGFR3 GFAP",
  "ASTRO FGFR3 COL5A3",
  # "ASTRO FGFR3 RGS20",
  # "MICRO TYROBP LILRB5",
  # "VLMC NOSTRIN SLC6A13",
  # "ENDO NOSTRIN SLC7A5",
  "OLIGO OPALIN MYRF",
  "OPC PDGFRA LRRC4C")

plist <-
  cross2(selected_cellTypes,
         sample_rep)


plots = list()
for (i in 1:length(plist)) {
  sample_ = plist[[i]][[2]]
  celltype_ = plist[[i]][[1]]
  pz <- filter(spotAll_long,
               sample == sample_,
               cellType == celltype_)
  density_pz = kde2d.weighted(
    pz$imagecol,
    pz$imagerow,
    w = pz$probability,
    n = 100,
    # h = 20,
    # lims = c(
    #   min(imagecol),
    #   max(imagecol),
    #   min(imagerow),
    #   max(imagerow)
    # ) 
  ) # row is x ,col is y  =.=#
  colnames(density_pz$z) <- density_pz$y
  rownames(density_pz$z) <- density_pz$x
  density_pz_ <-
    density_pz$z %>% as.data.frame %>%
    rownames_to_column("x") %>%
    pivot_longer(cols = !x,
                 names_to = "y",
                 values_to = "z") %>%
    mutate(y = as.numeric(y),
           x = as.numeric(x))
  
  plots[[i]] <-
    ggplot(pz,
           aes(x = imagecol,
               y = imagerow)) +
    images[[sample_]](1)+
    geom_raster(data = density_pz_,
                mapping = aes(
                  x = x,
                  y = y,
                  fill = z,
                  alpha = z
                ),interpolate = T) +
    scale_fill_gradientn(colors = ColorBar(20)) +
    scale_alpha(range = c(0, 1)) +
    new_scale("alpha") +
    new_scale_fill() +
    geom_point(#shape = 21,
      #stroke = 0.3,
      #color = "black",
      aes(color = probability,
          size = probability,
          alpha = probability
      )) +
    scale_alpha(range = c(0, 1)) +
    scale_size_continuous(range = c(0, 2)) +
    scale_color_gradientn(colors = ColorBar(20)) +
    outline[[sample_]]() +
    coord_cartesian(
      expand = FALSE,
      # ylim = c(images_lim[sample_, "y_max"], images_lim[sample_, "y_min"]),
      # xlim = c(images_lim[sample_, "x_min"], images_lim[sample_, "x_max"])
      ylim = c(images_lim[sample_, "up"], images_lim[sample_, "down"]),
      xlim = c(images_lim[sample_, "left"], images_lim[sample_, "right"])
    ) +
    # shape = 21,
    # stroke = 0.25,
    # colour = "black") +
    # scale_colour_gradient2(low = "white",mid = "yellow", high = "red")+    # scale_color_gradientn(colors = ColorBar(100)) +
    
    # scale_fill_manual(values = ColorBar(20)) +
    # ggsci::scale_fill_aaas()+
    # scale_alpha_continuous() +
    # theme_set(theme_bw(base_size = 10)) +
    theme_cowplot() +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(),
      panel.grid.minor = element_line(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "mm"),
      plot.title  = element_text(hjust  = 0.5),
      axis.title = element_blank(),
      text = element_text(family = "ArialMT")
    ) +
    guides(fill = guide_legend("Probability"),
           size = guide_legend("Probability")) +
    Seurat::SpatialTheme()
}


ggsave(
  file.path(activedir, "Fig5F.png"),
  wrap_plots(plots,
             ncol = length(selected_cellTypes),
             byrow = T),
  type = "cairo",
  height = length(sample_rep)  * 3,
  width = length(selected_cellTypes) * 3,
  limitsize = F,
  dpi = 150
)


# for (sample in samples) {
#   for (ct in selected_ct) {
#     # values_name = "proportion"
#
#
#
#
#   }
# }
# FigS6F -- Non-neural Deconvolution plot for each tissue ---------------------------------------------
activedir <- "1-ANALY/FigS6"  %T>%  dir.create()

selected_cellTypes <- c(
  #non-neural
  # "ASTRO FGFR3 GFAP",
  # "ASTRO FGFR3 COL5A3",
  "ASTRO FGFR3 RGS20",
  "MICRO TYROBP LILRB5",
  "VLMC NOSTRIN SLC6A13",
  "ENDO NOSTRIN SLC7A5"
  # "OLIGO OPALIN MYRF",
  # "OPC PDGFRA LRRC4C"
)

plist <-
  cross2(selected_cellTypes,
         sample_rep)

plots = list()
for (i in 1:length(plist)) {
  sample_ = plist[[i]][[2]]
  celltype_ = plist[[i]][[1]]
  pz <- filter(spotAll_long,
               sample == sample_,
               cellType == celltype_)
  density_pz = kde2d.weighted(
    pz$imagecol,
    pz$imagerow,
    w = pz$probability,
    n = 100,
    # h = 20,
    # lims = c(
    #   min(imagecol),
    #   max(imagecol),
    #   min(imagerow),
    #   max(imagerow)
    # ) 
  ) # row is x ,col is y  =.=#
  colnames(density_pz$z) <- density_pz$y
  rownames(density_pz$z) <- density_pz$x
  density_pz_ <-
    density_pz$z %>% as.data.frame %>%
    rownames_to_column("x") %>%
    pivot_longer(cols = !x,
                 names_to = "y",
                 values_to = "z") %>%
    mutate(y = as.numeric(y),
           x = as.numeric(x))
  
  plots[[i]] <-
    ggplot(pz,
           aes(x = imagecol,
               y = imagerow)) +
    images[[sample_]](1)+
    geom_raster(data = density_pz_,
                mapping = aes(
                  x = x,
                  y = y,
                  fill = z,
                  alpha = z
                ),interpolate = T) +
    scale_fill_gradientn(colors = ColorBar(20)) +
    scale_alpha(range = c(0, 1)) +
    new_scale("alpha") +
    new_scale_fill() +
    geom_point(#shape = 21,
      #stroke = 0.3,
      #color = "black",
      aes(color = probability,
          size = probability,
          alpha = probability
      )) +
    scale_alpha(range = c(0, 1)) +
    scale_size_continuous(range = c(0, 2)) +
    scale_color_gradientn(colors = ColorBar(20)) +
    outline[[sample_]]() +
    coord_cartesian(
      expand = FALSE,
      # ylim = c(images_lim[sample_, "y_max"], images_lim[sample_, "y_min"]),
      # xlim = c(images_lim[sample_, "x_min"], images_lim[sample_, "x_max"])
      ylim = c(images_lim[sample_, "up"], images_lim[sample_, "down"]),
      xlim = c(images_lim[sample_, "left"], images_lim[sample_, "right"])
    ) +
    # shape = 21,
    # stroke = 0.25,
    # colour = "black") +
    # scale_colour_gradient2(low = "white",mid = "yellow", high = "red")+    # scale_color_gradientn(colors = ColorBar(100)) +
    
    # scale_fill_manual(values = ColorBar(20)) +
    # ggsci::scale_fill_aaas()+
    # scale_alpha_continuous() +
    # theme_set(theme_bw(base_size = 10)) +
    theme_cowplot() +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(),
      panel.grid.minor = element_line(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "mm"),
      plot.title  = element_text(hjust  = 0.5),
      axis.title = element_blank(),
      text = element_text(family = "ArialMT")
    ) +
    guides(fill = guide_legend("Probability"),
           size = guide_legend("Probability")) +
    Seurat::SpatialTheme()
}

ggsave(
  file.path(activedir, "FigS6F.png"),
  wrap_plots(plots,
             ncol = length(selected_cellTypes),
             byrow = T),
  type = "cairo",
  height = length(sample_rep)  * 3,
  width = length(selected_cellTypes) * 3,
  limitsize = F,
  dpi = 150
)


# for (sample in samples) {
#   for (ct in selected_ct) {
#     # values_name = "proportion"
#
#
#
#
#   }
# }

# Fig6B -- RORB Umap depth ---------------------------------------------
activedir <- "1-ANALY/Fig6"  %T>%  dir.create()

{
  {
    FeaturePlot(
      object = EXC_6region,
      features = "depth_m",
      split.by = "sample",
      reduction = "tsne",
      ncol = 3,
      # combine = T,
      pt.size = 0.01
    ) & scale_color_gradientn(
      name = "Depth",
      colours = ColorBar(11),
      limits = c(1, 7),
      breaks = 1:7,
      labels = Layers,
      guide = guide_colorbar(reverse = T,
                             ticks.colour = "black")
    ) & theme_nothing()  %+replace%
      theme(plot.title = element_text(family = "ArialMT"))
  } + theme(legend.position = "right")
} %>%
  ggsave(
    file.path("1-ANALY/Fig6", "Fig6B.pdf"),
    .,
    height = 3,
    width = 6 * 3.2,
    limitsize = F
  )


# Fig6B -- RORB Umap depth ---------------------------------------------
activedir <- "1-ANALY/Fig6"  %T>%  dir.create()

{
  {
    FeaturePlot(
      object = EXC_6region,
      features = "depth_m",
      split.by = "sample",
      reduction = "tsne",
      ncol = 3,
      # combine = T,
      pt.size = 0.01
    ) & scale_color_gradientn(
      name = "Depth",
      colours = ColorBar(11),
      limits = c(1, 7),
      breaks = 1:7,
      labels = Layers,
      guide = guide_colorbar(reverse = T,
                             ticks.colour = "black")
    ) & theme_nothing()  %+replace%
      theme(plot.title = element_text(family = "ArialMT"))
  } + theme(legend.position = "right")
} %>%
  ggsave(
    file.path("1-ANALY/Fig6", "Fig6B.pdf"),
    .,
    height = 3,
    width = 6 * 3.2,
    limitsize = F
  )

# Fig6B-2 -- RORB Umap depth ---------------------------------------------
activedir <- "1-ANALY/Fig6"  %T>%  dir.create()

{
  {
    FeaturePlot(
      object = EXC_6region,
      features = "aver.depth",
      split.by = "sample",
      reduction = "tsne",
      ncol = 3,
      # combine = T,
      pt.size = 0.01
    ) & scale_color_gradientn(
      name = "Depth",
      colours =  ColorBar(11),
      limits = c(1, 7),
      breaks = 1:7,
      labels = Layers,
      guide = guide_colorbar(reverse = T,
                             ticks.colour = "black")
    ) & theme_nothing()  %+replace%
      theme(plot.title = element_text(family = "ArialMT"))
  } + theme(legend.position = "right")
} %>%
  ggsave(
    file.path("1-ANALY/Fig6", "Fig6B-2.pdf"),
    .,
    height = 3,
    width = 6 * 3.2,
    limitsize = F
  )


# Fig6C -- RORB Umap ---------------------------------------------
activedir <- "1-ANALY/Fig6"  %T>%  dir.create()

{
  DimPlot(
    object = EXC_6region,
    reduction = "tsne",
    cols = color_subclass[EXC_6region$level4 %>% unique()],
    pt.size = 0.01,
    group.by = "level4",
    # label = T,
    # label.box = T,
    # label.color = "white",
    # label.size = 7,
    repel = T
  ) + theme_nothing()
} %>%
  ggsave(
    file.path("1-ANALY/Fig6", "Fig6C.png"),
    .,
    height = 5,
    width = 5,
    limitsize = F
  )

# Fig6D -- RORB Rose ------------------------------------
{
  activedir <- "1-ANALY/Fig6"  %T>%  dir.create()
  
  cellType_stat_rorb <-
    cellType_stat2 %>%  filter(level3 == "L3-L5_IT_RORB",
                               sample %in% sample_rep) %>%
    mutate(ymax = if_else((m + sd) > 7, 7 , m + sd),
           ymin = if_else((m - sd) < 3, 2, m - sd))
  
  orderofplot <- cellType_stat_rorb$sample %>% unique %>% sort
  
  tmpdata <-
    cellType_stat_rorb %>% filter(sample == orderofplot[1])
  
  tmpdata$sample = ""
  
  tmp <- rbind(cellType_stat_rorb, tmpdata)
  
  tmp$sample <- factor(tmp$sample, levels = c(orderofplot, ""))
}

# tmp %<>% arrange(sample)

{
  ggplot(
    tmp,
    aes(
      x = sample,
      y = m,
      group = cellType,
      colour = cellType,
      fill = cellType,
      ymax = ymax,
      ymin = ymin
    )
  ) +
    geom_hline(
      size = 0.2,
      color = "grey",
      alpha = 1,
      yintercept = 3:6
    ) +
    # geom_ribbon(alpha = 0.6, colour = NA) +
    geom_hline(size = 0.5,
               color = "black",
               # alpha = 0.5,
               yintercept = 4) +
    geom_line(size = 0.5) +
    geom_point(size = 3) +
    theme(legend.position = "none") +
    coord_polar("x", start = pi / 6) +
    facet_wrap( ~ cellType, nrow = 2) +
    theme(
      text = element_text(family = "ArialMT"),
      axis.title.x = element_text(),
      legend.position = "none",
      panel.background = element_blank(),
      axis.text.x = element_text(hjust = 1),
      panel.grid.major.x = element_line(colour = "grey90"),
      panel.grid.major.y = element_blank(),
      strip.background = element_blank(),
      strip.switch.pad.wrap = unit(0, units = "mm"),
      panel.border = element_blank()
    ) +
    scale_y_reverse(
      limits = c(6, 2),
      breaks = 6:3,
      labels = c("L6", "L5", "L4", "L3")
    ) +
    scale_colour_manual(values = color_cluster) +
    scale_fill_manual(values = color_cluster) +
    scale_x_discrete(expand = c(0, 0))
  
} %>% ggsave(
  file.path(activedir, "Fig6D.pdf"),
  .,
  height = 4,
  width = 12,
  limitsize = FALSE
)

# Fig6E -- RORB Umap depth ---------------------------------------------
activedir <- "1-ANALY/Fig6"  %T>%  dir.create()

EXC_6region$selected <-
  purrr::map(EXC_6region$hicat_cluster_merge_anno,
             ~ ifelse(
               . %in% c(
                 "EXC RORB ALDH1A1",
                 "EXC RORB TRABD2A",
                 "EXC RORB GRIK1",
                 "EXC RORB COL22A1",
                 # "EXC RORB COL5A2",
                 "EXC RORB GABRG1"
                 # "EXC RORB FOLH1"
               ),
               .,
               " "
             )) %>% unlist()

Idents(EXC_6region) <- "selected"


{
  {
    FeaturePlot(
      object = EXC_6region %>%
        subset(level4 == "L3-L5_IT_RORB"),
      features = "depth_m",
      split.by = "sample",
      reduction = "tsne",
      label = T,
      label.size = 2,
      # ncol = 7,
      # combine = T,
      pt.size = 0.01,
      repel = T
    ) & scale_color_gradientn(
      name = "Depth",
      colours = ColorBar(11),
      limits = c(1, 7),
      breaks = 1:7,
      labels = Layers,
      guide = guide_colorbar(reverse = T,
                             ticks.colour = "black")
    ) & theme_nothing()  %+replace%
      theme(plot.title = element_text(family = "ArialMT"))
  } + theme(legend.position = "right")
} %>%
  ggsave(
    file.path("1-ANALY/Fig6", "Fig6E.pdf"),
    .,
    height = 3,
    width = 6 * 3.2,
    limitsize = F
  )


# Fig6F -- RORB Deconvolution plot for each tissue ---------------------------------------------
activedir <- "1-ANALY/Fig6"  %T>%  dir.create()

selected_cellTypes <- c("EXC RORB ALDH1A1",
                        "EXC RORB TRABD2A",
                        "EXC RORB GRIK1",
                        # "EXC RORB COL22A1",
                        # "EXC RORB COL5A2",
                        "EXC RORB GABRG1")
# "EXC RORB FOLH1")

plist <-
  cross2(selected_cellTypes,
         sample_rep[1:4])

plots = list()
for (i in 1:length(plist)) {
  sample_ = plist[[i]][[2]]
  celltype_ = plist[[i]][[1]]
  pz <- filter(spotAll_long,
               sample == sample_,
               cellType == celltype_)
  density_pz = kde2d.weighted(
    pz$imagecol,
    pz$imagerow,
    w = pz$probability,
    n = 100,
    # h = 20,
    # lims = c(
    #   min(imagecol),
    #   max(imagecol),
    #   min(imagerow),
    #   max(imagerow)
    # ) 
  ) # row is x ,col is y  =.=#
  colnames(density_pz$z) <- density_pz$y
  rownames(density_pz$z) <- density_pz$x
  density_pz_ <-
    density_pz$z %>% as.data.frame %>%
    rownames_to_column("x") %>%
    pivot_longer(cols = !x,
                 names_to = "y",
                 values_to = "z") %>%
    mutate(y = as.numeric(y),
           x = as.numeric(x))
  
  plots[[i]] <-
    ggplot(pz,
           aes(x = imagecol,
               y = imagerow)) +
    images[[sample_]](1)+
    geom_raster(data = density_pz_,
                mapping = aes(
                  x = x,
                  y = y,
                  fill = z,
                  alpha = z
                ),interpolate = T) +
    scale_fill_gradientn(colors = ColorBar(20)) +
    scale_alpha(range = c(0, 1)) +
    new_scale("alpha") +
    new_scale_fill() +
    geom_point(#shape = 21,
      #stroke = 0.3,
      #color = "black",
      aes(color = probability,
          size = probability,
          alpha = probability
      )) +
    scale_alpha(range = c(0, 1)) +
    scale_size_continuous(range = c(0, 2)) +
    scale_color_gradientn(colors = ColorBar(20)) +
    outline[[sample_]]() +
    coord_cartesian(
      expand = FALSE,
      # ylim = c(images_lim[sample_, "y_max"], images_lim[sample_, "y_min"]),
      # xlim = c(images_lim[sample_, "x_min"], images_lim[sample_, "x_max"])
      ylim = c(images_lim[sample_, "up"], images_lim[sample_, "down"]),
      xlim = c(images_lim[sample_, "left"], images_lim[sample_, "right"])
    ) +
    # shape = 21,
    # stroke = 0.25,
    # colour = "black") +
    # scale_colour_gradient2(low = "white",mid = "yellow", high = "red")+    # scale_color_gradientn(colors = ColorBar(100)) +
    
    # scale_fill_manual(values = ColorBar(20)) +
    # ggsci::scale_fill_aaas()+
    # scale_alpha_continuous() +
    # theme_set(theme_bw(base_size = 10)) +
    theme_cowplot() +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(),
      panel.grid.minor = element_line(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "mm"),
      plot.title  = element_text(hjust  = 0.5),
      axis.title = element_blank(),
      text = element_text(family = "ArialMT")
    ) +
    guides(fill = guide_legend("Probability"),
           size = guide_legend("Probability")) +
    Seurat::SpatialTheme()
}

ggsave(
  file.path(activedir, "Fig6F.png"),
  wrap_plots(plots,
             ncol = length(selected_cellTypes),
             byrow = T),
  type = "cairo",
  height = length(sample_rep[1:4])  * 3,
  width =  length(selected_cellTypes) * 3,
  limitsize = F,
  dpi = 150
)
# Fig6G RORB cell types boxplot -----------------------
activedir <- "1-ANALY/Fig6"  %T>%  dir.create()

Rorb_celltype <- c("EXC RORB ALDH1A1",
                   "EXC RORB TRABD2A",
                   "EXC RORB GRIK1",
                   # "EXC RORB COL5A2",
                   "EXC RORB GABRG1")
# "EXC RORB FOLH1")

tmp <-
  FetchData(st.combined.sct,
            vars = c("sample", Rorb_celltype, "Layer2")) %>%
  pivot_longer(
    cols = !c(sample, Layer2),
    names_to = "celltype",
    values_to = "probability"
  ) %>%
  group_by(sample, celltype, Layer2) %>%
  summarise(aver_prob = mean(probability))

for (ct in Rorb_celltype) {
  tmp2 <- tmp %>% filter(celltype == ct)
  tmp2$CN <-
    tmp2$Layer2 %>% factor(levels = sort(Layers)) %>% as.numeric
  sgs = c("L4", "L5") #spots groups to compare
  # each group compare to the median of the rest
  signs <- map(sgs, function(sg) {
    comps <-
      tmp2 %>% mutate(split = Layer2 == sg) %>% group_by(split) %>% group_split()
    t <-
      wilcox.test(comps[[1]]$aver_prob, comps[[2]]$aver_prob)[["p.value"]]
  }) %>% unlist
  
  signs_star <-
    signs %>% sign2star() %>% as.character()
  comb <-
    paste(format(signs, digits = 3), signs_star,  sep = "\n")
  print(ct)
  print(signs)
  signLab = data.frame(
    star_x = sgs,
    star_y = rep(max(tmp2$aver_prob) * 1.1, length(sgs)),
    star = signs_star
  )
  # p_x = 0,
  # p_y = rep(0,length(sgs)),
  # p = format(signs,digits = 6))
  tmp2_med <-
    tmp2 %>% group_by(CN) %>% summarise(median = median(aver_prob))
  {
    ggboxplot(tmp2,
              x = "Layer2",
              y = "aver_prob",
              color = "Layer2")  +
      geom_smooth(
        data = tmp2,
        mapping = aes(x = CN),
        se = 0,
        color = "grey"
      ) +
      geom_point(mapping = aes(color  = Layer2)) +
      stat_compare_means(
        method = "anova",
        label.x = 2.5,
        label.y = max(tmp2$aver_prob) * 1.2
      ) +      # Add global p-value
      geom_text(
        data = signLab,
        mapping = aes(
          x = star_x,
          y = star_y,
          label = star,
          color = star_x
        ),
        size = 5
      ) +
      # geom_text(data = signLab,
      #           mapping = aes(x = p_x,
      #                         y = p_y,
      #                         label = p,
      #                         color = star_x),
      #           size = 5) +
      scale_color_manual(values = color_layer) +
      # stat_compare_means(label = "p.signif",
      #                    method = "wilcox",
      #                    ref.group = ".all.")  +
      labs(y = "Average Probability Per Spot", x = "Cluster") +
      guides(color = F)
    } %>%
    ggsave(
      str_c(activedir, "/Fig6G_", ct, ".pdf"),
      .,
      height = 3.5,
      width = 3.5,
      limitsize = F
    )
}

# Fig6G -- RORB celltype boxplot ------------------------------------
activedir <-
  "1-ANALY/Fig6"  %T>%  dir.create()

tmp <-
  spotAll_long %>%
  filter(
    cellType %in% c(
      "EXC RORB ALDH1A1",
      "EXC RORB TRABD2A",
      "EXC RORB GRIK1",
      # "EXC RORB COL22A1",
      # "EXC RORB COL5A2",
      "EXC RORB GABRG1"
      # "EXC RORB FOLH1"
    ),
    sample.x %in% sample_rep
  ) %>%
  #     group_by(sample.x, Layer2, cellType) %>%
  #     slice_sample(
  #       n =  spotAll %>%
  #         group_by(sample.x, Layer2) %>%
  #         summarise(Layer_n = n(), ) %>%
  #         ungroup() %$% max(Layer_n) ,
  #       replace = T
  #     )
  # }%>%
  group_by(sample.x, cellType, Layer2) %>%
  summarise(sum = mean(probability))

tmp$Layer2 <-
  factor(tmp$Layer2, levels = sort(unique(tmp$Layer2)))

tmp$LayerN <- as.numeric(tmp$Layer2)

tmp$cellType %<>% factor(levels =  c(
  "EXC RORB ALDH1A1",
  "EXC RORB TRABD2A",
  "EXC RORB GRIK1",
  # "EXC RORB COL22A1",
  # "EXC RORB COL5A2",
  "EXC RORB GABRG1"
  # "EXC RORB FOLH1"
))

tmp_plot <-
  ggboxplot(tmp,
            x = "Layer2",
            y = "sum",
            color  = "Layer2")  +
  geom_smooth(
    data = tmp,
    mapping = aes(x = LayerN, y = sum),
    se = 0,
    color = "grey"
  ) +
  geom_point(mapping = aes(color  = Layer2)) +
  scale_color_manual(values = color_layer) +
  stat_compare_means(method = "anova", label.y = 0.25) +      # Add global p-value
  stat_compare_means(label = "p.signif",
                     method = "wilcox",
                     ref.group = ".all.")        +
  facet_wrap( ~ cellType,
              ncol = 2) +
  labs(y = "Average Probability Per Spot", x = "Layer")

ggsave(
  file.path(activedir, "Fig6G.pdf"),
  tmp_plot,
  height = 5,
  width = 4,
  limitsize = FALSE
)
# theme(plot.subtitle = element_text(hjust = 0.5))


# FigS2 -- S2 level2 concentration ------------------------------
activedir <-
  "1-ANALY/Fig concentration"  %T>%  dir.create()

sd.stat <-
  cellType_stat %>% group_by(level2) %>% summarise(
    median = median(sd),
    mean = mean(sd),
    sd = sd(sd),
    min = mean - sd,
    max = mean + sd
  )

sd.stat$level2 %<>% factor(levels = sd.stat %>% group_by(level2) %>% arrange(mean) %$% level2)

# cellType_stat$level2 %>% factor(levels = sd.stat %>% group_by(level2) %>% arrange(mean) %$% level2) %>% as.numeric()

{
  ggplot(sd.stat,
         aes(
           x = level2,
           y = mean,
           ymin = min,
           ymax = max,
           color = level2,
           fill = level2
         )) +
    # geom_jitter(size = 0.01, alpha = 1) +
    geom_bar(stat = "identity", color = "black") +
    geom_errorbar(color = "black", width = 0.5) +
    geom_pointrange(color = "black") +
    geom_jitter(data = cellType_stat,
                mapping = aes(x = level2,
                              y = sd)) +
    # geom_violin(color = "black") +
    # geom_boxplot(color = "black") +
    scale_fill_manual(values = color_ct) +
    # coord_polar("x", clip = "off", start = -pi / 8) +
    cowplot::theme_minimal_grid() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, hjust = 1)) +
    labs(x = "", y = "s.d of dissection of classes ")
} %>%
  ggsave(file.path(activedir, "main_class.pdf"),
         .,
         width = 4,
         height = 4)

# FigS2 -- level4 concentration ---------------------------------
activedir <-
  "1-ANALY/Fig concentration"  %T>%  dir.create()

cellType_stat$level4 %<>% factor(
  levels = cellType_stat %>% group_by(level4) %>% summarise(median = median(sd)) %>% arrange(median) %$% level4
)

{
  ggplot(cellType_stat,
         aes(
           x = level4,
           y = sd,
           fill = level4,
           color = level4
         )) +
    geom_jitter(alpha = 1, size = 1) +
    geom_violin(color = "black") +
    geom_boxplot(color = "black") +
    # ggsci::scale_fill_npg(alpha = 0.5) +
    scale_fill_manual(values = color_subclass) +
    # cowplot::theme_cowplot() +
    # facet_grid(vars(level2)) +
    labs(x = "", y = "S.d of layer dissection") +
    # coord_polar("x")+
    scale_y_continuous(
      breaks = c(0, 0.5, 1, 1.5, 2, 2.5),
      limits = c(0, 2.5),
      expand = c(0, 0)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    cowplot::theme_cowplot() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 270,
        hjust = 0,
        vjust = 0.5,
        family = "ArialMT"
      )
    ) +
    labs(x = "", y = "S.d of layer distribution")
} %>%
  ggsave(
    file.path(activedir, "subclass_distribution.pdf"),
    height = 7,
    width = 12,
    .
  )


# FigS3 -- INH CGE concentration ---------------------------------
activedir <-
  "1-ANALY/Fig concentration"  %T>%  dir.create()

{
  ggplot(
    cellType_stat %>%
      filter(str_detect(level3, "CGE")) %>%
      mutate(cellType = factor(cellType, levels = ct.order)),
    aes(
      x = cellType,
      y = sd,
      fill = cellType,
      color = cellType
    )
  ) +
    # geom_hline(yintercept = 1.5) +
    geom_violin(color = "black") +
    geom_boxplot(color = "black") +
    geom_jitter(alpha = 0.5, size = 1) +
    # ggsci::scale_fill_npg(alpha = 0.5) +
    scale_colour_manual(values = color_celltype) +
    scale_fill_manual(values = color_celltype) +
    # cowplot::theme_cowplot() +
    # facet_grid(vars(level2)) +
    labs(x = "", y = "S.d of layer dissection") +
    # coord_polar("x")+
    scale_y_continuous(breaks = c(1, 1.5, 2, 2.5),
                       # limits = c(0, 2.5),
                       expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    cowplot::theme_cowplot() +
    theme(
      panel.grid.major.y = element_line(colour = "grey50"),
      legend.position = "none",
      axis.text.x = element_text(
        angle = 270,
        hjust = 0,
        vjust = 0.5,
        family = "ArialMT"
      )
    ) +
    labs(x = "", y = "S.d of layer distribution") +
    coord_fixed(ratio = 1)
  } %>%
  ggsave(
    file.path(activedir, "INH_cge_distribution.pdf"),
    .,
    height = 5,
    width = 12
  )

# FigS3 -- INH MGE concentration ---------------------------------
activedir <-
  "1-ANALY/Fig concentration"  %T>%  dir.create()

{
  ggplot(
    cellType_stat2 %>%
      filter(str_detect(level3, "MGE")) %>%
      mutate(cellType = factor(cellType, levels = ct.order)),
    aes(
      x = cellType,
      y = sd,
      fill = cellType,
      color = cellType
    )
  ) +
    # geom_hline(yintercept = 1.5) +
    geom_violin(color = "black") +
    geom_boxplot(color = "black") +
    geom_jitter(alpha = 0.5, size = 1) +
    # ggsci::scale_fill_npg(alpha = 0.5) +
    scale_colour_manual(values = color_celltype) +
    scale_fill_manual(values = color_celltype) +
    # cowplot::theme_cowplot() +
    # facet_grid(vars(level2)) +
    labs(x = "", y = "S.d of layer dissection") +
    # coord_polar("x")+
    scale_y_continuous(breaks = c(1, 1.5, 2, 2.5),
                       # limits = c(0, 2.5),
                       expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    cowplot::theme_cowplot() +
    theme(
      panel.grid.major.y = element_line(colour = "grey50"),
      legend.position = "none",
      axis.text.x = element_text(
        angle = 270,
        hjust = 0,
        vjust = 0.5,
        family = "ArialMT"
      )
    ) +
    labs(x = "", y = "S.d of layer distribution") +
    coord_fixed(ratio = 1)
  } %>%
  ggsave(
    file.path(activedir, "INH_mge_distribution.pdf"),
    .,
    height = 5,
    width = 12
  )

# Fig4 -- EXC concentration ---------------------------------
activedir <-
  "1-ANALY/Fig concentration"  %T>%  dir.create()

{
  ggplot(
    cellType_stat %>%
      filter(str_detect(level2, "EXC")) %>%
      mutate(cellType = factor(cellType, levels = ct.order)),
    aes(
      x = cellType,
      y = sd,
      fill = cellType,
      color = cellType
    )
  ) +
    # geom_hline(yintercept = 1.5) +
    geom_violin(color = "black") +
    geom_boxplot(color = "black") +
    geom_jitter(alpha = 0.5, size = 1) +
    # ggsci::scale_fill_npg(alpha = 0.5) +
    scale_colour_manual(values = color_celltype) +
    scale_fill_manual(values = color_celltype) +
    # cowplot::theme_cowplot() +
    # facet_grid(vars(level2)) +
    labs(x = "", y = "S.d of layer dissection") +
    # coord_polar("x")+
    scale_y_continuous(breaks = c(1, 1.5, 2, 2.5),
                       # limits = c(0, 2.5),
                       expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    cowplot::theme_cowplot() +
    theme(
      panel.grid.major.y = element_line(colour = "grey50"),
      legend.position = "none",
      axis.text.x = element_text(
        angle = 270,
        hjust = 0,
        vjust = 0.5,
        family = "ArialMT"
      )
    ) +
    labs(x = "", y = "S.d of layer distribution") +
    coord_fixed(ratio = 1)
  } %>%
  ggsave(
    file.path(activedir, "EXC_distribution.pdf"),
    .,
    height = 5,
    width = 12
  )

# FigS5 -- Non-neural concentration ---------------------------------
activedir <-
  "1-ANALY/Fig concentration"  %T>%  dir.create()

{
  ggplot(
    cellType_stat %>% filter(str_detect(level1, "Non-neuron")) %>% mutate(cellType = factor(cellType, levels = ct.order)),
    aes(
      x = cellType,
      y = sd,
      fill = cellType,
      color = cellType
    )
  ) +
    # geom_hline(yintercept = 1.5) +
    geom_violin(color = "black") +
    geom_boxplot(color = "black") +
    geom_jitter(alpha = 0.5, size = 1) +
    # ggsci::scale_fill_npg(alpha = 0.5) +
    scale_colour_manual(values = color_celltype) +
    scale_fill_manual(values = color_celltype) +
    # cowplot::theme_cowplot() +
    # facet_grid(vars(level2)) +
    labs(x = "", y = "S.d of layer dissection") +
    # coord_polar("x")+
    scale_y_continuous(breaks = c(1, 1.5, 2, 2.5),
                       # limits = c(0, 2.5),
                       expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    cowplot::theme_cowplot() +
    theme(
      panel.grid.major.y = element_line(colour = "grey50"),
      legend.position = "none",
      axis.text.x = element_text(
        angle = 270,
        hjust = 0,
        vjust = 0.5,
        family = "ArialMT"
      )
    ) +
    labs(x = "", y = "S.d of layer distribution") +
    coord_fixed(ratio = 1)
  } %>%
  ggsave(
    file.path(activedir, "Non_neural_distribution.pdf"),
    .,
    height = 5,
    width = 12
  )

# FigS5F -- Non-neural Deconvolution plot for each tissue ---------------------------------------------
activedir <- "1-ANALY/FigS5"  %T>%  dir.create()

selected_cellTypes <- c(#non-neural
  # "ASTRO FGFR3 GFAP",
  # "ASTRO FGFR3 COL5A3",
  "ASTRO FGFR3 RGS20",
  "MICRO TYROBP LILRB5",
  "VLMC NOSTRIN SLC6A13",
  "ENDO NOSTRIN SLC7A5")
# "OLIGO OPALIN MYRF",
# "OPC PDGFRA LRRC4C"


plist <-
  cross2(selected_cellTypes,
         sample_rep)
plots = list()

for (i in 1:length(plist)) {
  sample_ = plist[[i]][[2]]
  celltype_ = plist[[i]][[1]]
  pz <-
    filter(spotAll_long,
           sample == sample_,
           cellType == celltype_) %>%
    arrange(probability) %>%
    mutate(zs = (probability - mean(probability)) / sd(probability))
  
  plots[[i]] <-
    ggplot(pz,
           aes(
             x = imagecol,
             y = imagerow,
             color = zs,
             size = zs
           )) +
    images[[sample_]](0.5) +
    geom_point() +
    outline[[sample_]]() +
    # shape = 21,
    # size = 1.5,
    # stroke = 0.01,
    # colour = "black") +
    scale_size_continuous(range = c(0.5, 2)) +
    # scale_alpha(range = c(0.2, 0.5)) +
    coord_cartesian(
      expand = FALSE,
      # ylim = c(images_lim[sample_, "y_max"], images_lim[sample_, "y_min"]),
      # xlim = c(images_lim[sample_, "x_min"], images_lim[sample_, "x_max"])
      ylim = c(images_lim[sample_, "up"], images_lim[sample_, "down"]),
      xlim = c(images_lim[sample_, "left"], images_lim[sample_, "right"])
    ) +
    # scale_fill_gradientn(colors = ColorBar(100)) +
    scale_colour_gradientn(colors = ColorBar(10)[1:9]) +
    # scale_alpha_continuous() +
    # theme_set(theme_bw(base_size = 10)) +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(),
      panel.grid.minor = element_line(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "mm"),
      plot.title  = element_text(hjust  = 0.5),
      axis.title = element_blank(),
      text = element_text(family = "ArialMT")
    ) +
    guides(fill = guide_legend("Probability"),
           size = guide_legend("Probability"))
  
}


ggsave(
  file.path(activedir, "FigS5F.png"),
  wrap_plots(plots,
             ncol = length(selected_cellTypes),
             byrow = T),
  type = "cairo",
  height = length(sample_rep)  * 3,
  width = length(selected_cellTypes) * 3,
  limitsize = F,
  dpi = 150
)



# L6b deconvolution ---------------------------------------------
activedir <- "1-ANALY/Fig 7"  %T>%  dir.create()

#HE的话只能画代表的6个，其他的还没勾线

cellType.meta %>% filter(str_detect(level4, "L6B"))
selected_ct <- c(
  "EXC FEZF2 FSHR",
  "EXC THEMIS ITGB8",
  "EXC THEMIS LINC00343",
  "EXC FEZF2 CHST8",
  "EXC FEZF2 TSPAN18",
  "EXC FEZF2 NREP",
  "EXC LINC00507 LINC01331",
  "EXC RORB ALDH1A1",
  "EXC RORB COL5A2",
  "EXC RORB GABRG1",
  "EXC RORB TRABD2A",
  "EXC RORB GRIK1",
  "EXC RORB FOLH1",
  "EXC RORB COL22A1",
  #non-neural
  "ASTRO FGFR3 GFAP",
  "ASTRO FGFR3 COL5A3",
  "ASTRO FGFR3 RGS20",
  "MICRO TYROBP LILRB5",
  "VLMC NOSTRIN SLC6A13",
  "ENDO NOSTRIN SLC7A5",
  "OLIGO OPALIN MYRF",
  "OPC PDGFRA LRRC4C"
)


for (sample in samples) {
  for (ct in selected_ct) {
    # values_name = "proportion"
    
    image <-
      st.combined.sct@images[[sample]]@image
    
    plot <-
      spotAll_long %>%
      filter(sample == .env$sample, cellType == ct) %>%
      ggplot(
        .,
        aes(
          x = imagecol,
          y = imagerow,
          fill = probability,
          alpha = probability,
          size = probability
        )
      ) +
      annotation_custom(
        images[[sample]],
        xmin = -Inf,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf
      ) +
      geom_point(shape = 21,
                 size = 1.5,
                 # stroke = 0.01,
                 colour = "black") +
      annotation_custom(
        outline[[sample]],
        xmin = -Inf,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf
      ) +
      ylim(nrow(image), 0) +
      xlim(0, ncol(image)) +
      scale_size_continuous(range = c(0, 3)) +
      scale_alpha(range = c(0.2, 1.5)) +
      coord_cartesian(expand = FALSE) +
      scale_fill_gradientn(colors = ColorBar(100)) +
      scale_colour_gradientn(colors = ColorBar(100)) +
      # labs(str_c(sample, "", ct)) +
      theme_set(theme_bw(base_size = 10)) +
      theme(
        legend.position = "none",
        panel.grid.major = element_line(),
        panel.grid.minor = element_line(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        plot.title  = element_text(hjust  = 0.5),
        axis.title = element_blank(),
        text = element_text(family = "ArialMT")
      ) +
      guides(fill = guide_legend("Proportion"),
             # alpha = guide_legend()
             size = guide_legend("Proportion"))
    
    ggsave(
      file.path(activedir, str_replace_all(str_c(
        sample, "_", ct, ".png"
      ),
      " ",
      "_")),
      plot,
      type = "cairo",
      height = unit(nrow(image) / 150, "in"),
      width = unit(ncol(image) / 150, "in"),
      limitsize = F,
      dpi = 300
    )
  }
}


# Section 7 L6b and WM clustering ---------------------------

# DeepLayer <-
#   st.combined.sct %>%
#   subset(Layer2 %in% c("L6", "WM")) %>%
#   RunPCA(verbose = FALSE) %>%
#   RunUMAP(
#     n.neighbors = 50,
#     min.dist = 0.05,
#     reduction = "pca",
#     dims = 1:20,
#     verbose = FALSE
#   ) %>%
#   RunTSNE() %>%
#   FindNeighbors(reduction = "pca",
#                 dims = 1:50) %>%
#   FindClusters(resolution = 0.1)

# change the cluster number following a order of depth for plotting

oldCL2New <- c(
  "2" = 1,
  "0" = 2,
  "4" = 3,
  "1" = 4,
  "3" = 5
)
DeepLayer$seurat_clusters2 <-
  factor(oldCL2New[as.character(DeepLayer$seurat_clusters)], levels = 1:5)
Idents(DeepLayer) <- DeepLayer$seurat_clusters2
cl_L6WM <-
  Idents(DeepLayer) %>% unique() %>% sort()
color_L6WM <-
  c("#008B45FF",
    "#3B4992FF",
    "#008280FF",
    "#EE0000FF",
    "#631879FF") %>%
  set_names(cl_L6WM)

# Fig7A -- L6 WM clustering  ------------------------------------
activedir <- "1-ANALY/Fig7"  %T>%  dir.create()


# for (sample in samples) {
#
#   SpatialPlot2(DeepLayer,
#                images = sample,
#                ncol = 4,
#                cells.highlight = CellsByIdentities(DeepLayer),
#                facet.highlight = T,
#                crop = F,
#   ) %>%
#     ggsave(
#       file.path("1-ANALY/Fig7",str_c("Fig7A_", sample, ".pdf")),
#       .,
#       height = 12,
#       width = 16,
#       limitsize = F
#     )
# }

SpatialPlot2(DeepLayer,
             ncol = 4,
             cols = color_L6WM,
             # cells.highlight = CellsByIdentities(DeepLayer),
             # facet.highlight = T,
             crop = F,
) %>%
  ggsave(
    file.path(activedir, "Fig7A.pdf"),
    .,
    height = 12,
    width = 16,
    limitsize = F
  )

# Fig7B -- umap of L6 and WM --------------------

activedir <- "1-ANALY/Fig7"  %T>%  dir.create()

umap_L6WM <- FetchData(DeepLayer,
                       vars = c("UMAP_1", "UMAP_2",
                                "seurat_clusters"))


label.loc.L6WM <- umap_L6WM %>%
  group_by(seurat_clusters) %>%
  summarise(x = mean(UMAP_1),
            y = mean(UMAP_2))

{
  ggplot(data = umap_L6WM,
         mapping = aes(x = UMAP_1, y = UMAP_2,
                       colour = seurat_clusters)) +
    geom_point(size = 0.1) +
    theme_nothing() +
    geom_label(data = label.loc.L6WM,
               mapping = aes(x = x,
                             y = y,
                             label = seurat_clusters)) +
    scale_colour_manual(values = color_L6WM)
  } %>%
  ggsave(
    file.path(activedir, "Fig7B.pdf"),
    .,
    height = 5,
    width = 5,
    limitsize = F
  )



# Fig7C -- L6 WM original Layer Bar ------------------------------------

activedir <- "1-ANALY/Fig7"  %T>%  dir.create()

{
  ggplot(data =  DeepLayer@meta.data,
         mapping = aes(x = seurat_clusters, fill = Layer2)) +
    geom_bar(position = "fill")   +
    labs(x = "Cluster", y = "Percentage") +
    theme_cowplot() +
    scale_fill_manual(name = "Layer", values = color_layer[unique(DeepLayer@meta.data$Layer2)])
} %>%
  ggsave(
    file.path(activedir, "Fig7C.pdf"),
    .,
    height = 5,
    width = 5,
    limitsize = F
  )

#! Fig 7D traditional Marker ----------------------------------------------

# SpatialFeaturePlot(
#   DeepLayer,
#   features = c("CCN2", "MBP", "THEMIS","DRD1"),
#   crop = F,
#   stroke = 0
# ) %>%
#   ggsave(
#     file.path(activedir, "Fig7D.pdf"),
#     .,
#     height = 15,
#     width = 55,
#     limitsize = F
#   )

# Fig7D -- DE dotplot ------------------------------------

DeGs_L6bWM <-
  FindAllMarkers(DeepLayer,
                 # ident.1 = "0",
                 # ident.2 = c("1", "2"),
                 # features = cl0_markers,
                 # logfc.threshold = 1,
                 # min.pct = 0.5,
                 only.pos = T)


cl0_markers <- c(
  # "NR4A2",
  # "C22orf23",
  # "DIRAS2",
  "DIRAS3",
  "ZFHX3",
  "CCN2",
  "RPL29",
  "MGST1",
  "GLIS2",
  "SULF1",
  "TMEM275",
  "SEMA3E",
  "HS3ST4",
  "LRP1B",
  "NXPH4",
  "MGST1",
  # "CPLX3",
  "DRD1",
  "THEMIS",
  "MBP"
)

# marker <- FetchData(DeepLayer, vars = c(cl0_markers,"seurat_clusters","sample")) %>% pivot_longer(cols = cl0_markers, names_to = "gene", values_to = "exp") %>% group_by(seurat_clusters,gene,sample) %>% summarise(exp_m = mean(exp))
#
# ggplot(marker, aes(x = seurat_clusters, y = exp_m)) +
#   geom_point() +
#   facet_wrap(~gene)

{
  DotPlot(
    DeepLayer,
    features = cl0_markers,
    dot.min = 0,
    dot.scale = 10
  ) +
    labs(x = "Cluster", y = "Gene") +
    theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ))  +
    scale_color_gradientn(colours = c("yellow", "orange", "red"))
} %>%
  ggsave(
    file.path(activedir, "Fig7D.pdf"),
    .,
    height = 4,
    width = 7,
    limitsize = F,
    # dpi = 100
  )

# Fig 7E marker on each slides -----------------------
st.combined.sct@active.assay <- "SCT"

for (gene in cl0_markers) {
  {
    SpatialPlot2(
      st.combined.sct,
      features = gene,
      images = sample_rep2,
      cols = ColorBar(10),
      ncol = 5,
      crop = F
    )
  } %>%
    ggsave(
      file.path(activedir, str_c("Fig7E_", gene, ".png")),
      .,
      height = 3,
      width = 15,
      limitsize = F,
      dpi = 100
    )
}

# # Fig 7E(right) marker boxplot -----------------------
# activedir <- "1-ANALY/Fig7"  %T>%  dir.create()
# cl0_markers5 <- c("MBP",
#                   "THEMIS",
#                   "CCN2",
#                   "SEMA3E",
#                   "MGST1")
#
# DeepLayer@active.assay <- "SCT"
# tmp <- FetchData(DeepLayer,
#                  vars = c(cl0_markers5, "sample", "seurat_clusters")) %>%
#   pivot_longer(
#     cols = !c(seurat_clusters, sample),
#     names_to = "gene",
#     values_to = "exp"
#   ) %>%
#   group_by(gene, seurat_clusters, sample) %>%
#   summarise(aver_exp = mean(exp)) %>% ungroup %>% group_split(gene)
#
# for (tmp2 in tmp) {
#   gene = tmp2$gene[1]
#   {
#     ggboxplot(tmp2,
#               x = "seurat_clusters", y = "aver_exp")  +
#       geom_smooth(
#         data = tmp2,
#         mapping = aes(x = seurat_clusters, y = aver_exp),
#         se = 0,
#         color = "grey"
#       ) +
#       geom_point(mapping = aes(color  = seurat_clusters)) +
#
#       stat_compare_means(method = "anova", label.y = 0.25) +      # Add global p-value
#       stat_compare_means(label = "p.signif",
#                          method = "wilcox",
#                          ref.group = ".all.")        +
#       labs(y = "Expression Per Spot", x = "Layer")
#   } %>%
#     ggsave(
#       str_c("1-ANALY/Fig7/", "Fig7ER_", gene, ".pdf"),
#       .,
#       height = 3.5,
#       width = 3.5,
#       limitsize = F
#     )
#
# }



# Fig 7F cell type of each cluster ---------------------------------
activedir <- "1-ANALY/Fig7"  %T>%  dir.create()
{
  ggplot(
    DotPlot(DeepLayer,
            features = ct.order) %$% data,
    aes(
      x = id,
      y = features.plot,
      colour = avg.exp.scaled,
      size = avg.exp.scaled
    )
  ) +
    geom_point() +
    theme_cowplot() +
    labs(x = "Cluster", y = "cellType") +
    theme(text = element_text(family = "ArialMT"),
          axis.text.x = element_text(# angle = 90,
            hjust = 1,
            vjust = 0.5)) +
    scale_y_discrete(position = "right") +
    scale_color_gradientn(colours = ColorBar(10), breaks = c(-2:2)) +
    scale_size_continuous(breaks = c(-2:2)) +
    guides(size = guide_legend("Scaled Probability"),
           colour = guide_colorbar("Scaled Probability"))
} %>%
  ggsave(
    file.path(activedir, "Fig7F.pdf"),
    .,
    height = 30,
    width = 6,
    limitsize = F,
    # dpi = 100
  )


# Fig 7G cell types on each slides -----------------------
activedir <- "1-ANALY/Fig7"  %T>%  dir.create()

cl0_celltype <- c(
  "ASTRO FGFR3 COL5A3",
  "EXC FEZF2 TSPAN18",
  "ASTRO FGFR3 RGS20",
  "OPC PDGFRA MYT1",
  "OPC PDGFRA PRDM16",
  "OPC PDGFRA LRRC4C",
  "EXC LINC00507 VAT1L",
  "EXC FEZF2 S100Z",
  "EXC FEZF2 IL4R",
  "EXC FEZF2 SLC15A5",
  "EXC FEZF2 SCUBE1"
)


for (ct in cl0_celltype) {
  {
    SpatialPlot2(
      DeepLayer,
      features = ct,
      images = sample_rep2,
      cols = ColorBar(10),
      legend = F,
      ncol = 5,
      crop = F
    ) & theme(panel.border = element_blank())
  } %>%
    ggsave(
      file.path(activedir, str_c("Fig7G_", ct, ".png")),
      .,
      height = 3,
      width = 15,
      limitsize = F,
      dpi = 100
    )
}


# FigSAB -- L6 WM clustering  ------------------------------------
SpatialPlot2(
  DeepLayer,
  images = c("V1_2", "V1_3", "M1", "FPPFC_2",
             "S1_2", "DLPFC_2"),
  ncol = 3,
  cols = color_L6WM,
  # cells.highlight = CellsByIdentities(DeepLayer),
  # facet.highlight = T,
  crop = F,
) %>%
  ggsave(
    file.path(activedir, "FigS9A.pdf"),
    .,
    height = 8,
    width = 12 ,
    limitsize = F
  )


# Fig S9B l6b cell types boxplot -----------------------
activedir <- "1-ANALY/FigS9"  %T>%  dir.create()

cl0_celltype <- c(
  "ASTRO FGFR3 COL5A3",
  "EXC FEZF2 TSPAN18",
  "ASTRO FGFR3 RGS20",
  "OPC PDGFRA MYT1",
  "OPC PDGFRA PRDM16",
  "OPC PDGFRA LRRC4C",
  "EXC LINC00507 VAT1L",
  "EXC FEZF2 S100Z",
  "EXC FEZF2 IL4R",
  "EXC FEZF2 SLC15A5",
  "EXC FEZF2 SCUBE1"
)

tmp <-
  FetchData(DeepLayer,
            vars = c("sample", cl0_celltype, "seurat_clusters", "seurat_clusters2")) %>%
  pivot_longer(
    cols = !c(sample, seurat_clusters, seurat_clusters2),
    names_to = "celltype",
    values_to = "probability"
  ) %>%
  group_by(sample, celltype, seurat_clusters, seurat_clusters2) %>%
  summarise(aver_prob = mean(probability))

for (ct in cl0_celltype) {
  tmp2 <- tmp %>% filter(celltype == ct)
  tmp2$CN <-
    tmp2$seurat_clusters2 %>% as.numeric()
  sgs = c(1, 2, 4) %>% sort()#spots groups to compare
  # each group compare to the median of the rest
  # signs <- map(sgs, function(sg) {
  #   comps <-
  #     tmp2 %>% mutate(split = seurat_clusters2 == sg) %>% group_by(split) %>% group_split()
  #   t <-
  #     wilcox.test(comps[[1]]$aver_prob, comps[[2]]$aver_prob)[["p.value"]]
  # }) %>% unlist %>% p.adjust(method = "fdr")%>% sign2star() %>% as.character()
  # signLab = data.frame( x = sgs,
  #                       y = rep(max(tmp2$aver_prob)*1.1,length(sgs)),
  #                       lab = signs)
  # tmp2_med <- tmp2 %>% group_by(CN) %>% summarise(median = median(aver_prob))
  {
    ggboxplot(tmp2,
              x = "seurat_clusters2",
              y = "aver_prob",
              color = "seurat_clusters2")  +
      geom_smooth(
        data = tmp2,
        mapping = aes(x = CN),
        se = 0,
        color = "grey"
      ) +
      geom_jitter(mapping = aes(color  = seurat_clusters2)) +
      scale_color_manual(values = color_L6WM) +
      stat_compare_means(method = "anova",
                         label.x = 2.5,
                         label.y = 0.2) +      # Add global p-value
      geom_signif(
        comparisons = list(c(1, 2), c(2, 4)),
        test = "wilcox.test",
        # map_signif_level = T,
        step_increase = 0.15,
        textsize = 4
      ) +
      # geom_text(data = signLab,
      #           mapping = aes(x = x,
      #                         y = y,
      #                         label = lab),
      #           size = 10) +
      # stat_compare_means(label = "p.signif",
      #                    method = "wilcox",
      #                    ref.group = ".all.")  +
      labs(y = "Average Probability Per Spot", x = "Cluster") +
      guides(color = F)
  } %>%
    ggsave(
      str_c(activedir, "/Fig9B_", ct, "_withP.pdf"),
      .,
      height = 3.5,
      width = 3.5,
      limitsize = F
    )
}

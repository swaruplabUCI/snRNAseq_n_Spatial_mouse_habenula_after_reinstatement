
custom_vln <- function(
  seurat_obj,
  features,
  group.by,
  split.by=NULL,
  groups = NULL,
  selected_split = NULL,
  plot.margin = margin(0,0,0,0, "cm"),
  slot = 'data',
  assay = NULL,
  raster_dpi=200,
  add_boxplot = TRUE,
  pt.size = 0,
  add.noise = TRUE,
  line.size=NULL,
  adjust=1,
  quantiles=c(0.5),
  stat_method = 'wilcox.test',
  comparisons = NULL,
  pval_y_adjust=0.7,
  ref_group = '.all.',
  group_color_df = NULL,
  split_colors = NULL,
  add_colorbar = TRUE,
  plot_ymin = 0
){

  if(is.null(assay)){assay <- seurat_obj@active.assay}

  seurat_meta <- seurat_obj@meta.data

  if(class(seurat_meta[[group.by]]) != 'factor'){
    seurat_meta[[group.by]] <- as.factor(seurat_meta[[group.by]])
  }

  # get list of groups from seurat object
  if(is.null(groups)){
      groups <- as.character(unique(seurat_obj@meta.data[[group.by]]))
      groups <- groups[order(groups)]
  } else{
    if(!all(groups %in% seurat_obj@meta.data[[group.by]])){
      stop('Some groups not present in seurat_obj@meta.data[[group.by]]')
    }
  }

  if(!is.null(selected_split)){
    if(!all(selected_split %in% seurat_obj@meta.data[[split.by]])){
      stop('Some selected_split not present in seurat_obj@meta.data[[split.by]] ')
    }
  }

  # number of individual plots to make
  n_plots <- length(features) * length(groups)


  # colors for groups
  if(is.null(group_color_df)){

    factor_df <- data.frame(
      level = 1:length(levels(seurat_meta[[group.by]])),
      group_name = levels(seurat_meta[[group.by]])
    )

    p <- seurat_meta %>%
      ggplot(aes_string(x=1, y=1, color=group.by)) +
      geom_point()
    g <- ggplot_build(p)
    g_df <- g$data[[1]]
    group_color_df <- dplyr::select(g_df, c(colour, group)) %>% distinct() %>% arrange(group)
    group_color_df$group <- factor_df$group_name

  } else if(is.character(group_color_df)){
    # check if it's a named character vec:
    if(is.null(names(group_color_df))){
      stop("color_df must be a dataframe with a group column and a colour column, or color_df can be a named character vector where the values are the colors and the names are the corresponding groups.")
    }
    group_color_df <- data.frame(
      colour = as.character(group_color_df),
      group = names(group_color_df)
    )
  }

  # only keep groups that we are using
  group_color_df <- subset(group_color_df, group %in% groups)
  group_color_df$var <- 1

    # color scheme
  group_colors <- group_color_df$colour
  names(group_colors) <- as.character(group_color_df$group)

  # makr color bar in ggplot
  if(add_colorbar){
    colorbar <- group_color_df %>%
      ggplot(aes(x=group, y=var, fill=group)) +
      geom_tile() +
      scale_fill_manual(values=group_color_df$colour) +
      NoLegend() +
      theme(
        plot.title=element_blank(),
        axis.line=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        plot.margin=margin(0,0,0,0)
      ) + RotatedAxis()
  }

  # initialize progress bar
  pb <- utils::txtProgressBar(min = 0, max = n_plots,style = 3, width = 50, char = "=")

  patch_list <- list()
  i <- 1
  for(feature in features){

    plot_df <- seurat_meta

    if(feature %in% colnames(plot_df)){
      plot_df$PlotFeature <- plot_df[[feature]]
    }
    else{
      plot_df$PlotFeature <- GetAssayData(seurat_obj, slot=slot, assay=assay)[feature,]
    }

    plot_range <- range(plot_df$PlotFeature)

    if(add.noise){
      noise <- rnorm(n = length(x = plot_df$PlotFeature)) / 100000
      plot_df$PlotFeature <- plot_df$PlotFeature + noise
    }

    # subset selected groups to compare only:
    if(!is.null(selected_split)){
      plot_df <- plot_df[plot_df[[split.by]] %in% selected_split,]
    }

    plot_list <- list()

    for(cur_group in groups){
      cur_df <- plot_df[as.character(plot_df[[group.by]]) == cur_group,]

      if(is.null(split.by)){
        p <- cur_df %>% ggplot(aes_string(x=group.by, y='PlotFeature', fill=group.by)) +
        scale_fill_manual(values=group_colors, na.value = 'grey90')
      } else{
        p <- cur_df %>% ggplot(aes_string(x=group.by, y='PlotFeature', fill=split.by))
        if(!is.null(split_colors)){
          p <- p + scale_fill_manual(values=split_colors)
        }
      }

      # add violin
      p <- p + geom_violin(
          trim=FALSE, adjust=adjust,
          scale='width', draw_quantiles = quantiles,
          color='black', lwd=0.5
        )

      if(add_boxplot){
        p <- p + geom_boxplot(width=0.1, fill='white', outlier.shape=NA)
      }

      if(pt.size > 0){
        p <- p +  ggrastr::rasterise(ggbeeswarm::geom_beeswarm(size=0.5), dpi=raster_dpi)
      }

      # code for stacked vln plot theme from CellChat:
      p <- p + theme(text = element_text(size = 10)) + theme(axis.line = element_line(size=line.size)) +
        theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 8), axis.line.x = element_line(colour = 'black', size=line.size),axis.line.y = element_line(colour = 'black', size= line.size))
      p <- p + theme(plot.title= element_blank(), # legend.position = "none",
                     axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.y = element_text(size = rel(1), angle = 0),
                     axis.text.y = element_text(size = rel(1)),
                     plot.margin = plot.margin ) +
        theme(axis.text.y = element_text(size = 8))
      p <- p + theme(element_line(size=line.size))

      if(length(plot_list) > 0){
        p <- p + theme(
          legend.position = "none",
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank()
        )
      }


      # y limits:
      ymax <- ceiling(plot_range[2]); ymin <- floor(plot_range[1])
      p <- p + expand_limits(y=c(ymin, ymax)) +
        scale_y_continuous(expand=c(0,0), breaks = c(ymax), limits=c(plot_ymin, NA))

      # stats:
      if(!is.null(split.by)){
        if(is.null(comparisons)){
          p <- p + stat_compare_means(method=stat_method, aes(label = ..p.signif..), label.y=ymax*pval_y_adjust)
        } else{
          p <- p + stat_compare_means(
            method=stat_method,
            label='p.signif',
            comparisons=comparisons,
            label.y=ymax-1.5,
            ref.group=ref_group
          )
        }
      }

      p <- p + ylab(feature)

      # different settings for the last gene:
      if(!add_colorbar & feature == features[length(features)]){
        p <- p + theme(axis.text.x = element_text(), axis.ticks.x = element_line()) +
          RotatedAxis()
      }

      plot_list[[cur_group]] <- p

      # update progress bar
      setTxtProgressBar(pb, i)
      i <- i+1

    }
    patch_list[[feature]] <- wrap_plots(plot_list, ncol=length(groups))

  }

  # close progress bar
  close(pb)

  if(length(patch_list) == 1){
    return(patch_list[[1]])
  }

  plot_heights <- rep(50 / length(features), length(features))

  if(add_colorbar){
    plot_heights <- c(plot_heights, 1)
    patch_list <- c(patch_list, list(colorbar))
  }

  wrap_plots(patch_list, ncol=1) + plot_layout(heights=plot_heights, guides='collect')
}

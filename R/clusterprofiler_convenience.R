#' Dotplot for compareClusterResult object
#'
#' Dotplot for compareClusterResult object
#' Provides more control over the plot than the method from clusterProfiler
#'
#'
#' @param compareClusterResult compareClusterResult object
#' @param top_n numeric value for how many top gene sets to keep from each cluster
#' @param facet_rgx character value for regex to str_extract from the 'Cluster' column for facetting
#' @param facet_lvls character value for levels for facetting variable. Should correspond to 'facet_rgx'
#' @param color character for variable for dot color aesthetic
#' @param size character for variable for dot size aesthetic
#' @param xlab x axis title
#' @param ylab y axis title
#' @param fontsize font size
#' @param title plot title
#' @param cluster_lvls numeric vector to rearrange the order of the cluster levels. Use 'print_clust_lvls' to figure out the original order of the cluster levels.
#' @param print_clust_lvls Don't plot anything. Just print the cluster levels after dropping those with no results (significant results)
#' @return GGplot object if print_clust_lvls set to FALSE
#' @import dplyr stringr ggplot2 tidyr tibble
#' @importFrom DOSE theme_dose
#' @importFrom rlang .data
#' @export
dotplot_compareClusterResult <- function(compareClusterResult,
                                         top_n=5,
                                         facet_rgx="\\S",
                                         facet_lvls="foo",
                                         color="p.adjust",
                                         size="GeneRatio",
                                         xlab="",
                                         ylab="",
                                         fontsize=12,
                                         title="",
                                         cluster_lvls=1:length(levels(droplevels(as_tibble(compareClusterResult)$Cluster))),
                                         print_clust_lvls=FALSE){

  if (print_clust_lvls) {
    all_nonmissing_levels <- levels(forcats::fct_drop(as_tibble(compareClusterResult)$Cluster))
    message("Just printing the non-missing cluster levels for user to edit 'cluster_lvls' param.")
    message(paste(paste0(1:length(all_nonmissing_levels), "\t", all_nonmissing_levels), collapse="\n"))
  } else{
    stopifnot(is(compareClusterResult, "compareClusterResult"))
    if(!all(stringr::str_detect(facet_lvls, facet_rgx))) stop("facet_rgx does not match facet_lvls")
    if(!all(1:length(levels(forcats::fct_drop(as_tibble(compareClusterResult)$Cluster))) == sort(cluster_lvls))) stop("Number of cluster levels and provided cluster level indices do not match.")

    num_genes_annot <- as_tibble(compareClusterResult)[, c("Cluster","GeneRatio")] %>%
      dplyr::mutate(num_genes=as.numeric(str_extract(.data$GeneRatio, "\\d+$"))) %>%
      dplyr::select(-.data$GeneRatio) %>% unique()
    num_genes_annot <- rlang::set_names(num_genes_annot$num_genes, num_genes_annot$Cluster)

    df <- as_tibble(compareClusterResult) %>%
      dplyr::filter(.data$Description %in% (as_tibble(compareClusterResult) %>%
                                        dplyr::group_by(.data$Cluster) %>%
                                        dplyr::arrange(.data$pvalue) %>%
                                        dplyr::filter(dplyr::row_number() <= top_n) %>% # get the top n genesets for each cluster
                                        dplyr::ungroup() %>% dplyr::pull(.data$Description) %>%
                                        unique())) %>%
      dplyr::mutate(facet_var=factor(str_extract(.data$Cluster, facet_rgx), levels=facet_lvls),
                    GeneRatioNum=as.numeric(str_extract(.data$GeneRatio, "^\\d+")),
                    GeneRatioDenom=as.numeric(str_extract(.data$GeneRatio, "\\d+$")),
                    GeneRatio=round(.data$GeneRatioNum/.data$GeneRatioDenom, 3),
                    Cluster=factor(.data$Cluster, levels=levels(forcats::fct_drop(as_tibble(compareClusterResult)$Cluster))[cluster_lvls])) %>%
      dplyr::group_by(.data$Description) %>%
      dplyr::mutate(count=n()) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(.data$count, .data$Cluster) %>% # Sort by the # of occurences of each geneset
      dplyr::mutate(Description=forcats::fct_relevel(.data$Description, unique(.data$Description))) # relevel the Description by the # of occurences of each geneset

    plot <- ggplot(df, aes_string(x="Cluster", y="Description", color=color, size=size)) +
      geom_point() +
      scale_color_continuous(low = "red", high = "blue", guide = guide_colorbar(reverse = TRUE)) +
      scale_x_discrete(labels = function(x) str_replace(paste0(x, " \n(", num_genes_annot[x], ")"), "[\\._]vs[\\._]" , "\nvs\n")) +
      scale_y_discrete(labels = function(x) str_trunc(x, 50)) +
      xlab(xlab) + ylab(ylab) + ggtitle(title) +
      DOSE::theme_dose(fontsize)

    if (facet_rgx!="\\S") {
      plot <- plot + facet_wrap(~facet_var, nrow=1, scales = "free_x")
    }

    return(plot)

  }
}

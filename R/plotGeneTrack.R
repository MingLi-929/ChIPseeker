#' add gene track to plot 
#' 
#' @param plot ggplot2 object, gene track will be added to it.
#' @param txdb txdb object, providing gene annotation.
#' @param chr chromosome id.
#' @param start_pos start coordiante of windows.
#' @param end_pos end coordiante of windows.
#' @param xlab x lab.
#' @param gene_symbol_oi gene of interest to show in the form of gene symbol.
#' @param gene_id_oi gene of interest to show in the form of gene id
#' @param sample_n_gene sample N gene to show if there are too many genes. 
#'  This parameter is mutually exclusive with the gene_symbol_oi and gene_id_oi parameter.
#' @param palette palette, default "Set3".
#' @param position position of gene track be added to("bottom", "top").
#' @param height the relative height of gene track(default: 0.5).
#' @param OrgDb OrgDb for change gene id to gene symbol.
#' @param show_legend show legend or not.
#' @return ggplot object
#' @importFrom IRanges IRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicFeatures genes
#' @importFrom gggenes geom_gene_arrow
#' @importFrom gggenes theme_genes
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 scale_fill_brewer
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_text
#' @importFrom aplot insert_bottom
#' @importFrom aplot insert_top
#' @importFrom rlang sym
#' @importFrom rlang check_installed
#' @export 
addGenePlot <- function(plot, txdb, chr, start_pos, end_pos, xlab = "",
                        gene_symbol_oi = NULL, gene_id_oi = NULL, sample_n_gene = NULL,
                        palette = NULL, position = "bottom", height = 1,
                        OrgDb = NULL, show_legend = FALSE) {
    
    position <- match.arg(position, c("bottom", "top"))
    add_func <- switch(position,
        bottom = aplot::insert_bottom,
        top = aplot::insert_top
    )

    if(sum(c(is.null(gene_symbol_oi), is.null(gene_id_oi), is.null(sample_n_gene))) < 2){
        stop("gene_symbol_oi, gene_id_oi and sample_n_gene are mutually exclusive...")
    }

    # Get genes in the region
    win <- GRanges(seqnames = chr,
                   ranges = IRanges::IRanges(start = start_pos, end = end_pos))
    
    # Get genes related to region
    all_genes <- suppressMessages(GenomicFeatures::genes(txdb))
    gene_df <- data.frame(subsetByOverlaps(x = all_genes, ranges = win, type = "any"))
    gene_df$gene_id <- factor(gene_df$gene_id, levels = unique(gene_df$gene_id))
    
    # Process gene data
    gene_df <- gene_df[,c("seqnames", "gene_id", "start", "end", "strand")]
    colnames(gene_df)[1] <- "chromosome"
    gene_df$forward <- ifelse(gene_df$strand=="+", TRUE, FALSE)

    gene_df$start <- pmax(gene_df$start, start_pos)
    gene_df$end <- pmin(gene_df$end, end_pos)  

    # Convert gene IDs
    if(!is.null(OrgDb)){
        # check clusterProfiler install or not 
        rlang::check_installed('clusterProfiler', reason = 'For coverting gene ids.')

        changeid <- clusterProfiler::bitr(geneID = gene_df$gene_id, 
                                          fromType = "ENTREZID", toType = "SYMBOL",
                                          OrgDb = OrgDb)
        colnames(changeid)[1] <- c("gene_id") 
        gene_df <- merge(gene_df, changeid, all.x = TRUE)
        gene_df <- gene_df[, c(2:ncol(gene_df))]
        colnames(gene_df)[ncol(gene_df)] <- "gene_symbol"
        legend_label <- "gene_symbol"
        
    }else{

        legend_label <- "gene_id"
    }

    # Subset gene_df
    if(!is.null(gene_symbol_oi)){
        if(!"gene_symbol" %in% colnames(gene_df)){
            stop("Please change gene id to gene symbol...")
        }

        gene_df <- gene_df[gene_df$gene_symbol %in% gene_symbol_oi,]

    }

    if(!is.null(gene_id_oi)){
        gene_df <- gene_df[gene_df$gene_id %in% gene_id,]
    }

    if(!is.null(sample_n_gene)){
        n_gene <- min(sample_n_gene, nrow(gene_df))
        index <- sample(rownames(gene_df), n_gene)
        gene_df <- gene_df[index,]
    }


    p <- ggplot(gene_df, aes(xmin = start, xmax = end, y = !!rlang::sym(legend_label), 
                             forward=forward, fill = !!rlang::sym(legend_label))) + 
            geom_gene_arrow() +
            xlim(c(start_pos, end_pos)) +
            theme_genes() +
            labs(fill = gsub("_", " ",legend_label), x = xlab) +
            theme(panel.grid.minor = element_blank(),
                  axis.line = element_line(colour = "black"),
                  panel.border = element_blank(),
                  plot.title = element_text(hjust = 0.5))
    
    if(!is.null(palette)){
        p <- p + scale_fill_brewer(palette = palette)
    }

    if (!show_legend) {
        p <- p + theme(legend.position = "none")
    }

    add_plot <- add_func(plot, p, height = height)
    
    return(add_plot)
}
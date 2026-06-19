#' @func: get_named_list
get_named_list <- function(tab, val_col, name_col){
    tab <- as.data.table(tab)
    value <- tab[, get(val_col)]
    names(value) <- tab[, get(name_col)]
    value
}
#' @func: 
cal_name_coverage <- function(named_list, return=FALSE){
    new <- named_list[names(named_list) != ""]
    if(return){
        new
    }else{
        length(new)/length(named_list)
    }
}

#' Calculate UniProt Mapping Rate
#'
#' This function computes the rate of UniProt identifiers within a named list.
#' It assumes that the names of the list elements are UniProt IDs (or "None" if absent)
#' and calculates:
#'   - The total number of elements (num_gene)
#'   - The number of elements with a valid UniProt ID (num_uniprot)
#'   - The rate of valid UniProt IDs (rate_uniprot = num_uniprot / num_gene)
#'
#' The function uses a regular expression (`regex_name`) to detect valid UniProt IDs.
#'
#' @param named_list A named list whose names are expected to be UniProt IDs.
#'
#' @return A named vector with three elements:
#'   - `num_uniprot`: Number of elements with a valid UniProt ID.
#'   - `num_gene`: Total number of elements.
#'   - `rate_uniprot`: The ratio of valid UniProt IDs to total elements.
#'
#' @examples
#' my_list <- list(a=1, b=2, c=3)
#' names(my_list) <- c("P12345", "None", "Q8N158")
#' cal_uniprot_rate(my_list)
cal_uniprot_rate <- function(named_list) {
  # Define the regular expression for detecting valid UniProt IDs.
  # This regex matches two patterns:
  #   1. IDs starting with O, P, or Q followed by a digit, then three alphanumeric characters, and a digit.
  #   2. IDs starting with A-N or R-Z followed by a digit, then one or two blocks of (an uppercase letter and two alphanumeric characters and a digit).
  regex_name <- "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
  
  uniport_id <- names(named_list) |> sort()
  det_uniprotid <- any(map_lgl(unique(uniport_id), ~ str_detect(..1, regex_name)))
  # If no valid UniProt IDs are detected, log an error.
  if (det_uniprotid == FALSE) {
    flog.error(glue("The list is not named by Uniport ID"))
  } else {
    tot_num <- length(unique(named_list))
    none_gene_num <- names(named_list) == "None"
    annotated_num <- length(unique(named_list[!none_gene_num]))
    annotated_rate <- annotated_num / tot_num
    uni_num <- round(length(uniport_id[uniport_id != "None"]))
    # Compute the UniProt mapping rate.
    uni_rate <- uni_num / tot_num
    # Return a named vector containing the number of valid UniProt IDs, total number, and the rate.
    as.vector(c(uni_num, round(uni_rate, 4),
                tot_num, annotated_num, round(annotated_rate, 4))) |>
      `names<-`(c("num_uniprot", "rate_uniprot", "num_tc_gene", "num_anno_gene", "rate_annotated"))
  }
}

#' Extract transfer gene id to ENTREZID-kegg paird ids
#'
#' @dependency:
#'  - get_named_list()
#' @param gene_groupList A nested list whose first-order list is a named list, with key-value pair as flyid-bettle_id
#' @return A new named & neted list, with key-value pair as ENTREZID-kegg.
convertID_grouplist <- function(gene_groupList){
    gene_convert_entrez <- map(gene_groupList, ~ cal_name_coverage(..1, return=TRUE) |> names()) |>
        map(bitr, fromType = "FLYBASE", toType = "ENTREZID", OrgDb = org.Dm.eg.db)
    gene_entrez_list <- map(gene_convert_entrez, function(x) x$ENTREZID)
    # Convert ENTREZ IDs to KEGG IDs
    gene_convert_kegg <- gene_entrez_list |>
        map(function(x) bitr_kegg(x, fromType = "ncbi-geneid", toType = "kegg", organism="dme"))
    gene_convert_kegg |>
        map(function(x) get_named_list(x, val_col = "ncbi-geneid", name_col = "kegg"))
}
#' @func:go enrichment
#keytype ENTREZID
EnrichGo <- function(geneL, bgenes=NULL, keytype = "ENTREZID",
          ont_key = "ALL", simplifyGO=F,db=org.Dm.eg.db) {
  if(simplifyGO & ont_key %in% c("BP","CC","MF")){
    enrichGO(geneL, db, keyType = keytype, ont = ont_key, universe = bgenes,
           pAdjustMethod = "fdr", qvalueCutoff = 0.2, pvalueCutoff = 0.05,
           minGSSize = 10, maxGSSize = 500, readable = T, pool = FALSE) |>
    clusterProfiler::simplify()
  }else{
    enrichGO(geneL, db, keyType = keytype, ont = ont_key, universe = bgenes,
            pAdjustMethod = "fdr", qvalueCutoff = 0.2, pvalueCutoff = 0.05,
            minGSSize = 10, maxGSSize = 500, readable = T, pool = FALSE)
  }
}
EnrichKEGG <- function(geneL, bgenes=NULL){
  enrichKEGG(gene=geneL, keyType = "kegg", universe = bgenes,
    organism = 'dme', pAdjustMethod  = "fdr", pvalueCutoff = 0.05)
}
#' @func: GSEA
gsea_enri <- function(ranked_gene ){
  gmtfile <- 'Reference/c5.all.v7.1.entrez.gmt'
  #gmtfile <- system.file("extdata", "c5.all.v7.1.entrez.gmt", package="clusterProfiler")
  c5 <- ReadGmt(gmtfile)
  egmt2 = ranked_gene[  !is.na(names(.)) & (names(.) !='NULL')]  |> 
    GSEA(.,TERM2GENE=c5,pvalueCutoff=0.1,verbose=T,pAdjustMethod = "fdr") 
  return(egmt2)
}
#' @func: pathways analysis in one
# fname_gsea_nofix=NULL,
BatchEnrich <- function(rank_gene, bg_go, bg_kegg,
          fname_go_nofix=NULL, ontology="ALL", simplifyGO=T,plot_if=T){
  if(!is.null(fname_go_nofix) ){
    if(simplifyGO && ontology %in% c("BP","CC","MF")){
        ego = EnrichGo(rank_gene, bgenes = bg_go, ont_key = ontology,simplifyGO=T)
    }else{
        ego = EnrichGo(rank_gene, bgenes = bg_go)
    }
    # kegg
    ekegg =  EnrichKEGG(names(rank_gene), bgenes = bg_kegg)
    if(plot_if){
      library(openxlsx)
      sheets = list("GO"=ego,'kegg'=ekegg) 
      write.xlsx(sheets,file =  paste0(fname_go_nofix,".xlsx"))
      if(nrow(ego) == 0 | nrow(ekegg) == 0 ){
        print(paste0("No terms enriched, the gene number is:",length(rank_gene)))
        }else{
        goplot = ego |>
          enrichplot::dotplot(showCategory = 15) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=20))
        keggplot =  ekegg |>  
          barplot(showCategory = 15) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=20))
        # output plot
        suppressMessages(library("gridExtra"))
        pdf(width =16,height=10 ,paste0(fname_go_nofix,'.pdf'))
        par(mar=c(5.1,4.1,4.1,2.1))
        print(grid.arrange(goplot, keggplot, ncol = 2))
        dev.off()
      }
    }else{
      library(openxlsx)
      sheets = list("GO"=ego,'kegg'=ekegg) 
      write.xlsx(sheets,file =  paste0(fname_go_nofix,".xlsx"))
    }
    print('GO KEGG done!')
    }
      # if(!is.null(fname_gsea_nofix)){
  #   print('-------- GSEA ----- ')
  #   emgt2 = gsea_enri(rank_gene)
  #   if(plot_if){
  #     # Make plots.
  #     plot_list = list()
  #     foreach::foreach(i = 1:10) %do% {
  #       p = emgt2 |> gseaplot(.,geneSetID=i,title=.$Description[i],pvalue_table=T)
  #       plot_list[[i]] = p}
  #     pdf(paste0(fname_gsea_nofix,".pdf"),width = 10)
  #     par(mar=c(5.1,4.1,4.1,2.1))
  #     foreach::foreach(i = 1:length(plot_list)) %do% {
  #       print(plot_list[[i]])
  #     }
  #     dev.off()
  #   }else{
  #     emgt2 |> as.data.frame |> dplyr::filter(qvalues<0.3) |>
  #       write.csv( paste0(fname_gsea_nofix,".csv"))
  #   }
  #   print('GSEA done!')}
}

# mapped_ids <- c(D6WBW6="TC000230", D6WBD3="TC000282", D6WBC2="TC000288", D6WB99 = "TC000302")
gprofile_go_kegg <- function(mapped_ids, bg_list, fname_nofix=NULL,
                              organism = "tcastaneum",
                            top=15, maxgsetsize = 500, mingsize=10){
  # none_id <- names(mapped_ids) == "None"
  # mapped_ids <- mapped_ids[! none_id]
  enrichment_results <- gost(query = mapped_ids,
                            organism = organism,
                            sources = c("GO", "KEGG"),
                            custom_bg = bg_list,
                            correction_method = "fdr")
  # Extract the result table
  result_table <- as.data.table(enrichment_results$result)
  nrow_res <- nrow(result_table)
  nrow_res <- nrow(result_table)
  if( nrow_res > 0){
    result_filter <- result_table[(term_size > mingsize) & (term_size < maxgsetsize)]
    # ensure pathway number is not 0 
    result_table <- if(nrow(result_filter) > 0){
      result_filter
      }else{
        result_table
      }
    # Calculate the percentage of query genes in each pathway
    result_table[, percentage := (intersection_size / query_size) * 100]
    # print(head(result_table[order(-percentage), head(.SD, top)]))
    # Create the dot plot
    dotplot <- ggplot(result_table[order(-percentage), head(.SD, top)],
                  aes( x = percentage,   # x-axis: percentage of query genes
                  y = reorder(term_name, percentage),      # y-axis: pathway name (ordered by percentage)
                  size = intersection_size,                # Dot size: number of genes in the pathway
                  color = -log10(p_value),
                  shape = source )) +
      geom_point(alpha = 0.8) +                # Add dots with transparency
      scale_color_gradient(low = "blue", high = "red") + 
      scale_y_discrete(labels=function(x) str_wrap(x, width=30)) +
      # facet_grid(~source, space="free") +
      labs(
        x = "Percentage of Query Genes (%)",   # x-axis label
        y = "",                    
        size = "Intersection Size",            # Legend for dot size
        color = "-log10(p-value)"              # Legend for dot color
      ) +
      theme_minimal() +                        # Use a minimal theme
      theme(
        text = element_text(size = 16)   # Adjust axis title size
      )
    if(is.null(fname_nofix)){
      dotplot
    }else{
      ggsave(paste0(fname_nofix,'.pdf'), dotplot,
            height = ifelse(nrow_res < 10, 5+round(0.6*sqrt(nrow_res+0.6), 1), 9),
            width = ifelse(nrow_res < 10, 6+round(0.6*sqrt(nrow_res+0.6), 1), 9))
      openxlsx::write.xlsx(result_table,file =  paste0(fname_nofix,".xlsx"))
    }
  }else{
    # Open a connection to "test.txt" in append mode
    con <- file(paste0(fname_nofix,'_empty.txt'), open = "a")
    # Write the current date and time as a character string
    writeLines(as.character(Sys.time()), con)
    writeLines(glue("total input gene number is {length(unique(mapped_ids))}"), con)
    # Close the connection
    close(con)
  }
}

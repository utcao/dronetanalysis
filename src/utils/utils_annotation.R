#' Annotate FlyBase gene IDs (FBgn) with gene symbols, names, and cross-references
#'
#' Uses the local org.Dm.eg.db Bioconductor database to map FBgn IDs
#' to human-readable annotations. Handles one-to-many mappings by keeping
#' the first match.
#'
#' @param fbgn_ids Character vector of FlyBase gene IDs (e.g., "FBgn0013275")
#' @param columns Character vector of annotation columns to retrieve.
#'   Available: "SYMBOL", "GENENAME", "ENTREZID", "ENSEMBL", "GO", etc.
#'   Run `columns(org.Dm.eg.db)` for full list.
#' @return data.frame with FLYBASE column plus requested annotation columns.
#'   One row per input FBgn ID (duplicates removed, NAs for unmapped IDs).
#' @examples
#' \dontrun{
#' annotations <- annotate_fbgn(c("FBgn0013275", "FBgn0267001"))
#' }
annotate_fbgn <- function(fbgn_ids,
                           columns = c("SYMBOL", "GENENAME", "ENTREZID", "ENSEMBL")) {
    require(org.Dm.eg.db)

    # Query org.Dm.eg.db
    annot <- AnnotationDbi::select(
        org.Dm.eg.db,
        keys = unique(fbgn_ids),
        keytype = "FLYBASE",
        columns = columns
    )

    # Handle one-to-many mappings: keep first match per FBgn ID
    annot <- annot[!duplicated(annot$FLYBASE), ]

    return(annot)
}

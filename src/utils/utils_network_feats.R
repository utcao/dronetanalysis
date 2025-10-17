suppressPackageStartupMessages({
    library(purrr)
    library(data.table)
    library(glue)
    library(stringr)
    library(ggplot2)
    library(tibble)  # for rownames_to_column / column_to_rownames
    library(futile.logger)
})

#' Render plot for deciding soft-thresholding power for WGCNA
#' 
#' This function generates and saves a plot to help decide the appropriate soft-thresholding power
#' for constructing a scale-free network using WGCNA. The plot includes the scale-free topology fit index
#' and mean connectivity as functions of the soft-thresholding power.
#'
#' @param sft Output from WGCNA::pickSoftThreshold() function containing fit indices
#' @param output_file Full path to output PDF file
#' @param powers Vector of powers used in soft threshold analysis (optional, extracted from sft if not provided)
#' @return NULL (saves plot to file)
#' @export

sft_plot <- function(sft, output_file, powers = c(1:20)){
    # Extract powers from sft if not provided
    if(is.null(powers)) {
        powers <- sft$fitIndices[,1]
    }
    
    # Create PDF with 2x2 layout
    pdf(file = output_file, width = 12, height = 10)
    par(mfrow = c(2,2), mar = c(5, 5, 4, 2) + 0.1) 
    cex1 = 1.0

    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
        xlab="Soft Threshold (power)",
        ylab="Scale Free Topology Model Fit,signed R^2",
        type="n", main = "Scale independence")
        text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
            labels=powers, cex=cex1, col="red")
        abline(h=0.9, col="red")

        # Mean connectivity as a function of the soft-thresholding power
        plot(sft$fitIndices[,1], sft$fitIndices[,5],
            xlab="Soft Threshold (power)",
            ylab="Mean Connectivity", type="n",
            main = "Mean connectivity")
        text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

        plot(sft$fitIndices[,1],sft$fitIndices[,6],
            xlab="Soft Threshold (power)",
            ylab="Median Connectivity",type="n",
            main = "Median connectivity")
        text(sft$fitIndices[,1], sft$fitIndices[,6], labels=powers,cex=cex1,col="red")

        plot(sft$fitIndices[,1],sft$fitIndices[,7],
            xlab="Soft Threshold (power)",
            ylab="Max Connectivity",type="n",
            main = "Max connectivity")
        text(sft$fitIndices[,1], sft$fitIndices[,7], labels=powers,cex=cex1,col="red")
        dev.off()
}
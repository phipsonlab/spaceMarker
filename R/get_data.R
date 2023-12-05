
#' Read 10x Xenium data.

#' @description
#' This function will read the 10x Xenium data and return a list of matrices
#' and vectors. The coordinates will be in microns.
#'
#' @param path The directory containing the matrix.mtx,
#' genes.tsv (or features.tsv), and barcodes.tsv files provided by 10X.
#' @param mtx_name The folder name under path containing the matrix.mtx,
#' genes.tsv (or features.tsv), and barcodes.tsv files.
#' @param trans_name The file name under path containing the transcript
#' information files. Default is "transcript_info.csv.gz"
#' @param cells_name The file name under path containing the cell information
#' files. Default is "cell_info.csv.gz"
#'
#' @importFrom data.table fread
#' @importFrom Seurat Read10X
#'
#' @return A named list with the following components
#' \item{\code{cm}  }{ The count matrix. Each row refers to a gene,
#' and each column refers to a cell}
#' \item{\code{cm_neg}  }{ The count matrix for negative control
#' probe and codeword genes. Each row refers to a probe/codeword gene,
#' and each column refers to a cell.}
#' \item{\code{trans_info}  }{A matrix contains the transcript information
#' for each transcript. All the columns read from the input \code{trans_name}
#' path is kept.}
#' \item{\code{cell_info}  }{A matrix contains the cell information
#' for each cell. The coordinates of each cell is multiplied with the input
#' \code{pixel_size}. All the columns read from the input \code{cells_name}
#' path is kept.}
#' \item{\code{zero_cells}  }{A vector giving the names of cell with zero count
#' on every cell. }
#' \item{\code{probe}  }{A vector giving the names of negative control probe
#' genes.}
#' \item{\code{codeword}  }{A vector giving the names of negative control
#' codeword genes.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' path <- 'path/to/data/directory'
#' data = get_data(path, mtx_name="cell_feature_matrix",
#'                 trans_name="transcripts.csv.gz",
#'                 cells_name="cells.csv.gz")
#'                 }

get_data<-function(path,mtx_name, trans_name="transcript_info.csv.gz",
                    cells_name="cell_info.csv.gz"){

    transcript_info <- as.data.frame(fread(paste(path, trans_name,sep="")))
    cell_info <- as.data.frame(fread(paste(path,cells_name,sep="")))

    data <- Read10X(data.dir = paste(path,mtx_name, sep=""))

    cm <- as.matrix(data$`Gene Expression`)
    r_codeword <- as.matrix(data$`Negative Control Codeword`)

    r_probe <- as.matrix(data$`Negative Control Probe`)
    # merge negative control genes and real genes
    cm_neg <- as.data.frame(rbind(r_probe, r_codeword))
    zero_cells <- colnames(cm)[colSums(cm)==0]

    transcript_info$x_location <- as.numeric(transcript_info$x_location)
    transcript_info$y_location <- as.numeric(transcript_info$y_location)

    cell_info$x_centroid <- as.numeric(cell_info$x_centroid)
    cell_info$y_centroid <- as.numeric(cell_info$y_centroid)

    return (list(cm = cm, cm_neg=cm_neg, zero_cells = zero_cells,
                    trans_info=transcript_info, cell_info=cell_info,
                    probe = row.names(r_probe),
                    codeword=row.names(r_codeword)))

}

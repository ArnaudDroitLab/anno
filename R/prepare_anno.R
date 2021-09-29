#' Download and clean ref and prepare anno
#'
#' The goal of this function is to download the reference fasta file for a
#' specific release of Ensembl or Gencode. The reference is then cleaned. We
#' keep only the transcript id and we remove the transcript version by default.
#' It is also possible to add ERCC92 sequences.
#'
#' After calling this function, a <prefix>.raw_ref.fa.gz file will be
#' downloaded (if not already present) in the current working directory that
#' corresponds to the raw reference file. There will also be a clean version of
#' the transcriptome in the <prefix>.fa.gz format that will be different
#' based on the parameter used to call the function. There will be a
#' <prefix>.info that will contains metadata about the file download and the
#' parameters used. Finally, there will be a <prefix>.csv that contains the
#' annotation formated correctly for the rnaseq packages.
#'
#' The <prefix>.info file contains the following columns:
#'    * prefix: The prefix of the file. Must match filename (i.e.: prefix of
#'              Hs.Gencode38.csv is Hs.Gencode38).
#'    * org: The organism name (i.e.: Homo sapiens)
#'    * db: Database where the annotation was downloaded.
#'    * release: The version of the database.
#'    * rnaseq_pkg_version: The rnaseq package version, if the rnaseq package
#'    was used to download the annotation.
#'    * download_date: The date the annotation was downloaded.
#'    * download_url: The URL that was used to download the annotation.
#'    * md5_raw_ref: md5sum of the raw transcriptome file.
#'    * md5_clean_ref: md5sum of the cleaned transcriptome.
#'    * md5_anno: md5sum of the annotation file.
#'
#' @param prefix The prefix to be used for the files that will be produced.
#' @param org The organism name. Currently accepted:
#'                * Homo sapiens (Ensembl and Gencode)
#'                * Mus musculus (Ensembl and Gencode)
#'                * Macaca mulatta (Ensembl only)
#'                * Rattus norvegicus (Ensembl only)
#' @param db The database to use: Ensembl or Gencode
#' @param release The version of the database to use. Must be greater than 100
#' for Ensembl, 35 for Gencode Homo sapiens and 25 for Gencode Mus musculus.
#' @param removeTxVersion Remove tx version? Default: TRUE.
#' @param ERCC92 Add ERCC92 sequence to reference and to anno? Default: TRUE
#' @param force_download Re-download raw reference if it is already present?
#' Default: FALSE
#'
#' @return Invisibly returns a  \code{list} including the reference
#' transcriptome as a the\code{DNAStringSet} object, the annotation, and the
#' infos (metadata).
#'
#' @examples
#' \dontrun{
#'   prepare_anno("Hs.Ensembl103", org = "Homo sapiens", db = "Ensembl",
#'                release = 103)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings writeXStringSet
#' @importFrom stringr str_extract
#' @importFrom stringr str_replace
#' @importFrom stringr str_detect
#' @importFrom readr write_csv
#' @importFrom tools md5sum
#' @importFrom AnnotationDbi mapIds
#' @importFrom dplyr mutate
#' @importFrom tidyr separate
#' @importFrom utils download.file
#' @importFrom utils packageVersion
#' @importFrom GenomeInfoDb genomeStyles
#' @importFrom methods is
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#' @import org.Mmu.eg.db
#'
#' @export
prepare_anno <- function(prefix, org, db, release,
                         removeTxVersion = TRUE, ERCC92 = FALSE,
                         force_download = FALSE) {

  # Validate params
  stopifnot(is.character(prefix))
  stopifnot(length(prefix) == 1)
  stopifnot(nchar(prefix) > 0)

  stopifnot(org %in% c("Homo sapiens", "Mus musculus", "Macaca mulatta",
                       "Rattus norvegicus"))

  stopifnot(db %in% c("Ensembl", "Gencode"))
  if (db == "Gencode") {
    stopifnot(org %in% c("Homo sapiens", "Mus musculus"))
  }

  stopifnot(is.numeric(release))
  if (db == "Ensembl") {
    stopifnot(release >= 100)
  }
  if (db == "Gencode") {
    if (org == "Homo sapiens") {
      stopifnot(release >= 35)
    }
    if (org == "Mus musculus") {
      stopifnot(release >= 25)
    }
  }

  stopifnot(is.logical(removeTxVersion))
  stopifnot(is.logical(ERCC92))
  stopifnot(is.logical(force_download))

  # Download anno
  raw_ref_infos <- get_filename_and_url(org, db, release)
  raw_ref_filename <- paste0(prefix, ".raw_ref.fa.gz")
  if (!file.exists(raw_ref_filename) | force_download) {
    download.file(raw_ref_infos$url, destfile = raw_ref_filename,
                  method = "curl", extra = "-L")
  }

  # Import and clean
  ref_fasta <- Biostrings::readDNAStringSet(raw_ref_filename)
  if (db == "Gencode") {
    ref_fasta <- ref_fasta[!stringr::str_detect(names(ref_fasta), "PAR_Y")]
  }
  if (db == "Ensembl") {
    if (org == "Macaca mulatta") {
      chromosomes <- names(ref_fasta) %>% str_extract("Mmul_[0-9]*:[^:]*") %>% str_extract("[^:]*$")
      std_chr <- c(1:20, "X", "Y", "MT")
    } else {
      chromosomes <- names(ref_fasta) %>% str_extract("chromosome:[^:]*:[^:]*") %>% str_extract("[^:]*$")
      std_chr <- GenomeInfoDb::genomeStyles(org)$NCBI
    }
    ref_fasta <- ref_fasta[chromosomes %in% std_chr]
  }
  ref_fasta <- ref_fasta[width(ref_fasta) != 0]
  anno <- extract_anno(ref_fasta, org, db, removeTxVersion)
  if (removeTxVersion) {
    names(ref_fasta) <- stringr::str_extract(names(ref_fasta),
                                             "^ENS[^\\.]*")
  } else {
    names(ref_fasta) <- stringr::str_extract(names(ref_fasta),
                                             "^ENS[^ ]*")
  }

  # Add ERCC92?
  if (ERCC92) {
    anno <- add_ercc92_anno(anno)
    ref_fasta <- add_ercc92_fasta(ref_fasta)
  }

  # Save results
  output_ref_fasta <- paste0(prefix, ".fa.gz")
  Biostrings::writeXStringSet(ref_fasta, output_ref_fasta, compress = TRUE)

  output_anno <- paste0(prefix, ".csv")
  readr::write_csv(anno, output_anno)

  # Infos
  output_info <- paste0(prefix, ".info")
  info <- data.frame(prefix = prefix,
                     org = org,
                     db = db,
                     release = release,
                     ERCC92 = ERCC92,
                     rnaseq_pkg_version = packageVersion("rnaseq"),
                     download_date = as.Date(Sys.Date(), format = "%B %d %Y"),
                     download_url = raw_ref_infos$url,
                     md5_raw_ref = tools::md5sum(raw_ref_filename),
                     md5_clean_ref = tools::md5sum(output_ref_fasta),
                     md5_anno = tools::md5sum(output_anno))
  readr::write_csv(info, output_info)

  list(ref_fasta = ref_fasta, anno = anno, info = info)
}

get_filename_and_url <- function(org, db, release) {
  stopifnot(org %in% c("Homo sapiens", "Mus musculus", "Macaca mulatta",
                       "Rattus norvegicus"))

  stopifnot(db %in% c("Ensembl", "Gencode"))
  if (db == "Gencode") {
    stopifnot(org %in% c("Homo sapiens", "Mus musculus"))
  }

  stopifnot(is.numeric(release))
  if (db == "Ensembl") {
    stopifnot(release >= 100)
  }
  if (db == "Gencode") {
    if (org == "Homo sapiens") {
      stopifnot(release >= 35)
    }
    if (org == "Mus musculus") {
      stopifnot(release >= 25)
    }
  }

  filename <- ""
  url <- ""
  base_url_ensembl <- "http://ftp.ensembl.org/pub/release-"
  base_url_gencode <- "http://ftp.ebi.ac.uk/pub/databases/gencode"

  if (org == "Homo sapiens") {
    if (db == "Ensembl") {
      filename <- "Homo_sapiens.GRCh38.cdna.all.fa.gz"
      url <- paste0(base_url_ensembl, release,
                    "/fasta/homo_sapiens/cdna/", filename)
    } else {
      filename <- paste0("gencode.v", release, ".transcripts.fa.gz")
      url <- paste0(base_url_gencode, "/Gencode_human/release_",
                    release, "/", filename)
    }
  }
  if (org == "Mus musculus") {
    if (db == "Ensembl") {
      if (release < 103) {
        filename <- "Mus_musculus.GRCm38.cdna.all.fa.gz"
      } else {
        filename <- "Mus_musculus.GRCm39.cdna.all.fa.gz"
      }
      url <- paste0(base_url_ensembl, release,
                    "/fasta/mus_musculus/cdna/", filename)
    } else {
      filename <- paste0("gencode.vM", release, ".transcripts.fa.gz")
      url <- paste0(base_url_gencode, "/Gencode_mouse/release_M",
                    release, "/", filename)
    }
  }
  if (org == "Macaca mulatta") {
    filename <- "Macaca_mulatta.Mmul_10.cdna.all.fa.gz"
    url <- paste0(base_url_ensembl, release, "/fasta/macaca_mulatta/cdna/",
                  filename)
  }
  if (org == "Rattus norvegicus") {
    filename <- "Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz"
    url <- paste0(base_url_ensembl, release, "/fasta/rattus_norvegicus/cdna/",
                  filename)
  }
  list(filename = filename, url = url)
}

extract_anno <- function(raw_ref, org, db, removeTxVersion) {
  stopifnot(is(raw_ref, "DNAStringSet"))
  stopifnot(db %in% c("Gencode", "Ensembl"))
  if (db == "Gencode") {
    col_names <- c("id",
                   "ensembl_gene",
                   "havana_gene",
                   "havana_transcript",
                   "transcript_name",
                   "symbol",
                   "length",
                   "transcript_type",
                   "filler")

    full_name <- NULL
    df <- data.frame(full_name = names(raw_ref)) %>%
      tidyr::separate(full_name, into = col_names, sep = "\\|")
    if (removeTxVersion) {
      df <- dplyr::mutate(df, id = stringr::str_replace(id, "\\..*$", ""),
                          ensembl_gene = stringr::str_replace(ensembl_gene, "\\..*$", ""))
    }
    if (org == "Homo sapiens") org_db <- org.Hs.eg.db
    if (org == "Mus musculus") org_db <- org.Mm.eg.db
    if (org == "Rattus norvegicus") org_db <- org.Rn.eg.db
    if (org == "Macaca mulatta") org_db <- org.Mmu.eg.db
    df$entrez_id <- AnnotationDbi::mapIds(org_db,
                                          keys = df$ensembl_gene,
                                          keytype = "ENSEMBL",
                                          column = "ENTREZID")
    dplyr::select(df, id, ensembl_gene, symbol, entrez_id, transcript_type)
  } else if (db == "Ensembl") {
    id <- stringr::str_extract(names(raw_ref), "^ENS[^ ]*")
    ensembl_gene <- stringr::str_extract(names(raw_ref), "gene:[^ ]*") %>%
      stringr::str_replace("gene:", "")
    if (removeTxVersion) {
      id <- stringr::str_extract(id, "^ENS[^\\.]*")
      ensembl_gene <- stringr::str_extract(ensembl_gene, "^ENS[^\\.]*")
    }
    symbol <- stringr::str_extract(names(raw_ref), "gene_symbol:[^ ]*") %>%
      stringr::str_replace("gene_symbol:", "")
    transcript_type <- stringr::str_extract(names(raw_ref),
                                            "transcript_biotype:[^ ]*") %>%
      stringr::str_replace("transcript_biotype:", "")
    if (org == "Homo sapiens") org_db <- org.Hs.eg.db
    if (org == "Mus musculus") org_db <- org.Mm.eg.db
    if (org == "Rattus norvegicus") org_db <- org.Rn.eg.db
    if (org == "Macaca mulatta") org_db <- org.Mmu.eg.db
    entrez_id <- AnnotationDbi::mapIds(org_db,
                                       keys = ensembl_gene,
                                       keytype = "ENSEMBL",
                                       column = "ENTREZID")
    data.frame(id = id,
               ensembl_gene = ensembl_gene,
               symbol = symbol,
               entrez_id = entrez_id,
               transcript_type = transcript_type)
  } else {
    stop("Invalid db value")
  }
}

add_ercc92_anno <- function(anno) {
  ercc92_filename <- system.file("extdata/ERCC92.fa", package = "rnaseq")
  ercc92_fasta <- Biostrings::readDNAStringSet(ercc92_filename)

  ercc92_anno <- data.frame(id = names(ercc92_fasta),
                            ensembl_gene = names(ercc92_fasta),
                            symbol = names(ercc92_fasta),
                            entrez_id = NA,
                            transcript_type = "Spike-in")
  rbind(anno, ercc92_anno)
}

add_ercc92_fasta <- function(ref_fasta) {
  ercc92_filename <- system.file("extdata/ERCC92.fa", package = "rnaseq")
  ercc92_fasta <- Biostrings::readDNAStringSet(ercc92_filename)
  c(ref_fasta, ercc92_fasta)
}

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

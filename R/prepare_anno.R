#' Download and clean ref and prepare anno
#'
#' The goal of this function is to download the reference fasta file for a
#' specific release of Ensembl or Gencode. The reference is then cleaned. We
#' keep only the transcript id and we remove the transcript version by default.
#' It is also possible to add ERCC92 sequences. Reference files without the
#' alternative chromosomes and only with the protein coding are also generated.
#'
#' #' After calling this function, a <prefix>.raw_ref.fa.gz file will be
#' downloaded (if not already present) to the current working directory that
#' corresponds to the raw reference file. There will also be a clean version,
#' without alternative chromosomes in the format <prefix>.no_alt_chr.fa.gz.
#' A <prefix>.protein_coding.fa.gz file is also generated, containing only the
#' protein_coding genes.
#' Finally, for all 3 fa.gz files, a <prefix>.csv file is created. The csv
#' file contains the annotation formated correctly for the rnaseq packages.
#' Finally, a <prefix>.info file is created. This file contains metadata about
#' every file and the parameters used.
#'
#' The <prefix>.info file contains the following columns:
#'    * prefix: The prefix of the file. Must match filename (i.e.: prefix of
#'              Hs.Gencode38.csv is Hs.Gencode38).
#'    * org: The organism name (i.e.: Homo sapiens)
#'    * db: Database where the annotation was downloaded.
#'    * release: The version of the database.
#'    * ERCC92: The value of the ERCC92 argument.
#'    * anno_pkg_version: The anno package version.
#'    * download_date: The date the annotation was downloaded.
#'    * download_url: The URL that was used to download the annotation.
#'    * A md5sum for every file generated, one column per file.
#'
#' @param org The organism name. Currently supported:
#'                * Homo sapiens (Ensembl and Gencode)
#'                * Mus musculus (Ensembl and Gencode)
#'                * Macaca mulatta (Ensembl only)
#'                * Rattus norvegicus (Ensembl only)
#'                * Bos taurus (Ensembl only)
#'                * Mesocricetus auratus (Ensembl only, annotation unsufficient to do the no_alt_chr)
#' This can also be any organism in the ensembl database, but they are not
#' officially supported and can result in bugs. If a fasta is provided, orgnaism
#' can be anything (containing only letters, numbers, spaces, _, -).
#' Default: NA
#' @param db The database to use: Ensembl or Gencode. Default: "Ensembl"
#' @param release The version of the database to use. Must be greater than 100
#' for Ensembl, 35 for Gencode Homo sapiens and 25 for Gencode Mus musculus.
#' Default: NA
#' @param annotation_version The version of the annotation database to use for
#' Ensembl database. Will be ignored if using Gencode database. Must be lower
#' or equal to Ensembl release.
#' Default: NA
#' @param fasta Path to the fasta file, if it has already been downloaded. For
#' now only Ensembl fasta are supported.
#' Default: NA
#' @param ERCC92 Add ERCC92 sequence to reference and to anno?
#' Default: FALSE
#' @param standard_chr Unsupported as of now
#' @param force_download Re-download raw reference if it is already present?
#' Default: FALSE
#' @param gtf Download the annotation corresponding to the fasta in gtf format?
#' Default: FALSE
#' @param outdir Directory in which to save the files.
#' Default: "."
#'
#' @return Returns a \code{list} including every information in the
#' <prefix>.info file.
#'
#' @examples
#' \dontrun{
#'   prepare_anno(org = "Homo sapiens", db = "Ensembl",
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
#' @importFrom XML getHTMLLinks
#' @importFrom AnnotationHub query
#' @importFrom AnnotationHub AnnotationHub
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#'
#' @export

prepare_anno <- function(org = NA, db = "Ensembl", release = NA, annotation_version = NA, fasta = NA, ERCC92 = FALSE, standard_chr = NA,
                         force_download = FALSE, gtf = FALSE, outdir = ".") {

  # Validate params
  new_params <- validate_params(org, db, release, annotation_version, fasta, ERCC92, standard_chr, force_download, gtf, outdir)
  release = new_params[["release"]]
  annotation_version = new_params[["annotation_version"]]
  valid_org = new_params[["valid_org"]]
  supported_org = new_params[["supported_org"]]
  org = new_params[["org"]]

  # Download raw ref file
  prefix <- get_prefix(org, db, release, annotation_version, outdir)

  if (is.na(fasta)) {
    raw_ref_infos <- get_filename_and_url(org, db, release, valid_org)
    raw_ref_filename <- paste0(prefix, ".raw_ref.fa.gz")
    if (!file.exists(raw_ref_filename) | force_download) {
      download.file(raw_ref_infos$url, destfile = raw_ref_filename,
                    method = "curl", extra = "-L")
    }
  } else {raw_ref_filename = fasta}

  gtf_filename <- paste0(prefix, ".gtf.gz")
  if (gtf & (!file.exists(gtf_filename) | force_download)) {
    gtf_url <- get_gtf_link(org, db, release)
    download.file(gtf_url, destfile = gtf_filename,
                  method = "curl", extra = "-L")
  }
  raw_ref_fasta <- Biostrings::readDNAStringSet(raw_ref_filename)
  if (any(BiocGenerics::width(raw_ref_fasta) == 0)) {
    stop("0 length transcripts in fasta file")
  }

  raw_ref_anno <- extract_anno_bis(raw_ref_fasta, org, db, annotation_version, supported_org)

  # Get cleaned ref
  anno <- raw_ref_anno
  ref_fasta <- raw_ref_fasta
  clean_names <- clean_ref(ref_fasta, db)
  names(ref_fasta) <- clean_names

  if (ERCC92) {
    ref_fasta <- add_ercc92_fasta(ref_fasta)
    anno <- add_ercc92_anno(anno)
  }

  save_anno_resuts(ERCC92, prefix, "cleaned_ref", ref_fasta, anno)

  # Get ref without alternative chromosomes
  if (supported_org) {
    ref_fasta <- no_alt_chr(raw_ref_fasta, db, org)
    names(ref_fasta) <- clean_ref(ref_fasta, db)
    anno <- anno %>% dplyr::filter(id %in% names(ref_fasta))

    if (ERCC92) {
      # anno <- add_ercc92_anno(anno)
      ref_fasta <- add_ercc92_fasta(ref_fasta)
    }

    save_anno_resuts(ERCC92, prefix, "no_alt_chr", ref_fasta, anno)
  }


  # Get ref with only protein_coding genes
  if (!supported_org) {
    ref_fasta <- raw_ref_fasta
    names(ref_fasta) <- clean_names
  }
  ref_fasta <- ref_fasta[BiocGenerics::width(ref_fasta) != 0]

  i <- anno$transcript_type == "protein_coding"
  anno <- anno[i,]
  ref_fasta <- ref_fasta[i]

  if (ERCC92) {
    anno <- add_ercc92_anno(anno)
    ref_fasta <- add_ercc92_fasta(ref_fasta)
  }
  save_anno_resuts(ERCC92, prefix, "protein_coding", ref_fasta, anno)

  info <- save_info(prefix, org, db, release, ERCC92, if (is.na(fasta)) raw_ref_infos$url else "NA", gtf, annotation_version, raw_ref = raw_ref_filename)
  info
}

validate_params <- function(org, db, release, annotation_version, fasta, ERCC92, standard_chr,
                            force_download, gtf, outdir) {

  if (!db %in% c("Ensembl", "Gencode")) {
    stop("db must be either Ensembl or Gencode")
  }

  is_supported = FALSE
  is_valid_ensembl_org = FALSE
  if (!is.character(org) & !is.na(org)) {
    stop("org must be NA or a supported organism")
  }
  org = tolower(org)
  if (tolower(org) %in% c("homo sapiens", "human", "grch38", "hg38", "hs", "h. sapiens", "homo_sapiens")) {
    org = "homo_sapiens"
    is_supported = TRUE
    is_valid_ensembl_org = TRUE
  } else if (tolower(org) %in% c("mus musculus", "mouse", "grcm38", "mm10", "mm", "grcm39", "mm39", "m. musculus", "mus_musculus")) {
    org = "mus_musculus"
    is_supported = TRUE
    is_valid_ensembl_org = TRUE
  } else if (tolower(org) %in% c("macaca mulatta", "rhesus monkey", "mmul_10", "rhemac10", "mmu", "m. mulatta", "macaca_mulatta")) {
    org = "macaca_mulatta"
    is_supported = TRUE
    is_valid_ensembl_org = TRUE
  } else if (tolower(org) %in% c("rattus norvegicus", "rat", "rrnor_6.0", "rn6", "rn", "r. norvegicus", "rattus_norvegicus")) {
    org = "rattus_norvegicus"
    is_supported = TRUE
    is_valid_ensembl_org = TRUE
  } else if (tolower(org) %in% c("bos taurus", "cow", "bos_taurus_umd_3.1.1", "bt", "b. taurus", "bos_taurus")) {
    org = "bos_taurus"
    is_supported = TRUE
    is_valid_ensembl_org = TRUE
  } else if (tolower(org) %in% c("mesocricetus auratus", "hamster", "mesaur1.0", "ma", "m. auratus", "mesocricetus_auratus")) {
    org = "mesocricetus_auratus"
    is_supported = TRUE
    is_valid_ensembl_org = TRUE
  } else if (!is.na(org)) {
    if (db == "Ensembl" & is.na(fasta)) {
      print(paste0(org, " is an unsupported organism, checking if it is a valid Ensembl organism"))
      is_valid_ensembl_org <- fetch_ensembl_organism(org, db, release)
      if (!is_valid_ensembl_org) {
        stop(paste0("fasta is not provided, and organism ", org, " is not a valid ensembl organism"))
      } else {
        print(paste0("fasta is not provided, but organism ", org, " is a valid ensembl organism, genome will be downloaded but ENTREZID will be skipped"))
      }
    } else if (db == "Ensembl" & !is.na(fasta)) {
      print(paste0(org, " is an unsupported organism, fasta has been provided and genome will not be downloaded"))
    }
  }
  if (is.na(org)) {
    print("organism is NA, a fasta must be provided and ENTREZID will be skipped")
  } else {
    org = gsub("_", " ", stringr::str_to_sentence(org))
  }

  if (db == "Gencode" & !org %in% c("Homo sapiens", "Mus musculus")) {
    stop("Gencode database only accepts Homo sapiens and Mus musculus organisms")
  }

  if (!is_supported & !is.na(org) & db == "Ensembl" & is.na(fasta)) {
    is_valid_ensembl_org <- fetch_ensembl_organism(org, db, release)
    if (!is_valid_ensembl_org) {
      stop(paste0("fasta is not provided, but organism ", org, " is not a valid ensembl organism"))
    }
  }

  if (is.na(release) & is_supported) {
    release <- fetch_latest_release(org, db)
    print(paste("using latest release :", db, release), sep = " ")
  } else if (is.na(release) & is.na(org)) {
    print("No organism specified, skipping EntrezID search")
  } else if (is_valid_ensembl_org) {
    if (!is.numeric(release) & db == "Ensembl") {
      print("Organism not officially supported, trying to get the latest Ensembl release for this organism")
      release <- fetch_latest_release(org, db)
      print(paste("using latest release :", db, release), sep = " ")

      # stop("release must be a number")
    }
    if (db == "Ensembl" & release < 100) {
      stop("release must be >= 100 for Ensembl")
    }
    if (db == "Gencode") {
      if (org == "Homo sapiens" & release < 35) {
        stop("release must be >= 35 for Homo sapiens in Gencode")
      }
      if (org == "Mus musculus" & release < 25) {
        stop("release must be >= 25 for Mus musculus in Gencode")
      }
    }
  }

  if (db == "Ensembl" & is_supported) {
    latest_annot <- fetch_latest_annotation(org)
    if (is.na(annotation_version)) {
      annotation_version = latest_annot
      print(paste0("Using latest annotation version : ", annotation_version))
    } else if (!is.numeric(annotation_version) & !is.na(annotation_version)) {
      stop("annotation_version must be a number")
    } else if (annotation_version > latest_annot & !is.na(annotation_version)) {
      annotation_version = latest_annot
      print(paste0("Annotation version not supported yet, using latest annotation version : "), annotation_version)
    }
  }

  if (!is.character(fasta) & !is.na(fasta)) {
    stop("fasta must be NA or a string")
  } else if (!is.na(fasta)) {
    if (sum(unlist(lapply(c("fa", "fasta"), function(x) x %in% stringr::str_split_1(fasta, "\\.")))) == 0) {
      stop("fasta must end with .fasta, .fa, .fa.gz or .fasta.gz")
    }
  }

  if (is.na(fasta) & is.na(org)) {
    stop("Must specify at least one of org and fasta")
  }

  if (!is.logical(ERCC92) | is.na(ERCC92)) {
    stop("ERCC92 must be either TRUE or FALSE")
  }
  if (!is.logical(force_download) | is.na(force_download)) {
    stop("force_download must be either TRUE or FALSE")
  }
  if (!is.logical(gtf) | is.na(gtf)) {
    stop("gtf must be either TRUE or FALSE")
  }
  if (is.character(outdir)) {
    if (!dir.exists(outdir)) {
      stop("outdir must be a valid directory")
    }
  } else {
    stop("outdir must be a string and a valid directory")
  }
  return(list(org = org, db = db, release = release, annotation_version = annotation_version,
              fasta = fasta, ERCC92 = ERCC92, standard_chr = standard_chr,
              force_download = force_download, gtf = gtf, outdir = outdir, supported_org = is_supported, valid_org = is_valid_ensembl_org))
}

clean_ref <- function(ref_fasta, db) {
  if (db == "Gencode") {
    clean_names <- stringr::str_extract(names(ref_fasta), "^ENS[^\\|]*")
    clean_names <- stringr::str_remove(clean_names, "\\.[0-9]*")
  } else {
    clean_names <- stringr::str_extract(names(ref_fasta), "^ENS[^\\.]*")
  }
  # Stop if duplicated transcripts ID
  if (length(clean_names) != length(unique(clean_names))){
    stop("Duplicated transcripts in fasta")
  }
  return(clean_names)
}

no_alt_chr <- function(ref_fasta, db, org) {
  if (db == "Gencode") {
    ref_fasta <- ref_fasta[!stringr::str_detect(names(ref_fasta), "PAR_Y")]
  }
  if (db == "Ensembl" & org != "Mesocricetus auratus") {
    if (org == "Macaca mulatta") {
      chromosomes <- names(ref_fasta) %>% str_extract("Mmul_[0-9]*:[^:]*") %>% str_extract("[^:]*$")
      std_chr <- c(1:20, "X", "Y", "MT")
    }
    else if (org == "Bos taurus"){
      chromosomes <- names(ref_fasta) %>% str_extract("ARS-UCD1.[0-9]*:[^:]*") %>% str_extract("[^:]*$")
      std_chr <- c(1:29, "X", "Y", "MT")
    }
    else if (org == "Rattus norvegicus" & db == "Ensembl" & release >= 105) {
      chromosomes <- names(ref_fasta) %>% str_extract("primary_assembly:[^:]*:[^:]*") %>% str_extract("[^:]*$")
      std_chr <- GenomeInfoDb::genomeStyles(org)$NCBI
    }
    else {
      chromosomes <- names(ref_fasta) %>% str_extract("chromosome:[^:]*:[^:]*") %>% str_extract("[^:]*$")
      std_chr <- GenomeInfoDb::genomeStyles(org)$NCBI
    }
    ref_fasta <- ref_fasta[chromosomes %in% std_chr]
  }
  return(ref_fasta)
}

get_prefix <- function(org, db, release, annotation_version, outdir){
  if (is.na(org)) org = "unknown"
  if (is.na(release)) release = "unknown"
  if (is.na(annotation_version)) annotation_version = "unknown"
  organism = str_replace(org, " ", "_")
  prefix <- paste0(organism, ".", db, release)
  if (db == "Ensembl") {prefix <- paste(prefix, annotation_version, sep = "_")}
  paste(outdir, prefix, sep = "/")
}

get_filename_and_url <- function(org, db, release, is_valid) {
  if (!is_valid) {stop("Not a valid organism")}
  # stopifnot(org %in% c("Homo sapiens", "Mus musculus", "Macaca mulatta",
  #                      "Rattus norvegicus", "Bos taurus", "Mesocricetus auratus"))

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
    if (release >= 105) {
      filename <- "Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz"
    } else {
      filename <- "Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz"
    }
    url <- paste0(base_url_ensembl, release, "/fasta/rattus_norvegicus/cdna/",
                  filename)
  }
  if (org == "Bos taurus") {
    filename <- "Bos_taurus.ARS-UCD1.2.cdna.all.fa.gz"
    url <- paste0(base_url_ensembl, release, "/fasta/bos_taurus/cdna/",
                  filename)
  }
  if (org == "Mesocricetus auratus") {
    filename <- "Mesocricetus_auratus.MesAur1.0.cdna.all.fa.gz"
    url <- paste0(base_url_ensembl, release, "/fasta/mesocricetus_auratus/cdna/",
                  filename)
  } else {
    # general download with grep "cdna.all.fa.gz"
    low_org <- gsub(" ", "_", stringr::str_to_lower(org))
    temp <- XML::getHTMLLinks(paste0(base_url_ensembl, release, "/fasta/", low_org, "/cdna/"))
    filename <- grep("cdna.all.fa.gz", temp, value = TRUE)[1]
    url <- paste0(base_url_ensembl, release, "/fasta/", low_org, "/cdna/", filename)
  }
  list(filename = filename, url = url)
}

extract_anno_bis <- function(raw_ref, org, db, annotation_version, valid_org) {
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
      tidyr::separate(full_name, into = col_names, sep = "\\|") %>%
      dplyr::mutate(id = stringr::str_replace(id, "\\.[0-9]*", ""),
                    ensembl_gene = stringr::str_replace(ensembl_gene, "\\..*$", ""))
    if (org == "Homo sapiens") org_db <- org.Hs.eg.db
    if (org == "Mus musculus") org_db <- org.Mm.eg.db
    df$entrez_id <- AnnotationDbi::mapIds(org_db,
                                          keys = df$ensembl_gene,
                                          keytype = "ENSEMBL",
                                          column = "ENTREZID")
    dplyr::select(df, id, ensembl_gene, symbol, entrez_id, transcript_type)
  } else if (db == "Ensembl") {
    df <- parse_fasta_names(raw_ref)
    df$ids <- stringr::str_replace(df$ids, "\\.[0-9]*", "")
    if (valid_org) {
      hub <- AnnotationHub::AnnotationHub()
      all_db <- AnnotationHub::query(hub, c(org, "EnsDb", annotation_version))
      if (length(all_db) == 0){
        stop("Unsupported annotation version")
      }
      db_name <- names(all_db@.db_uid[length(all_db)])
      ensdb <- hub[[db_name]]
    }


    annot = data.frame(id = df$ids,
                       ensembl_gene = df$gene,
                       symbol = df$gene_symbol,
                       entrez_id = if (valid_org) mapIds(ensdb, keys = df$ids, keytype = "TXID", column = "ENTREZID") else NA,
                       transcript_type = df$transcript_biotype)
  } else {
    stop("Invalid db value")
  }

}

extract_anno <- function(raw_ref, org, db, annotation_version) {
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
      tidyr::separate(full_name, into = col_names, sep = "\\|") %>%
      dplyr::mutate(id = stringr::str_replace(id, "\\.[0-9]*", ""),
                    ensembl_gene = stringr::str_replace(ensembl_gene, "\\..*$", ""))
    if (org == "Homo sapiens") org_db <- org.Hs.eg.db
    if (org == "Mus musculus") org_db <- org.Mm.eg.db
    df$entrez_id <- AnnotationDbi::mapIds(org_db,
                                          keys = df$ensembl_gene,
                                          keytype = "ENSEMBL",
                                          column = "ENTREZID")
    dplyr::select(df, id, ensembl_gene, symbol, entrez_id, transcript_type)
  } else if (db == "Ensembl") {
    id <- stringr::str_extract(names(raw_ref), "^ENS[^ ]*") %>%
      stringr::str_extract("^ENS[^\\.]*")
    transcript_type <- stringr::str_extract(names(raw_ref), "transcript_biotype:[^ ]*") %>%
      stringr::str_replace("transcript_biotype:", "")
    ensembl_gene <- stringr::str_extract(names(raw_ref), "gene:[^ ]*") %>%
      stringr::str_replace("gene:", "") %>%
      stringr::str_extract("^ENS[^\\.]*")
    symbol <- stringr::str_extract(names(raw_ref), "gene_symbol:[^ ]*") %>%
      stringr::str_replace("gene_symbol:", "")

    hub <- AnnotationHub::AnnotationHub()
    all_db <- AnnotationHub::query(hub, c(org, "EnsDb", annotation_version))
    if (length(all_db) == 0){
      stop("Unsupported annotation version")
    }
    db_name <- names(all_db@.db_uid[length(all_db)])
    ensdb <- hub[[db_name]]

    annot = data.frame(id = id,
                       ensembl_gene = ensembl_gene,
                       symbol = symbol,
                       entrez_id = mapIds(ensdb, keys = id, keytype = "TXID", column = "ENTREZID"),
                       transcript_type = transcript_type)
  } else {
    stop("Invalid db value")
  }
}

add_ercc92_anno <- function(anno) {
  ercc92_filename <- system.file("extdata/ERCC92.fa", package = "anno")
  ercc92_fasta <- Biostrings::readDNAStringSet(ercc92_filename)

  ercc92_anno <- data.frame(id = names(ercc92_fasta),
                            ensembl_gene = names(ercc92_fasta),
                            symbol = names(ercc92_fasta),
                            entrez_id = NA,
                            transcript_type = "Spike-in")
  rbind(anno, ercc92_anno)
}

add_ercc92_fasta <- function(ref_fasta) {
  ercc92_filename <- system.file("extdata/ERCC92.fa", package = "anno")
  ercc92_fasta <- Biostrings::readDNAStringSet(ercc92_filename)
  c(ref_fasta, ercc92_fasta)
}

save_anno_resuts <- function(ERCC92, prefix, ref_type, ref_fasta, anno){
  if (ERCC92){
    output_ref_fasta <- paste(prefix, "ERCC92", ref_type, "fa.gz", sep = ".")
    output_anno <- paste(prefix, "ERCC92", ref_type, "csv", sep = ".")
  }
  else{
    output_ref_fasta <- paste(prefix, ref_type, "fa.gz", sep = ".")
    output_anno <- paste(prefix, ref_type, "csv", sep = ".")
  }

  if (all(anno$id == names(ref_fasta)) & nrow(anno) == length(ref_fasta)) {
    Biostrings::writeXStringSet(ref_fasta, output_ref_fasta, compress = TRUE)
    csv <- readr::write_csv(anno, output_anno)
  } else {
    stop("Fasta and annotation IDs do not match")
  }
}

save_info <- function(prefix, org, db, release, ERCC92, url, gtf, annotation_version, raw_ref){
  if (ERCC92){
    prefix <- paste(prefix, "ERCC92", sep = ".")
  }

  info <- data.frame(prefix = prefix,
                     org = org,
                     db = db,
                     release = release,
                     annotation_version = ifelse(db == "Ensembl", annotation_version, NA),
                     ERCC92 = ERCC92,
                     anno_pkg_version = packageVersion("anno"),
                     download_date = as.Date(Sys.Date(), format = "%B %d %Y"),
                     download_url = url,
                     md5_raw_ref = tools::md5sum(raw_ref),
                     md5_cleaned_raw_ref = tools::md5sum(paste0(prefix, ".cleaned_ref.fa.gz")),
                     md5_no_alt_chr_ref = tools::md5sum(paste0(prefix, ".no_alt_chr.fa.gz")),
                     md5_protein_coding_ref = tools::md5sum(paste0(prefix, ".protein_coding.fa.gz")),
                     md5_cleaned_raw_anno = tools::md5sum(paste0(prefix, ".cleaned_ref.csv")),
                     md5_no_alt_chr_anno = tools::md5sum(paste0(prefix, ".no_alt_chr.csv")),
                     md5_protein_coding_anno = tools::md5sum(paste0(prefix, ".protein_coding.csv")))
  if (gtf) {
    info$md5_gtf <- tools::md5sum(paste0(prefix,".gtf.gz"))
  }
  readr::write_csv(info, paste0(prefix, ".info"))
  info
}

# TODO créer une fonction pour télécharger le gtf
get_gtf_link <- function(org, db, release){

  stopifnot(org %in% c("Homo sapiens", "Mus musculus", "Macaca mulatta",
                       "Rattus norvegicus", "Bos taurus"))

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
  url_ensembl <- paste0("http://ftp.ensembl.org/pub/release-", release, "/gtf/",
                        gsub(" ", "_", stringr::str_to_lower(org)), "/")
  url_gencode <- "http://ftp.ebi.ac.uk/pub/databases/gencode"

  if (org == "Homo sapiens") {
    if (db == "Ensembl") {
      filename <- paste0("Homo_sapiens.GRCh38.", release, ".gtf.gz")
      url <- paste0(url_ensembl, filename)
    } else {
      filename <- paste0("gencode.v", release, ".annotation.gtf.gz")
      url <- paste0(url_gencode, "/Gencode_human/release_",
                    release, "/", filename)
    }
  }
  if (org == "Mus musculus") {
    if (db == "Ensembl") {
      if (release < 103) {
        filename <- paste0("Mus_musculus.GRCm38.", release, ".gtf.gz")
      } else {
        filename <- paste0("Mus_musculus.GRCm39.", release, ".gtf.gz")
      }
      url <- paste0(url_ensembl, filename)
    } else {
      filename <- paste0("gencode.vM", release, ".annotation.gtf.gz")
      url <- paste0(url_gencode, "/Gencode_mouse/release_M",
                    release, "/", filename)
    }
  }
  if (org == "Macaca mulatta") {
    filename <- paste0("Macaca_mulatta.Mmul_10.", release, ".gtf.gz")
    url <- paste0(url_ensembl, filename)
  }
  if (org == "Rattus norvegicus") {
    if (release >= 105) {
      filename <-  paste0("Rattus_norvegicus.mRatBN7.2.", release, ".gtf.gz")
    } else {
      filename <-  paste0("Rattus_norvegicus.Rnor_6.0.", release, ".gtf.gz")
    }
    url <- paste0(url_ensembl, filename)
  }
  if (org == "Bos taurus") {
    filename <- paste0("Bos_taurus.ARS-UCD1.2.", release, ".gtf.gz")
    url <- paste0(url_ensembl, filename)
  }
  if (org == "Mesocricetus auratus") {
    filename <- paste0("Mesocricetus_auratus.MesAur1.0.", release, ".gtf.gz")
    url <- paste0(url_ensembl, filename)
  }
  url
}

fetch_latest_release <- function(org, db){
  if (db == "Gencode") {
    if (org == "Homo sapiens") {
      versions = XML::getHTMLLinks("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
                                   xpQuery = "//a/@href[contains(., 'release_')]")
      versions = as.integer(gsub("[^0-9]*", "", versions))
    } else {
      versions = XML::getHTMLLinks("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/",
                                   xpQuery = "//a/@href[contains(., 'release_M')]")
      versions = as.integer(gsub("[^0-9]*", "", versions))
    }
    version = max(versions)
  } else {
    versions = XML::getHTMLLinks("http://ftp.ensembl.org/pub/", xpQuery = "//a/@href[contains(., 'release-')]")
    versions = as.integer(gsub("[^0-9]*", "", versions))
    version = max(versions)
    retry = TRUE
    while(retry) {
      if (version<100) {
        stop(paste0("could not find a valid Ensembl version for organism "), org)
      }
      temp = XML::getHTMLLinks(paste0("http://ftp.ensembl.org/pub/release-", version, "/"))
      if (length(grep("fasta", temp))>0) {
        temp2 <- XML::getHTMLLinks(paste0("http://ftp.ensembl.org/pub/release-", version, "/fasta/"))
        if (length(grep(gsub(" ", "_", stringr::str_to_lower(org)), temp2))>0) {
          retry = FALSE
        } else {
          version = version-1
        }
      } else {
        version = version-1
      }
    }
  }
  version
}

fetch_latest_annotation <- function(org){
  hub <- AnnotationHub::AnnotationHub()
  all_db <- AnnotationHub::query(hub, c(org, "EnsDb"))
  versions <- stringr::str_extract(all_db$title, "Ensembl [0-9]*") %>% stringr::str_remove("Ensembl ") %>% as.numeric
  max(versions)
}

fetch_ensembl_organism <- function(org, db, release){
  if (is.na(release)) {
    release = fetch_latest_release(org, db)
  }
  query_org = gsub(" ", "_", stringr::str_to_lower(org))
  is_found <- length(grep(query_org, XML::getHTMLLinks(paste0("http://ftp.ensembl.org/pub/release-", release, "/fasta/")))) > 0
  return(is_found)
}

read_fasta_to_biostrings <- function(path){
  if (!is.character(path)) {
    stop("Path to the fasta file must be a string")
  }
  if (!c("fa", "fasta") %in% stringr::str_split_1(path, "\\.")) {
    stop("File name must end with .fasta, .fa, .fa.gz or .fasta.gz")
  }
  fasta_biostrings = Biostrings::readDNAStringSet(path)
}

parse_fasta_names <- function(fasta_biostrings){
  temp <- lapply(fasta_biostrings@ranges@NAMES, function(name) {
    colnames <- stringr::str_extract_all(name, "\\s\\w+:") %>% unlist() %>% stringr::str_remove_all("[\\s:]")
    spaced <- stringr::str_split_1(name, " ")
    values = c(spaced[2])
    temp_values = lapply(colnames, function(x){
      search_name = paste0(x, ":");
      if (search_name == "description:") {
        stringr::str_split_1(name, search_name)[2]
      } else {
        grep(search_name, spaced, value = TRUE) %>% gsub(pattern = search_name, replacement = "", x = .)}
      }) %>% unlist()
    values = c(values, temp_values)
    colnames = c("type", colnames)
    ids = rep(spaced[1], times = length(values))
    data.frame(ids = ids, colnames = colnames, values = values)
  })
  parsed = dplyr::bind_rows(temp, .id = "column_label") %>% tidyr::pivot_wider(., names_from = colnames, values_from = values)

}

# anno_from_fasta <- function(fasta_path, anno_path) {
#   fasta_biostrings <- read_fasta_to_biostrings(fasta_path)
#   fasta_df <- parse_fasta_names(fasta_biostrings)
#   anno_df <- format_anno(fasta_df) # add entrez id + rename col
#   readr::write_csv(anno_df, anno_path)
#
# }


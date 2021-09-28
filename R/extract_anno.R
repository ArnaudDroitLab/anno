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

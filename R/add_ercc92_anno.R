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

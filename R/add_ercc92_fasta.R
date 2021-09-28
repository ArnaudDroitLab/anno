add_ercc92_fasta <- function(ref_fasta) {
  ercc92_filename <- system.file("extdata/ERCC92.fa", package = "rnaseq")
  ercc92_fasta <- Biostrings::readDNAStringSet(ercc92_filename)
  c(ref_fasta, ercc92_fasta)
}

output_verification <- function(org, db, release, ERCC92 = FALSE, prefix) {
  outdir = system.file(paste0("extdata/tests/", prefix), package = "anno")
  prepare_anno(org, db, release, ERCC92, force_download = FALSE, gtf = FALSE,
               outdir = outdir)
  expect_equal(
    unname(tools::md5sum(system.file(paste0("extdata/tests/", prefix, "/", prefix, ".cleaned_ref.fa.gz"), package = "anno"))),
    unname(tools::md5sum(system.file(paste0("extdata/tests/reference/", prefix, ".cleaned_ref.fa.gz"), package = "anno"))))
  expect_equal(
    unname(tools::md5sum(system.file(paste0("extdata/tests/", prefix, "/", prefix, ".cleaned_ref.csv"), package = "anno"))),
    unname(tools::md5sum(system.file(paste0("extdata/tests/reference/", prefix, ".cleaned_ref.csv"), package = "anno"))))
  expect_equal(
    unname(tools::md5sum(system.file(paste0("extdata/tests/", prefix, "/", prefix, ".no_alt_chr.fa.gz"), package = "anno"))),
    unname(tools::md5sum(system.file(paste0("extdata/tests/reference/", prefix, ".no_alt_chr.fa.gz"), package = "anno"))))
  expect_equal(
    unname(tools::md5sum(system.file(paste0("extdata/tests/", prefix, "/", prefix, ".no_alt_chr.csv"), package = "anno"))),
    unname(tools::md5sum(system.file(paste0("extdata/tests/reference/", prefix, ".no_alt_chr.csv"), package = "anno"))))
  expect_equal(
    unname(tools::md5sum(system.file(paste0("extdata/tests/", prefix, "/", prefix, ".protein_coding.fa.gz"), package = "anno"))),
    unname(tools::md5sum(system.file(paste0("extdata/tests/reference/", prefix, ".protein_coding.fa.gz"), package = "anno"))))
  expect_equal(
    unname(tools::md5sum(system.file(paste0("extdata/tests/", prefix, "/", prefix, ".protein_coding.csv"), package = "anno"))),
    unname(tools::md5sum(system.file(paste0("extdata/tests/reference/", prefix, ".protein_coding.csv"), package = "anno"))))
  files <- list.files(system.file(paste0("extdata/tests/", prefix, "/"), package = "anno"))
  files <- grep("raw_ref", files, invert = TRUE, value = TRUE)
  file.remove(paste(outdir, files, sep = "/"))

}


test_that("Arguments checking works", {
  valid_org <- "Mus musculus"
  msg <- "org must be a string and a supported organism"
  expect_error(prepare_anno(org = NA), msg)
  msg <- "Abracadabra is an unsupported organism"
  expect_error(prepare_anno(org = "Abracadabra"), msg)
  msg <- "db must be either Ensembl or Gencode"
  expect_error(prepare_anno(org = valid_org, db = NA), msg)
  expect_error(prepare_anno(org = valid_org, db = "Abracadabra"), msg)
  msg <- "Gencode database only accepts Homo Sapiens and Mus musculus organisms"
  expect_error(prepare_anno(org = "Bos taurus", db = "Gencode"), msg)
  msg <- "release must be a number"
  expect_error(prepare_anno(org = valid_org, release = "Abracadabra"), msg)
  msg <- "release must be >= 100 for Ensembl"
  expect_error(prepare_anno(org = valid_org, db = "Ensembl", release = 0), msg)
  msg <- "release must be >= 35 for Homo sapiens in Gencode"
  expect_error(prepare_anno(org = "Homo sapiens", db = "Gencode", release = 0), msg)
  msg <- "release must be >= 25 for Mus musculus in Gencode"
  expect_error(prepare_anno(org = "Mus musculus", db = "Gencode", release = 0), msg)
  msg <- "ERCC92 must be either TRUE or FALSE"
  expect_error(prepare_anno(org = valid_org, ERCC92 = NA), msg)
  expect_error(prepare_anno(org = valid_org, ERCC92 = "Abracadabra"), msg)
  msg <- "force_download must be either TRUE or FALSE"
  expect_error(prepare_anno(org = valid_org, force_download = NA), msg)
  expect_error(prepare_anno(org = valid_org, force_download = "Abracadabra"), msg)
  msg <- "gtf must be either TRUE or FALSE"
  expect_error(prepare_anno(org = valid_org, gtf = NA), msg)
  expect_error(prepare_anno(org = valid_org, gtf = "Abracadabra"), msg)
  msg <- "outdir must be a string and a valid directory"
  expect_error(prepare_anno(org = valid_org, outdir = NA), msg)
  expect_error(prepare_anno(org = valid_org, outdir = 123), msg)
  msg <- "outdir must be a valid directory"
  expect_error(prepare_anno(org = valid_org, outdir = "Abracadabra"), msg)

})

# Mus musculus Ensembl 10000 will be used with download = FALSE because it is a
# manually generated file with 0 length transcripts in it
test_that("0 length transcripts checking works", {
 msg <- "0 length transcripts in fasta file"
 expect_error(prepare_anno(org = "Mus musculus", db = "Ensembl", release = 10000,
                           outdir = system.file("extdata/tests/error", package = "anno")), msg)
})

# release too high verify error

# Mus musculus Ensembl 10001 will be used with download = FALSE because it is a
# manually generated file with duplicated transcripts in it
test_that("Duplicated transcripts checking works", {
  msg <- "Duplicated transcripts in fasta"
  expect_error(prepare_anno(org = "Mus musculus", db = "Ensembl", release = 10001,
                            outdir = system.file("extdata/tests/error", package = "anno")), msg)
})

# verify getting last release

test_that("Hs Ensembl 108 works", {
  output_verification("Homo sapiens", "Ensembl", 108, FALSE, "Hs.Ensembl108")
})

test_that("Mm Ensembl 102 works", {
  output_verification("Mus musculus", "Ensembl", 102, FALSE, "Mm.Ensembl102")
})

test_that("Mm Ensembl 108 works", {
  output_verification("Mus musculus", "Ensembl", 108, FALSE, "Mm.Ensembl108")
})

test_that("Mmu Ensembl 108 works", {
  output_verification("Macaca mulatta", "Ensembl", 108, FALSE, "Mmu.Ensembl108")
})

test_that("Rn Ensembl 102 works", {
  output_verification("Rattus norvegicus", "Ensembl", 102, FALSE, "Rn.Ensembl102")
})

test_that("Rn Ensembl 108 works", {
  output_verification("Rattus norvegicus", "Ensembl", 108, FALSE, "Rn.Ensembl108")
})

test_that("Bt Ensembl 108 works", {
  output_verification("Bos taurus", "Ensembl", 108, FALSE, "Bt.Ensembl108")
})

test_that("Mm Gencode 31 works", {
  output_verification("Mus musculus", "Gencode", 31, FALSE, "Mm.Gencode31")
})

test_that("Hs Gencode 42 works", {
  output_verification("Homo sapiens", "Gencode", 42, FALSE, "Hs.Gencode42")
})

# R compression does not give the same md5sum for ERCC92 and human 108 than bash compression
# test_that("Hs Ensembl 108 with ERCC92 works", {
#   output_verification("Homo sapiens", "Ensembl", 108, TRUE, "Hs.Ensembl108.ERCC92")
# })

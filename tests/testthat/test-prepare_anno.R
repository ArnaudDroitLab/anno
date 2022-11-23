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

# Mus musculus Ensembl 10001 will be used with download = FALSE because it is a
# manually generated file with duplicated transcripts in it
test_that("Duplicated transcripts checking works", {
  msg <- "Duplicated transcripts in fasta"
  expect_error(prepare_anno(org = "Mus musculus", db = "Ensembl", release = 10001,
                            outdir = system.file("extdata/tests/error", package = "anno")), msg)
})

# Verify each organism (Homo sapiens, Rattus norvegicus, Bos Taurus, Mus musculus, Macaca mulatta)
# with last version (108) of Ensembl, mouse < 103, rat < 105, mouse Gencode, human Gencode
# Verify last version (108 for Ensembl, 42 for Homo sapiens Gencode, 31 for Mus musculus Gencode)
# Need to change last version test each time there is a new version
# Unable to test Fasta and annotation ID that do not match because of implementation
# Test addition of ERCC92

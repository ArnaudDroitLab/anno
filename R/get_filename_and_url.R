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

### Code to generate source files for a MgDb-object with data from the
### Silva v128
library(DECIPHER)
library(Biostrings)
library(metagenomeFeatures) ## Needed for the make_mgdb_sqlite function
library(tidyr)
library(tidyverse)
library(stringr)
library(R.utils)
library(data.table)
library(digest)

## Database URL
db_root_url <- "https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/"
taxa_url <- paste0(db_root_url, "taxonomy/taxmap_slv_ssu_ref_128.txt.gz")
taxmap_url <- paste0(db_root_url, "taxonomy/tax_slv_ssu_128.txt")
aligned_seq_url <- paste0(db_root_url, "SILVA_128_SSURef_tax_silva_full_align_trunc.fasta.gz")
seq_url <- paste0(db_root_url, "SILVA_128_SSURef_tax_silva.fasta.gz")
tree_url <- paste0(db_root_url, "taxonomy/tax_slv_ssu_128.tre")
rnacentral_url <- "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/id_mapping.tsv.gz"
rnacentral_md5_url <- "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/md5/md5.tsv.gz"
## Downloaded files
taxa_file <- tempfile()
taxagz_file <- tempfile()
seq_file <- tempfile()
seqgz_file <- tempfile()
taxmap_file <- tempfile()
rnacentralgz_file <- tempfile()
rnacentral_file <- tempfile()
rnacentral_md5_file <- tempfile()
rnacentralgz_md5_file <- tempfile()
tree_file <- "../extdata/silva128.1.tre"

## MD5 check sums from initial download
taxa_md5 <- "766e219360430c7117fd6d76a4525b04"
taxmap_md5 <- "fe5ff3ce9681edc8cf523a4706e76243"
seq_md5 <- "e6feea458c5f1194192c0710c0f40938"
tree_md5 <- "81fe22e0f57a478bda5cea67857baade"

## MgDb database files name
db_file <- "../extdata/silva128.1.sqlite"
metadata_file <- "../extdata/silva128.1_metadata.RDS"

### Download database files ####################################################
download_db <- function(url, file_name, md5){
    ## Downloade file and check to make sure MD5 checksum matches checksum for
    ## previously downloaded version

    download.file(url,file_name)
    new_md5 <- tools::md5sum(file_name)
    if (md5 != new_md5) warning("checksum does not match downloaded file.")
}

## Taxa Data
download_db(taxa_url, taxagz_file, taxa_md5)
gunzip(filename = taxagz_file, destname = taxa_file)

## Taxonomy mapping Data
download_db(taxmap_url, taxmap_file, taxmap_md5)

## Seq Data
download_db(seq_url, seqgz_file, seq_md5)
gunzip(filename = seqgz_file, destname = seq_file)

## Tree Data
download_db(tree_url, tree_file, tree_md5)

##Directly downloading RNAcentral mapping file
download.file(rnacentral_url, rnacentralgz_file)
gunzip(filename = rnacentralgz_file, destname = rnacentral_file)

##Downloading RNAcentral id to md5 mapping file
download.file(rnacentral_md5_url, rnacentralgz_md5_file)
gunzip(filename = rnacentralgz_md5_file, destname = rnacentral_md5_file)

### Parse silva taxonomy
parse_silva <- function(taxa_file, taxmap_file){
    taxa <- fread(taxa_file, sep = "\t", header = TRUE , stringsAsFactors = FALSE)
    taxamap <- fread(taxmap_file, sep = "\t", stringsAsFactors = FALSE, header = FALSE)
    tpath <- taxa$path

    leaf <- lapply(taxamap$V1, function(n) {
        last_node <- str_split(n, ";")
        last_node[[1]][length(last_node[[1]]) - 1]
    })
    taxamap$leaf <- leaf
    taxa_sub <- unique(data.frame(node = unlist(leaf), level = taxamap$V3))
    taxa_list <- taxa_sub$level
    names(taxa_list) <- taxa_sub$node
    req <- c("domain", "phylum", "class", "order", "family", "genus")
    hierarchy <- lapply(tpath, function(path) {
        last_node <- str_split(path, ";")
        nodes <- last_node[[1]][1:length(last_node[[1]]) - 1]
        level <- lapply(nodes, function(node) {
            as.character(taxa_list[[node]][1])
        })
        names(nodes) <- level
        out <- nodes[req]
        out
    })
    df_hierarchy <- data.frame(matrix(unlist(hierarchy),
                                      nrow = length(hierarchy),
                                      byrow = T))
    taxa$domain <- df_hierarchy$X1
    taxa$phylum <- df_hierarchy$X2
    taxa$class <- df_hierarchy$X3
    taxa$order <- df_hierarchy$X4
    taxa$family <- df_hierarchy$X5
    taxa$genus <- df_hierarchy$X6

    taxa$path <- NULL
    taxa$taxid <- NULL
    colnames(taxa) <- c("Accession", "start", "stop","Species", "Kingdom", "Phylum", "Class", "Ord", "Family", "Genus" )
    taxa$Keys <- paste0(taxa$Accession,".", taxa$start,".", taxa$stop)
    taxa_tbl_new <- taxa[,c(11, 1, 2, 3, 4, 5, 6 , 7, 8, 9, 10)]
    ## Return as a data.frame
    data.frame(taxa_tbl_new)
}

taxa_tbl <- parse_silva(taxa_file, taxmap_file)

rna_seqs <- Biostrings::readRNAStringSet(seq_file)
seqs <- Biostrings::DNAStringSet(rna_seqs)
seq_names_final <- vapply(strsplit(names(seqs)," ", fixed = TRUE),
                          `[`, 1, FUN.VALUE = character(1))
names(seqs) <- seq_names_final


add_rnacentral_mapping <- function(rnacentral_md5_file, rnacentral_file, taxa_tbl){
    md5mapping <- fread(rnacentral_md5_file, sep = "\t", header = FALSE )
    id_map <- fread(rnacentral_file, sep = "\t", header = FALSE )
    new_idmap <- id_map[,c("V1","V4")]
    dedup_idmap <- subset(new_idmap,!duplicated(new_idmap$V1))
    dedup_idmap$digest <- md5mapping$V2[match(dedup_idmap$V1,md5mapping$V1)]
    seqsdigest <- sapply(as.character(seqs), digest, algo = "md5", serialize = F)
    seqsdigest_tbl <- as.data.frame(seqsdigest)
    seqsdigest_tbl$Keys <- names(seqsdigest)
    colnames(seqsdigest_tbl) <- c("md5digest", "Keys")
    seqsdigest_tbl$RNAcentralID <- dedup_idmap$V1[match(seqsdigest_tbl$md5digest, dedup_idmap$digest)]
    seqsdigest_tbl$NCBItaxonID <- dedup_idmap$V4[match(seqsdigest_tbl$md5digest, dedup_idmap$digest)]
    taxa_tbl$RNAcentralID <- seqsdigest_tbl$RNAcentralID[match(taxa_tbl$Keys, seqsdigest_tbl$Keys)]
    taxa_tbl$NCBItaxonID <- seqsdigest_tbl$NCBItaxonID[match(taxa_tbl$Keys, seqsdigest_tbl$Keys)]
    ## Return as a data.frame
    data.frame(taxa_tbl)
}

taxa_tbl_final <- add_rnacentral_mapping(rnacentral_md5_file, rnacentral_file, taxa_tbl)
# save(taxa_tbl_final, file = "taxa_tbl_final.Rdata")
# save(seqs, file = "seqs.Rdata")
metagenomeFeatures:::make_mgdb_sqlite(db_name = "silva",
                                     db_file = db_file,
                                     taxa_tbl = taxa_tbl_final,
                                     seqs = seqs)

# rnacentral_tbl <- fread(rnacentral_file, sep = '\t', header = FALSE)


### Database Metadata ##########################################################
metadata <- list(ACCESSION_DATE = date(),
                 URL = "https://www.arb-silva.de/fileadmin/silva_databases/release_128/",
                 DB_TYPE_NAME = "SILVA",
                 DB_VERSION = "128 99% OTUS",
                 DB_TYPE_VALUE = "MgDb",
                 DB_SCHEMA_VERSION = "2.0")

saveRDS(metadata, file = metadata_file)

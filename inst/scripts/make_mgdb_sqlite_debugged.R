## Using Debugged Seqs2DB function
make_mgdb_sqlite <- function(db_name,db_file,taxa_tbl,seqs){
   if (!(is.character(db_name) & length(db_name) == 1)) {
        stop("db_name must be a character string")
    }
    if (!(is.character(db_file) & length(db_file) == 1)) {
        stop("db_name must be a character string")
    }
    if (file.exists(db_file)) {
        warning("db_file exists adding sequence and taxa data will be added to the database")
    }
    if (!is.data.frame(taxa_tbl)) {
        stop("taxa_tbl must be a data frame")
    }
    if (is.character(seqs)) {
        if (file.exists(seqs)) {
            seqs <- Biostrings::readDNAStringSet(seqs)
        }
        else {
            stop("seqs is a character string but no file exists, check filename")
        }
    } else if (!is(seqs, "DNAStringSet")) {
        stop("seqs must be either a DNAStringSet class object or path to a fasta file")
    }
    taxa_keys <- taxa_tbl$Keys
    seq_keys <- names(seqs)
    if (length(taxa_keys) != length(seq_keys)) {
        stop("taxa_tbl$Keys and names(seqs) must match")
    }
    if (sum(taxa_keys %in% seq_keys) != length(taxa_keys)) {
        stop("taxa_tbl$Keys and names(seqs) must match")
    }
    if (sum(seq_keys %in% taxa_keys) != length(seq_keys)) {
        stop("taxa_tbl$Keys and names(seqs) must match")
    }
    if (length(unique(taxa_keys)) != length(taxa_keys)) {
        stop("taxa_tbl$Keys must be unique")
    }
    taxa_tbl$Keys <- as.character(taxa_tbl$Keys)
    taxa_tbl <- taxa_tbl[match(names(seqs), taxa_tbl$Keys), ]
    rownames(taxa_tbl) <- 1:nrow(taxa_tbl)
    db_conn <- RSQLite::dbConnect(RSQLite::SQLite(), db_file)
    # DECIPHER::
    Seqs2DB(seqs = seqs, type = "DNAStringSet", dbFile = db_conn, identifier = "MgDb")
    db_seqs <- RSQLite::dbReadTable(db_conn, "Seqs")
    taxa_tbl$row_names <- seq(1:nrow(taxa_tbl))
    db_merge_table <- merge(db_seqs, taxa_tbl, by = "row_names")
    RSQLite::dbWriteTable(db_conn, "Seqs", db_merge_table, overwrite = TRUE)
    RSQLite::dbDisconnect(db_conn)
}


## Standin with bug fix for dbWriteTables function call
Seqs2DB <- function(seqs, type, dbFile, identifier, tblName = "Seqs", chunkSize = 1e+07,
    replaceTbl = FALSE, fields = c(accession = "ACCESSION", rank = "ORGANISM"),
    processors = 1, verbose = TRUE, ...){
    time.1 <- Sys.time()
    if (length(seqs) == 0)
        stop("seqs is zero length.")
    SEQTYPES <- c("FASTA", "FASTQ", "GenBank", "XStringSet",
        "DNAStringSet", "RNAStringSet", "AAStringSet", "BStringSet",
        "QualityScaledXStringSet", "QualityScaledDNAStringSet",
        "QualityScaledRNAStringSet", "QualityScaledAAStringSet",
        "QualityScaledBStringSet")
    type <- pmatch(type, SEQTYPES)
    if (is.na(type))
        stop("Invalid seqs type.")
    if (type == -1)
        stop("Ambiguous seqs type.")
    if (length(type) > 1)
        stop("type must be length 1.")
    if (!is.character(identifier))
        stop("Identifier must be a character!")
    if (!is.character(tblName))
        stop("tblName must be a character string.")
    if (substr(tblName, 1, 1) == "_")
        stop("Invalid tblName.")
    if (tblName == "taxa")
        stop("taxa is a reserved tblName.")
    if (!is.numeric(chunkSize))
        stop("chunkSize must be a numeric.")
    if (floor(chunkSize) != chunkSize)
        stop("chunkSize must be a whole number.")
    if (chunkSize <= 0)
        stop("chunkSize must be greater than zero.")
    if (!is.logical(replaceTbl))
        stop("replaceTbl must be a logical.")
    if (!is.logical(verbose))
        stop("verbose must be a logical.")
    if (!is.null(processors) && !is.numeric(processors))
        stop("processors must be a numeric.")
    if (!is.null(processors) && floor(processors) != processors)
        stop("processors must be a whole number.")
    if (!is.null(processors) && processors < 1)
        stop("processors must be at least 1.")
    if (is.null(processors)) {
        processors <- detectCores()
    } else {
        processors <- as.integer(processors)
    }
    if (!is(seqs, "XStringSet") && !is(seqs, "QualityScaledXStringSet") &&
        type > 3)
        stop("seqs must be an XStringSet or QualityScaledXStringSet.")
    if (type == 3 && length(fields) > 0) {
        if (!is.character(fields))
            stop("fields must be a character vector.")
        if (any(is.na(fields)))
            stop("fields must be a character vector.")
        if (is.null(names(fields)) || any(is.na(names(fields))))
            stop("fields must be a named character vector.")
        if (any(names(fields) == ""))
            stop("The names of fields cannot be empty.")
        if (any(fields == ""))
            stop("fields cannot be empty.")
        if (any(duplicated(names(fields))))
            stop("The names of fields must be unique.")
        if (any(duplicated(fields)))
            stop("fields must be unique.")
        if (any(fields %in% "DEFINITION"))
            stop("fields contains a field that is already imported.")
        if (any(names(fields) %in% c("description", "identifier",
            "row_names")))
            stop("The names of fields contains a reserved name.")
        if (any(grepl("[^A-Za-z]", fields)))
            stop("fields may only contain letters.")
        if (any(nchar(fields) > 12))
            stop("fields may be at most 12 characters wide.")
        fields <- c(description = "DEFINITION", fields)
    } else {
        fields <- c(description = "DEFINITION")
    }
    driver = dbDriver("SQLite")
    numSeq <- integer(1)
    if (is.character(dbFile)) {
        dbConn = dbConnect(driver, dbFile)
        on.exit(dbDisconnect(dbConn))
    } else {
        dbConn = dbFile
        if (!inherits(dbConn, "SQLiteConnection"))
            stop("'dbFile' must be a character string or connection.")
        if (!dbIsValid(dbConn))
            stop("The connection has expired.")
    }
    result <- dbListTables(dbConn)
    w <- which(result == tblName)
    if (length(w) == 1 && length(which(result == paste("_", tblName,
        sep = ""))) != 1)
        stop("Table is corrupted")
    f <- character(0)
    if (length(w) == 1 && !replaceTbl) {
        searchExpression <- paste("select max(row_names) from ",
            tblName, sep = "")
        numSeq <- as.integer(dbGetQuery(dbConn, searchExpression))
        searchExpression <- paste("select max(row_names) from _",
            tblName, sep = "")
        if (as.integer(dbGetQuery(dbConn, searchExpression)) !=
            numSeq)
            stop("Table is corrupted.")
        f <- dbListFields(dbConn, paste("_", tblName, sep = ""))
        colIDs <- c("row_names", "sequence", "quality")
        fts <- c("INTEGER PRIMARY KEY ASC", "BLOB", "BLOB")
        for (i in 1:length(colIDs)) {
            if (is.na(match(colIDs[i], f))) {
                expression1 <- paste("alter table _", tblName,
                  " add column ", colIDs[i], " ", fts[i], sep = "")
                if (verbose)
                  cat("Expression:  ", expression1, "\n", sep = "")
                rs <- dbSendStatement(dbConn, expression1)
                dbClearResult(rs)
            }
        }
        f <- dbListFields(dbConn, tblName)
        if (type != 3) {
            colIDs <- c("row_names", "identifier", "description")
            fts <- c("INTEGER PRIMARY KEY ASC", "TEXT", "TEXT")
        }
        else {
            colIDs <- c("row_names", "identifier", "description",
                names(fields))
            fts <- c("INTEGER PRIMARY KEY ASC", rep("TEXT", length(fields) +
                2))
        }
        for (i in 1:length(colIDs)) {
            if (is.na(match(colIDs[i], f))) {
                expression1 <- paste("alter table ", tblName,
                  " add column ", colIDs[i], " ", fts[i], sep = "")
                if (verbose)
                  cat("Expression:  ", expression1, "\n", sep = "")
                rs <- dbSendStatement(dbConn, expression1)
                dbClearResult(rs)
            }
        }
    } else {
        numSeq <- 0
        replaceTbl <- TRUE
    }
    if (verbose) {
        it <- 0L
        priorWidth <- 0L
        .cat <- function(text, ...) {
            width <- nchar(text)
            if (width < priorWidth) {
                text <- paste(text, paste(rep(" ", priorWidth -
                  width), collapse = ""), sep = "")
            }
            priorWidth <<- width
            cat(text, sep = "")
            flush.console()
        }
    }
    if (type == 1) {
        if (is.character(seqs)) {
            if (substr(seqs, 1, 7) == "http://" || substr(seqs,
                1, 8) == "https://" || substr(seqs, 1, 7) ==
                "ftps://" || substr(seqs, 1, 6) == "ftp://") {
                con <- gzcon(url(seqs, "rb"))
            }
            else {
                con <- gzfile(seqs, "rb")
            }
            on.exit(close(con))
        }
        else if (inherits(seqs, "connection")) {
            con <- seqs
            if (!isOpen(con)) {
                open(con, type = "rb")
                on.exit(close(con))
            }
        }
        else {
            stop("seqs must be a file path or connection.")
        }
        newSeqs <- 0
        buffer <- ""
        enter <- TRUE
        newline <- FALSE
        while (enter) {
            if (verbose) {
                it <- it + 1L
                .cat(paste("\rReading FASTA file chunk ", it,
                  sep = ""))
            }
            r <- readChar(con = con, nchars = chunkSize, useBytes = TRUE)
            if (length(r) == 0L) {
                if (buffer == "") {
                  break
                }
                else {
                  enter <- FALSE
                }
            }
            else if (nchar(r) < chunkSize) {
                r2 <- readChar(con = con, nchars = chunkSize -
                  nchar(r), useBytes = TRUE)
                while (length(r2) > 0 && nchar(r) < chunkSize) {
                  r <- paste(r, r2, sep = "")
                  r2 <- readChar(con = con, nchars = chunkSize -
                    nchar(r), useBytes = TRUE)
                }
                if (nchar(r) < chunkSize)
                  enter <- FALSE
            }
            if (buffer != "")
                r <- paste(buffer, r, sep = ifelse(newline, "\n",
                  ""))
            newline <- substr(r, nchar(r), nchar(r)) == "\n"
            r <- strsplit(r, "\n", useBytes = TRUE, fixed = TRUE)[[1]]
            r <- gsub("\r", "", r, useBytes = TRUE, fixed = TRUE)
            descriptions <- which(substr(r, 1L, 1L) == ">")
            if (length(descriptions) == 0)
                stop("No FASTA records found.")
            dp <- descriptions + 1L
            dm <- descriptions - 1L
            end <- c(dm[-1], length(r))
            if (enter) {
                buffer <- paste(r[descriptions[length(descriptions)]:length(r)],
                  collapse = "\n")
                length(descriptions) <- length(descriptions) -
                  1L
            }
            numF <- length(descriptions)
            if (numF == 0)
                next
            myData <- data.frame(row_names = seq(from = (numSeq +
                1), to = (numSeq + length(descriptions))), identifier = identifier)
            myData$description <- substr(r[descriptions], 2L,
                nchar(r[descriptions]))
            sequence <- .Call("collapse", r, dp[seq_len(numF)],
                end[seq_len(numF)], PACKAGE = "DECIPHER")
            sequence <- gsub(" ", "", sequence, useBytes = TRUE,
                fixed = TRUE)
            myData_ <- data.frame(row_names = seq(from = (numSeq +
                1), to = (numSeq + length(descriptions))), sequence = I(Codec(sequence,
                processors = processors, ...)), quality = I(rep(list(raw()),
                length(sequence))))
            if (length(f) > 0) {
                for (i in 1:length(f)) {
                  if (is.na(match(f[i], names(myData)))) {
                    d <- data.frame(rep(NA, numF))
                    names(d) <- f[i]
                    myData <- data.frame(myData, d)
                  }
                }
            }
            numSeq <- numSeq + length(descriptions)
            newSeqs <- newSeqs + length(descriptions)
            if (replaceTbl) {
                ft <- c(row_names = "INTEGER PRIMARY KEY ASC",
                  identifier = "TEXT", description = "TEXT")
                ft_ <- c(row_names = "INTEGER PRIMARY KEY ASC",
                  sequence = "BLOB", quality = "BLOB")
            }
            else {
                ft <- ft_ <- NULL
            }
            dbWriteTable(dbConn, tblName, myData, row.names = FALSE,
                overwrite = replaceTbl, append = !replaceTbl,
                field.types = ft)
            dbWriteTable(dbConn, paste("_", tblName, sep = ""),
                myData_, row.names = FALSE, overwrite = replaceTbl,
                append = !replaceTbl, field.types = ft_)
            replaceTbl <- FALSE
        }
    } else if (type == 2) {
        if (is.character(seqs)) {
            if (substr(seqs, 1, 7) == "http://" || substr(seqs,
                1, 8) == "https://" || substr(seqs, 1, 7) ==
                "ftps://" || substr(seqs, 1, 6) == "ftp://") {
                con <- gzcon(url(seqs, "rb"))
            }
            else {
                con <- gzfile(seqs, "rb")
            }
            on.exit(close(con))
        }
        else if (inherits(seqs, "connection")) {
            con <- seqs
            if (!isOpen(con)) {
                open(con, type = "rb")
                on.exit(close(con))
            }
        }
        else {
            stop("seqs must be a file path or connection.")
        }
        newSeqs <- 0
        buffer <- ""
        enter <- TRUE
        newline <- FALSE
        while (enter) {
            if (verbose) {
                it <- it + 1L
                .cat(paste("\rReading FASTQ file chunk ", it,
                  sep = ""))
            }
            r <- readChar(con = con, nchars = chunkSize, useBytes = TRUE)
            if (length(r) == 0L) {
                if (buffer == "") {
                  break
                }
                else {
                  enter <- FALSE
                }
            }
            else if (nchar(r) < chunkSize) {
                r2 <- readChar(con = con, nchars = chunkSize -
                  nchar(r), useBytes = TRUE)
                while (length(r2) > 0 && nchar(r) < chunkSize) {
                  r <- paste(r, r2, sep = "")
                  r2 <- readChar(con = con, nchars = chunkSize -
                    nchar(r), useBytes = TRUE)
                }
                if (nchar(r) < chunkSize)
                  enter <- FALSE
            }
            if (buffer != "")
                r <- paste(buffer, r, sep = ifelse(newline, "\n",
                  ""))
            newline <- substr(r, nchar(r), nchar(r)) == "\n"
            r <- strsplit(r, "\n", useBytes = TRUE, fixed = TRUE)[[1]]
            r <- gsub("\r", "", r, useBytes = TRUE, fixed = TRUE)
            descriptions <- which(substr(r, 1L, 1L) == "@")
            if (length(descriptions) == 0)
                stop("No FASTQ records found.")
            descriptions <- descriptions[which((descriptions -
                descriptions[1])%%4 == 0)]
            if (enter) {
                buffer <- paste(r[descriptions[length(descriptions)]:length(r)],
                  collapse = "\n")
                length(descriptions) <- length(descriptions) -
                  1L
            }
            numF <- length(descriptions)
            if (numF == 0)
                next
            myData <- data.frame(row_names = seq(from = (numSeq +
                1), to = (numSeq + length(descriptions))), identifier = identifier)
            myData$description <- substr(r[descriptions], 2L,
                nchar(r[descriptions]))
            myData_ <- data.frame(row_names = seq(from = (numSeq +
                1), to = (numSeq + length(descriptions))), sequence = I(Codec(r[descriptions +
                1L], processors = processors, ...)), quality = I(Codec(r[descriptions +
                3L], compression = c("qbit", "gzip"), processors = processors)))
            if (length(f) > 0) {
                for (i in 1:length(f)) {
                  if (is.na(match(f[i], names(myData)))) {
                    d <- data.frame(rep(NA, numF))
                    names(d) <- f[i]
                    myData <- data.frame(myData, d)
                  }
                }
            }
            numSeq <- numSeq + length(descriptions)
            newSeqs <- newSeqs + length(descriptions)
            if (replaceTbl) {
                ft <- c(row_names = "INTEGER PRIMARY KEY ASC",
                  identifier = "TEXT", description = "TEXT")
                ft_ <- c(row_names = "INTEGER PRIMARY KEY ASC",
                  sequence = "BLOB", quality = "BLOB")
            }
            else {
                ft <- ft_ <- NULL
            }
            dbWriteTable(dbConn, tblName, myData, row.names = FALSE,
                overwrite = replaceTbl, append = !replaceTbl,
                field.types = ft)
            dbWriteTable(dbConn, paste("_", tblName, sep = ""),
                myData_, row.names = FALSE, overwrite = replaceTbl,
                append = !replaceTbl, field.types = ft_)
            replaceTbl <- FALSE
        }
    } else if (type == 3) {
        if (is.character(seqs)) {
            if (substr(seqs, 1, 7) == "http://" || substr(seqs,
                1, 8) == "https://" || substr(seqs, 1, 7) ==
                "ftps://" || substr(seqs, 1, 6) == "ftp://") {
                con <- gzcon(url(seqs, "rb"))
            } else {
                con <- gzfile(seqs, "rb")
            }
            on.exit(close(con))
        } else if (inherits(seqs, "connection")) {
            con <- seqs
            if (!isOpen(con)) {
                open(con, type = "rb")
                on.exit(close(con))
            }
        } else {
            stop("seqs must be a file path or connection.")
        }
        newSeqs <- 0
        buffer <- ""
        enter <- TRUE
        newline <- FALSE
        while (enter) {
            if (verbose) {
                it <- it + 1L
                .cat(paste("\rReading GenBank file chunk ", it,
                  sep = ""))
            }
            r <- readChar(con = con, nchars = chunkSize, useBytes = TRUE)
            if (length(r) == 0L) {
                if (buffer == "") {
                  break
                } else {
                  enter <- FALSE
                }
            } else if (nchar(r) < chunkSize) {
                r2 <- readChar(con = con, nchars = chunkSize -
                  nchar(r), useBytes = TRUE)
                while (length(r2) > 0 && nchar(r) < chunkSize) {
                  r <- paste(r, r2, sep = "")
                  r2 <- readChar(con = con, nchars = chunkSize -
                    nchar(r), useBytes = TRUE)
                }
                if (nchar(r) < chunkSize)
                  enter <- FALSE
            }
            if (buffer != "")
                r <- paste(buffer, r, sep = ifelse(newline, "\n",
                  ""))
            newline <- substr(r, nchar(r), nchar(r)) == "\n"
            r <- strsplit(r, "\n", useBytes = TRUE, fixed = TRUE)[[1]]
            r <- gsub("\r", "", r, useBytes = TRUE, fixed = TRUE)
            descriptions <- which(substr(r, 1L, 5L) == "LOCUS")
            if (length(descriptions) == 0)
                stop("No GENBANK records found.")
            seq_start <- which(substr(r, 1L, 6L) == "ORIGIN")
            seq_end <- which(substr(r, 1L, 2L) == "//") - 1L
            temp <- integer(length(seq_end))
            for (i in seq_along(seq_end)) {
                w <- which(seq_start > descriptions[i] & seq_start <=
                  seq_end[i])
                if (length(w) == 0) {
                  temp[i] <- seq_end[i] + 1L
                } else if (length(w) == 1) {
                  temp[i] <- seq_start[w] + 1L
                } else {
                  stop("More than one ORIGIN found in GenBank record.")
                }
            }
            seq_start <- temp
            if (enter) {
                buffer <- paste(r[descriptions[length(descriptions)]:length(r)],
                  collapse = "\n")
                length(descriptions) <- length(descriptions) -
                  1L
                length(seq_start) <- length(seq_end) <- length(descriptions)
            }
            numF <- length(descriptions)
            if (numF == 0)
                next
            x <- .Call("extractFields", r, fields, descriptions,
                seq_start, PACKAGE = "DECIPHER")
            names(x) <- names(fields)
            myData <- data.frame(row_names = seq(from = (numSeq +
                1), to = (numSeq + length(descriptions))), identifier = identifier,
                x)
            sequence <- .Call("collapse", substr(r, 11, nchar(r)),
                seq_start[seq_len(numF)], seq_end[seq_len(numF)],
                PACKAGE = "DECIPHER")
            sequence <- gsub(" ", "", sequence, useBytes = TRUE,
                fixed = TRUE)
            myData_ <- data.frame(row_names = seq(from = (numSeq +
                1), to = (numSeq + length(descriptions))), sequence = I(Codec(sequence,
                processors = processors, ...)), quality = I(rep(list(raw()),
                length(sequence))))
            if (length(f) > 0) {
                for (i in 1:length(f)) {
                  if (is.na(match(f[i], names(myData)))) {
                    d <- data.frame(rep(NA, numF))
                    names(d) <- f[i]
                    myData <- data.frame(myData, d)
                  }
                }
            }
            numSeq <- numSeq + length(descriptions)
            newSeqs <- newSeqs + length(descriptions)
            if (replaceTbl) {
                ft <- lapply(1:dim(myData)[2], function(x) return("TEXT"))
                ft[[1]] <- "INTEGER PRIMARY KEY ASC"
                names(ft) <- names(myData)
                ft_ <- c(row_names = "INTEGER PRIMARY KEY ASC",
                  sequence = "BLOB", quality = "BLOB")
            } else {
                ft <- ft_ <- NULL
            }
            dbWriteTable(dbConn, tblName, myData, row.names = FALSE,
                overwrite = replaceTbl, append = !replaceTbl,
                field.types = ft)
            dbWriteTable(dbConn, paste("_", tblName, sep = ""),
                myData_, row.names = FALSE, overwrite = replaceTbl,
                append = !replaceTbl, field.types = ft_)
            replaceTbl <- FALSE
        }
    } else {
        newSeqs <- length(seqs)
        if (length(newSeqs) == 0)
            stop("No sequences in seqs.")
        if (verbose) {
            cat("Adding", newSeqs, "sequences to the database.")
            flush.console()
        }
        if (is.null(names(seqs)))
            names(seqs) <- (numSeq + 1):(numSeq + newSeqs)
        myData <- data.frame(row_names = seq(from = (numSeq +
            1), to = (numSeq + newSeqs)), identifier = identifier,
            description = names(seqs))
        if (type > 8) {
            quality <- Codec(as.character(quality(seqs)), compression = c("qbit",
                "gzip"), processors = processors)
        }
        else {
            quality <- rep(list(raw()), length(seqs))
        }
        myData_ <- data.frame(row_names = seq(from = (numSeq +
            1), to = (numSeq + newSeqs)), sequence = I(Codec(as.character(seqs),
            processors = processors)), quality = I(quality))
        if (length(f) > 0) {
            for (i in 1:length(f)) {
                if (is.na(match(f[i], names(myData)))) {
                  d <- data.frame(rep(NA, newSeqs))
                  names(d) <- f[i]
                  myData <- data.frame(myData, d)
                }
            }
        }
        if (replaceTbl) {
            ft <- c(row_names = "INTEGER PRIMARY KEY ASC",
                identifier = "TEXT", description = "TEXT")
            ft_ <- c(row_names = "INTEGER PRIMARY KEY ASC",
                sequence = "BLOB", quality = "BLOB")
        } else {
            ft <- ft_ <- NULL
        }
        dbWriteTable(dbConn, tblName, myData, row.names = FALSE,
            overwrite = replaceTbl, append = !replaceTbl, field.types = ft)
        dbWriteTable(dbConn, paste("_", tblName, sep = ""), myData_,
            row.names = FALSE, overwrite = replaceTbl, append = !replaceTbl,
            field.types = ft_)
    }
    searchExpression <- paste("select count(*) from ", tblName,
        sep = "")
    numSeq <- as.integer(dbGetQuery(dbConn, searchExpression))
    if (verbose) {
        time.2 <- Sys.time()
        cat("\n")
        if (newSeqs != numSeq)
            cat("\nAdded ", newSeqs, " new sequence", ifelse(newSeqs >
                1, "s", ""), " to table ", tblName, ".", sep = "")
        cat("\n", numSeq, " total sequence", ifelse(numSeq >
            1, "s", ""), " in table ", tblName, ".", "\n", sep = "")
        print(round(difftime(time.2, time.1, units = "secs"),
            digits = 2))
        cat("\n")
        flush.console()
    }
    invisible(numSeq)
}

#!/usr/bin/env Rscript
#- scalpelR reference pages of the documentation: one page per exported function, from the package help files (man/*.Rd)
#- usage: Rscript scalpelR_reference.R <scalpelR source directory> <output directory (docs/scalpelR)>
suppressPackageStartupMessages(library(tools, quietly = TRUE))

args = commandArgs(TRUE)
PKG = args[1]
OUT = args[2]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

#- exported functions & their help file...
exports = sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", readLines(file.path(PKG, "NAMESPACE")), value = TRUE))
version = sub("^Version: *", "", grep("^Version:", readLines(file.path(PKG, "DESCRIPTION")), value = TRUE))
pages = c()
for(f in sort(exports)){
    rd = parse_Rd(file.path(PKG, "man", paste0(f, ".Rd")))
    title = trimws(paste(unlist(rd[sapply(rd, attr, "Rd_tag") == "\\title"]), collapse = ""))
    #- body of the help page (without the page header & title), embedded in the MyST page
    html = capture.output(Rd2HTML(rd, out = "", no_links = TRUE))
    html = html[(grep("<body", html)[1] + 1):(grep("</body>", html)[1] - 1)]
    #- page header (first table: name & "R Documentation") and title...
    h = grep("<table", html)[1]
    if(!is.na(h) && grepl("R Documentation", paste(html[h:(h + 3)], collapse = ""))) html = html[-(h:grep("</table>", html)[1])]
    html = html[!grepl("^<h2|^<div class=\"container\">|^</div>$", html)]
    writeLines(c(paste0("# ", f), "", title, "", "```{raw} html", html, "```"), file.path(OUT, paste0(f, ".md")))
    pages = c(pages, f)
    message("  ", f)
}

#- index page: table of the functions & toctree...
titles = sapply(pages, function(f) sub("^\\s+", "", readLines(file.path(OUT, paste0(f, ".md")))[3]))
writeLines(c(
    "# scalpelR reference", "",
    paste0("Functions of the R package [scalpelR](https://github.com/Franzx7/scalpelR) (version ", version, ") for the downstream analysis of the SCALPEL results."), "",
    "| Function | Description |", "| --- | --- |",
    paste0("| [`", pages, "`](", pages, ".md) | ", titles, " |"), "",
    "```{toctree}", ":hidden:", ":maxdepth: 1", "", pages, "```"), file.path(OUT, "index.md"))

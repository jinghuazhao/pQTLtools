utils::globalVariables(c("SNP","p"))

#' QTL lookup with LD-based proxy selection
#'
#' This function performs QTL lookup by integrating MR results with GWAS summary
#' statistics and identifying proxy SNPs based on linkage disequilibrium (LD).
#'
#' It supports:
#' - 1000 Genomes LD reference (ieugwasr API)
#' - Local PLINK reference panels (recommended for large-scale runs)
#'
#' Loci are classified as:
#' - same_locus (r^2 >= threshold)
#' - independent_locus (r^2 < threshold)
#'
#' @param d Directory containing GWAS files.
#' @param dat Data frame of MR results. Must contain:
#'   protein, id, pqtl, qtl, p_qtl, file_gwas.
#' @param panel LD reference type: "1000Genomes" or "local".
#' @param p_threshold P-value cutoff for selecting candidate SNPs.
#' @param r2_threshold LD r^2 threshold for locus classification.
#' @param pop 1000 Genomes population (default "EUR").
#' @param plink_bin Path to PLINK executable (required for local mode).
#' @param max_snps Maximum SNPs used in LD computation.
#' @param xlsx Optional Excel output path.
#' @param verbose Logical; show progress messages.
#'
#' @return Updated dat with columns:
#' proxy, p_proxy, rsq, classification.
#'
#' @export
#'
#'
qtl_lookup <- function(
    d,
    dat,
    panel = c("1000Genomes", "local"),
    p_threshold = 1e-3,
    r2_threshold = 0.8,
    pop = "EUR",
    plink_bin = NULL,
    max_snps = 500,
    xlsx = NULL,
    verbose = TRUE
) {
    panel <- match.arg(panel)

    required <- c("protein","id","pqtl","qtl","p_qtl","file_gwas")
    miss <- setdiff(required, names(dat))
    if (length(miss) > 0) {
        stop("Missing required columns: ", paste(miss, collapse = ", "))
    }

    if (panel == "local" && is.null(plink_bin)) {
        stop("plink_bin required for local panel mode")
    }

    for (col in c("proxy","p_proxy","rsq","classification")) {
        if (!col %in% names(dat)) dat[[col]] <- NA
    }

    strip_alleles <- function(x) gsub("_[A-Z]+$", "", x)

    compute_ld <- function(snps, bfile = NULL) {
        if (panel == "1000Genomes") {
            if (length(snps) > max_snps) {
                stop("Too many SNPs (> max_snps) for 1000G mode")
            }
            ieugwasr::ld_matrix(variants = snps, pop = pop)
        } else {
            ieugwasr::ld_matrix_local(
                variants = snps,
                bfile = bfile,
                plink_bin = plink_bin,
                with_alleles = TRUE
            )
        }
    }

    ld_cache <- new.env(parent = emptyenv())

    for (i in seq_len(nrow(dat))) {

        z <- dat[i, , drop = FALSE]
        gwas_file <- file.path(d, basename(z$file_gwas))

        if (verbose) {
            message("[", i, "/", nrow(dat), "] ", basename(gwas_file))
        }

        if (!file.exists(gwas_file)) next

        gwas <- tryCatch(
            data.table::fread(gwas_file, select = c("SNP","p")),
            error = function(e) NULL
        )
        if (is.null(gwas)) next

        hits <- gwas[!is.na(SNP) & !is.na(p) & p <= p_threshold, ]
        if (nrow(hits) == 0) next

        if (nrow(hits) > max_snps) {
            hits <- hits[order(p)][seq_len(max_snps)]
        }

        panel_snps <- unique(c(z$pqtl, hits$SNP))
        panel_snps <- panel_snps[!is.na(panel_snps)]
        if (length(panel_snps) < 2) next

        cache_key <- paste(sort(panel_snps), collapse = ":")

        if (exists(cache_key, envir = ld_cache)) {
            ld_mat <- get(cache_key, envir = ld_cache)
        } else {
            ld_mat <- tryCatch(
                compute_ld(panel_snps, bfile = z$bfile),
                error = function(e) NULL
            )
            assign(cache_key, ld_mat, envir = ld_cache)
        }

        if (is.null(ld_mat) || !is.matrix(ld_mat)) next

        colnames(ld_mat) <- strip_alleles(colnames(ld_mat))
        rownames(ld_mat) <- strip_alleles(rownames(ld_mat))

        pqtl <- z$pqtl
        if (!(pqtl %in% colnames(ld_mat))) next

        r2_mat <- ld_mat^2
        snps <- intersect(hits$SNP, colnames(r2_mat))
        if (length(snps) == 0) next

        proxy_df <- data.frame(
            proxy = snps,
            rsq = as.numeric(r2_mat[pqtl, snps])
        )

        proxy_df <- merge(
            proxy_df,
            gwas[, c("SNP","p")],
            by.x = "proxy",
            by.y = "SNP",
            all.x = TRUE
        )

        names(proxy_df)[3] <- "p_proxy"
        proxy_df <- proxy_df[!is.na(proxy_df$rsq), ]
        if (nrow(proxy_df) == 0) next

        same <- proxy_df[proxy_df$rsq >= r2_threshold, ]
        ind  <- proxy_df[proxy_df$rsq < r2_threshold, ]

        chosen <- NULL
        class <- NA_character_

        if (nrow(same) > 0) {
            chosen <- same[order(-same$rsq, same$p_proxy), ][1, ]
            class <- "same_locus"
        } else if (nrow(ind) > 0) {
            chosen <- ind[order(ind$p_proxy), ][1, ]
            class <- "independent_locus"
        }

        if (!is.null(chosen)) {
            dat$proxy[i] <- chosen$proxy
            dat$p_proxy[i] <- chosen$p_proxy
            dat$rsq[i] <- chosen$rsq
            dat$classification[i] <- class
        }
    }

    if (!is.null(xlsx)) {
        wb <- openxlsx::createWorkbook()

        style <- openxlsx::createStyle(
            textDecoration = "bold",
            fontColour = "#FFFFFF",
            fgFill = "#4F80BD"
        )

        out <- dat[, c(
            "protein","id","pqtl","qtl","p_qtl",
            "proxy","p_proxy","rsq","classification"
        )]

        openxlsx::addWorksheet(wb, "proxies")
        openxlsx::writeDataTable(wb, "proxies", out, headerStyle = style)
        openxlsx::freezePane(wb, "proxies", firstActiveRow = 2)
        openxlsx::saveWorkbook(wb, xlsx, overwrite = TRUE)
    }

    dat
}

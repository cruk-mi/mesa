# Offline stand-ins for biomaRt::useMart() and biomaRt::getBM(), so the
# plotGeneHeatmap() tests need no network (#123). Lookups are served from
# fixtures/biomart-genes.csv, which make-biomart-genes.R in that folder
# regenerates.

.cachedBiomartGenes <- function() {
    utils::read.csv(testthat::test_path("fixtures", "biomart-genes.csv"),
        stringsAsFactors = FALSE)
}

localMockBiomart <- function(env = parent.frame()) {
    genes <- .cachedBiomartGenes()
    testthat::local_mocked_bindings(
        useMart = function(biomart, dataset, host, ...) {
            methods::new("Mart", biomart = biomart, dataset = dataset,
                host = host)
        },
        getBM = function(attributes, filters, values, mart, ...) {
            atHost <- vapply(genes$host, grepl, logical(1), x = mart@host)
            keep <- atHost & genes$filter == filters & genes$value %in% values
            hits <- genes[keep, ]
            if (nrow(hits) == 0) {
                stop("No cached gene for ", values, " at ", mart@host,
                    "; add it to fixtures/make-biomart-genes.R.")
            }
            hits[, attributes]
        },
        .package = "biomaRt",
        .env = env
    )
}

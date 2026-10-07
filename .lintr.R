# Lint configuration for mesa.
#
# Scope: R/ only. The rules below come from the "R / Bioconductor coding
# standards" section of CLAUDE.md, which govern package source. They are
# deliberately NOT applied to tests/ or vignettes/:
#
#   * tests/testthat.R legitimately calls library(testthat)/library(mesa) --
#     the CLAUDE.md rule is about library() inside R/, not inside tests.
#   * tests/ and vignettes/ carry 167 lines over 80 characters. Reflowing them
#     is a cosmetic churn that CLAUDE.md explicitly discourages ("do not
#     rewrite working code to match a style preference", "keep diffs small").
#
# The linter set is listed EXPLICITLY (defaults = list()) rather than built on
# lintr::default_linters(), because two of those defaults are wrong here:
#
#   * object_name_linter() defaults to snake_case and would flag every exported
#     function (makeQset, calculateDMRs, plotRegionsHeatmap, ...). mesa follows
#     the qsea/Bioconductor camelCase convention instead.
#   * cyclocomp_linter() would flag the large functions in R/qseaExtra.R and
#     R/PCA.R. That is a refactoring question, not a style gate.
#
# object_usage_linter() is omitted too: it has to load the package, which would
# drag the whole Bioconductor dependency tree into what is meant to be a
# dependency-free, sub-minute job. `R CMD check`'s "checking R code for
# possible problems" already reports undefined globals.
#
# indentation_linter(indent = 4L) is omitted: CLAUDE.md asks for 4-space
# indentation, but lintr also enforces continuation/hanging indents and reports
# 198 violations across R/ today. Enabling it is a separate piece of work.
#
# Reference: https://contributions.bioconductor.org/r-code.html
linters <- lintr::linters_with_defaults(
    defaults = list(),
    # "Line length <= 80 characters for R code"
    line_length_linter = lintr::line_length_linter(80L),
    # "No tabs"
    whitespace_linter = lintr::whitespace_linter(),
    # "Use `<-` for assignment, never `=` at the top level"
    assignment_linter = lintr::assignment_linter(),
    # "Never use T / F as substitutes for TRUE / FALSE"
    T_and_F_symbol_linter = lintr::T_and_F_symbol_linter(),
    # "Avoid 1:n loops" -- note this catches the 1:length(x) and 1:nrow(x)
    # forms specifically; a bare 1:n is not statically detectable.
    seq_linter = lintr::seq_linter(),
    semicolon_linter = lintr::semicolon_linter(),
    # Covers ":::", "<<-" and "->>", all three named in CLAUDE.md
    undesirable_operator_linter = lintr::undesirable_operator_linter(),
    # "Do not use library() or require() inside package code"
    undesirable_function_linter = lintr::undesirable_function_linter(c(
        library = "use the DESCRIPTION Imports field and pkg::fn() calls",
        require = "use the DESCRIPTION Imports field and pkg::fn() calls",
        attach = "modifies the search path"
    ))
)

exclusions <- list("tests", "vignettes", "data-raw", "docs", "man")

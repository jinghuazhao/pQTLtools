# Auxiliary files

File/Directory   | Description
-----------------|------------------------
Bioconductor/    | Bioconductor-related
eQTL-Catalogue/  | eQTL-catalog-related
PPI/             | PPI-related
STRING/          | StringDB-related
snakemake/       | snakemake showcase
scripts/         | Bash scripts[^scripts]
tests/           | SCALLOP-INF-related
turboman/        | turboman reference data
UniProt/         | UniProt-related
README.md        | This file
REFERENCES.bib   | BibTeX bibliography
nature-genetics.csl| CSL-style file

[^scripts]: **Two Bash scripts**

    1. [docs.sh](scripts/docs.sh). To build pkgdown-style website on HPC and interact with GitHub.
    2 [cran.sh](scripts/cran.sh). A chain to build, install, and check (--as-cran) the package.

    The first script is run first to establish inst/doc/{lz,stack}.html on pQTLtools/articles/{lz,stack}.html, followed by validation through the second script. Generally, this implies that a separate website can be used.

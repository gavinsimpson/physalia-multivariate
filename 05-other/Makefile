all: slides purl examples
examples: coca prc varpart lda

slides: slides.Rmd slides.css
	Rscript -e "rmarkdown::render(\"slides.Rmd\")"

purl: slides.Rmd
	Rscript -e "knitr::purl(\"slides.Rmd\")"

coca: cocorrespondence-analysis.Rmd
	Rscript -e "rmarkdown::render(\"cocorrespondence-analysis.Rmd\")"

prc: prc.Rmd
	Rscript -e "rmarkdown::render(\"prc.Rmd\")"

varpart: variation-partitioning.Rmd
	Rscript -e "rmarkdown::render(\"variation-partitioning.Rmd\")"

lda: linear-discriminants.Rmd
	Rscript -e "rmarkdown::render(\"linear-discriminants.Rmd\")"

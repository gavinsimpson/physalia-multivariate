all: slides purl coca prc dbrda permanova

slides: slides.Rmd slides.css
	Rscript -e "rmarkdown::render(\"slides.Rmd\")"

purl: slides.Rmd
	Rscript -e "knitr::purl(\"slides.Rmd\")"

coca: cocorrespondence-analysis.Rmd
	Rscript -e "rmarkdown::render(\"cocorrespondence-analysis.Rmd\")"

prc: prc.Rmd
	Rscript -e "rmarkdown::render(\"prc.Rmd\")"

dbrda: dbrda.Rmd
	Rscript -e "rmarkdown::render(\"dbrda.Rmd\")"

permanova: permanova.Rmd
	Rscript -e "rmarkdown::render(\"permanova.Rmd\")"

all: slides purl ponds

slides: slides.Rmd slides.css
	Rscript -e "rmarkdown::render(\"slides.Rmd\")"

purl: slides.Rmd
	Rscript -e "knitr::purl(\"slides.Rmd\")"

ponds: constrained-ordination.Rmd
	Rscript -e "rmarkdown::render(\"constrained-ordination.Rmd\")"


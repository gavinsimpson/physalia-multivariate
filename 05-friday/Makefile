all: slides purl

slides: slides.Rmd slides.css
	Rscript -e "rmarkdown::render(\"slides.Rmd\")"

purl: slides.Rmd
	Rscript -e "knitr::purl(\"slides.Rmd\")"


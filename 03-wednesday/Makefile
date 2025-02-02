all: slides purl ponds meadows

slides: slides.Rmd slides.css
	Rscript -e 'library("rmarkdown"); render("slides.Rmd")'


purl: slides.Rmd
	Rscript -e "knitr::purl(\"slides.Rmd\")"

ponds: constrained-ordination.Rmd
	Rscript -e "rmarkdown::render(\"constrained-ordination.Rmd\")"

meadows: spring-meadows.Rmd
	Rscript -e "rmarkdown::render(\"spring-meadows.Rmd\")"


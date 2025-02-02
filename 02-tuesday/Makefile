all: slides birds ponds dune

slides: slides.Rmd slides.css
	Rscript -e 'library("rmarkdown"); render("slides.Rmd")'

ponds: ponds.Rmd
	Rscript -e "rmarkdown::render(\"ponds.Rmd\")"

birds: birds.Rmd
	Rscript -e "rmarkdown::render(\"birds.Rmd\")"

dune: dune-meadows.Rmd
	Rscript -e "rmarkdown::render(\"dune-meadows.Rmd\")"
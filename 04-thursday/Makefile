all: slides purl random permutations

slides: slides.Rmd slides.css
	Rscript -e "rmarkdown::render(\"slides.Rmd\")"

purl: slides.Rmd
	Rscript -e "knitr::purl(\"slides.Rmd\")"

random: randomised-complete-block.Rmd
	Rscript -e "rmarkdown::render(\"randomised-complete-block.Rmd\")"

permutations: permutation-tests-solutions.Rmd
	Rscript -e "rmarkdown::render(\"permutation-tests-solutions.Rmd\")"

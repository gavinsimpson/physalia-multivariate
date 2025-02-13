all: slides purl random permutations solutions

slides: slides.Rmd slides.css
	Rscript -e "rmarkdown::render(\"slides.Rmd\")"

purl: slides.Rmd
	Rscript -e "knitr::purl(\"slides.Rmd\")"

random: randomised-complete-block.Rmd
	Rscript -e "rmarkdown::render(\"randomised-complete-block.Rmd\")"

solutions: permutation-tests-solutions.Rmd
	Rscript -e "rmarkdown::render(\"permutation-tests-solutions.Rmd\")"

permutations: permutation-tests.Rmd
	Rscript -e "rmarkdown::render(\"permutation-tests.Rmd\")"

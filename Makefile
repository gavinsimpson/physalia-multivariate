.PHONY : all
all:
	cd ./01-monday && make $@
	cd ./02-tuesday && make $@
	cd ./03-wednesday && make $@
	cd ./04-thursday && make $@

.PHONY : slides
slides:
	cd ./01-monday && make $@
	cd ./02-tuesday && make $@
	cd ./03-wednesday && make $@
	cd ./04-thursday && make $@

.PHONY : purl
purl:
	cd ./01-monday && make $@
	cd ./02-tuesday && make $@
	cd ./03-wednesday && make $@
	cd ./04-thursday && make $@


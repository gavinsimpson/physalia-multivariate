## ----setup-options, echo = FALSE, results = "hide", message = FALSE-----------
options(htmltools.dir.version = FALSE)
knitr::opts_chunk$set(cache = TRUE, dev = 'svg', echo = TRUE, message = FALSE, warning = FALSE,
                      fig.height=6, fig.width = 1.777777*6)
library("vegan")
library("ggplot2")
data(varespec)
data(varechem)

## plot defaults
theme_set(theme_minimal(base_size = 16, base_family = 'Fira Sans'))


## ----permustats-1, results = "hide"-------------------------------------------
cca1 <- cca(varespec ~ ., data = varechem)
pstat <- permustats(anova(cca1))
summary(pstat)


## ----permustats-1, echo = FALSE-----------------------------------------------
cca1 <- cca(varespec ~ ., data = varechem)
pstat <- permustats(anova(cca1))
summary(pstat)


## ----permustats-2, fig.width = 6, fig.height = 6------------------------------
densityplot(pstat)


## ----cca-anova----------------------------------------------------------------
set.seed(42)
(perm <- anova(cca1))


## ----cyclic-shift-figure, echo = FALSE----------------------------------------
knitr::include_graphics("./resources/cyclic-shifts-figure.png")


## ----shuffle-time-series, echo = TRUE-----------------------------------------
shuffle(10, control = how(within = Within(type = "series")))


## ----set-up-toroidal----------------------------------------------------------
set.seed(4)
h <- how(within = Within(type = "grid",
                         ncol = 3, nrow = 3))
perm <- shuffle(9, control = h)
matrix(perm, ncol = 3)


## ----toroidal-shifts-figure, echo = FALSE-------------------------------------
knitr::include_graphics("./resources/Toroidal_coord.png")


## ----sketch-1, echo = FALSE, out.width = "75%", fig.align = "center"----------
knitr::include_graphics("./resources/permutation-designs-sketch-1.png")


## ----sketch-2, echo = FALSE, out.width = "75%", fig.align = "center"----------
knitr::include_graphics("./resources/permutation-designs-sketch-2.png")


## ----sketch-3, echo = FALSE, out.width = "75%", fig.align = "center"----------
knitr::include_graphics("./resources/permutation-designs-sketch-3.png")


## ----cyclic-shift-mirror-figure, echo = FALSE---------------------------------
knitr::include_graphics("./resources/cyclic-shifts-with-mirror-figure.svg")


## -----------------------------------------------------------------------------
plt <- gl(3, 10)
h <- how(within = Within(type = "series"), plots = Plots(strata = plt))


## ----helper-funs--------------------------------------------------------------
args(Within)
args(Plots)


## ----how-args-----------------------------------------------------------------
args(how)


## ----ts-perm-example1---------------------------------------------------------
plt <- gl(3, 10)
h <- how(within = Within(type = "series"),
         plots = Plots(strata = plt))
set.seed(4)
p <- shuffle(30, control = h)
do.call("rbind", split(p, plt)) ## look at perms in context


## ----ts-perm-example2---------------------------------------------------------
plt <- gl(3, 10)
h <- how(within = Within(type = "series", constant = TRUE),
         plots = Plots(strata = plt))
set.seed(4)
p <- shuffle(30, control = h)
do.call("rbind", split(p, plt)) ## look at perms in context


## ----worked-example-devel-1---------------------------------------------------
## load vegan
library("vegan")

## load the data
spp <- read.csv("data/ohraz-spp.csv", header = TRUE, row.names = 1)
env <- read.csv("data/ohraz-env.csv", header = TRUE, row.names = 1)
molinia <- spp[, 1]
spp <- spp[, -1]

## Year as numeric
env <- transform(env, year = as.numeric(as.character(year)))


## ----worked-example-devel-2---------------------------------------------------
c1 <- rda(spp ~ year + year:mowing + year:fertilizer + year:removal + Condition(plotid), data = env)
(h <- how(within = Within(type = "free"), plots = Plots(strata = env$plotid, type = "none")))


## ----worked-example-devel-2a--------------------------------------------------
set.seed(42)
anova(c1, permutations = h, model = "reduced")


## ----worked-example-devel-2b--------------------------------------------------
set.seed(24)
anova(c1, permutations = h, model = "reduced", by = "axis")


## ----load-crayfish------------------------------------------------------------
## load data
crayfish <- head(read.csv("data/crayfish-spp.csv")[, -1], -1)
design <- read.csv("data/crayfish-design.csv", skip = 1)[, -1]

## fixup the names
names(crayfish) <- gsub("\\.", "", names(crayfish))
names(design) <- c("Watershed", "Stream", "Reach", "Run",
                   "Stream.Nested", "ReachNested", "Run.Nested")


## ----crayfish-unconstrained---------------------------------------------------
m.pca <- rda(crayfish)
summary(eigenvals(m.pca))


## ----crayfish-pca-plot, fig.show = "hide", collapse = TRUE--------------------
layout(matrix(1:2, ncol = 2))
biplot(m.pca, type = c("text", "points"), scaling = "species")
set.seed(23)
ev.pca <- envfit(m.pca ~ Watershed, data = design, scaling = "species")
plot(ev.pca, labels = levels(design$Watershed), add = FALSE)
layout(1)


## ----crayfish-pca-plot, fig.show = "hold", out.width = "75%", echo = FALSE, fig.height = 5----
layout(matrix(1:2, ncol = 2))
biplot(m.pca, type = c("text", "points"), scaling = "species")
set.seed(23)
ev.pca <- envfit(m.pca ~ Watershed, data = design, scaling = "species")
plot(ev.pca, labels = levels(design$Watershed), add = FALSE)
layout(1)


## ----crayfish-watershed-------------------------------------------------------
m.ws <- rda(crayfish ~ Watershed, data = design)
m.ws


## ----crayfish-watershed-2-----------------------------------------------------
summary(eigenvals(m.ws, constrained = TRUE))


## ----crayfish-watershed-3-----------------------------------------------------
set.seed(1)
ctrl <- how(nperm = 499, within = Within(type = "none"),
            plots = with(design, Plots(strata = Stream, type = "free")))
(sig.ws <- anova(m.ws, permutations = ctrl))


## ----crayfish-stream----------------------------------------------------------
m.str <- rda(crayfish ~ Stream + Condition(Watershed), data = design)
m.str


## ----crayfish-stream-2--------------------------------------------------------
summary(eigenvals(m.str, constrained = TRUE))


## ----crayfish-stream-3--------------------------------------------------------
set.seed(1)
ctrl <- how(nperm = 499, within = Within(type = "none"),
            plots = with(design, Plots(strata = Reach, type = "free")),
            blocks = with(design, Watershed))
(sig.str <- anova(m.str, permutations = ctrl))


## ----crayfish-reach-----------------------------------------------------------
(m.re <- rda(crayfish ~ Reach + Condition(Stream), data = design))


## ----crayfish-reach-2---------------------------------------------------------
set.seed(1)
ctrl <- how(nperm = 499, within = Within(type = "none"),
            plots = with(design, Plots(strata = Run, type = "free")),
            blocks = with(design, Stream))
(sig.re <- anova(m.re, permutations = ctrl))


## ----crayfish-run-------------------------------------------------------------
(m.run <- rda(crayfish ~ Run + Condition(Reach), data = design))


## ----crayfish-run-2-----------------------------------------------------------
set.seed(1)
ctrl <- how(nperm = 499, within = Within(type = "free"),
            blocks = with(design, Reach))
(sig.run <- anova(m.run, permutations = ctrl))


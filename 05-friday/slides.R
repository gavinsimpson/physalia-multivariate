## ----setup-options, echo = FALSE, results = "hide", message = FALSE-----------
options(htmltools.dir.version = FALSE)
knitr::opts_chunk$set(cache = TRUE, dev = 'svg', echo = TRUE, message = FALSE, warning = FALSE,
                      fig.height=6, fig.width = 1.777777*6)
library("vegan")
library("ggplot2")
library("dplyr")
library("tibble")
library("janitor")
library("readr")
data(varespec)
data(varechem)

## plot defaults
theme_set(theme_minimal(base_size = 16, base_family = 'Fira Sans'))


## ----load-cocorresp-data------------------------------------------------------
library("cocorresp")
data(beetles)
## log transform the beetle data
beetles <- log1p(beetles)
data(plants)


## ----fit-symcoca, message = TRUE----------------------------------------------
bp.sym <- coca(beetles ~ ., data = plants, method = "symmetric")
bp.sym


## ----out.width = "80%", fig.align = "center"----------------------------------
screeplot(bp.sym)


## ----plot-symcoca, fig.show = "hold", out.width = "80%", fig.align = "center"----
layout(matrix(1:2, ncol = 2))
biplot(bp.sym, which = "y1", main = "Beetles")
biplot(bp.sym, which = "y2", main = "Plants")
layout(1)


## ----eval = FALSE-------------------------------------------------------------
# rda(comm ~ F:t + Condition(t), data = df)


## -----------------------------------------------------------------------------
data(pyrifos)
dim(pyrifos)

ditch <- gl(12, 1, length = 132)
week <- gl(11, 12, labels = c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
dose <- factor(rep(c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44, 0.9, 0, 6), 11))


## -----------------------------------------------------------------------------
mod <- prc(pyrifos, dose, week)


## -----------------------------------------------------------------------------
mod


## -----------------------------------------------------------------------------
ctrl <- how(plots = Plots(strata = ditch,type = "free"),
    within = Within(type = "series"), nperm = 99)
anova(mod, permutations = ctrl, first = TRUE)


## ----out.width = "90%", fig.align = "center"----------------------------------
plot(mod, species = FALSE, legpos = "topright")


## ----fig.keep = "none"--------------------------------------------------------
plot(mod, species = FALSE, legpos = "topright")
logabu <- colSums(pyrifos)
scrs <- scores(mod, display = "species", choices = 1)
linestack(scrs[logabu > 150, , drop = FALSE]); axis(side = 2)


## ----echo = FALSE-------------------------------------------------------------
plot(mod, species = FALSE, legpos = "topright")


## ----fig.width = 3, fig.height = 7, echo = FALSE------------------------------
logabu <- colSums(pyrifos)
scrs <- scores(mod, display = "species", choices = 1, scaling = "symmetric")
linestack(scrs[logabu > 150, , drop = FALSE])
axis(side = 2)


## -----------------------------------------------------------------------------
df <- data.frame(y1 = c(0, 0, 1),
    y2 = c(4, 1, 0),
    y3 = c(8, 1, 0))
D <- as.matrix(vegdist(df))
D
D[upper.tri(D)] # unroll


## -----------------------------------------------------------------------------
crossprod(1:3, 1:3)

sum(1:3 * 1:3)


## -----------------------------------------------------------------------------
data(varespec, varechem, package = "vegan")
D_veg <- vegdist(varespec) # Bray-Curtis
D_env <- vegdist(scale(varechem), "euclidean")
mantel(D_veg, D_env, permutations = how(nperm = 999))


## -----------------------------------------------------------------------------
data(mite, mite.xy, package = "vegan")
D_mite <- vegdist(mite) # Bray-Curtis
D_xy <- vegdist(mite.xy, "euclidean")
mantel(D_mite, D_xy, permutations = how(nperm = 999), method = "spearman")


## ----goodness-----------------------------------------------------------------
upr <- cca(varespec ~ ., data = varechem)
lwr <- cca(varespec ~ 1, data = varechem)
set.seed(1)

mods <- ordistep(lwr, scope = formula(upr), trace = 0)
head(goodness(mods))


## ----inertcomp----------------------------------------------------------------
head(inertcomp(mods, proportional = TRUE))


## ----spenvcor-----------------------------------------------------------------
spenvcor(mods)


## ----intersetcor--------------------------------------------------------------
intersetcor(mods)


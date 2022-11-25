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


## ----permanova-idea-plot, echo = FALSE----------------------------------------
data(varespec)
     
## Bray-Curtis distances between samples
dis <- vegdist(varespec)
     
## First 16 sites grazed, remaining 8 sites ungrazed
groups <- factor(c(rep(1,16), rep(2,8)), labels = c("grazed","ungrazed"))
     
## Calculate multivariate dispersions
mod <- betadisper(dis, groups)
plot(mod)


## ----adonis2-by-terms---------------------------------------------------------
data(dune, dune.env)
adonis2(dune ~ Management*A1, data = dune.env, by = "terms")


## ----adonis2-by-terms-flipped-------------------------------------------------
data(dune, dune.env)
adonis2(dune ~ A1 * Management, data = dune.env, by = "terms")


## ----adonis2-by-margin--------------------------------------------------------
data(dune, dune.env)
adonis2(dune ~ Management * A1, data = dune.env, by = "margin")


## ----adonis2-margin-2---------------------------------------------------------
adonis2(dune ~ Management + A1, data = dune.env, by = "margin")


## ----permanova-idea-plot, echo = FALSE----------------------------------------
data(varespec)
     
## Bray-Curtis distances between samples
dis <- vegdist(varespec)
     
## First 16 sites grazed, remaining 8 sites ungrazed
groups <- factor(c(rep(1,16), rep(2,8)), labels = c("grazed","ungrazed"))
     
## Calculate multivariate dispersions
mod <- betadisper(dis, groups)
plot(mod)


## ----permdisp-----------------------------------------------------------------
data(varespec)
dis <- vegdist(varespec) # Bray-Curtis distances
## First 16 sites grazed, remaining 8 sites ungrazed
groups <- factor(c(rep(1,16), rep(2,8)),
                 labels = c("grazed","ungrazed"))

mod <- betadisper(dis, groups)
mod


## ----permdisp-plot, fig.height = 6, fig.width = 6-----------------------------
boxplot(mod)


## ----permdisp-anova-----------------------------------------------------------
set.seed(25)
permutest(mod)


## ----permdisp-plot-it, fig.width = 6, fig.height = 6--------------------------
plot(mod)


## ----permdisp-anova-2---------------------------------------------------------
set.seed(4)
permutest(mod, pairwise = TRUE)


## ----adonis2-by-margin-as-db-rda----------------------------------------------
data(dune, dune.env)
dune_dbrda <- dbrda(dune ~ Management * A1, data = dune.env,
    method = "bray")


## -----------------------------------------------------------------------------
spp <- read_csv(url("https://bit.ly/ohraz-spp")) %>%
    rename(label = "...1") %>%
    janitor::clean_names()

molinia <- spp %>%
    select(label:molicaer)

spp <- spp %>%
    select(-molicaer) %>%
    column_to_rownames("label")

env <- read_csv(url("https://bit.ly/ohraz-env")) %>%
    rename(label = "...1") %>%
    mutate(across(c(mowing:removal, plotid), ~ factor(.x))) %>%
    column_to_rownames("label")


## ----worked-example-devel-2-with-dbrda----------------------------------------
ohraz_dbrda <- dbrda(spp ~ year +
    year:mowing + year:fertilizer + year:removal +
    Condition(plotid), data = env, method = "bray", add = "lingoes")
h <- how(within = Within(type = "free"),
    plots = Plots(strata = env$plotid, type = "none"))
set.seed(42)
anova(ohraz_dbrda, permutations = h, model = "reduced")


## ----load-cocorresp-data------------------------------------------------------
library("cocorresp")
data(beetles)
## log transform the beetle data
beetles <- log1p(beetles)
data(plants)


## ----fit-symcoca, message = TRUE----------------------------------------------
bp.sym <- coca(beetles ~ ., data = plants, method = "symmetric")
bp.sym


## ---- out.width = "80%", fig.align = "center"---------------------------------
screeplot(bp.sym)


## ----plot-symcoca, fig.show = "hold", out.width = "80%", fig.align = "center"----
layout(matrix(1:2, ncol = 2))
biplot(bp.sym, which = "y1", main = "Beetles")
biplot(bp.sym, which = "y2", main = "Plants")
layout(1)


## ---- eval = FALSE------------------------------------------------------------
## rda(comm ~ F:t + Condition(t), data = df)


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


## ---- out.width = "90%", fig.align = "center"---------------------------------
plot(mod, species = FALSE, legpos = "topright")


## ---- fig.keep = "none"-------------------------------------------------------
plot(mod, species = FALSE, legpos = "topright")
logabu <- colSums(pyrifos)
scrs <- scores(mod, display = "species", choices = 1)
linestack(scrs[logabu > 150, , drop = FALSE]); axis(side = 2)


## ---- echo = FALSE------------------------------------------------------------
plot(mod, species = FALSE, legpos = "topright")


## ---- fig.width = 3, fig.height = 7, echo = FALSE-----------------------------
logabu <- colSums(pyrifos)
scrs <- scores(mod, display = "species", choices = 1, scaling = "symmetric")
linestack(scrs[logabu > 150, , drop = FALSE])
axis(side = 2)


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


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
adonis2(dune ~ A1*Management, data = dune.env, by = "terms")


## ----adonis2-by-margin--------------------------------------------------------
data(dune, dune.env)
adonis2(dune ~ Management*A1, data = dune.env, by = "margin")


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
h <- how(within = Within(type = "none"),
    plots = Plots(strata = env$plotid, type = "free"))
set.seed(42)
anova(ohraz_dbrda, permutations = h, model = "reduced")


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


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


## ----cca-model----------------------------------------------------------------
cca1 <- cca(varespec ~ ., data = varechem)
cca1


## ----rda-model----------------------------------------------------------------
rda1 <- rda(varespec ~ ., data = varechem)
rda1


## ----eigenvals----------------------------------------------------------------
eigenvals(cca1)


## ----your-turn-fit-cca-1------------------------------------------------------
library("vegan")
data(varechem, varespec)


## ----your-turn-fit-cca-2------------------------------------------------------
mycca1 <- cca(varespec ~ N + P + K, data = varechem)
mycca1


## ----your-turn-fit-cca-3------------------------------------------------------
ev <- eigenvals(mycca1, model = "constrained")
head(ev)
length(ev)


## ----scores-------------------------------------------------------------------
str(scores(cca1, choices = 1:4, display = c("species","sites")), max = 1)
head(scores(cca1, choices = 1:2, display = "sites"))


## ----scaling-example, results = "hide"----------------------------------------
scores(cca1, choices = 1:2, display = "species", scaling = 3)


## ----your-turn-cca-4----------------------------------------------------------
scrs <- scores(mycca1, display = "sites", choices = c(2,3),
               scaling = "sites", hill = TRUE)
head(scrs)


## ----partial-ordination-1-----------------------------------------------------
pcca <- cca(X = varespec,
            Y = varechem[, "Ca", drop = FALSE],
            Z = varechem[, "pH", drop = FALSE])


## ----partial-ordination-2-----------------------------------------------------
pcca <- cca(varespec ~ Ca + Condition(pH), data = varechem) ## easier!


## ----partial-ordination-3-----------------------------------------------------
pcca <- cca(varespec ~ Ca + Condition(pH), data = varechem) ## easier!
pcca


## ----triplot-1, fig.height = 6, fig.width = 6---------------------------------
plot(cca1)


## ----lc-vs-wa-scores-1, echo = TRUE, fig.keep = "none"------------------------
# example from Design Decision vignette
data(dune, dune.env, package = "vegan")
ord <- cca(dune ~ Moisture, data = dune.env)
plot(ord, display = "lc", type = "points")


## ----lc-vs-wa-scores-1, echo = FALSE------------------------------------------
# example from Design Decision vignette
data(dune, dune.env, package = "vegan")
ord <- cca(dune ~ Moisture, data = dune.env)
plot(ord, display = "lc", type = "points")


## ----lc-vs-wa-scores-2, echo = TRUE, fig.keep = "none"------------------------
# example from Design Decision vignette
data(dune, dune.env, package = "vegan")
ord <- cca(dune ~ Moisture, data = dune.env)
plot(ord, display = "wa", type = "points")
ordispider(ord, col = "red")
text(ord, display = "cn", col = "blue")


## ----lc-vs-wa-scores-2, echo = FALSE------------------------------------------
# example from Design Decision vignette
data(dune, dune.env, package = "vegan")
ord <- cca(dune ~ Moisture, data = dune.env)
plot(ord, display = "wa", type = "points")
ordispider(ord, col = "red")
text(ord, display = "cn", col = "blue")


## ----cca-model-build1---------------------------------------------------------
vare.cca <- cca(varespec ~ Al + P*(K + Baresoil), data = varechem)
vare.cca


## ----vif-cca1-----------------------------------------------------------------
vif.cca(cca1)


## ----stepwise-1---------------------------------------------------------------
upr <- cca(varespec ~ ., data = varechem)
lwr <- cca(varespec ~ 1, data = varechem)
set.seed(1)
mods <- ordistep(lwr, scope = formula(upr), trace = 0)


## ----stepwise-cca-------------------------------------------------------------
mods


## ----stepwise-anova-----------------------------------------------------------
mods$anova


## ----stepwise-reverse, results = "hide"---------------------------------------
mods2 <- step(upr, scope = list(lower = formula(lwr), upper = formula(upr)), trace = 0,
              test = "perm")
mods2


## ----stepwise-reverse---------------------------------------------------------
mods2 <- step(upr, scope = list(lower = formula(lwr), upper = formula(upr)), trace = 0,
              test = "perm")
mods2


## ----rsq-cca1-----------------------------------------------------------------
RsquareAdj(cca1)


## ----stopping-rules-----------------------------------------------------------
ordiR2step(lwr, upr, trace = FALSE)


## ----permustats-1, results = "hide"-------------------------------------------
pstat <- permustats(anova(cca1))
summary(pstat)


## ----permustats-1, echo = FALSE-----------------------------------------------
pstat <- permustats(anova(cca1))
summary(pstat)


## ----permustats-2, fig.width = 6, fig.height = 6------------------------------
densityplot(pstat)


## ----cca-anova----------------------------------------------------------------
set.seed(42)
(perm <- anova(cca1))


## ----anova-args---------------------------------------------------------------
args(anova.cca)


## ----anova-by-axis------------------------------------------------------------
set.seed(1)
anova(mods, by = "axis")


## ----anova-by-term------------------------------------------------------------
set.seed(5)
anova(mods, by = "terms")


## ----anova-by-margin----------------------------------------------------------
set.seed(10)
anova(mods, by = "margin")


## ----meadows-setup------------------------------------------------------------
# load vegan, dplyr & readr
library("vegan"); library("dplyr"); library("readr")

# load the data
spp <- read_csv("https://bit.ly/meadows-species") %>%
    rename("sample_id" = "...1") %>%
    tibble::column_to_rownames("sample_id")
env <- read_csv("https://bit.ly/meadows-env") %>%
    rename("sample_id" = "...1")


## ----meadows-cca-full---------------------------------------------------------
m1 <- cca(spp ~ ., data = env)
set.seed(32)
anova(m1)


## ----meadows-cca-full-triplot, fig.show = "hide"------------------------------
plot(m1)


## ----meadows-cca-full-triplot, fig.height = 6, fig.width = 6, echo = FALSE----
plot(m1)


## ----meadows-cca-stepwise-----------------------------------------------------
set.seed(67)
lwr <- cca(spp ~ 1, data = env)
( m2 <- ordistep(lwr, scope = formula(m1), trace = FALSE) )


## ----meadows-cca-reduced-triplot, fig.show = "hide"---------------------------
plot(m2)


## ----meadows-cca-reduced-triplot, fig.height = 6, fig.width = 6, echo = FALSE----
plot(m2)


## ----meadows-cca-anova--------------------------------------------------------
m2$anova


## ----meadows-rda--------------------------------------------------------------
spph <- decostand(spp, method = "hellinger")
m3 <- rda(spph ~ ., data = env)
lwr <- rda(spph ~ 1, data = env)
m4 <- ordistep(lwr, scope = formula(m3),
               trace = FALSE)


## ----meadows-rda-print--------------------------------------------------------
m4


## ----meadows-rda-reduced-triplot, fig.show = "hide"---------------------------
plot(m4)


## ----meadows-rda-reduced-triplot, fig.height = 6, fig.width = 6, echo = FALSE----
plot(m4)


## ----meadows-rda-adjrsquare---------------------------------------------------
m5 <- ordiR2step(lwr, scope = formula(m3), trace = FALSE)
m5$anova


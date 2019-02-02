# gibby

Gibbs sampler experiments.

In theory, you should be able to install by invoking:

`devtools::install_github("maj-tki/gibby", build_opts = c("--no-resave-data", "--no-manual"), force = TRUE)`

However, you may need to do `devtools::build_vignettes("gibby")` independently to build the vignettes - I do not know why.

See: https://community.rstudio.com/t/vignettes-suddenly-stopped-installing/18391 and https://github.com/r-lib/devtools/issues/1896 for more on that issue.

When I am developing, prior to pushing to github, I use `devtools::build()`.

I tend to use `utils::browseVignettes()` to render vignettes in a browser (e.g. chrome) but you can also use  `RShowDoc("gibby_intro", package= "gibby" )`. The rstudio viewer does not render the equations.

There is negligible help doc, but the `gibby_intro` vignette spells out what you need to know. Minimal implementation and no tests so use at your own risk.


# gibby

In theory, you should be able to just do:

`devtools::install_github("maj-tki/gibby", build_opts = c("--no-resave-data", "--no-manual"), force = TRUE)`

However, you may need to do `devtools::build_vignettes("gibby")` independently - I do not know why.

See: https://community.rstudio.com/t/vignettes-suddenly-stopped-installing/18391 and https://github.com/r-lib/devtools/issues/1896.

Gibbs sampler experiments. When I am developing, prior to pushing to github, I use `devtools::build()`.

I tend to use `utils::browseVignettes()` to render vignettes in a browser (e.g. chrome) but you can also use  `RShowDoc("gibby_intro", package= "gibby" )`. The rstudio viewer does not render the equations.




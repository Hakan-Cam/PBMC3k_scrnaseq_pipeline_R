dir.create("data/raw", recursive = TRUE)
dir.create("data/processed", recursive = TRUE)
dir.create("scripts")
dir.create("figures")
dir.create("results")



install.packages("renv")
renv::init()

install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork", "Matrix"))

# DoubletFinder isn't on CRAN — needs remotes
install.packages("remotes")
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")

install.packages("usethis")
usethis::use_git()

system("git add -A")
system('git commit -m "Phase 0 complete: env, folders, packages, git"')

getwd()
.libPaths()
installed.packages()[c("Seurat","dplyr","ggplot2","patchwork","Matrix","DoubletFinder"), "Version"]

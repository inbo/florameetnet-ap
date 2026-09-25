rstudioapi::restartSession()

library(targets)
stopifnot(
  getwd() == file.path(
    rprojroot::find_root(rprojroot::is_git_root),
    "source/florameetnet"
  )
)

tar_visnetwork(physics = TRUE)
#tar_invalidate(florabank_data) #use this when florabank data need to be queried
tar_make(use_crew = TRUE, as_job = TRUE)

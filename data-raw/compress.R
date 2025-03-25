# Meta-analysis of graduated compression stockings (GCS) for prevention of DVT
# As featured in https://doi.org/10.3310/hta9490

event.gcs <- c(15, 0, 11, 4, 7, 8, 5, 0, 7)
n.gcs <- c(97, 8, 50, 110, 65, 25, 126, 104, 80)
event.control <- c(37, 5, 23, 16, 7, 8, 17, 4, 16)
n.control <- c(103, 10, 48, 110, 32, 25, 126, 92, 81)
study <- c(
  "Allan",
  "Barnes",
  "Holford",
  "Inada",
  "Muir",
  "Rosengarten",
  "Shirai",
  "Turner",
  "Turple"
)

compress <- data.frame(study, event.gcs, n.gcs, event.control, n.control)

usethis::use_data(compress, overwrite = TRUE)

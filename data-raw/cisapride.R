# Meta-analysis of Cisapride for non-ulcer dyspepsia
# As featured in Hartung and Knapp 2001: Statist. Med. 2001; 20:3875–3889

event.cisa <- c(15, 12, 29, 42, 14, 44, 14, 29, 10, 17, 38, 19, 21)
n.cisa <- c(16, 16, 34, 56, 22, 54, 17, 58, 14, 26, 44, 29, 38)
event.plac <- c(9,  1, 18, 31,  6, 17,  7, 23,  3,  6, 12, 22, 19)
n.plac <- c(16, 16, 34, 56, 22, 55, 15, 58, 15, 27, 45, 30, 38)
study <- c("Creytens",
           "Milo",
           "Francois and De Nutte",
           "Deruyttere et al.",
           "Hannon",
           "Roesch",
           "De Nutte et al.",
           "Hausken and Bestad",
           "Chung",
           "Van Outryve et al.",
           "Al-Quorain et al.",
           "Kellow et al.",
           "Yeoh et al."
)

cisapride <- data.frame(study, event.cisa, n.cisa, event.plac, n.plac)

usethis::use_data(cisapride, overwrite = TRUE)

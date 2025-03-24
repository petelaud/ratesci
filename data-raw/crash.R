#CRASH meta-analysis updated with final results- Lancet
#Roberts I, Yates D, Sandercock P et al. Effect of intravenous corticosteroids
# on death within 14 days in 10008 adults with clinically significant head
# injury (mrc crash trial): randomised placebo-controlled trial.
#Lancet 2004; 364(9442): 1321–1328. DOI:10.1016/S0140-6736(04)17188-2.
# URL http://www.ncbi.nlm.nih.gov/pubmed/15474134.

event.steroid <- c(9,16,16,26,35,
                   114,8,44,34,33,
                   1,4,13,19,38,
                   0, 1052)
length(event.steroid)
event.control <- c(13,22,16,13,36,
                   38,9,47,7,21,
                   0,4,5,21,49,
                   0, 893)
n.steroid <- c(17,55,67,49,81,
               201,50,81,72,68,
               5,12,98,133,175,
               30, 4985)
n.control <- c(18,55,28,27,83,
               74,50,80,16,62,
               5,12,54,136,195,
               30, 4979)
study <- c(
  "Ransohoff 1972",
  "Alexander 1972",
  "Faupel 1976",
  "Cooper 1979",
  "Hernesniemi 1979",
  "Pitts 1980",
  "Saul 1981",
  "Braakman 1983",
  "Giannotta 1984",
  "Dearden 1986",
  "Chacon 1987",
  "Zagara 1987",
  "Stubbs 1989",
  "Gaab 1994",
  "Grumme 1995",
  "Zarate 1995",
  "MRC CRASH trial 2004"
)

crash <- data.frame(study, event.steroid, n.steroid, event.control, n.control)

usethis::use_data(crash, overwrite = TRUE)



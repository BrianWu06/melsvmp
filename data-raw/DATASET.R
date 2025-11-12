## code to prepare `DATASET` dataset goes here

riesby <- read.table("RIESBY.DAT.txt", na.strings = ".", 
                     col.names = c("id", "hamd", "intcpt", "week", "endog", "endweek"))

posmood <- read.table("moods_example_augmented.dat", header = FALSE, 
                      col.names = c("id", "posmood", "negmood", "t1", "t2", "t3", "t4", "w1", 
                                    "w2", "w3", "w4", "w5", "w6", "other_bs", "other_ws", 
                                    "genderf", "age15", "tirbor", "frustr"))

usethis::use_data(riesby, posmood, overwrite = TRUE)

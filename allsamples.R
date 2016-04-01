inv <- expand.grid(Colony=c(1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 19, 20, 25, 26, 31, 32, 39, 40, 43, 44, 51, 52, 53, 54, 
                            57, 58, 65, 66, 69, 70, 71, 72, 77, 78, 79, 80, 83, 84, 91, 92, 109, 110, 111, 112, 119, 
                            120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 137, 138, 201, 202,
                            205, 206, 207, 208, 211, 212, 215, 216, 223, 224, 227, 228, 233, 234, 237, 238, 239, 240),
                   Date=as.Date(c("2015-08-11", "2015-09-15", "2015-10-01", "2015-10-21", "2015-11-04", "2015-12-04",
                                  "2015-12-18", "2016-01-20", "2016-02-11")))
head(inv)
dim(inv)
# Remove Reef 42 from August and September
inv <- inv[-which(inv$Colony > 200 & inv$Date < as.Date("2015-09-30")), ]
# Remove dead colonies
inv <- inv[-which(inv$Colony %in% c(1, 2, 5, 6, 7, 12, 39, 91, 92, 111, 128, 129)), ]
inv <- inv[-which(inv$Colony==127 & inv$Date > as.Date("2015-09-30")), ]

inv$Reef <- droplevels(cut(inv$Colony, breaks = c(0, 50, 100, 150, 200, 250), labels = c("HIMB", "25", "44", "", "42")))

inv <- inv[, c(1,3,2)]
head(inv)

inv.s <- split(inv, f=inv$Colony)
inv.sr <- sample(inv.s)
inv.sr

inv.sr.df <- ldply(inv.sr, data.frame)[, -1]
inv.sr.df

write.csv(inv.sr.df, file="plate_order.csv")


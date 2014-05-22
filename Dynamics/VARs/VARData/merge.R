
##################################################

# Merges data

##################################################

fnames <- c("rgdp.txt", "unemp.txt", "prices.txt",
	    "wages.txt", "m1.txt")

data <- read.table(fnames[1], header=TRUE)

for(i in 2:length(fnames)){

    toMerge <- read.table(fnames[i], header=TRUE)

    data <- merge(data, toMerge, by.x="date", by.y="date")
}

write.table(data, file="data.txt")



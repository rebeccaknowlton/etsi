
#run all the functions in the master file
#set the working directory to where the txt files are

setwd("...")
setting <- 1

outputfile=c()
for(u in 1:10) {
	 outputfile = rbind(outputfile,read.table(paste("etsi.outputfile", setting, "_030724",u,".txt", sep=""), header = T))
}

each.rows = 3
results <- matrix(NA, nrow = each.rows, ncol = ncol(outputfile))

for (i in 1:each.rows) {
  # want idx 1, 4, 7, ..., 28 
  idx <- seq(from = i, to = 10*each.rows, by = each.rows)
  results[i,] <- colMeans(outputfile[idx,])
}

results <- as.data.frame(results)
colnames(results) <- colnames(outputfile)
print(results)


#run all the functions in the master file
#set the working directory to where the txt files are

setwd("C:/Users/rkkno/Documents/University of Texas at Austin/etsi/output files/Setting 1")
setting <- 1

output.summary=c()
output.all = c()
for(u in 1:10) {
	 output.summary = rbind(output.summary,read.table(paste("etsi.output.summary", setting, "_032824",u,".txt", sep=""), header = T))
	 output.all = rbind(output.all, read.table(paste("etsi.output.all",setting, "_032824",u,".txt",sep=""), header = T))
}

each.rows = 3
results <- matrix(NA, nrow = each.rows, ncol = ncol(output.summary))

for (i in 1:each.rows) {
  # want idx 1, 4, 7, ..., 28 
  idx <- seq(from = i, to = 10*each.rows, by = each.rows)
  results[i,] <- colMeans(output.summary[idx,])
}

results <- as.data.frame(results)
colnames(results) <- colnames(output.summary)
rownames(results) <- rownames(output.summary)[1:3]
print(results)

output.all


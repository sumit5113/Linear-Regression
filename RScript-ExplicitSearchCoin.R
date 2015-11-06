# Create an empty list to store the results
num_of_flips=500
sample_size=10000
prb_head=.50
prb_tail=1-prb_head

results <- list()
for(i in 1:sample_size) {
  coinTosses   <- cumsum(sample(c(-1,1), num_of_flips, replace = TRUE,prob=c(prb_head,prb_tail))) 
  results[i] <- coinTosses[length(coinTosses)]
}

# Unlist the list and create a histogram. Set a title and set the color and breaks
hist(unlist(results),xlab='count of the total numebr head and tail', main = "Histogram of all sample outcoms of coin flips",col = "lightblue")

print("heads and tail count, total number of times appeared, probability of that outcome")

for(i in -num_of_flips:num_of_flips){
  y=sum(results==i)
  print(c(i,y,y/sample_size))
}


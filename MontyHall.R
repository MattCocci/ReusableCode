
library(e1071)
doors = c("A","B","C")
nsim=1000000

sim.mat = matrix(NA, nrow=nsim, ncol=5)
winning.loc = rdiscrete(nsim,doors, prob=rep(1,length(doors)))
first.choice = rdiscrete(nsim,doors, prob=rep(1,length(doors)))
sim.mat[,1] = winning.loc
sim.mat[,2] = first.choice

for(i in 1:length(winning.loc)){
if(winning.loc[i] != first.choice[i]){
## Contestant has selected the goat. Therefore the host must show the other goat and not the winner.
sim.mat[i,3] = doors[ !(doors %in% c(winning.loc[i], first.choice[i])) ]
} else {
rand.goat.options = doors[ !(doors %in% c(winning.loc[i], first.choice[i])) ]
## Contestant has selected the winner and the host then randomly select one of the other two.
rand.goat = rdiscrete(1, rand.goat.options , prob=rep(1,length(rand.goat.options)) )
sim.mat[i,3] = rand.goat
}
## Second choice changes
sim.mat[i,5] = doors[ !(doors %in% c(sim.mat[i,3], first.choice[i])) ]
}
## Second choice is the same as the first choice
sim.mat[,4] = sim.mat[,2]
## names = TRUE.VAL, FIRST.CHOICE, GOAT.LOC.SHOWN, SECOND.CHOICE.STAY
sim.val.stay = as.numeric( (sim.mat[,4]==sim.mat[,1])*1 )
sim.val.change = as.numeric( (sim.mat[,5]==sim.mat[,1])*1 )
mean(sim.val.stay)
mean(sim.val.change)

## Probability by iteration
covergence.prob.change = cumsum(sim.val.change)/cumsum(rep(1,length(sim.val.change)))
covergence.prob.stay = cumsum(sim.val.stay)/cumsum(rep(1,length(sim.val.stay)))

par(mfrow=c(1,2))
plot(head(covergence.prob.stay, 100), pch=16, cex=.5,
ylab="Probability",
main="Probability of Winning\nKeeping the Same Door")
plot(head(covergence.prob.change, 100), pch=16, cex=.5, ylab="Probability",
main="Probability of winning\nChanging Doors")

###########################################################################
## Using R version 3.0.0												
## Using igraph version 0.6.5-1 -WARNING do not use an older version!! 	
## Authors: Despina Stasi, Sonja Petrovic, Elizabeth Gross.				
##                                                                      
## Code for example in Section 4.1. in http://arxiv.org/abs/1401.4896   
##                                                                      
###########################################################################
###########################################################################

## Load library and code:
library("igraph")
source("p1walk.R") 

###########################################################################
## THE GRAPH
d.OHT = graph.empty(8)
b.OHT = graph(c(1,2, 2,3, 1,4, 4,5, 5,6, 5,7, 7,8), d=FALSE)

Plot.Mixed.Graph(d.OHT, b.OHT)

###########################################################################
## Calculate the MLE:
mle = Get.MLE(d.OHT, b.OHT, maxiter=30000, tol = 1e-04)
#write.table(mle, file="mle-OHT.txt", col.names=F, row.names=F)


###########################################################################
## RUN a chain with 500K steps for p-value estimate:
oldMLE=mle  #use line below instead if you do not wish to re-compute the MLE:
#oldMLE = matrix(scan(file = "mle-OHT.txt"), ncol=4,byrow=TRUE)

runtime500K = system.time( OHT.500Kwalk <- Estimate.p.Value.for.Testing( d.OHT,b.OHT,steps.for.walk=500000, mleMatr=oldMLE) )
 
# p-value estimates after 50K burn-in steps:  
gof.OHT.500K = OHT.500Kwalk[[3]] #these are the chi-square values
pvalues.OHT.500K.50Kburnin = Estimate.p.Value.From.GoFs(gof.OHT.500K , 50000)
pvalues.OHT.500K.50Kburnin[[1]] #gets new p-value estimate
 
 ## GRAPHICS: p-values and GoF histograms 
hist(c(gof.OHT.500K[1],gof.OHT.500K[50001:450000]), main="OHT graph",sub="450K steps (after 50K burn-in steps)",  xlab="Goodness of Fit Statistic", ylab="Frequency", xlim=c(22,36))
abline(v= gof.OHT.500K[1], col="red")   

plot(pvalues.OHT.500K.50Kburnin[[2]], main="OHT",sub="450K steps (after 50K burn-in steps)",  xlab="Length of Walk", ylab="p-value estimates", ylim=c(0,1))




###########################################################################
## Fiber enumeration tests: 
t.15K = system.time(OHT.15K <- Enumerate.Fiber(d.OHT,b.OHT,numsteps=15000))

length(OHT.15K[[1]]) #this is the number of graphs discovered

OHT.15K[[4]] #this is the tv-distance

# Various plots from sampling:
graph.counts.OHT.15K = as.numeric(OHT.15K[[3]])
barplot(graph.counts.OHT.15K , main = "OHT - Graphs", sub="15,000 steps", xlab="Distinct Graphs", ylab="Frequency in walk")

empty.move.counts.OHT.15K = as.numeric(OHT.15K[[5]]) 
barplot(empty.move.counts.OHT.15K, main = "OHT - Empty Moves", sub="15,000 steps", xlab="Distinct Graphs", ylab="Number of empty moves per graph in walk")

distinct.visits.counts.OHT.15K = as.numeric(OHT.15K[[3]])-as.numeric(OHT.15K[[5]]) 
barplot(distinct.visits.counts.OHT.15K, main = "OHT - Number of distinct visits per graph", sub="15,000 steps", xlab="Distinct Graphs", ylab="Number of visits")


## Fiber enumeration for 50,000 steps:
t.50K = system.time(OHT.50K <- Enumerate.Fiber(D.OHT,B.OHT,numsteps=50000))
length(OHT.50K[[1]]) #591  graphs in the fiber
OHT.50K[[4]] #tv-distance

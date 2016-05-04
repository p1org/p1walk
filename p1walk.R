###########################################################################
## Using R version 3.0.0												##
## Using igraph version 0.6.5-1 -WARNING do not use an older version!! 	##
## Authors: Despina Stasi, Sonja Petrovic, Elizabeth Gross.				##
###########################################################################


 ########################################################################
# Estimate.p.Value														#
# Estimate the percentage of networks in the fiber of G= (gdir,gbdir)	#
# that are further from the MLE than G.									#
# Input: directed graph gdir, undirected graph gbidir
#	Optional input:
#		- steps.for.walk
#		- coin:  a fair coin by default. 
#		c[1]=P(directed move); 	c[2]=P(bidirected move); c[3]=P(mixed move).

 ########################################################################
Estimate.p.Value<-function(gdir, gbidir, steps.for.walk=100, coin=c(1/3,1/3,1/3), mle.maxiter = 10000, mle.tol = 1e-03){
	#Error Checking
	if(!is.simple(as.undirected(gdir,mode=c("each")))){
		stop("Reciprocated edges in directed graph or gdir not simple.")
	}
	if(!is.simple(gbidir)){
		stop("gbidir must be a simple graph.")
	}
	if(!is.directed(gdir)){		
		stop("gdir must be a directed graph.")
	}
	if(is.directed(gbidir)){		
		stop("gbidir must be an undirected graph.")
	}
	#Ensure graphs have same number of vertices
	nd = vcount(gdir)
	nb = vcount(gbidir)
	if (nd>nb){
		gbidir = add.vertices(gbidir,nd-nb)	
	}
	else if (nd<nb){
		gdir = add.vertices(gdir,nb-nd)	
	}
	
	mleMatr = Get.MLE(gdir,gbidir, maxiter = mle.maxiter, tol = mle.tol)
	obs.gf = Get.GoF.Statistic(gdir, gbidir, mleMatr)
	if (is.nan(obs.gf)){
		print("NaN error in calculation of GF statistic.")
	}
	if (obs.gf== Inf){print("Error: Infinite GF statistic for this network.")}
	next.network = list(gdir,gbidir)
	count = 0
	for(i in 1: steps.for.walk){
		next.network = Get.Next.Network(next.network[[1]],next.network[[2]], coin)	
		new.gf= Get.GoF.Statistic(next.network[[1]], next.network[[2]], mleMatr)
		# If the GoF statistic for new network is larger than the GoF statistic
		# for the observed network. Note that a badly estimated MLE can give 
		# significant errors in the p-value.
		if (new.gf>=obs.gf){
			count = count +1
		}
	}
	return (count/steps.for.walk)
}
 ########################################################################
# Implements the IPS algorithm for fitting the probability parameters	#
# of the p1 model based on Holland and Leinhardt (1981) page 40			#
# Input:																#
# - network: an nxn adjacency matrix for a directed graph. Used to		#
# extract sufficient statistics: indegrees, outdegrees, number of		#  
# bidirectional edges													#
# - maxiter: maximal number of iterations								#
# - tol: tolerance for declaring convergence (based on the				#
# ell-infinity norm of the difference between observed and fitted		#
# row and columns sums)												#
# Output:															#
# - fit: 4x(n choose 2) vector of estimated probabilities				#
#       (4 for each dyad)												#
# - parameters: alpha, beta, rho (rho[1] - rho[n]), rhoconst,			#
#       theta (see Holland and Leihartdt (1981), page 40).				#
 ########################################################################
p1.ips.general <- function(network, maxiter = 3000, tol=1e-6, alpha = array(0,dim = c(1, length(outdeg))), beta = array(0, dim = c(1, length(outdeg))), rho= array(0, dim = c(1, length(outdeg))), theta = 0, rhoconst=0){
	# outdegrees and indegrees are the row and column sums of the observed network
	outdeg = apply(network, 1, sum)
	indeg = apply(network, 2, sum)
	degree.sum = sum(indeg)
			
	v = length(indeg) #number of vertices
	vchoose2 = v*(v-1)/2	
	if (length(outdeg)!=v){ print("error: outdeg, indegree vectors must be of same dimension")} 
	M=0 			#number of bidirected edges
	for (i in 1:(v-1)){ 
		for (j in i:v) { 
			M=M+network[i,j]*network[j,i] 
		} 
	}	
	#initialize m_ij, a_ij, n_ij from equations (22), (23), (24) on page 40, Holland and Leinhardt 1981 p1 paper - available from JASA
	k=array(0, dim=c(v,v))
	m=array(0, dim=c(v,v))
	a=array(0, dim=c(v,v))
	n=array(0, dim=c(v,v))
	
	for (i in 1:v){
		for(j in 1:v){
			if (i!=j){
				k[i,j] = 1+ exp(theta + alpha[i]+beta[j])+ exp(theta+alpha[j]+beta[i])+ exp(rhoconst+rho[i]+rho[j]+2*theta+alpha[i]+beta[j]+alpha[j]+beta[i])
				n[i,j] = 1/k[i,j] # probability of null edge
				m[i,j] = exp(rhoconst +rho[i]+rho[j]+2*theta + alpha[i] + alpha[j]+beta[i]+beta[j])*n[i,j]
				a[i,j] = exp(theta+alpha[i]+beta[j])*n[i,j]
			}
		}
	}
	
	F = array(0,dim=c(1,v))
	G = array(0,dim=c(1,v))
	H = array(0,dim=c(1,v))
	L = array(0,dim=c(1,v))
	R = array(0,dim=c(v,v))
	
	num.iter=0
	converge=0
	while(!converge && num.iter < maxiter){
		############ Row step ##############
		mi = apply(m, 1, sum)
		ai = apply(a, 1, sum)
		#update F
		### f = apply(f, 1:length(outdeg), function(x) x/(mi[i]+ai[i]))
		for (i in 1:v){
			if (outdeg[i]!=0){
				F[i]=outdeg[i]/(mi[i]+ai[i]) 
			}
		}
		#update K
		K= (v*(v-1) - degree.sum) /(sum(ai)+sum(n))	
		for (i in 1:v){
			for(j in 1:v){
				if (i!=j){
					m[i,j] = m[i,j]*sqrt(F[i]*F[j])
					a[i,j] = a[i,j]*sqrt(F[i]*K)
					n[i,j] = n[i,j]*K
				}
			}
		}
		
		########### Column step ###########
		mj = apply(m, 2, sum)
		aj = apply(a, 2, sum)
		#update G
		### f = apply(f, 1:length(outdeg), function(x) x/(mi[i]+ai[i]))
		for (j in 1:v){
			if (indeg[j]!=0){
				G[j] = indeg[j]/(mj[j]+aj[j]) 
			}				
		}
		#update K
		K= (v*(v-1) - degree.sum) /(sum(aj)+sum(n))	
		for (i in 1:v){
			for(j in 1:v){
				if (i!=j){
					m[i,j] = m[i,j]*sqrt(G[i]*G[j])
					a[i,j] = a[i,j]*sqrt(G[j]*K)
					n[i,j] = n[i,j]*K
				}
			}
		}
		
		########## Mutual step ##########
		if(M!=0){
			m.total = sum(m)
			# Update H
			H=M/(m.total/2)
			# Update L
			L = (vchoose2-M) /(vchoose2-m.total/2)	
			for (i in 1:v){
				for(j in 1:v){
					if (i!=j){
						m[i,j] = m[i,j]*H
						a[i,j] = a[i,j]*L
						n[i,j] = n[i,j]*L
					}
				}
			}
		}
		
		########## Normalizing Step ##########
		for (i in 1:v){
			for(j in 1:v){
				if (i!=j){
					R[i,j] = m[i,j] + a[i,j] + a[j,i] + n[i,j]
					r = 1/R[i,j]
					m[i,j] = m[i,j]*r
					a[i,j] = a[i,j]*r
					n[i,j] = n[i,j]*r
				}
			}
		}		
		######### Estimate Convergence #########
		ai = apply(a, 1, sum)
		mi = apply(m, 1, sum)
		aj = apply(a, 2, sum)
		mj = apply(m, 2, sum)
		pi = ai+mi			# estimate of indegree probabilities
		pj = aj+mj 			# estimate of outdegree probabilites
		Mest = sum(m)/2 		# estimate of number of bidirected edges
		maxDifference = max(max(abs(pi - outdeg)), max(abs(pj - indeg)), abs(Mest-M))
		if (!is.nan(maxDifference)){
			converge = (maxDifference < tol )
		}
		else {
			print("Error: NaN issue in convergence constant calculation." )
			break
		}
		num.iter = num.iter +1
  	}
  	if (num.iter == maxiter){
  		print("Warning: no convergence detected while calculating MLE. Try increasing maxiter.")
  	}
  	fit = numeric(4 * vchoose2)
  	k=1
  	for(i in 1:(v-1)){
    	for(j in (i+1):v){
      		fit[k:(k+3)] = c( n[i,j], a[i,j], a[j,i],m[i,j])
      		k = k+4
    	}
    }
	fit
	out=list()
	out[[1]]=fit
	out[[2]]=n
	out[[3]]=a
	out[[4]]=m
	out[[5]]=Mest
    #out
    return(out)
}

 #######################################################################
# Returns the MLE in the form of an {(n choose 2) X 4} matrix, 		    #
#   with each row representing the vector 							    #
#   (pij(0,0), pij(1,0), pij(0,1), pij(1,1)) with (i,j)th row 		    #
#   appearing in lexicographic order.								    #
 #######################################################################
Get.MLE<-function(gdir, gbidir, maxiter=3000, tol = 1e-06,alpha = array(0, dim = c(1, length(outdeg))), beta = array(0, dim = c(1, length(outdeg))), rho = array(0, dim = c(1, length(outdeg))), theta = 0, rhoconst = 0){
	nd = vcount(gdir)
	nb = vcount(gbidir)
	if (nd>nb){
		gbidir = add.vertices(gbidir,nd-nb)	
	}
	else if (nd<nb){
		gdir = add.vertices(gdir,nb-nd)	
	}

	adjMatr = get.adjacency(gbidir)+get.adjacency(gdir)
	out = p1.ips.general(adjMatr, maxiter, tol)
	mleMatr = t(matrix(out[[1]], nrow = 4))
	return (mleMatr)
}
 ########################################################################
# Get.GoF.Statistic													#
# Estimates difference between current network and MLE					#
# Input:																#
# - gdir:    directed graph 						                  	#
# - gbidir:  bidirected graph 						                  	#
# - mleMatr: the mle or extended mle in a 4x(n choose 2) matrix	format. 	#
 ########################################################################
Get.GoF.Statistic<- function(gdir, gbidir, mleMatr){
	confMatr = Get.Configuration.Matrix(gdir,gbidir)
	diff = confMatr-mleMatr
	gf=0
	for (i in 1:nrow(mleMatr)){
		for (j in 1:ncol(mleMatr)){
			if (!is.nan(diff[i,j])){
				if(diff[i,j]^2 !=0){
#					gf = gf + (diff[i,j]^2)/(mleMatr[i,j]^2)
					gf = gf + (diff[i,j]^2)/(mleMatr[i,j])
					if (gf==Inf) {return(Inf)}
				}
			}
			else print("Warning: NaN occurence in Get.GoF.Statistic")
		}
	}
	return(gf)
}
 #######################################################################
# Returns a 4x(n choose 2) matrix where each row is an indicator        #
# vector for the state of the dyad (i, j) in the order				    #
# i---j, i-->j, i<--j, i<->j										    #
 #######################################################################
Get.Configuration.Matrix<-function(gdir,gbidir){
	num.vertices = max(vcount(gdir), vcount(gbidir))
	nrows = num.vertices*(num.vertices-1)/2
	x = matrix(data=0, nrow = nrows , ncol=4)	
	gdir.vector = as.vector(t(get.edgelist(gdir)))
	gbidir.vector = as.vector(t(get.edgelist(gbidir)))
	if(ecount(gdir)!=0){
		for(k in seq(1,2*ecount(gdir),2)){
			if(gdir.vector[k]<gdir.vector[k+1]){
				i=gdir.vector[k]
				j=gdir.vector[k+1]
				c =2
			}
			else {
				i=gdir.vector[k+1]
				j=gdir.vector[k]
				c = 3
			}
			index = num.vertices*(i-1)-(i-1)*(i)/2+j-i# index of edge (i,j) in configuration matrix
			x[index,c]=1
		}
	}
	if(ecount(gbidir)!=0){
		for(k in seq(1,2*ecount(gbidir),2)){
			if(gbidir.vector[k]<gbidir.vector[k+1]){
				i=gbidir.vector[k]
				j=gbidir.vector[k+1]	
			}
			else {
				i=gbidir.vector[k+1]
				j=gbidir.vector[k]	
			}
			index = num.vertices*(i-1)-(i-1)*(i)/2+j-i # index of edge (i,j) in configuration matrix
			x[index,4]=1
		}
	}
	x[,1]=matrix(1,nrow = nrows, ncol=1)-x[,2]-x[,3]-x[,4]
	return(x)
}
 #######################################################################
# Write.Walk.To.File: 												#
# performs a walk and saves the consecutive networks 	    				#
# in string format to file.	Prints each network in a new line as a     #
# sequence of integers separated by commas. First two integers 	   		#
# signify first edge etc. 										    #
 #######################################################################
Write.Walk.To.File<-function(gdir,gbidir, steps=20,coin=c(1/3,1/3,1/3), filename = "walk.txt"){
	write("====================", filename)
	num.cols = 2*max(ecount(gdir),ecount(gbidir)) #to pass to the write function so that all entries are in one row.
	network = list(gdir,gbidir)
	for (i in 1:steps){
		stepnum = sprintf("step %d",i)
		write(stepnum, filename,append=TRUE)
		write("====================", filename,append=TRUE)
		write("Directed Graph", filename, append=TRUE)
		write(t(get.edgelist(network[[1]])), filename, append=TRUE,ncolumns=num.cols, sep = ", ")
		write("Bidirected Graph", filename, append=TRUE)
		write(t(get.edgelist(network[[2]])), filename, append=TRUE,ncolumns=num.cols, sep = ", ")
		write("====================", filename, append=TRUE)
		network = Get.Next.Network(network[[1]],network[[2]], coin)
	}
}
 #######################################################################
# Write.Network.To.File												    #
# Prints a network in a new line as a sequence of integers separated    #
# by commas. First two integers signify first edge etc. 			    #
 #######################################################################
Write.Network.To.File<-function(gdir,gbidir, filename = "walk.txt"){
	num.cols = 2*max(ecount(gdir),ecount(gbidir)) #to pass to the write function so that all entries are in one row.
	write("Directed Graph", filename, append=TRUE)
	write(t(get.edgelist(gdir)), filename, append=TRUE,ncolumns=num.cols, sep = ", ")
	write("Bidirected Graph", filename, append=TRUE)
	write(t(get.edgelist(gbidir)), filename, append=TRUE,ncolumns=num.cols, sep = ", ")
}
 #######################################################################
# Save.Walk.Plots														#
# This function performs a walk and saves the plots of the 				#
# consecutive networks to file. Currently the best way to view the 		#
# walk (in a mac) is to select all generated files and open in preview 	#
# and use the down arrow to view each one. Should possibly be replaced 	#
# by an animation function sometime in the future.						#
 #######################################################################
Save.Walk.Plots<-function(gdir,gbidir, steps=20,coin=c(1/3,1/3,1/3)){
	network = list(gdir,gbidir)
	for (i in 1:steps){
		network = Get.Next.Network(network[[1]],network[[2]], coin)
		filename = sprintf("FiberWalk%d.png",i)
		png(filename,width=800, height=600,bg="white")
		Plot.Mixed.Graph(network[[1]],network[[2]])	
		dev.off()
	}
}
 #######################################################################
# Plot.Walk 															#
# Performs a walk on the fiber and plots the consecutive networks. 		#
# It does not store consecutive networks.								#
 #######################################################################
Plot.Walk<-function(gdir,gbidir, steps=20,coin=c(1/3,1/3,1/3)){	network = list(gdir,gbidir)
	# Should be replaced by an animation function sometime in the future.
	for (i in 1:steps){
		network = Get.Next.Network(network[[1]],network[[2]], coin)
		Plot.Mixed.Graph(network[[1]],network[[2]])	
	}
}
 ########################################################################
# Plot.Mixed.Graph													#
# Given gdir:   the directed part of a mixed graph,						#
#   and gbidir: the bidirected part of a mixed graph 					#							
# plot the graph.														#
# For the purposes of the p1 model view the undirected edges as 			#
# bidirected.														#
 ########################################################################
Plot.Mixed.Graph<- function(gdir,gbidir){
	# Prints the mixed graph, with the vertices in the shape of a circle
	gmixed = graph(t(get.edgelist(gdir)),n=max(vcount(gdir),vcount(gbidir)))
	E(gmixed)$arrow.mode = 2
	if(ecount(gbidir)!=0){
		g.became.dir = as.directed(gbidir,mode = "arbitrary")
		gmixed = add.edges(gmixed,as.vector(t(get.edgelist(g.became.dir))))
		undirected_edge_idxs = c((ecount(gdir)+1):(ecount(gdir)+ecount(gbidir)))
		E(gmixed)[undirected_edge_idxs]$arrow.mode = 0
	}
	plot(gmixed,layout=layout.circle,vertex.shape="none")
}
 #######################################################################
# Get.Next.Network
#	Input:
#		- d: a directed graph
#		- b: a bidirected graph
#	Optional input:
#		- coin:  a fair coin by default. 
#		c[1]=P(directed move); 	c[2]=P(bidirected move); c[3]=P(mixed move).
# Given a mixed graph G=(d,b) returns new mixed graph G' in the p1 fiber
# F(G) with reciprocation after applying a random Graver basis element 
# The move could be 		#
#	only directed, or only bidirected, or a composite of the two.		#
#   coin is optional input; by default it's "fair": 					#
#	c[1]=P(directed move); 	c[2]=P(bidirected move); c[3]=P(mixed move).#
 #######################################################################
Get.Next.Network <- function(d,b,coin=c(1/3,1/3,1/3)){
  	markov.move = Get.Move(d,b,coin)
	if (!ecount(markov.move[[1]])==0 || !ecount(markov.move[[3]])==0){
		#d minus directed.to.be.removed plus directed.to.be.added:
		new.directed.graph = graph.union(graph.difference(d,markov.move[[1]]),markov.move[[2]])
		##Possible speed up: can we use the following instead of union to speed up? d<-add.edges(d,c(1,2))
		#b minus bidirected.to.be.removed plus bidirected.to.be.added
		new.bidirected.graph = graph.union(graph.difference(b,markov.move[[3]]),markov.move[[4]])
	} 
	else{
		#empty move, graphs unchanged
		new.directed.graph = d
		new.bidirected.graph = b
	}
	return(list(new.directed.graph,new.bidirected.graph))		

}
 #######################################################################
# Get.Move																#
# 	Given a mixed graph G=(d,b)											#
#		with d: directed graph, b: bidirected graph						#
#	returns a random move, not necessarily primitive that is 			#
#	applicable to the observed network G that is guaranteed to move		#
#	to a network in the fiber, including itself. The move could be 		#
#	only directed, or only bidirected, or a composite of the two.		#
#   coin is optional input; by default it's "fair": 					#
#	c[1]=P(directed move); 	c[2]=P(bidirected move); c[3]=P(mixed move).#
 #######################################################################
Get.Move <- function(d,b,coin=c(1/3,1/3,1/3)){
  	# if the coin options do not sum up to 1 exit.
  	if (! sum(coin)==1) {
    	stop("invalid coin")
  	}
	#Generate a random real number between 0 and 1.
  	coin.value = runif(1)
  	# Now just see where coin.value is in relation to the coin vector (a,b,c):
  	# first coin option is : coin \in [0.0, a]:
  	if (coin.value <= coin[1]) {
    	dir.move = Get.Directed.Move(d,b)
    	return(list(dir.move[[1]],dir.move[[2]], graph.empty(vcount(b),directed=FALSE), graph.empty(vcount(b),directed=FALSE)))
  	}
  	# second coin option is: coin \in (a,a+b]:
  	else if (coin[1]<coin.value && coin.value <= coin[1]+coin[2]) {
    	bidir.move = Get.Bidirected.Move(d,b)
    	return(list(graph.empty(vcount(d)),graph.empty(vcount(d)),bidir.move[[1]],bidir.move[[2]]))
    }
    # third coin option is : coin \in (a+b,a+b+c]:
    else if (coin[2]<coin.value) {
    	return(Get.Mixed.Move(d,b))
    }
}
 #######################################################################
# Get.Directed.Move													#
# 	Given a mixed graph G=(d,b)										#
#		with d: directed graph, b: bidirected graph					#
#	returns a random move consisting of directed edges only				#
#	applicable to the observed network G that is guaranteed to move		#
#	to a network in the fiber, including itself. 						#
 #######################################################################
Get.Directed.Move <- function(d,b){
	dir.piece=Get.Directed.Piece(d)
	if (is.null(dir.piece[[1]])){ 
		return(list(graph.empty(vcount(d)),graph.empty(vcount(d))))
	}
	else{
    # Finally, check :
    g.add = graph(dir.piece[[2]])
    g.remove = graph(dir.piece[[1]])
	#(1) edges.to.add makes a simple graph, and has no bidirected edges.
	# This can happen if more than one partitions.
	if (!is.simple(as.undirected(g.add,mode="each"))) #instead of: if (!is.simple(g.add)) 
    	return(list(graph.empty(vcount(d)),graph.empty(vcount(d))))
    #(2) edges.to.add does not intersect d - edges.to.remove in any direction [i.e. no conflicts created!]:
    if (!ecount(graph.intersection(as.undirected(graph.difference(d,g.remove)),as.undirected(g.add)))==0) 
    	return(list(graph.empty(vcount(d)),graph.empty(vcount(d)))) 
    #(3) unordered edges.to.add does not intersect B:
    if (!is.null(b)) {
		if (!ecount(graph.intersection(as.undirected(g.add),b))==0){
        	return(list(graph.empty(vcount(d)),graph.empty(vcount(d))))
        }
    }
    return (list(g.remove,g.add))
  }
}
 #######################################################################
# Get.Bidirected.Piece												#
# Given a mixed graph G=(d,b) d:directed, b: bidirected, returns			#
# an applicable bidirected only Markov move in the form of a list		#
# (g.remove, g.add) where g.remove is a graph containing the edges to 	#
# remove and g.add is a graph containing the edges to add.				#
# The move may be empty.
 #######################################################################
Get.Bidirected.Move <- function(d, b) {
	bidir.piece = Get.Bidirected.Piece(b)
	if (is.null(bidir.piece[[1]])) 
		return(list(graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
	else {
		# Finally, check :
		g.add = graph(bidir.piece[[2]], directed = FALSE)
		g.remove = graph(bidir.piece[[1]], directed = FALSE)
		#(1) edges.to.add makes a simple graph. This can happen if more than one partitions.
		if (!is.simple(g.add))
			return(list(graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
		#(2) edges.to.add does not intersect b-edges.to.remove [i.e. no conflicts created!]:
		# and edges.to.add does not intersect d-edges.to.remove [i.e. no conflicts created!]:
		if (!ecount(graph.intersection(graph.difference(b, g.remove), g.add)) == 0) 
			return(list(graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
		#(3) neither order of edges.to.add intersects D:
		if (!is.null(d)) {
			if (!ecount(graph.intersection(as.directed(g.add), d)) == 0) 
				return(list(graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
		}
		return(list(g.remove, g.add))
	}
}
 #######################################################################
# Get.Mixed.Move													    #
# Given d: directed part of a mixed graph							    #
#   and b: bidirected part of a mixed graph,  						    #
# returns an applicable composite move consisting of a directed piece    #
# and a bidirected piece. Either or both of the pieces can be empty.	    #
 #######################################################################
Get.Mixed.Move <- function(d, b) {
    dir.piece = Get.Directed.Piece(d)
    if (is.null(dir.piece[[1]])){
    	g.add.dir = graph.empty(vcount(d))
    	g.remove.dir = graph.empty(vcount(d))
    }
	else{	
	    g.add.dir = graph(dir.piece[[2]])
    	g.remove.dir = graph(dir.piece[[1]])
    }
	# get bidir.piece too and then:
	bidir.piece=Get.Bidirected.Piece(b)
    if (is.null(bidir.piece[[1]])){
    	g.add.bidir = graph.empty(vcount(b), directed=FALSE)
    	g.remove.bidir = graph.empty(vcount(b), directed=FALSE)   	
    } 
	else{
		g.add.bidir = graph(bidir.piece[[2]],directed=FALSE)
    	g.remove.bidir = graph(bidir.piece[[1]],directed=FALSE)
    }
	# carry out checks:
	#(1) edges.to.add makes a simple graph. This can happen if more than one partitions. We also check that the directed edges to be added do not create new reciprocal edges
	if ( (!is.simple(as.undirected(g.add.dir,mode="each")))   ||  (!is.simple(g.add.bidir)) || 
				(!ecount(graph.intersection(g.add.bidir,as.undirected(g.add.dir)))==0) )
		return(list(graph.empty(vcount(d)),graph.empty(vcount(d)),graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
	#(2) edges.to.add does not intersect b-edges.to.remove [i.e. no conflicts created!]:
	if ((!ecount(graph.intersection(graph.difference(b, g.remove.bidir), g.add.bidir)) == 0) || 
	 			!ecount(graph.intersection(as.undirected(graph.difference(d,g.remove.dir)),as.undirected(g.add.dir)))==0)
	 			######This line causing the bug!! will fix soon
	 	return( list(graph.empty(vcount(d)), graph.empty(vcount(d)), graph.empty(vcount(b),directed=FALSE),graph.empty(vcount(b), directed=FALSE)) )
    #(3i) neither order of g.add.bidir intersects d-g.remove.dir:
    if (!is.null(d)) {
    	if (!ecount(graph.intersection(as.directed(g.add.bidir), graph.difference(d,g.remove.dir))) == 0)
          return(list(graph.empty(vcount(d)),graph.empty(vcount(d)),graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
		}
    #(3ii) unordered g.add.dir does not intersect b-g.remove.bidir:
	if (!is.null(b)) {
      if (!ecount(graph.intersection(as.undirected(g.add.dir), graph.difference(b,g.remove.bidir)))==0)
        return(list(graph.empty(vcount(d)),graph.empty(vcount(d)),graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
      }
    return(list(g.remove.dir, g.add.dir,g.remove.bidir, g.add.bidir)) 
}
#######################################################################
#######################################################################
Get.Directed.Piece <- function(d){
  # d = directed part of G.
  # pick a random subset E of edges of d and randomly shuffle it
  # (i.e., E = random sample from d of random size):
  if (ecount(d)==2) {# avoid unwanted behaviour of sample function
    random.subset.of.d=get.edges(d,1:2) 
    subset.size=2
  }
  else if (ecount(d)>2){
    subset.size = sample(2:ecount(d) ,1) #this is a random integer
    random.edge.indices = sample(1:(ecount(d)),subset.size)
    random.subset.of.d = get.edges(d,random.edge.indices)
  }
  else return(NULL)
  # randomly partition E,
  # and for every part E_i, call Bipartite.Walk(E_i)
  # and merge the edges.to.add_i from each of the partitions into a big set edges.to.add
  number.of.partitions = sample(1:(floor(subset.size/2)), 1)
  # initialize where to store the pieces of the walk:
  edges.to.add = c()
  edges.to.remove = c()
  more.edges = c()
  num.edges.left =subset.size
  s=1 #index
  while(num.edges.left>1) {
    if (num.edges.left==2) k=2 #avoid unwanted behaviour of sample function
    else k = sample(2:num.edges.left,1) #size of current part.
    if (num.edges.left-k == 1) k=k+1 #E's assumption on not leaving out that last edge hanging. 
    more.edges=Bipartite.Walk(random.subset.of.d[s:(s+k-1),])#stupid comma to get the rows!
    if (is.null(more.edges)) return(NULL)
    else edges.to.add = c(edges.to.add,more.edges ) 
    num.edges.left=num.edges.left-k
    s=s+k
  }
  # edges.to.remove has to be in the same format as edges.to.add, so do this:
  if ( !is.null(edges.to.add)) as.vector(t(random.subset.of.d)) -> edges.to.remove
  return(list(edges.to.remove,edges.to.add))
}
#######################################################################
#######################################################################
Get.Bidirected.Piece <- function(b) {
	## THIS function computes bidirected move ONLY without checks for conflicts.
	#this calls Bipartite.Walk but first checks if edges are a matching?
	# Randomly direct the entire bidirected graph and call Get.Directed.Piece

	if (ecount(b) < 2) 
		return(NULL)
	b.directed = as.arbitrary.directed(b)	
	# now cheat: just call directedPiece on this:
	return(Get.Directed.Piece(b.directed))
}
 #######################################################################
# Bipartite.Walk														#
# Given a (randomized) list of edges (edges.to.remove) return a list 	#
# of edges (edges.to.add) that complete an even closed walk by 			#
# connecting the endpoints of successive edges.							#
 #######################################################################
Bipartite.Walk <- function(edges.to.remove) {
	#connect head of (i+1)st edge to tail of ith edge to complete a walk:
	num.edges = nrow(edges.to.remove)
	edges.to.add = c()
	for (i in 1:(num.edges - 1)) {
		edges.to.add = c(edges.to.add, edges.to.remove[i + 1, 1], edges.to.remove[i, 2])
	}
	edges.to.add = c(edges.to.add, edges.to.remove[1, 1], edges.to.remove[num.edges, 
		2])
	# Ensure that edges.to.add form no loops or multiple edges
	if (!is.simple(graph(edges.to.add))) 
		return(NULL)
	return(edges.to.add)
}
#######################################################################
 #######################################################################
# as.arbitrary.directed													#
# Given an undirected graph b, return a directed graph containing an	#
# arbitrarily directed copy of each edge of b.							#
 #######################################################################
as.arbitrary.directed <- function(b) {
	# Create a directed graph out of the edges of b
	b.decr = graph(t(get.edges(b, 1:ecount(b))))
	# Pick a random integer from 0 to #edges in b
	num.edges.to.reverse = sample(0:ecount(b), 1) 
	# Direct the first num.edges.to.reverse edges in one way and the others the other way
	if (num.edges.to.reverse==0) {
		b.directed = b.decr
	} 
	else{
		random.edge.indices = sample(1:ecount(b), num.edges.to.reverse)
		b.subset.decr = graph(t(get.edges(b, random.edge.indices))) #get.edges and get.edgelist direct edges in a different order somehow!
		el = get.edgelist(b.subset.decr, names = FALSE) ## magically swap cols to reverse direction
		b.subset.incr = graph(rbind(el[, 2], el[, 1]))
		# make the directed graph out of:
		# (reversed.edges.of.b being directed in the decreasing order) union (remaining edges in incr.order):    
		b.directed = graph.union(graph.difference(b.decr, b.subset.decr), b.subset.incr)
	}
	return(b.directed)
}
########################################################################
Estimate.p.Value.for.Testing<-function(gdir, gbidir, steps.for.walk=100, mleMatr = NULL, coin=c(1/3,1/3,1/3), mle.maxiter = 10000, mle.tol = 1e-03){
	#Error Checking
	if(!is.simple(as.undirected(gdir,mode=c("each")))){
		stop("Reciprocated edges in directed graph or gdir not simple.")
	}
	if(!is.simple(gbidir)){
		stop("gbidir must be a simple graph.")
	}
	if(!is.directed(gdir)){		
		stop("gdir must be a directed graph.")
	}
	if(is.directed(gbidir)){		
		stop("gbidir must be an undirected graph.")
	}

	nd = vcount(gdir)
	nb = vcount(gbidir)
	if (nd>nb){
		gbidir = add.vertices(gbidir,nd-nb)	
	}
	else if (nd<nb){
		gdir = add.vertices(gdir,nb-nd)	
	}
	# if mleMatr was not given as argument use generate the MLE
	if (is.null(mleMatr)){
		print("Now estimating MLE.")
		mleMatr = Get.MLE(gdir,gbidir, maxiter = mle.maxiter, tol = mle.tol)
		print("MLE estimate completed.")
	}
	# Error Check: inputted mleMatr dimension
	else if ((dim(mleMatr)[[1]]!=nd*(nd-1)/2) || (dim(mleMatr)[[2]]!=4)){
		stop("mleMatr dimension is incorrect.")
	}
	obs.gf = Get.GoF.Statistic(gdir, gbidir, mleMatr)
	if (is.nan(obs.gf)){
		print("NaN error in calculation of GF statistic.")
	}
	if (obs.gf== Inf){print("Error: Infinite GF statistic for this network.")}
	next.network = list(gdir,gbidir)
	count = 0
	int.values=c() # To estimate convergence of count/i to p-value
	gof.values=c(obs.gf) # To record the  goodness of fit statistics for all networks in walk
	for(i in 1: steps.for.walk){
		next.network = Get.Next.Network(next.network[[1]],next.network[[2]], coin)	
		new.gf= Get.GoF.Statistic(next.network[[1]], next.network[[2]], mleMatr)
		# If the GoF statistic for new network is larger than the GoF statistic
		# for the observed network. Note that a badly estimated MLE can give 
		# significant errors in the p-value.
		if (new.gf>=obs.gf){
			count = count +1
		}
		int.values<-c(int.values,count/i)
		gof.values<-c(gof.values,new.gf)
	}
	return (list(count/steps.for.walk,int.values, gof.values, mleMatr))
}
#######################################################################
########################################################################
# Input: 
#		-gofs: list of goodness of fit statistics, with the first one the gof of the observed network
#       -burnsteps: the number of first entries we should ignore when estimating p-value of observed network
# Output:
#		- A list of 
#			- p-value estimate
#			- a list of p-value estimates for each step of the loop
#######################################################################
Estimate.p.Value.From.GoFs<-function(gofs, burnsteps){
	count = 0	# To estimate convergence of count/i to p-value
	p.values = c()
	for(i in (burnsteps +2):length(gofs)){
		# If the GoF statistic for new network is larger than the GoF statistic
 		# for the observed network. Note that a badly estimated MLE can give 
 		# significant errors in the p-value.
 		if (gofs[i]>=gofs[1]){
 			count = count +1
 		}
 		p.values = c(p.values, count/(i-(burnsteps+1)))
 	}
 	return (list(count/(length(gofs)-length(burnsteps)-1), p.values))
}
# Enumerates the fiber, and returns the list of graps visited, directed+bidirected parts,  a vector of counts for each graph, the Total Variation Distance of the walk, and a count of all empty moves made in each graph. #NOTE: Does not currently keep track of empty moves.
Enumerate.Fiber<-function(gdir, gbidir, numsteps=1000, coin = c(1/3,1/3,1/3)){
	counts=list(1)
	empty.move.counts=list(0)	
	network=list(gdir,gbidir)
	
	graphsD=list()
	graphsB=list()
	
####### REPLACE	####### 
#	graphsD[[1]]=gdir
#	graphsB[[1]]=gbidir
#######   WITH  ####### 
	graphsD[[1]]= as.numeric(t(get.edgelist(gdir)))
	graphsB[[1]]= as.numeric(t(get.edgelist(gbidir)))
####### END REPLACE

	numGraphs=1
	
	for (i in 1:numsteps){
		# In case there are errors note that i is the number of steps 
		on.exit(return(list(graphsD, graphsB, counts, TV<-(sum(abs(as.numeric(counts)/i-1/numGraphs)))/2, empty.move.counts, i)))

		flag = FALSE
		empty.move.flag=FALSE
		prev.network = network
		network = Get.Next.Network(network[[1]],network[[2]],coin)
		if (ecount(graph.difference(network[[1]],prev.network[[1]]))==0 && ecount(graph.difference(network[[2]],prev.network[[2]]))==0){
			# new network is same as previous network
			empty.move.flag=TRUE
		}
		for(j in 1:numGraphs){
####### REPLACE	####### 
#			if(ecount(graph.difference(network[[1]],graphsD[[j]]))==0 &&ecount(graph.difference(network[[2]],graphsB[[j]]))==0){
#######   WITH  ####### 
			if (ecount(graph.difference(network[[1]],graph(graphsD[[j]], n=vcount(gdir), directed=TRUE)))==0 && ecount(graph.difference(network[[2]],graph(graphsB[[j]], n=vcount(gbidir), directed=FALSE)))==0){
####### END REPLACE
				# New network was encountered before
				counts[[j]]=counts[[j]]+1
				flag = TRUE
				if (empty.move.flag==TRUE){
					empty.move.counts[[j]]=empty.move.counts[[j]]+1
				}
			} 
		}
		if (!flag){
			# Encountered new graph
			counts = append(counts, list(1))
			empty.move.counts = append(empty.move.counts,list(0))
			#old code - delete: counts[[numGraphs+1]]=counts[[numGraphs+1]]+1
			numGraphs = numGraphs+1
####### REPLACE	####### 
#			graphsD[[numGraphs]]=network[[1]]
#			graphsB[[numGraphs]]=network[[2]]
#######   WITH  ####### 
			graphsD[[numGraphs]]=as.numeric(t(get.edgelist(network[[1]])))
			graphsB[[numGraphs]]=as.numeric(t(get.edgelist(network[[2]])))			
####### END REPLACE
		}
	}
	
	#calculate the TV distance
	TV=(sum(abs(as.numeric(counts)/numsteps-1/numGraphs)))/2
	return (list(graphsD, graphsB, counts, TV, empty.move.counts))
}

Write.Graphs.to.File<-function(graphs, filename){
	for (i in 1:length(graphs)){
		if (is.igraph(graphs[[i]])){		
			num.cols = 2*ecount(graphs[[i]]) #to pass to the write function so that all entries are in one row.
			write(t(get.edgelist(graphs[[i]])), filename, append=TRUE,ncolumns=num.cols, sep = ", ")	
		}
		else{
			num.cols=2*length(graphs[[i]])
			if (length(graphs[[i]])!=0){
				write(graphs[[i]], filename, append=TRUE,ncolumns=num.cols, sep = ", ")			
			}
			else { write("\n",filename, append=TRUE)}
		}
	}
}

Get.GoF.Statistics.From.File<-function(dir.graphs.filename,bidir.graphs.filename,mleMatr){
	print("Not implemented yet.")
}

# Enumerates the fiber for use with large fibers: writes to file the list of graps visited, directed+bidirected parts,  a vector of counts for each graph, the Total Variation Distance of the walk, and a count of all empty moves made in each graph. #NOTE: Does not currently keep track of empty moves.
Enumerate.Fiber.to.File<-function(gdir, gbidir, numsteps=1000, coin = c(1/3,1/3,1/3), filename.extension){
	# METHOD NOT COMPLETE
	print("Buggy code: don't know what is wrong with this method yet!")
	counts=list(1)
	empty.move.counts=list(0)	
	network=list(gdir,gbidir)
	numGraphs=1
	num.cols.d = 2*ecount(gdir)
	num.cols.b = 2*ecount(gbidir)
	write(t(get.edgelist(gdir)), paste(filename.extension, "dir.graphs.txt", sep="."), append=FALSE,ncolumns=num.cols.d, sep = ", ")	
	write(t(get.edgelist(gbidir)), paste(filename.extension, "bidir.graphs.txt", sep="."), append=FALSE,ncolumns=num.cols.b, sep = ", ")	
		
	for (i in 1:numsteps){
#		on.exit(print(paste("Number of steps: ", numsteps, " ------- Number of Graphs discovered: ", numGraphs, "\n total variation distance: ",(sum(abs(counts/numsteps-1/numGraphs)))/2, "\n counts: \n", counts, "\n empty move counts: \n", empty.move.counts)))
		flag = FALSE
		empty.move.flag=FALSE
		prev.network = network
		network = Get.Next.Network(network[[1]],network[[2]],coin)
		if (ecount(graph.difference(network[[1]],prev.network[[1]]))==0 && ecount(graph.difference(network[[2]],prev.network[[2]]))==0){
			empty.move.flag=TRUE 			# new network is same as previous network
		}		
		# OPEN CONNECTION  
		conD = file(paste(filename.extension, "dir.graphs.txt", sep="."), open = "r")
		conB = file(paste(filename.extension, "bidir.graphs.txt", sep="."), open ="r")
		graph.index = 0
		# (WHILE FILE HAS MORE GRAPHS TO READ: READ GRAPH, COMPARE IT TO CURRENT. REST AS BEFORE)
		while ( (length(str.dir.graph <- readLines(conD, n = 1, warn = FALSE)) > 0) && 
		(length(str.bidir.graph <- readLines(conB, n = 1, warn = FALSE)) > 0) && flag==FALSE) {
	#		dir.graphs = readLines(con = paste(filename.extension, "dir.graphs.txt", n = 50)		
	#		bidir.graphs = readLines(con = paste(filename.extension, "bidir.graphs.txt", n = 50)
			graph.index = graph.index+1

			dir.what = unlist(strsplit(str.dir.graph,split=','))
			bidir.what = unlist(strsplit(str.bidir.graph,split=','))
			dir.graph = as.numeric(dir.what)
			bidir.graph = as.numeric(bidir.what)
			d = graph(  dir.graph, n = vcount(gdir), directed = TRUE)
			b = graph(bidir.graph, n = vcount(gbidir), directed = FALSE)

			if (ecount(graph.difference(network[[1]],d)) == 0 && ecount(graph.difference(network[[2]],b)) == 0){
				# Current network was visited before
				counts[[graph.index]]=counts[[graph.index]]+1
				flag = TRUE
				if (empty.move.flag==TRUE){
					empty.move.counts[[graph.index]]=empty.move.counts[[graph.index]]+1
				}
			} 
		}
		# CLOSE CONNECTION  
		close(conD)
		close(conB)
		if (!flag){
			# Encountered new graph
			counts = append(counts, list(1))
			empty.move.counts = append(empty.move.counts,list(0))
			numGraphs = numGraphs+1
			# Write new graph to files 
			### CAUTION: writing to the file line by line can be very time consuming. 
			### Rewrite this to only write to file every so often (every 1000 graphs?)
			write(t(get.edgelist(network[[1]])), paste(filename.extension, "dir.graphs.txt", sep="."), append=TRUE,ncolumns=num.cols.d, sep = ", ")	
			write(t(get.edgelist(network[[2]])), paste(filename.extension, "bidir.graphs.txt", sep="."), append=TRUE,ncolumns=num.cols.b, sep = ", ")	
			}	
		
		# Logging information. Time consuming, but avoids losing all info if the program is interrupted early
		
		# Write counts to file
		write(counts, paste(filename.extension, "distinct.graph.counts.txt", sep="."), append=FALSE, ncolumns=length(counts), sep = ", ")	
		# Write empty move counts to file
		write(empty.move.counts, paste(filename.extension, "empty.move.counts.txt", sep="."), append=FALSE, ncolumns=length(counts), sep = ", ")	

	}
	# calculate the TV distance
	TV=(sum(abs(as.numeric(counts)/numsteps-1/numGraphs)))/2
	write(TV, paste(filename.extension, "tv.distance.txt", sep="."), append=FALSE, ncolumns=length(counts), sep = ", ")	
	return (numGraphs)
}

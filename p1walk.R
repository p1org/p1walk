###########################################################################
## Using R version 3.0.0                                                 ##
## Using igraph version 0.6.5-1 -WARNING do not use an older version!!   ##
## Authors: Despina Stasi, Sonja Petrovic, Elizabeth Gross.              ##
###########################################################################
###########################################################################

Estimate.p.Value<-function(gdir, gbidir, steps.for.walk=100, coin=c(1/3,1/3,1/3), mle.maxiter = 10000, mle.tol = 1e-03){
########################################################################
# 	Input: 
# 		- gdir:    directed graph (represents the directed-only edges in the network)		                  	
# 		- gbidir:  undirected graph (represents the reciprocated edges in the network)			                  	
#	Optional input:
#		- steps.for.walk, integer indicating the length of the walk used for the estimation.
#		- coin:  a numeric vector indicating the probability of each of three types of moves. c[1]=P(directed move); c[2]=P(bidirected move); c[3]=P(mixed move).
#		- mle.maxiter: integer indicating the maximum number of iterations to be used in the IPS algorithm.
#		- mle.tol:	numeric indicating the smallest error to be passed to the IPS algorithm.
#	Output:
#		- numeric: an estimate of the p-value, i.e. the percentage of networks in the fiber of G= (gdir,gbdir) that are further from the MLE than G.
########################################################################
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
		if (new.gf>=obs.gf){
			count = count +1
		}
	}
	return (count/steps.for.walk)
}
#######################################################################
#######################################################################
p1.ips.general <- function(network, maxiter = 3000, tol=1e-6, alpha = array(0,dim = c(1, length(outdeg))), beta = array(0, dim = c(1, length(outdeg))), rho= array(0, dim = c(1, length(outdeg))), theta = 0, rhoconst=0){
########################################################################
# Implements the IPS algorithm for fitting the probability parameters of the p1 model based on Holland and Leinhardt (1981) pp40			
# Input:																
# 		- network: an nxn adjacency matrix for a directed graph. Used to extract sufficient statistics: indegrees, outdegrees, number of bidirectional edges
# 		- maxiter: maximal number of iterations								
# 		- tol: tolerance for declaring convergence (based on the ell-infinity norm of the difference between observed and fitted row and columns sums)	
# Optional Input:
#		- alpha,  beta, rho, theta, rhoconst: initial values for the parameters of the p1 model
# Output:	
#		- a list containing:															
# 			- fit: 4x(n choose 2) vector of estimated probabilities (4 for each dyad)								
#			- n: estimate of null edges
#			- a: estimate of directed edges
#			- m: estimate of reciprocated edges
#			- Mest: estimate of number of reciprocated edges
########################################################################
	outdeg = apply(network, 1, sum)
	indeg = apply(network, 2, sum)
	degree.sum = sum(indeg)
			
	v = length(indeg)
	vchoose2 = v*(v-1)/2	
	if (length(outdeg)!=v){ print("error: outdeg, indegree vectors must be of same dimension")} 
	M=0
	for (i in 1:(v-1)){ 
		for (j in i:v) { 
			M=M+network[i,j]*network[j,i] 
		} 
	}	
	k=array(0, dim=c(v,v))
	m=array(0, dim=c(v,v))
	a=array(0, dim=c(v,v))
	n=array(0, dim=c(v,v))
	
	for (i in 1:v){
		for(j in 1:v){
			if (i!=j){
				k[i,j] = 1+ exp(theta + alpha[i]+beta[j])+ exp(theta+alpha[j]+beta[i])+ exp(rhoconst+rho[i]+rho[j]+2*theta+alpha[i]+beta[j]+alpha[j]+beta[i])
				n[i,j] = 1/k[i,j]
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
		pi = ai+mi
		pj = aj+mj 				
		Mest = sum(m)/2 		
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
	out=list()
	out[[1]]=fit
	out[[2]]=n
	out[[3]]=a
	out[[4]]=m
	out[[5]]=Mest
    return(out)
}
#######################################################################
#######################################################################
Get.MLE<-function(gdir, gbidir, maxiter=3000, tol = 1e-06,alpha = array(0, dim = c(1, length(outdeg))), beta = array(0, dim = c(1, length(outdeg))), rho = array(0, dim = c(1, length(outdeg))), theta = 0, rhoconst = 0){
#######################################################################
# Returns the MLE in the form of an {(n choose 2) X 4} matrix, with each row representing the vector (pij(0,0), pij(1,0), pij(0,1), pij(1,1)) with (i,j)th row appearing in lexicographic order.								    
# Input: 
# 		- gdir:    directed graph (represents the directed-only edges in the network)	
# 		- gbidir:  undirected graph (represents the reciprocated edges in the network)	
# Optional Input:
#		- maxiter: maximum number of iterations to be passed to IPS algorith, 
#		- alpha,  beta, rho, theta, rhoconst: initial values for the parameters of the p1 model
#######################################################################
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
########################################################################
Get.GoF.Statistic<- function(gdir, gbidir, mleMatr){
########################################################################
# Estimates difference between current network and MLE					
# Input:																
# 		- gdir:    directed graph (represents the directed-only edges in the network)	
# 		- gbidir:  undirected graph (represents the reciprocated edges in the network)	
# 		- mleMatr: the mle or extended mle in a 4x(n choose 2) matrix format.
# Output:
#		- numeric: the goodness-of-fit statistic for the given network	
########################################################################
	confMatr = Get.Configuration.Matrix(gdir,gbidir)
	diff = confMatr-mleMatr
	gf=0
	for (i in 1:nrow(mleMatr)){
		for (j in 1:ncol(mleMatr)){
			if (!is.nan(diff[i,j])){
				if(diff[i,j]^2 !=0){
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
#######################################################################
Get.Configuration.Matrix<-function(gdir,gbidir){
#######################################################################
# Input: 
# 		- gdir:    directed graph (represents the directed-only edges in the network)	
# 		- gbidir:  undirected graph (represents the reciprocated edges in the network) 
# Output:
#		- The configuration matrix of the network: a 4x(n choose 2) matrix where each row is an indicator vector for the state of the dyad (i, j) in the order i---j, i-->j, i<--j, i<->j
#######################################################################
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
			index = num.vertices*(i-1)-(i-1)*(i)/2+j-i
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
			index = num.vertices*(i-1)-(i-1)*(i)/2+j-i
			x[index,4]=1
		}
	}
	x[,1]=matrix(1,nrow = nrows, ncol=1)-x[,2]-x[,3]-x[,4]
	return(x)
}
#######################################################################
#######################################################################
Write.Walk.To.File<-function(gdir,gbidir, steps=20,coin=c(1/3,1/3,1/3), filename = "walk.txt"){
#######################################################################
# Input: 
# 		- gdir:    directed graph (represents the directed-only edges in the network)		
# 		- gbidir:  undirected graph (represents the reciprocated edges in the network) 
# Optional Input:
#		- filename: string with name of the file where the walk will be saved	
#		- coin:  a numeric vector indicating the probability of each of three types of moves. c[1]=P(directed move); c[2]=P(bidirected move); c[3]=P(mixed move).
# Performs a walk and saves the consecutive networks in string format to file.	    
#######################################################################
	write("====================", filename)
	num.cols = 2*max(ecount(gdir),ecount(gbidir)) 
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
#######################################################################
Write.Network.To.File<-function(gdir,gbidir, filename = "walk.txt"){
#######################################################################
# Input: 
# 		- gdir:    directed graph (represents the directed-only edges in the network)
# 		- gbidir:  undirected graph (represents the reciprocated edges in the network) 
# Optional Input:
#		- filename: string with name of file to save.	
# Prints a network in a new line as a sequence of integers separated by commas. First two integers signify first edge etc. 			    
#######################################################################
	num.cols = 2*max(ecount(gdir),ecount(gbidir))
	write("Directed Graph", filename, append=TRUE)
	write(t(get.edgelist(gdir)), filename, append=TRUE,ncolumns=num.cols, sep = ", ")
	write("Bidirected Graph", filename, append=TRUE)
	write(t(get.edgelist(gbidir)), filename, append=TRUE,ncolumns=num.cols, sep = ", ")
}
#######################################################################
#######################################################################
Save.Walk.Plots<-function(gdir,gbidir, steps=20,coin=c(1/3,1/3,1/3)){
#######################################################################
# Input: 
# 		- gdir:    directed graph (represents the directed-only edges in the network)	
# 		- gbidir:  undirected graph (represents the reciprocated edges in the network) 
# Optional Input:
#		- steps: 
#		- coin:  a numeric vector indicating the probability of each of three types of moves. c[1]=P(directed move); c[2]=P(bidirected move); c[3]=P(mixed move).
# Performs a walk and saves the plots of the consecutive networks to file. Currently the best way to view the walk (in a mac) is to select all generated files and open in preview and use the down arrow to view each one.
#######################################################################
	network = list(gdir,gbidir)
	for (i in 1:steps){
		network = Get.Next.Network(network[[1]],network[[2]], coin)
		filename = sprintf("FiberWalk%d.png",i)
		png(filename,width=800, height=600,bg="white")
		Plot.Mixed.Graph(network[[1]],network[[2]])	
		dev.off()
	}
}
########################################################################
########################################################################
Plot.Walk<-function(gdir,gbidir, steps=20,coin=c(1/3,1/3,1/3)){	
#######################################################################
# Input: 
# 		- gdir:    directed graph (represents the directed-only edges in the network)
# 		- gbidir:  undirected graph (represents the reciprocated edges in the network) 
# Optional Input:
#		- steps: 
#		- coin:  a numeric vector indicating the probability of each of three types of moves. c[1]=P(directed move); c[2]=P(bidirected move); c[3]=P(mixed move).
# Performs a walk on the fiber and plots the consecutive networks. 		
# It does not store consecutive networks.								
#######################################################################
	network = list(gdir,gbidir)
	for (i in 1:steps){
		network = Get.Next.Network(network[[1]],network[[2]], coin)
		Plot.Mixed.Graph(network[[1]],network[[2]])	
	}
}
########################################################################
########################################################################
Plot.Mixed.Graph<- function(gdir,gbidir){
########################################################################
# Input: 
# 		- gdir:    directed graph (represents the directed-only edges in the network)	
# 		- gbidir:  undirected graph (represents the reciprocated edges in the network) 
# Plots the graph.														
########################################################################
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
#######################################################################
Get.Next.Network <- function(d,b,coin=c(1/3,1/3,1/3)){
#######################################################################
# Input: 
# 		- d:    directed graph (represents the directed-only edges in the network)		
# 		- b:  undirected graph (represents the reciprocated edges in the network) 
# Optional input:
#		- coin:  a numeric vector indicating the probability of each of three types of moves. c[1]=P(directed move); c[2]=P(bidirected move); c[3]=P(mixed move).
# Output:	
#		- A directed network in the p1 fiber in the form:
# 			- directed graph (represents the directed-only edges in the network)		
# 			- undirected graph (represents the reciprocated edges in the network) 
# Applies a random move to the given netowrk and returns the new network in the fiber.		
#######################################################################
  	markov.move = Get.Move(d,b,coin)
	if (!ecount(markov.move[[1]])==0 || !ecount(markov.move[[3]])==0){
		new.directed.graph = graph.union(graph.difference(d,markov.move[[1]]),markov.move[[2]])
		new.bidirected.graph = graph.union(graph.difference(b,markov.move[[3]]),markov.move[[4]])
	} 
	else{
		new.directed.graph = d
		new.bidirected.graph = b
	}
	return(list(new.directed.graph,new.bidirected.graph))		

}
#######################################################################
#######################################################################
Get.Move <- function(d,b,coin=c(1/3,1/3,1/3)){
#######################################################################
# Input: 
# 		- gdir:    directed graph (represents the directed-only edges in the network)
# 		- gbidir:  undirected graph (represents the reciprocated edges in the network) 
# Optional Input:
#		- coin:  a numeric vector indicating the probability of each of three types of moves. c[1]=P(directed move); c[2]=P(bidirected move); c[3]=P(mixed move).
# Output:
#		- a random applicable move. # The move could be composed of directed-only edges, or reciprocated-only edges, or a composite of the two.		
#######################################################################
  	if (! sum(coin)==1) {
    	stop("invalid coin")
  	}
  	coin.value = runif(1)
  	if (coin.value <= coin[1]) {
    	dir.move = Get.Directed.Move(d,b)
    	return(list(dir.move[[1]],dir.move[[2]], graph.empty(vcount(b),directed=FALSE), graph.empty(vcount(b),directed=FALSE)))
  	}
  	else if (coin[1]<coin.value && coin.value <= coin[1]+coin[2]) {
    	bidir.move = Get.Bidirected.Move(d,b)
    	return(list(graph.empty(vcount(d)),graph.empty(vcount(d)),bidir.move[[1]],bidir.move[[2]]))
    }
    else if (coin[2]<coin.value) {
    	return(Get.Mixed.Move(d,b))
    }
}
#######################################################################
#######################################################################
Get.Directed.Move <- function(d,b){
#######################################################################
# Input: 
# 		- d:    directed graph (represents the directed-only edges in the network)				
# 		- b:  undirected graph (represents the reciprocated edges in the network) 
# Output:	
#		- a random applicable move consisting of directed-only edges only applicable to the observed network G that is guaranteed to move to a network in the fiber, including itself. 						
#######################################################################
	dir.piece=Get.Directed.Piece(d)
	if (is.null(dir.piece[[1]])){ 
		return(list(graph.empty(vcount(d)),graph.empty(vcount(d))))
	}
	else{
    g.add = graph(dir.piece[[2]])
    g.remove = graph(dir.piece[[1]])
	if (!is.simple(as.undirected(g.add,mode="each")))
    	return(list(graph.empty(vcount(d)),graph.empty(vcount(d))))
    if (!ecount(graph.intersection(as.undirected(graph.difference(d,g.remove)),as.undirected(g.add)))==0) 
    	return(list(graph.empty(vcount(d)),graph.empty(vcount(d)))) 
    if (!is.null(b)) {
		if (!ecount(graph.intersection(as.undirected(g.add),b))==0){
        	return(list(graph.empty(vcount(d)),graph.empty(vcount(d))))
        }
    }
    return (list(g.remove,g.add))
  }
}
#######################################################################
#######################################################################
Get.Bidirected.Move <- function(d, b) {
#######################################################################
# Input: 
# 		- d: directed graph (represents the directed-only edges in the network)	
# 		- b: undirected graph (represents the reciprocated edges in the network) 
# Output:
#		- a random applicable move consisting of reciprocated edges applicable to the observed network G that is guaranteed to move to a network in the fiber, including itself. 						
#######################################################################
	bidir.piece = Get.Bidirected.Piece(b)
	if (is.null(bidir.piece[[1]])) 
		return(list(graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
	else {
		g.add = graph(bidir.piece[[2]], directed = FALSE)
		g.remove = graph(bidir.piece[[1]], directed = FALSE)
		if (!is.simple(g.add))
			return(list(graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
		if (!ecount(graph.intersection(graph.difference(b, g.remove), g.add)) == 0) 
			return(list(graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
		if (!is.null(d)) {
			if (!ecount(graph.intersection(as.directed(g.add), d)) == 0) 
				return(list(graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
		}
		return(list(g.remove, g.add))
	}
}
#######################################################################
#######################################################################
Get.Mixed.Move <- function(d, b) {
#######################################################################
# Input: 
# 		- d: directed graph (represents the directed-only edges in the network)	
# 		- b: undirected graph (represents the reciprocated edges in the network) 
# Output: 
#		- an applicable composite move possibly consisting of both directed-only edges and reciprocated edges. 	
#######################################################################
    dir.piece = Get.Directed.Piece(d)
    if (is.null(dir.piece[[1]])){
    	g.add.dir = graph.empty(vcount(d))
    	g.remove.dir = graph.empty(vcount(d))
    }
	else{	
	    g.add.dir = graph(dir.piece[[2]])
    	g.remove.dir = graph(dir.piece[[1]])
    }
	bidir.piece=Get.Bidirected.Piece(b)
    if (is.null(bidir.piece[[1]])){
    	g.add.bidir = graph.empty(vcount(b), directed=FALSE)
    	g.remove.bidir = graph.empty(vcount(b), directed=FALSE)   	
    } 
	else{
		g.add.bidir = graph(bidir.piece[[2]],directed=FALSE)
    	g.remove.bidir = graph(bidir.piece[[1]],directed=FALSE)
    }
	if ( (!is.simple(as.undirected(g.add.dir,mode="each")))   ||  (!is.simple(g.add.bidir)) || 
				(!ecount(graph.intersection(g.add.bidir,as.undirected(g.add.dir)))==0) )
		return(list(graph.empty(vcount(d)),graph.empty(vcount(d)),graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
	if ((!ecount(graph.intersection(graph.difference(b, g.remove.bidir), g.add.bidir)) == 0) || 
	 			!ecount(graph.intersection(as.undirected(graph.difference(d,g.remove.dir)),as.undirected(g.add.dir)))==0)
	 	return( list(graph.empty(vcount(d)), graph.empty(vcount(d)), graph.empty(vcount(b),directed=FALSE),graph.empty(vcount(b), directed=FALSE)) )
    if (!is.null(d)) {
    	if (!ecount(graph.intersection(as.directed(g.add.bidir), graph.difference(d,g.remove.dir))) == 0)
          return(list(graph.empty(vcount(d)),graph.empty(vcount(d)),graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
		}
	if (!is.null(b)) {
      if (!ecount(graph.intersection(as.undirected(g.add.dir), graph.difference(b,g.remove.bidir)))==0)
        return(list(graph.empty(vcount(d)),graph.empty(vcount(d)),graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
      }
    return(list(g.remove.dir, g.add.dir,g.remove.bidir, g.add.bidir)) 
}
#######################################################################
#######################################################################
Get.Directed.Piece <- function(d){
#######################################################################
# Input: 
# 		- d: directed graph (represents the directed-only edges in the network)		
# Output: 
#		- a move on the directed-only part of the network. May not be applicable
#######################################################################
  if (ecount(d)==2) {
    random.subset.of.d=get.edges(d,1:2) 
    subset.size=2
  }
  else if (ecount(d)>2){
    subset.size = sample(2:ecount(d) ,1)
    random.edge.indices = sample(1:(ecount(d)),subset.size)
    random.subset.of.d = get.edges(d,random.edge.indices)
  }
  else return(NULL)
  number.of.partitions = sample(1:(floor(subset.size/2)), 1)
  edges.to.add = c()
  edges.to.remove = c()
  more.edges = c()
  num.edges.left =subset.size
  s=1
  while(num.edges.left>1) {
    if (num.edges.left==2) k=2
    else k = sample(2:num.edges.left,1)
    if (num.edges.left-k == 1) k=k+1
    more.edges=Bipartite.Walk(random.subset.of.d[s:(s+k-1),])
    if (is.null(more.edges)) return(NULL)
    else edges.to.add = c(edges.to.add,more.edges ) 
    num.edges.left=num.edges.left-k
    s=s+k
  }
  if ( !is.null(edges.to.add)) as.vector(t(random.subset.of.d)) -> edges.to.remove
  return(list(edges.to.remove,edges.to.add))
}
#######################################################################
#######################################################################
Get.Bidirected.Piece <- function(b) {
#######################################################################
# Input: 
# 		- b: undirected graph (represents the reciprocated edges in the network) 
# Output: 
#		- a move on the reciprocated edges of the network.	
#######################################################################
	if (ecount(b) < 2) 
		return(NULL)
	b.directed = as.arbitrary.directed(b)	
	return(Get.Directed.Piece(b.directed))
}
#######################################################################
#######################################################################
Bipartite.Walk <- function(edges.to.remove) {
#######################################################################
# Input: 
# 		- edges.to.remove: list of edges 
# Output: 
#  		- a list of edges (edges.to.add) that complete an even closed walk by connecting the endpoints of successive edges.							
#######################################################################
	num.edges = nrow(edges.to.remove)
	edges.to.add = c()
	for (i in 1:(num.edges - 1)) {
		edges.to.add = c(edges.to.add, edges.to.remove[i + 1, 1], edges.to.remove[i, 2])
	}
	edges.to.add = c(edges.to.add, edges.to.remove[1, 1], edges.to.remove[num.edges, 
		2])
	if (!is.simple(graph(edges.to.add))) 
		return(NULL)
	return(edges.to.add)
}
#######################################################################
#######################################################################
as.arbitrary.directed <- function(b) {
#######################################################################
# Input:
#		- b: an undirected graph
# Output:
#		- a directed graph containing an arbitrarily directed copy of each edge of b.	 
#######################################################################
	b.decr = graph(t(get.edges(b, 1:ecount(b))))
	num.edges.to.reverse = sample(0:ecount(b), 1) 
	if (num.edges.to.reverse==0) {
		b.directed = b.decr
	} 
	else{
		random.edge.indices = sample(1:ecount(b), num.edges.to.reverse)
		b.subset.decr = graph(t(get.edges(b, random.edge.indices)))
		el = get.edgelist(b.subset.decr, names = FALSE)
		b.subset.incr = graph(rbind(el[, 2], el[, 1]))
		b.directed = graph.union(graph.difference(b.decr, b.subset.decr), b.subset.incr)
	}
	return(b.directed)
}
#######################################################################
#######################################################################
Estimate.p.Value.for.Testing<-function(gdir, gbidir, steps.for.walk=100, mleMatr = NULL, coin=c(1/3,1/3,1/3), mle.maxiter = 10000, mle.tol = 1e-03){
########################################################################
# Input: 
#		-gdir:   directed graphs
#       -gbidir: bidirected-reciprocated 
# Output:
#		- numeric: p-value estimate
#		- numeric list: the p-value estimates for each step of the loop
#######################################################################
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
	if (is.null(mleMatr)){
		print("Now estimating MLE.")
		mleMatr = Get.MLE(gdir,gbidir, maxiter = mle.maxiter, tol = mle.tol)
		print("MLE estimate completed.")
	}
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
	int.values=c()
	gof.values=c(obs.gf)
	for(i in 1: steps.for.walk){
		next.network = Get.Next.Network(next.network[[1]],next.network[[2]], coin)	
		new.gf= Get.GoF.Statistic(next.network[[1]], next.network[[2]], mleMatr)
		if (new.gf>=obs.gf){
			count = count +1
		}
		int.values<-c(int.values,count/i)
		gof.values<-c(gof.values,new.gf)
	}
	return (list(count/steps.for.walk,int.values, gof.values, mleMatr))
}
########################################################################
########################################################################
Estimate.p.Value.From.GoFs<-function(gofs, burnsteps){
########################################################################
# Input: 
#		- gofs: list of goodness of fit statistics, with the first one the gof of the observed network
#       - burnsteps: the number of first entries we should ignore when estimating p-value of observed network
# Output:
#		- numeric: p-value estimate
#		- numeric list: the p-value estimates for each step of the loop
#######################################################################
	count = 0
	p.values = c()
	for(i in (burnsteps +2):length(gofs)){
 		if (gofs[i]>=gofs[1]){
 			count = count +1
 		}
 		p.values = c(p.values, count/(i-(burnsteps+1)))
 	}
 	return (list(count/(length(gofs)-length(burnsteps)-1), p.values))
}
#######################################################################
#######################################################################
Enumerate.Fiber<-function(gdir, gbidir, numsteps=1000, coin = c(1/3,1/3,1/3)){
#######################################################################
# Input: 
#		- gdir
#		- gbidir
# Optional Input:
#		- numsteps: the number of steps in the walk
#		- coin: a numeric vector of length 3 containing the probabilies of each of three types of moves to be passed to the getNextNetwork metod.
# Output:
# 		- graphsD: a list of directed graphs representing the directed-only part of the graphs encountered in the fiber
#		- graphsB: a list of undirected graphs representing the reciprocated part of the graphs encountered in the fiber
#		- counts: a numeric vector counting the number of times each graph in the fiber was encounterd during the talk
#		- TV: an estimate of the total variation distance of the walk
#		- empty.move.counts: a numeric vector of counts of empty moves made in each graph encountered in the fiber
#######################################################################
	counts=list(1)
	empty.move.counts=list(0)	
	network=list(gdir,gbidir)
	graphsD=list()
	graphsB=list()
	graphsD[[1]]= as.numeric(t(get.edgelist(gdir)))
	graphsB[[1]]= as.numeric(t(get.edgelist(gbidir)))

	numGraphs=1
	
	for (i in 1:numsteps){
		on.exit(return(list(graphsD, graphsB, counts, TV<-(sum(abs(as.numeric(counts)/i-1/numGraphs)))/2, empty.move.counts, i)))

		flag = FALSE
		empty.move.flag=FALSE
		prev.network = network
		network = Get.Next.Network(network[[1]],network[[2]],coin)
		if (ecount(graph.difference(network[[1]],prev.network[[1]]))==0 && ecount(graph.difference(network[[2]],prev.network[[2]]))==0){
			empty.move.flag=TRUE
		}
		for(j in 1:numGraphs){
			if (ecount(graph.difference(network[[1]],graph(graphsD[[j]], n=vcount(gdir), directed=TRUE)))==0 && ecount(graph.difference(network[[2]],graph(graphsB[[j]], n=vcount(gbidir), directed=FALSE)))==0){
				counts[[j]]=counts[[j]]+1
				flag = TRUE
				if (empty.move.flag==TRUE){
					empty.move.counts[[j]]=empty.move.counts[[j]]+1
				}
			} 
		}
		if (!flag){
			counts = append(counts, list(1))
			empty.move.counts = append(empty.move.counts,list(0))
			numGraphs = numGraphs+1
			graphsD[[numGraphs]]=as.numeric(t(get.edgelist(network[[1]])))
			graphsB[[numGraphs]]=as.numeric(t(get.edgelist(network[[2]])))			
		}
	}
	TV=(sum(abs(as.numeric(counts)/numsteps-1/numGraphs)))/2
	return (list(graphsD, graphsB, counts, TV, empty.move.counts))
}

#######################################################################
#######################################################################
Write.Graphs.to.File<-function(graphs, filename){
########################################################################
# Input: 
#		- graphs: a list of graphs
#		- filename: a string containing the name of the file where the graphs will be written
# Writes to file a list of graphs.
#######################################################################
	for (i in 1:length(graphs)){
		if (is.igraph(graphs[[i]])){		
			num.cols = 2*ecount(graphs[[i]])
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
#######################################################################
#######################################################################

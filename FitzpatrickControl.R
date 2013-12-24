#######################################
#FITZPATRICK CONTROL FOR PAIRWISE DATA
#######################################
#Getting Node weighted avergaes for pairwise trait where direction matters and data may be missing
#Yaniv Brandvain (w. help from Mike May) sept 27 2012
library(ape)
	getDescendants<-function(tree,node,curr=NULL){	#note modified from liam revell's blogpost
		#takes tree and node and gets descendents (including focal node)
	  if(is.null(curr)) curr<-vector()
	  daughters<-tree$edge[which(tree$edge[,1]==node),2]
	  curr<-c(curr,daughters)
	  w<-which(daughters>=length(tree$tip))
	  if(length(w)>0) for(i in 1:length(w)) 
	    curr<-getDescendants(tree,daughters[w[i]],curr)
	   if(length(curr)==0){return(node)}
	  return(c(node,curr))
	}
#
	GetNWA = function(dats,tmp.tree2){
		#Get Node weighted average of dats on the subtree tmp.tree2
		tmp = data.frame(node = c(1:length(tmp.tree2$tip)), dats = dats[tmp.tree2$tip] )   #Fills in the tip data
		tmp = rbind(tmp,data.frame(node= as.numeric(names(sort(branching.times(tmp.tree2)))), dats=NA))  #makes NA internal node data, orderd by node age
		for(i in (length(tmp.tree2$tip)+1):nrow(tmp)  ){	#loops acrossing missing node data and gets average vals for descendents
			desc = tmp.tree2$edge[tmp.tree2$edge[,1]%in%tmp$node[i],2]
			tmp$dats[i] = with(tmp,mean(dats[node%in%desc],na.rm=TRUE))
		}
		return(rev(tmp$dats)[1]) #returns the value for the focal (i.e. deepest node)
	}
#
	FormatRunNWA = function(this.data,tmp.tree){ #formating our data for nwa
		this = this.data[!is.na(this.data)]	#only take non-missing data
		if(length(this)<3) return(ifelse(length(this)==0,NA,mean(unlist(this),na.rm=TRUE)) ) #return NA if all is missing, mean if we have 1 or 2 obs
		GetNWA(this, drop.tip(tmp.tree, which(!tmp.tree$tip%in%names(this)) ) )   #otherwise send to node weighted averge function
	}
#
	NodeComp = function(NODE,tmp.tree,data.matrix){
		immediate.descend = tmp.tree$edge[tmp.tree$edge[,1]== NODE,2]	#the two descendents of NODE
		all.descend = lapply(immediate.descend, function(X){getDescendants(tmp.tree,X)}) #find all descendents for each immediate descendent node
		Bs  = tmp.tree$tip[ all.descend[[2]][all.descend[[2]] <= length(tmp.tree$tip)] ]  #all B get species 
		A.names = tmp.tree$tip[all.descend[[1]][all.descend[[1]] <= length(tmp.tree$tip)]]
		a.NWAs = sapply(A.names ,function(spA){ # loop across all descents of one immediate node, get nwa for each
			tmp.data = data.matrix[spA,Bs]
			to.use = !is.na(data.matrix[spA,Bs])
			tmp2 = tmp.data[to.use]
			names(tmp2) = names(tmp.data)[to.use]
			return(FormatRunNWA(tmp2,tmp.tree))
		})
		names(a.NWAs) = tmp.tree$tip[all.descend[[1]][all.descend[[1]] <= length(tmp.tree$tip)]]
		return(FormatRunNWA(a.NWAs,tmp.tree))
	}
#
	FitzControl = function(tmp.tree,data.matrix){
	 	tmp.tree = drop.tip(tmp.tree,which(!tmp.tree$tip.label%in%rownames(data.matrix))) #get rid of tips without any data
		int.nodes = unique(tmp.tree$edge[,1]) #get all internal nodes
		NWA = t(sapply(int.nodes,function(NODE){c(NODE, NodeComp(NODE,tmp.tree, data.matrix) )  } ) ) #node weighted avg (nwa) for each internal node
		data.frame( node = NWA[,1] , age = branching.times(tmp.tree), node.avg = NWA[,2])
	}
#
NodeWeightedAverage = function(tmp.tree,data.matrix){
	AXB = FitzControl(tmp.tree,data.matrix)
	#if(!is.ultrametric(tmp.tree)){warning("Tree is not ultrametric, don't believe ages, do believe node weighted averages")}
 	tmp.tree = drop.tip(tmp.tree,which(!tmp.tree$tip.label%in%rownames(data.matrix)))
	nodes = as.list(AXB$node);names(nodes)= nodes
	#NOW THE GOAL IS TO GET THE DISTANCE BETWEEN NODES
	combos = lapply(nodes,function(NODE){
		immediate.descendents = tmp.tree$edge[tmp.tree$edge[,1]==NODE,2]
		return( lapply(immediate.descendents,function(INTERNAL_NODE){
			all.descendants = getDescendants(tmp.tree,INTERNAL_NODE)
			tmp.tree$tip[all.descendants[  all.descendants <= length(tmp.tree$tip)]]
		})	)
	})
	return(list(AXB=AXB,combos=combos))
	}
#
#Example usage: NodeWeightedAverage(toad.tree,cross.matrix)

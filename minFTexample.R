#MINIMAL EXAMPLE FOR FT CORRECTION

library(ape)
source('~/FitzTurr/FitzpatrickControl.R', chdir = TRUE) # note make sure we are in the correct directory
my.tree = read.tree(text="(D:0.6,(C:0.4,(B:0.1,A:0.1):0.3):.2);") #CREATES A SIMPLE TREE



my.xmat = data.frame(matrix(c(NA,NA,NA,NA,.8,NA,NA,NA,.5,.2,NA,NA,.1, .3,0,NA),nrow=4,byrow=TRUE)) 
	rownames(my.xmat) = colnames(my.xmat) = LETTERS[1:4]
	# CREATES A SIMPLE CROSSING FILE. NOTE WE ASSUME THE MATRX IS SYMMETRIC AND USE ONLY THE LOWER DIAGNOL
	# AN UPDATED SCRIPT WHICH HANDKES ASSYMETRIC MATRICES WILL BE ADDED SHORTLY
	
NodeWeightedAverage(tmp.tree=my.tree,data.matrix=my.xmat)$AXB


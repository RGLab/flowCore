## last changed November 2, 2006
##   I abstracted out the code for checking whether or not the parent exists

##getParentSet = function(parent,ncells) {
## 	if(missing(parent)) {
   		## if the parent is missing assume that it is all cells
##		parentSet = rep(1,ncells)
## 	}
##	else {  	
##		if(class(parent) != "filterResult") {
##			stop("parent must be of class filterResult.")
##		}
##		else {
##			if(is.null(parent@subSet)) {
##				stop("The filterResult parent must have a subSet slot.")
##			}
##			else {
##				parentSet = parent@subSet
##			}  	
##		}
##	}
##	parentSet
##}
 

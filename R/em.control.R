em.control <- 
function(maxit=500,tol=1e-8,crit="relative",random.start=TRUE,classification=c("soft","hard")) {
    classification <- match.arg(classification)
	return(list(maxit=maxit,tol=tol,crit=crit,random.start=random.start,classification=classification))}

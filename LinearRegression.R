#this function accepts file to be use that contains data in the form of [X<separator>Y] data format, user of the function can
#pass the separator as used in the passed file.
#for example if data in the file contains 1.2,3.56 then separator can be ','
#if no separator is passed than white-space is taken as default
###ASSUMPTIONS:
#the function assumes that data are cleaned form and it must contain valid data such tat log(of y) is valid for all possible values.
#this function draws graph for log(y) and x, where first column data in the file is x and 2nd column data is y, whose log is yet to be calculated in the file
#here it is also assumed that 2nd column data is dependent variable and 1st column data is in-dependent cariable

linearRegression <- function(filePath,sep=""){

	#reading data forms , assumed cleaned data
	XYPoints = read.table(filePath, col.names=c("x", "y"),fill=FALSE,strip.white=TRUE,sep=sep);
	
	#getting y and x values
	y<-log(XYPoints[,2],10)
	x<-XYPoints[,1]
	
	myfittingLine <- lm(y ~ x) # log[10](y) and x are used for displaying 

	#required for displaying grid lines
	minX <- min(XYPoints[,1])
	minY <- min(XYPoints[,2])
	maxX <- max(XYPoints[,1])
	maxY <- max(XYPoints[,2])
	
	#required for displaying the intercept and slope details
	interceptCoord <-coef(myfittingLine)
	yIntercept <- interceptCoord['(Intercept)']
	xInterCord <- interceptCoord['x']
	
	# x and y axis label and main label
	ylab <- expression(paste({log[10](l)},"-axis"))
	xlab <- expression(paste(delta,"-axis"))
	mainlab <- expression({paste({log[10](l)}," vs ")*paste({paste(delta," Plot")})})
	
	sub <- paste("Slope: ",toString(xInterCord)," , Y-Intercept :",toString(yIntercept),collapse=" ")
	 
	#plotting only points of raw data
	plot(x,y, type = "p", xlab = xlab, ylab = ylab,main=mainlab, col="blue",sub=sub)

	
	
	#displaying the horizontal and vertical line for (0,0)
	if(minY<1)
	 abline(h = 1, v = (maxX+minX)/2, col = "green")
	
	#plotting grids on line
	abline(h = minY:maxY, v = minX:maxX, col = "lightgray", lty = 3)
	
	#plotting the linear regression line
	abline(myfittingLine,col="red")
	
	#displaying the intercept coordinate
	text(0,yIntercept, "Y-Intercept", col = 2)
 
	#summary of regression line is displayed below
	summary(myfittingLine)
 }

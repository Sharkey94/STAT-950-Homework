# contour plots -> three ways
library(Rcpp)
library(optimization)


# Himmelblau's function
himm <- function(x){(-1)*((x[1]**2 + x[2] - 11)**2 + (x[1] + x[2]**2 -7)**2)}

out_sa <- optim_sa(fun = himm, start = c(runif(2, min = -1, max = 1)),
                   trace = TRUE, lower = c(-4, -4) ,upper=c(4,4),
                   control = list(t0 = 1000, nlimit = 1500,r = 0.8))

# Examples for optimization results via 'Simulated Annealing' method.
plot(out_sa)
plot(out_sa, type = "contour")

library(ContourFunctions)
f1 <- function(x){(-1)*((x[1]**2 + x[2] - 11)**2 + (x[1] + x[2]**2 -7)**2)}
cf_func(f1, xlim = c(-5, 5), ylim = c(-5, 5))


# surface plot
f <- function(x, y) { (-1)*((x^2+y-11)^2 + (x+y^2-7)^2) }
x <- seq(-6, 6, length = 100)
y <- x
z <- outer(x, y, f)
library(plotly)
plot_ly(x = x, y = y, z = ~z) %>% add_surface()

# use Newton's method to obtain each of the optima
himmelblauR=function(x,der=0){
  part1=x[1]^2+x[2]-11
  part2=x[1]+x[2]^2-7
  value=-(part1 ^2+part2 ^2)
  if (der==0) return(value)
  der1=matrix(NA, 2 , 1 )
  der1 [1]=-2*(2*part1*x [1]+ part2 )
  der1 [2]=-2*( part1+2*part2*x[2] )
  if ( der==1) return ( list ( value=value , der1=der1 ) )
  der2=matrix (NA, 2 , 2 )
  der2 [1 ,1]=-(8*x [1]^2+4*part1+2)
  der2 [1 ,2]=-(4*x [1]+4*x [2])
  der2 [2 ,1]= der2[1,2]
  der2 [2 ,2]=-(8*x [2]^2+4*part2+2)
  return( list ( value=value , der1=der1 , der2=der2 ) )
}

newtonR=function(f,xInit,maxIt=20,relConvCrit =1e-10,.. ){
  p=length(xInit)
  results=matrix (NA, maxIt , p+2)
  colnames (results)=c("value" , paste( "x" , 1:p , sep="" )," Conv " )
  xCurrent=xInit
  for ( t in 1:maxIt ){
    evalF=f(xCurrent,der=2)
    results[t,"value"]=evalF$value
    results[ t , 1+( 1:p ) ]=xCurrent
    xNext=xCurrent-solve( evalF$der2 , evalF$der1 )
    Conv=sqrt(crossprod(xNext-xCurrent ) ) / (sqrt( crossprod( xCurrent) )+relConvCrit )
    results [ t , " Conv "]=Conv
    if ( Conv < relConvCrit ) break
    xCurrent=xNext
  }
  return ( list ( x=xNext , value=f(xNext) , convergence=(Conv < relConvCrit ) ,
                  results=results[ 1:t,] ) )
}

# four minima
# first: (3,2)
ans1=newtonR( himmelblauR , c(-2,-2) , relConvCrit = 1e-14)
ans1
himmelblauR(c(ans1$x) , 2 )

# second: (-3.779310,-3.283186) -> Do I need value=0 if converge=0?
ans2=newtonR(himmelblauR , c(-4,-4) , relConvCrit= 1e-14)
ans2
himmelblauR(c(ans2$x) , 2 )

# third: (-2.805118,3.131313) -> Do I need value=0 if converge=0?
ans3=newtonR(himmelblauR , c(-4,4) , relConvCrit= 1e-14)
ans3
himmelblauR(c(ans3$x) , 2 )

# fourth: (3.584428,-1.848126)
ans4=newtonR(himmelblauR , c(3,-2) , relConvCrit= 1e-14)
ans4
himmelblauR(c(ans4$x) , 2 )

# maxima: (-0.270845,-0.923039)
ans5=newtonR(himmelblauR , c(0,0) , relConvCrit= 1e-14)
ans5
himmelblauR(c(ans5$x) , 2 )

##################################################################################
# this program checks to see what the optima are
f <- function(x1,y1) {(-1)*((x1^2 + y1 - 11)^2 + (x1 + y1^2 - 7)^2)}
x <- seq(-4.5,4.5,by=.2)
y <- seq(-4.5,4.5,by=.2)
z <- outer(x,y,f)
persp(x,y,z,phi=-45,theta=45,col="yellow",shade=.65 ,ticktype="detailed")
f <- function(x) (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2

optim(c(-4,-4), f)$par
optim(c(2,-2), f)$par
optim(c(2,2), f)$par
optim(c(-4,4),f)$par



##############part c######
####part c

#part c
set.seed(9102021)

l1=function(x1) {-8*x1^4+168*x1^2+80/3*x1-18224/15}
l2=function(x2) {-8*x2^4+104*x1^2+272/3*x2-4368/5}

l1.p=function(x1) {-32*x1^3+336*x1+80/3}
l2.p=function(x2) {-32*x2^3+208*x1+272/3}

#l1.p2=function(x1) {-96*x1^2+336}
#l2.p2=function(x2) {-96*x2^2+208}

## newton method
## INITIAL VALUES
x1 = 2
x2 = -2
x1.old=x1
x1 = x1 - l1(x1)/l1.p(x1)
x2.old=x2
x2 = x2 - l2(x2)/l2.p(x2)

## MAIN
while (abs((x1-x1.old)/x1)<1e-10) {
  x1.old=x1
  x1 = x1 - l1(x1)/l1.p(x1)
  iter=iter+1
}
iter=0
while (abs((x2-x2.old)/x2)<1e-10) {
  x2.old=x2
  x2 = x2 - l2(x2)/l2.p(x2)
  iter=iter+1
}
x1
x2



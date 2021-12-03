###stat 950 hw1 q2
set.seed(08312021)
lambda=0.2
n=100
exp_rvs <- rexp(n=100,rate=lambda)
#mean(exp_rvs)

### likelihood
L=function(x,rvs) {x^n*exp(-x*sum(rvs))}
### loglikelihood
l=function(x,rvs){n*log(x)-x*sum(rvs)}
l.prime=function(x,rvs){n/x-sum(rvs)}
l.2prime=function(x){-n/x^2}

##bisection
## INITIAL VALUES
a = 0
b = 2
x = a+(b-a)/2
itr = 40
## MAIN
for (i in 1:itr){
  if (l.prime(a,exp_rvs)*l.prime(x,exp_rvs) < 0) {b = x}
  else {a = x}
  x = a+(b-a)/2
}

x		# FINAL ESTIMATE
1/x # FINAL ESTIMATE for mean ###4.91355


#95% asymptotic Confidence Interval for mean
ci.lower = x-(1.96*(sqrt(-1/l.2prime(x))))
ci.mean.upper=1/ci.lower
ci.mean.upper
ci.upper = x+(1.96*(sqrt(-1/l.2prime(x))))
ci.mean.lower=1/ci.upper
ci.mean.lower
###(4.11,6.11)

## newton method
## INITIAL VALUES
x = 0.001
## MAIN
for(i in 1:itr){x = x - l.prime(x,exp_rvs)/l.2prime(x)}

## OUTPUT
x		# FINAL ESTIMATE
1/x # FINAL ESTIMATE for mean  ###4.91355
#95% asymptotic Confidence Interval for mean
ci.lower = x-(1.96*(sqrt(-1/l.2prime(x))))
ci.mean.upper=1/ci.lower
ci.mean.upper
ci.upper = x+(1.96*(sqrt(-1/l.2prime(x))))
ci.mean.lower=1/ci.upper
ci.mean.lower
###(4.11,6.11)


###secant method
###Input:
0             #First initial guess at the desired root
2             #Second initial guess at the desired root
1e-6          #Error Tolerance
Func          #Pointer to function whose root is desired
40           #Maximum number of iterations
Answer        #Best estimate obtained for desired root

##Initialize:
X1 = 0.001
X2 = 0.01
Y1 = l(X1,exp_rvs )
Y2 = l(X2,exp_rvs )


for (i in 1:itr){
  Answer = X2 - Y2 * (X2 - X1) / (Y2 - Y1)
  Y3 = l(Answer,exp_rvs )
  #if ( abs(Y3) < 1e-6 || abs(Answer - X2) < 1e-6 ) break;

  X1 = X2
  Y1 = Y2
  X2 = Answer
  Y2 = Y3

}

Answer



sec <- function(x){
  
  x1 = x
  x2 = x/10 # I have changed the initial values so that x1 is the larger of the two initial estimates.
  
  while(abs(x2 - x1) > 0.0000001){ # This is how I'm trying to converge the points.
    
    # The secant function to determine a new x:
    
    x_new = x2 - l(x2,exp_rvs)*(x2 - x1)/(l(x2,exp_rvs) - l(x1,exp_rvs))
    
    
    if( abs(x2 - x_new) > abs(x_new - x1)){
      
      x2 = x_new
      
    }else{
      
      x1 = x_new
      
    }
    
  }
  
  m_dot = x_new
  m_dot
  
}
sec(0.01)

pre=seq(0.1,10,0.1)
res=l(pre,exp_rvs)
plot(pre,res)

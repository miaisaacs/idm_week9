# This script provides an example of alternative coding for the HIV model 
# Can use similar coding when there are many similar groups in the model

library(deSolve)

# Alternative coding for the HIV model
HIV5riskGrs_alt <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Alternative coding - align all variables and equations
    # initialize the state variables using the same order passed in 
    S = state[1:num_grp]
    I = state[num_grp + (1:num_grp)]
    A = state[num_grp * 2 + (1:num_grp)]
    cumA = state[num_grp * 3 + (1:num_grp)]
    
    ODES <- vector(length = 4*num_grp) # place holder for the equations: 4 state variables (S, I, A, cumA)
    
    # make sure the order of ode eqns. and state variables are the same!
    for (i in 1:num_grp){ 
      infection = sum(BETA[i,] * I) * S[i] # new infections
      ODES[i] = NU[i] - infection - mu * S[i]  #dS/dt
      ODES[num_grp + i] = infection - mu * I[i] - gamma * I[i] #dI/dt
      ODES[num_grp * 2 + i] = d * gamma * I[i] - mu * A[i] - m * A[i]  #dA/dt
      ODES[num_grp * 3 + i] =  d * gamma * I[i] # dcumA/d: cumulative incidence 
    }
    list(ODES)
  })
}

num_grp = 5 # number of risk groups
mu=.0312;  # in years
N=c(.06,.31,.52,.08,.03); # population in each group
NU=mu*N; # recruits to each group
gamma=.2; m=1; d=.3; # in years
beta=.0217; # scaling factor
BETA=matrix(c(rep(0,5),
              0,.65,2.15,12.9,21.5,
              0,2.15,7.17,43.1,71.8,
              0,12.9,43.1,258,431,
              0,21.5,71.8,431,718),5,5,byrow=T)*beta
I0=c(0,0,0,0,1e-5);
S0=N-I0;
A0=cumA0=rep(0,5);

# CAUTION: NOTE HERE WE CODE THE PARAMETER SET AS A LIST, SO VARIABLE NAMES ARE PRESERVED 
parameters=list(BETA=BETA,mu=mu,gamma=gamma,m=m,d=d,
             NU=NU) # use list here for vectors etc.

state=c(S=S0, I=I0, A=A0, cumA=cumA0) # only numeric variables are acceptable for state
times=0:100;

simHIV_alt=ode(y=state,times=times,func=HIV5riskGrs_alt,parms=parameters);


# SAMPLE CODE FOR PROCESSING MODEL OUTPUT FOR EACH GROUP AND AGGREGATE ALL GROUPS
# incidence per 1000 for each group
Incidence_alt=(simHIV_alt[-1,paste0('cumA', 1:num_grp)]-
             simHIV_alt[-length(times),paste0('cumA', 1:num_grp)])/
  matrix(N,100,5,byrow=T)*1000
# total incidence per 1000 population
totInci_alt=rowSums((simHIV_alt[-1,paste0('cumA', 1:num_grp)]-
                   simHIV_alt[-length(times),paste0('cumA', 1:num_grp)]))*1e3

# susceptibility for each group
S_alt=simHIV_alt[,paste0('S', 1:num_grp)]/matrix(N,101,5,byrow=T)*100
# susceptibility combining all groups
totS_alt=rowSums(simHIV_alt[,paste0('S', 1:num_grp)])*100

## plot results
par(mfrow=c(1,1),cex=1.2,mar=c(3,3,1,1),mgp=c(1.8,.5,0))
matplot(Incidence_alt[,2:5],type='l',lty=2:5,col='black',lwd=1.5,ylim=c(0,35),
        ylab='Incidence AIDS per year per 1000',xlab='Time (years)')
lines(totInci_alt,lwd=2)
legend('topright',leg=c('1-5','6-50','51-100','100+','total'),lty=c(2:5,1),
       col='black',lwd=c(rep(1.5,4),2),cex=.9,bty='n')


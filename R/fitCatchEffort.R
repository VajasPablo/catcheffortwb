fitCatchEffort <-
function(listData, n.chains=4, thin=10, n.iter=100000, n.burnin=5000)
{

    require(rjags)
    require(coda)

    ## Defines the model used in JAGS
    codeModel <- "model {

   ########### Priors
   for (i in 1:Nua) {
        density[i]~dunif(0,0.3)
        N[i]<-round(density[i]*area[mug[i]])
   }

   ## overdispersion in the catchability
   t0~dnorm(0,0.01)
   phi <- exp(t0)


   for (m in 1:nmonths) {
      for (i in 1:ntypes) {
         catchability[m,i]~dunif(0,1)
      }
    }


   ########### Likelihood
   for (i in 1:Nua) {

     ## Initial value for N
     Ncurrent[1,i] <- N[i]
     for (j in 1:nmus[i]) {

         ## Expected hunting pressure
         pcurrentm[idstart[i]+j-1] <- 1-exp(-catchability[month[idstart[i]+j-1],type[mu[idstart[i]+j-1]]]*nhunters[idstart[i]+j-1])

         ## Conversion to parameters of the beta distribution
         parama[idstart[i]+j-1] <- pcurrentm[idstart[i]+j-1]*phi
         paramb[idstart[i]+j-1] <- (1-pcurrentm[idstart[i]+j-1])*phi

         ## Actual hunting pressure
         pcurrent[idstart[i]+j-1] ~ dbeta(parama[idstart[i]+j-1],paramb[idstart[i]+j-1])

         ## Culling
         cull[idstart[i]+j-1] ~ dbin(pcurrent[idstart[i]+j-1], Ncurrent[j,i])

         ## Update N
         Ncurrent[j+1,i] <- Ncurrent[j,i]-cull[idstart[i]+j-1]
     }
   }
}"


    ## Starting values
    density <- rep(.2, listData$Nua)
    catchability <- matrix(0.0005, listData$nmonths, listData$ntypes)

    ## Fit the model
    mod <- jags.model(textConnection(codeModel), data=listData,
                        inits=list(density=density, catchability=catchability),
                        n.chains = n.chains)

    cod2 <- coda.samples(mod, variable.names=c("catchability","density", "phi"),
                         thin = thin, n.iter=n.iter, n.burnin = n.burnin)

    ## Return the results
    return(cod2)
}

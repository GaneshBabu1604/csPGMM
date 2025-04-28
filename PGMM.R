library(LaplacesDemon)
library(ggplot2)
library(caTools)
library(mclust)
library(pgmm)
library(fossil)
library(MASS)
library(sys)
library(MixSim)
library(IMIFA)

#Creating a Parsimonious Gaussian Mixture Model(UUU) Alternating Expectation - Conditional Maximization Algorithm:
pgmmAECM <- function(x,G = 1,Q = 2,epoch=100, convergence_value = 0.1,zstart=1,zlist=c(),seed=123456,sim_data = 1, sim_type='wellseparated_'){
  #Input:
  #x : Dataframe(N,P). Where n is the number of observations for creating the clusters and p is the features of each observation.
  #epoch : int. Number of iterations.
  #G: int. Number of clusters.
  #Q: int. Number of factors.
  #convergence_value: convergence score. Aitken accelaration is used to assess convergence.
  #zstart: int. kmeans starting values (1) or user defined starting values(2)
  #zlist: vector. User defined starting values. Applicable only when zstart = 2
  #seed: int. seed for kmeans clustering. Applicable only when zstart = 1
  
  #Output:
  #predicted_cluster: final cluster membership of observations.
  #initial_cluster: initial cluster membership of observations.
  #predicted_score: matrix of posterior probability of cluster membership.
  #log_likelihood_list: vector of log-likelihood.
  #tau: prior probability of mixture components.
  #mu: mean of mixture components.
  #Sigma: Parsimonious covariance matrix of mixture components.
  #lambda: loadings matrix.
  #psi: diagonal covariance matrix of error/specific factor.
  #converged: whether the model is converged or not.
  #bic: Bayesian information criterion.
  #model_time: time taken to fit the model.
  
  #computing the total number of training pixels
  start_time = Sys.time()
  N=nrow(x)
  #computing the number of features for training pixels
  P=ncol(x)
  #Creating a initial cluster with kmeans clustering
  if (zstart == 1){
    set.seed(seed)
    init_cluster = kmeans(x,G,nstart = 10)$cluster
  }
  else if(zstart==2){
    init_cluster = zlist
  }
  #Creating the initial data with initial clusters for maximum log-likelihood estimation
  initial_data=cbind(x,init_cluster)
  colnames(initial_data)=c(colnames(x),"init_cluster")
  #Computing the initial parameters
  tau = matrix(0,G,1)
  mu = matrix(0,G,P)
  lambda = array(0,c(P,Q,G))
  psi = array(0,c(P,P,G))
  sigma = array(0,c(P,P,G))
  beta = matrix(0,Q,P)
  theta = matrix(0,Q,Q)
  s = array(0,c(P,P,G))
  Log_Likelihood = c()
  predicted_cluster=c()
  z=matrix(0,N,G)
  z_hat=matrix(0,N,G)
  z_hat=unmap(initial_data[,"init_cluster"])
  covw_computed = covw(x,z_hat)
  #Initialization
  for (g in 1:G){
    #Computing the probability of belonging to each class
    tau[g,]=nrow(initial_data[initial_data[,'init_cluster']==g,])/N
    #Computing the mean of each class
    mu[g,]=covw_computed$mean[,g]
    #Computing the covariance matrix of each class
    s[,,g]=covw_computed$S[,,g]
    eigen_g=eigen(s[,,g])
    for (q in 1:Q){
      for (p in 1:P){
        if (eigen_g$vectors[p,1] < 0){
          lambda[p,q,g] = -exp((0.5*log(eigen_g$values[q])+log(-eigen_g$vectors[p,1])))
        }
        else{
          lambda[p,q,g] = exp((0.5*log(eigen_g$values[q])+log(eigen_g$vectors[p,1])))
        }
      }
    }
    diag(psi[,,g]) = abs(diag(s[,,g]-lambda[,,g]%*%t(lambda[,,g])))
    # Updating sigma_g based on latent factors
    sigma[,,g] = lambda[,,g]%*%t(lambda[,,g]) + psi[,,g]
  }
  # AECM Algorithm
  for (e in 0:epoch){
    if(e%%100 == 0){
      print(e)      
    }
    #Step 1
    #Conditional Maximization
    covw_computed = covw(x,z_hat)
    for (g in 1:G){
      #Computing the probability of belonging to each class
      tau[g,]=(sum(z_hat[,g]))/N
      #Computing the mean of each class
      mu[g,]=covw_computed$mean[,g]
    }
    #Creating a temporary variable
    z_old = z_hat
    #Expectation
    if(e!=0){
      #Estimating the probability of each training pixels belonging to each class based on the estimated parameters
      for (g in 1:G){
        z[,g]=log(tau[g,])+dmvn(as.matrix(x),mu[g,],sigma[,,g],log=TRUE)
      }
      #Computing Log-Sum-Exp
      z_max = apply(z,1,max)
      z_sum=log(apply(exp(z-z_max),1,sum))+z_max
      z_hat=exp(z-z_sum)
    }
    # Step 2:
    # Update s,beta,theta and Conditional Maximization
    covw_computed = cov_weighted(x,G,mu,z_hat)
    for (g in 1:G){
      s[,,g]=covw_computed[,,g]
      # Using woodbury identity for computing the Inverse of lambda%*%t(lambda) + psi
      psi_inv = solve(psi[,,g])
      beta = t(lambda[,,g]) %*% (psi_inv - (psi_inv %*% lambda[,,g] %*% solve(diag(Q)+t(lambda[,,g]) %*% psi_inv %*% lambda[,,g]) %*% t(lambda[,,g]) %*% psi_inv))
      theta = diag(Q) - (beta %*% lambda[,,g]) + (beta %*% s[,,g] %*% t(beta))
      # Estimating new lambda
      lambda[,,g] = s[,,g] %*% t(beta) %*% solve(theta)
      # Estimating new psi
      diag(psi[,,g]) = (diag(s[,,g]-lambda[,,g]%*%beta%*%s[,,g]))
      # Estimating new sigma based on latent factors
      sigma[,,g] = lambda[,,g]%*%t(lambda[,,g]) + psi[,,g]
    }
    # Expectation
    #Estimating the probability of each training pixels belonging to each class based on the estimated parameters
    for (g in 1:G){
      z[,g]=log(tau[g,])+dmvn(as.matrix(x),mu[g,],sigma[,,g],log = TRUE)
    }
    #Computing Log-Sum-Exp
    z_max = apply(z,1,max)
    z_sum=log(apply(exp(z-z_max),1,sum))+z_max
    z_hat=exp(z-z_sum)
    print(sum(z_sum))
    Log_Likelihood =c(Log_Likelihood,sum(z_sum))
    write_time = Sys.time()
    temp_out_saver = list(X = x, N=N,P = P, init_cluster = init_cluster, epochs_ran = e, z_hat = z_hat,Log_Likelihood = Log_Likelihood, tau = tau, mu = mu, s = s, sigma = sigma, lambda = lambda, psi = psi, beta = beta, theta = theta, start_time = start_time,write_time = write_time,G = G,Q = Q,epoch=epoch, convergence_value = convergence_value, sim_data, sim_type)
    if(e%%500==0){
      save_loc = paste(sim_type,as.character(sim_data),'_G',as.character(G),'_epoch_',as.character(e),'_final.RData',sep = '')
      print(save_loc)
      save(temp_out_saver,file = save_loc)
    }
    #if(e>3){
    #  l_k_1 = aitken_acceleration(Log_Likelihood[length(Log_Likelihood)-2],Log_Likelihood[length(Log_Likelihood)-1],Log_Likelihood[length(Log_Likelihood)])
    #  if(toString(l_k_1) != "invalid"){
    #    if (abs(l_k_1-Log_Likelihood[length(Log_Likelihood)]) < convergence_value){
    #      converged = 'converged'
    #      print('Cluster Solution Converged')
    #      break
    #    }          
    #  }
    #}
    if(e >= 3){
      log_lik_change = (Log_Likelihood[length(Log_Likelihood)]-Log_Likelihood[length(Log_Likelihood)-1])/(Log_Likelihood[length(Log_Likelihood)-1])
      if(abs(log_lik_change)< convergence_value){
        converged = 'Converged'
        print('Cluster solution converged')
        break
      }
    }
    if (e == epoch){
      converged = 'not converged'
      print('Cluster Solution not Converged')
      break
    }
  }
  #bic_model = 2*tail(Log_Likelihood,1) - PGMM_dfree(Q,P,G,'UUU') * log(N)
  bic_model = 'not_computed'
  final_prediction=data.frame(z_hat)
  colnames(final_prediction) = c(1:G)
  predicted_cluster = colnames(final_prediction)[max.col(final_prediction, ties.method = "first")]
  end_time = as.numeric(Sys.time() - start_time,units = "mins")
  return(list(predicted_cluster=predicted_cluster,initial_cluster = init_cluster, probability_score=z_hat,log_likelihood_list = Log_Likelihood,tau=tau,mu=mu,sigma=sigma,lambda=lambda,psi=psi,converged = converged,bic = bic_model,time = end_time))
}

cov_weighted = function(x,G=1,mu,z_hat){
  s = array(0,c(ncol(x),ncol(x),G))
  for(g in 1:G){
    x_temp = sweep(x,2,mu[g,],'-')
    x_temp_2 = sweep(x_temp,1,z_hat[,g],'*')
    s[,,g] = (t(x_temp_2) %*% as.matrix(x_temp))/sum(z_hat[,g])
  }
  return(s)
}

aitken_acceleration = function(l1,l2,l3){
  a_k = (l3 - l2)/(l2-l1)
  if(a_k<1.0){
    l_inf = l2 + ((l3-l2)/(1-a_k))
  }
  else{
    l_inf = "invalid"
  }
  return(l_inf)
}

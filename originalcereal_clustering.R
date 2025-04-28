source('PGMM.R')
source('constrained_PGMM.R')
source('consensus_Clustering.R')

# Data
Data = read.csv('Data.csv')
Data$X<-NULL

load('cols_list.RData')
n_cols_list = 1:length(cols_list)

G = 4
N = nrow(Data)
P = ncol(Data)
Q = 2
e = 25000

# Fitting PGMM and constrained-PGMM on different settings of r & d
pgmmAECM_implementer<-function(x, colname,G=4,Q=2,epoch = 10000, zstart = 1,zlist = c(), seed = 123456){
  start_time = Sys.time()
  new <- try(pgmmAECM(x,G = G, Q = Q, epoch = epoch,convergence_value=1e-7,zstart = zstart,zlist = zlist,sim_data = colname, sim_type='originalpuffedcereal_cPGMM_UUU_'),silent = T)
  if(is.list(new)==TRUE){
    prod_probability = new$probability_score*log(new$probability_score)
    prod_probability[is.na(prod_probability)] = 0
    weight <- -1*sum(colSums(prod_probability))
    weighted_B <- (1/weight)
    end_time = as.numeric(Sys.time() - start_time,units='mins')
    return(list(new$probability_score,new$converged,new$predicted_cluster,weighted_B,end_time,new$log_likelihood_list,new$initial_cluster,colname,new$tau, new$mu, new$lambda, new$psi, new$sigma)) 
  }
  else{
    return(new)
  }
}

# Parallel run
library(parallel)  
library(doParallel)
no_cores <- detectCores()
no_cores
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
clusterEvalQ(cl, {
  library(LaplacesDemon)
  library(ggplot2)
  library(caTools)
  library(mclust)
  library(pgmm)
  library(fossil)
  library(MASS)
  library(IMIFA)
  library(MixSim)
})
library(ggplot2)


print('simulation begins')
r_index = 1
for(r in c(25)){
      start_time = Sys.time()
      #pgmmAECM
      result <- foreach(m = n_cols_list) %dopar% pgmmAECM_implementer(Data[,cols_list[[m]][[1]]],m,G=G,Q=Q,epoch = e, zstart = 1)
      print('Models computed')
      end_time = as.numeric(Sys.time() - start_time,units = "mins")
      print(end_time)
      parallel_list = list(no_models = r,parallel_time = end_time)
      n_result = 1
      for(s in result){
        location = paste('originalpuffedcereal_cPGMM_r',as.character(r),'_d',as.character(20),'_',as.character(n_result),'_final.RData',sep = '')
        save(s,file = location)
        n_result=n_result+1
      }
      location = paste('originalpuffedcereal_cPGMM_r',as.character(r),'_d',as.character(20),'_creation_time_final.RData',sep = '')
      save(parallel_list,file=location)
}

stopCluster(cl)
cols_list[[20]][[1]]

# Artificial cereal generation
Data = read.csv('Mixture_Image.csv')
Data$X<-NULL
pred_cluster_list = list()
G = 4
N = nrow(Data)
z=matrix(0,N,G)
z_hat=matrix(0,N,G)
for(ss in c(1:25)){
  loc = paste('originalpuffedcereal_cPGMM_r25_d20_',as.character(ss),'.RData',sep = '')
  load(loc)
  params = s
  for (g in 1:G){
    z[,g]=log(s[[9]][g,])+dmvn(as.matrix(Data[,cols_list[[ss]][[1]]]),s[[10]][g,],s[[13]][,,g],log=TRUE)
  }
  z_max = apply(z,1,max)
  z_sum=log(apply(exp(z-z_max),1,sum))+z_max
  z_hat=exp(z-z_sum)
  final_prediction=data.frame(z_hat)
  colnames(final_prediction) = c(1:G)
  predicted_cluster = colnames(final_prediction)[max.col(final_prediction, ties.method = "first")]
  pred_cluster_list[[ss]] = predicted_cluster
}

table(pred_cluster_list[[25]])
table(temp)
temp = rep(0,43621)
temp[which(pred_cluster_list[[25]]==1)] = 2

pred_cluster_list[[25]] = temp

pred_cluster_dataframe = matrix(0,nrow = 43621,ncol = 25)
for(i in c(1:25)){
  pred_cluster_dataframe[,i] = as.numeric(pred_cluster_list[[i]])
}

# Function to find most frequent value(s) and how many times it occurs
max_freq_value <- function(x) {
  freq_table <- table(x)
  max_count <- max(freq_table)
  max_values <- names(freq_table)[freq_table == max_count]
  
  # If there's more than one value with same max frequency, list all
  list(values = max_values, count = max_count)
}

# Apply to each row and get results
results <- apply(pred_cluster_dataframe, 1, max_freq_value)

# Transform into a readable data frame
output <- data.frame(
  row = paste0("Row", seq_along(results)),
  value = sapply(results, function(x) paste(x$values, collapse = ", ")),
  count = sapply(results, function(x) x$count)
)

print(output)

save(pred_cluster_list,file = 'pred_cluster_list.RData')

pred_cluster_dataframe
write.csv(pred_cluster_dataframe,file = 'mixture_image_cPGMM_all_predicted_cluster.csv')
write.csv(output,file = 'mixture_image_cPGMM_predicted_cluster_count.csv')
length(intersect(which(pred_cluster_list[[2]]==2),which(pred_cluster_list[[1]]==4)))

write.csv(predicted_cluster,file = 'mixture_image_PGMM_predicted_cluster.csv')
uncertainty = 1- apply(z_hat,1,max)
write.csv(uncertainty,file = 'mixture_image_PGMM_uncertainty.csv')
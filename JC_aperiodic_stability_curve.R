
######################################
# define functons here 
fun_with_subjects <- function(n) { 
  # input n is number of subjects to run over
  
  if(n==250){
    lapply(1, iterate_over_subject, n=n) #1:1000
  }else{
    lapply(1:100, iterate_over_subject, n=n) #1:1000
  }
  
}


######################################
# define functons here 
iterate_over_random_files <- function(j, n, rand_index_of_files_all, group1, group2, k) {  #iterations for resampling purpose under a defined participant subset
  
  # loop over possible number of files 
  rand_index_of_files_all=sample(num_files)
  
  beta_coef = vector(mode = "double", length = num_files)
  p_value = vector(mode = "double", length = num_files)
  AIC = vector(mode = "double", length = num_files)
  BIC = vector(mode = "double", length = num_files)
  BF = vector(mode = "double", length = num_files)
  DeltaAIC = vector(mode = "double", length = num_files)
  number_of_files = vector(mode = "double", length = num_files)
  files_used = vector(mode = "list", length = num_files)
  number_of_subj = vector(mode = "double", length = num_files)
  iteration_subj = vector(mode = "double", length = num_files)
  
  
  for (i in 1: num_files){
    # loop over random selection of files
    mean_data =c()
    rand_index_of_files=rand_index_of_files_all[1:i]
    
    if (i !=1){mean_data=rowMeans(data[c(group1, group2),rand_index_of_files])} else{mean_data = data[c(group1, group2),rand_index_of_files]}
    
    if (i !=1){mean_MSE=rowMeans(MSE[c(group1, group2),rand_index_of_files])} else{mean_MSE = MSE[c(group1, group2),rand_index_of_files]}
    
    #take mean matrix and fit model
    group= rep(c('Group 1', 'Group 2'), each= n)
    alpha_power = scale(mean_data)
    MSE_scale = scale(mean_MSE)
    lm0= lm(alpha_power ~ MSE_scale)
    lm1= lm(alpha_power ~ group + MSE_scale)
    sum_lm1= summary(lm1)
    bf=exp(0.5*(BIC(lm0)-BIC(lm1)))
    
    #write data to pre-allocated vectors
    vec_ind = i
    beta_coef[vec_ind] = lm1$coefficients[2]
    p_value[vec_ind] = sum_lm1$coefficients[2,4]
    AIC[vec_ind] = AIC(lm1)
    BIC[vec_ind] = BIC(lm1)
    BF[vec_ind] = bf
    DeltaAIC[vec_ind] = AIC(lm0)-AIC(lm1)
    number_of_files[vec_ind] = i
    #files_used[vec_ind] = paste0(rand_index_of_files, collapse=" ")
    number_of_subj[vec_ind] = n
    iteration_subj[vec_ind] = k
  } 
  
  # TRY MAKING THIS A MATRIX CALL TO SPEED IT UP
  outtt=  na.omit(matrix(c(beta_coef,p_value,AIC,BIC,BF,DeltaAIC,number_of_files,number_of_subj,iteration_subj), nrow = 100))
  return(outtt)
  
}


iterate_over_subject <- function(k, n) {  #iterations for a certain number of participants subjects
  
  # select random n subject per group
  group1= sample(250)[1:n]
  group2= sample(251:500)[1:n]
  
  #pre-allocate for efficiency #AIW: switched to pre-allocating as vectors and then binding into a data.frame before writing to file, to avoid re-allocation issues (which slow down the loop more and more each iteration)
  
  outtt= mclapply(1:num_file_iter, iterate_over_random_files, n=n, rand_index_of_files_all=rand_index_of_files_all, group1=group1, group2=group2, k=k ) #1:1000
  
  output = array(unlist(outtt), dim= c(100,9,num_file_iter))
  output <- aperm(output, c(1, 3, 2))
  dim(output) <- c(prod(dim(output)[-3]), dim(output)[3])
  output = na.omit(output)
  colnames(output)= c('beta_coef', 'p_value', 'AIC', 'BIC', 'BF', 'DeltaAIC', 'number_of_files', 'number_of_sbuj', 'iteration_subj')
  output= as.data.frame(output)
  
  output$number_of_subj= n
  output$iteration_subj= k 
  file=paste('/export02/data/Sapphire_stability/new_results/0.2_aperiodic_participant_subset_noSubj_', n,'_no_of_iteration_',k, '.csv', sep = '')
  write.table(x = output, file = file, append = TRUE, sep = ",", col.names = FALSE)
  
}


##### quantify the effect sizes #####
#library(dplyr)
library(parallel)

cbbPalette= c("#4f89e0", "#f5ec6c",'#156605',"#76D7C4","#268236", '#4d3d87','#593d80', '#3b49e3')

num_files=100
num_subj= 500
num_file_iter=1000

data=read.csv('/export02/data/Sapphire_stability/SimulatedData/aperiodic/0.2_250_PSD/0.2_250_exponent_extracted.csv',header  = FALSE)
MSE=read.csv('/export02/data/Sapphire_stability/SimulatedData/aperiodic/0.2_250_PSD/0.2_250_MSE_extracted.csv',header = FALSE)

# HERE WE APPLY TO ARRAY OF SAMPLE SIZE FUNCTION X

subj_array=c(10,20,40,60,80,100,250)

lapply(subj_array, fun_with_subjects)


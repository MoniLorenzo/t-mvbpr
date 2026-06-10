################################################################################
######################## MULTI-VIEW BPR#########################################
################################################################################

###FUNCTION TO INCLUDE##########################


##########################
#2. Sample from categorical
#(ountput numeric from 1 to n_cat)
##########################
SampleMult1=function(Probs){
  K=length(Probs)
  return(sample.int(K, 1, prob = Probs))
}


##########################
#3. Sample/Update   Dir
##########################
SampleDir=function(a ){
  l=sum(!is.na(a) )
  r=a
  N=rgamma(l,shape = a)  
  
  r[1:l]= N/sum(N, na.rm = TRUE)  
  
  return(r)
}
 






















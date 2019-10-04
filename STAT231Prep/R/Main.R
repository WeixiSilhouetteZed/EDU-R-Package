# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
purpose <- function() {
  print("The purpose of this library is to provide essential
        and simple support for undergraduate statistical education.
        There are a lot of functions related to different functionalities
        of this library. A detailed list of function will be embedded
        in later released versions.")
}

#Part I: Data Generation on Statistical Model
#In this part, the library can provide easy support
#related to generating tables on different statistical model


StatVec<-function(model,size,paraList=c()){
  # This function can accept any default R statistical models to generate random
  # vector based on that. The user is required to include a parameter list for
  # the model specified.
  strRe=paste0(paraList,collapse = ",")
  eval(parse(text=paste("r",model,"(",size,",",strRe,")",sep="")))
}
StatDF<-function(modelList,colNum,colNameList=c(1:colNum),rowNum,paraList){
  # This function can accept a list of statistical model names in R statical model
  # to generate a dataframe of random column vectors. The user should note that
  # the paraList parameter in this case is a list of list.
  mll<-length(modelList)
  pll<-length(paraList)
  matAcc<-cbind(c())
  if (mll<colNum){
    for (i in c(1:mll)){
      model<-modelList[i]
      paraListPart<-paraList[i]
      matAcc<-cbind(matAcc,StatVec(model,rowNum,paraListPart))
    }
    for (i in c(1:(colNum-mll))){
      model<-modelList[length(modelList)]
      paraListPart<-paraList[length(modelList)]
      matAcc<-cbind(matAcc,StatVec(model,rowNum,paraListPart))
    }
  }
  else{
    for (i in c(1:colNum)){
      model<-modelList[i]
      paraListPart<-paraList[i]
      matAcc<-cbind(matAcc,StatVec(model,rowNum,paraListPart))    }
  }
  resultDF<-as.data.frame(matAcc)
  colnames(resultDF)<-colNameList
  resultDF
}
#########################################################################
#Part II: Graph Report Generation on created dataframe
#Intended Graph List: pdf,cdf,scatterplot


gen_graph<-function(vec,name="data"){
  # This function is used to generate a set of normality check graphics
  # based on a given list and the size is 2x2.
  library(MASS)
  par(mfrow=c(2,2))
  boxplot(vec, main=paste("Boxplot of",name))
  truehist(vec, main=paste("Histogram of",name))
  curve(dnorm(x,mean(vec),sd(vec)),add = TRUE,col='red')
  plot(ecdf(vec), main=paste("Boxplot of",name))
  curve(pnorm(x,mean(vec),sd(vec)),add = TRUE,col='red')
  qqnorm(vec, main=paste("QQ-plot of", name))
  qqline(vec)
}

#########################################################################
#Part III: Calculation Package

mle_summary<-function(xvec,thetaList,probdf,lk_lvl){
  # This function is used to evaluate MLE for a given pdf/pmf based on observed data,
  # the user is required to test whether the starting point is viable in the thetalist function.
  # A relative likelihood function graph together with the 15% likelihood line will be presented.
  fn_lst<-NULL
  for (x in xvec){
    fn<-function(theta){
      eval(parse(text=probdf))
    }
    fn_lst<-c(fn_lst,fn)
  }
  fn_sum_base<-function(theta){
    0
  }
  fn_sum<-function(theta){
    result<-0
    if (length(xvec)==0){
      result<-fn_sum_base(theta)
    }
    else{
      lambda<-function(fn){
        fn(theta)
      }
      result<-sum(sapply(fn_lst,lambda))
    }
    result<-(-1)*result
    result
  }
  mle_result<-nlm(fn_sum, theta <- thetaList, hessian=TRUE)
  print(mle_result)
  plot_range<-seq(mle_result$estimate-3*sd(xvec),mle_result$estimate+3*sd(xvec),0.01)
  invert_fn_sum<-function(theta){
    fn_sum(theta)/fn_sum(mle_result$estimate)
  }
  plot(plot_range,sapply(plot_range,invert_fn_sum),type="l",col="blue",
       xlab=expression(theta),ylab=expression(R(theta)),main="Relative Likelihood Function")
  abline(a=lk_lvl,b=0,col="red")
  expr<-function(){
    eval(parse(text=probdf))
  }
  text(mle_result$estimate,0.8,paste("hat(theta)=",
    parse(text=round(mle_result$estimate,digits = 3))))
  list(mle_result,invert_fn_sum)
}

lk_int_summary<-function(xvec,thetaList,probdf,lk_lvl,root="DEFAULT"){
  #This function can execute mle_summary with additional likelihood interval
  #information
  mle_info<-mle_summary(xvec,thetaList,probdf,lk_lvl)
  if (root=="DEFAULT"){
    root_lst<-list(c(mle_info[[1]]$estimate-3*sd(xvec),mle_info[[1]]$estimate)
                  ,c(mle_info[[1]]$estimate,mle_info[[1]]$estimate+3*sd(xvec)))
  }
  else{
    root_lst<-root
  }
  sol<-c(NULL,NULL)
  f<-function(x){
    mle_info[[2]](x)-lk_lvl
  }
  sol[1]<-uniroot(f,root_lst[[1]])
  sol[2]<-uniroot(f,root_lst[[2]])
  print(sol)
  print(paste0("The ",lk_lvl*100,"% likelihood interval is:"))
  print(paste0("[",round(sol[[1]],digits = 3),",",round(sol[[2]],digits = 3),"]"))
}









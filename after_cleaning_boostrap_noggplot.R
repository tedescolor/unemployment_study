# if (!require("pacman")) install.packages("pacman")
# pacman::p_load("quantreg","dplyr","doParallel","foreach",
#                "abind")

list.of.packages <- c("quantreg","foreach","Rlab","optimbase",
                      "SimDesign","zeallot", "survival","utils","dppmix",
                      "microbenchmark","reshape2","dplyr","gdata","numbers","pracma","doParallel","foreach",
                      "abind","icesTAF","fastDummies","rockchalk","Rcpp","abind")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos='https://lib.ugent.be/CRAN/')
lapply(list.of.packages, library, character.only = TRUE)

print("all packages but ggplot loaded")

#list.of.packages <- c("ggplot2", "Rcpp","survival","quantreg","dplyr","doParallel","foreach","abind" )
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)
#lapply(list.of.packages, library, character.only = TRUE)

#library(ggplot2, lib.loc="/data/leuven/339/vsc33962/packages")
#print("loaded ggplot2")
library(ggplot2)
#data = read.csv("/scratch/leuven/339/vsc33962/cleaned_data_my_version2.csv")
data = read.csv("cleaned_data_my_version2.csv")
print("loaded data")
data = data[1:100,]
#setwd("/scratch/leuven/339/vsc33962/output")
n = nrow(data)
#output setting
cut.year = 2013 
current.comparison.string = paste(cut.year,"unemployment_study_size",n,sep = "")  #string to append to every image name before saving
file.output = paste(current.comparison.string,"_output.txt",sep = "") #name file output

#data = data[sample.int(nrow(data), size = 1000),]
nu = 0.09
tauU = 0.61
step = 0.01
taus = seq(nu,tauU,step)
rq1 =  quantreg::crq( survival::Surv(utime1,censored1)~sex+age+
              driving+trasportation+knowledge_dutch+immigrant+
              education,
            data = data , taus = taus , method = 'PengHuang')
rq2 =  quantreg::crq( survival::Surv(utime2,censored2)~sex+age2+
              driving+trasportation+knowledge_dutch+immigrant+
              education,
            data = data,  taus = taus , method = 'PengHuang')

#saveRDS(rq1, file = "rq1.rds")
#saveRDS(rq2, file = "rq2.rds")
#rq1 = readRDS("rq1.rds")
#rq2 = readRDS("rq2.rds")
taus = tail(head(taus,-1),-1)
t1 = coef(rq1,taus = taus)
t2 = coef(rq2,taus = taus)

t = t1 - t2
sum(is.na(t))
saveRDS(t, file = paste(current.comparison.string,"_t.rds",sep = ""))
n = nrow(data)
td = as.data.frame(t(t))
td$taus = taus
tdd <- melt(td,id.var="taus")
p = ggplot(tdd, aes(x=taus,y=value,group=variable,colour=variable)) +
  geom_line(aes(lty=variable),size=1) 
saveRDS(p, file = paste(current.comparison.string,"_plot_t.rds",sep = ""))

# png(file=paste(current.comparison.string,"lines_comparison_t",".png",sep = ""),width=1200, height=750)
# matplot(taus,t(t),type = "l",col = c(1:nrow(t)), 
#         main = paste("components of (beta1 - beta2)(tau): before and after", cut.year),
#         ylab = '', sub = paste(n,"obs"),lwd = 2 )
# legend("topleft", legend=rownames(t), col=c(1:nrow(t)),lty = 1:nrow(t))
# dev.off()

## dfs as matrix, necessary for bootcrq2samples imput
#df1$region = as.integer(df1$region)
#df2$region = as.integer(df2$region)

x1 = data %>% dplyr::select(sex,age,
                     driving,trasportation,knowledge_dutch,immigrant,
                     education) %>% as.matrix(.)
x1 = cbind(1,x1)
x2 = data %>% dplyr::select(sex,age2,
                     driving,trasportation,knowledge_dutch,immigrant,
                     education)  %>% as.matrix(.)
x2 = cbind(1,x2)
y1 = data %>% dplyr::select(utime1)  %>% as.matrix(.)
y2 = data %>% dplyr::select(utime2)  %>% as.matrix(.)
c1 = data$censored1  %>% as.matrix(.)
c2 = data$censored2  %>% as.matrix(.)
colnames(x1) = c("intercept", tail(colnames(x1),-1))
colnames(x2) = c("intercept", tail(colnames(x2),-1))

R = 50 # bootstrap iterations

sample.size.bootstrap =  nrow(data)#size of boostrap data
my.abind = function(x,y){ #personal abind function (just abind with along=3)
  abind(x,y,along = 4)
}
boot.crq.fixed.sample = function(x, y, c, taus, method, ctype = "right", R = 1, 
                                 mboot, bmethod = "jack", sample, ...) {
  n <- length(y)
  p <- ncol(x)
  if (length(taus) > 1) 
    A <- array(0, dim = c(p, length(taus), R))
  else A <- matrix(0, p, R)
  for (i in 1:R) {
    if (bmethod == "jack") {
      s <- sample #sample(1:n, mboot)
      yb <- y[-s]
      xb <- x[-s, ]
      cb <- c[-s]
      w <- rep(1, n - mboot)
    }
    else if (bmethod == "xy-pair") {
      w <- table(sample) #table(sample(1:n, mboot, replace = TRUE))
      s <- as.numeric(names(w))
      w <- as.numeric(w)
      yb <- y[s]
      xb <- x[s, ]
      cb <- c[s]
    }
    else if (bmethod == "Bose") {
      w <- sample
      yb <- y
      xb <- x
      cb <- c
    }
    else stop("invalid bmethod for boot.crq")
    if (method == "Portnoy") 
      a <- quantreg::crq.fit.por(xb, yb, cb, weights = w, ctype = ctype, 
                                 ...)
    else if (method == "PengHuang") 
      a <- quantreg::crq.fit.pen(xb, yb, cb, weights = w, ctype = ctype, 
                                 ...)
    else stop("Invalid method for boot.crq")
    if ((i%%floor(R/10)) == 0 & n > 1e+05) 
      cat(paste("bootstrap roughly ", 100 * (i/R), 
                " percent complete\n"))
    if (length(taus) > 1) 
      A[, , i] <- coef(a, taus)
    else A[, i] <- coef(a, taus)
  }
  list(A = A, n = length(y), mboot = mboot, bmethod = bmethod)
}


# boot.crq from quantreg package for 2 dataset at the same time: the same bootstrap sample is used for both the data:
# the function works in the following way: create the bootstrap sample depending on "bmethod" and than call the funtion boot.crq.fixed.sample
method = "PengHuang"; ctype = "right";bmethod = "Bose"; n <- length(y1); p <- ncol(x1);
#A = array(0, dim = c(p, length(taus), 2,R))
#setup parallel backend to use many processors
cores=detectCores()
print(cores)
cl <- makeCluster(cores[1]) #not to overload your computer
registerDoParallel(cl)

pb = txtProgressBar(min = 0, max = R, initial = 0) 
res = foreach(i=1:R,.combine=my.abind, .packages="foreach") %dopar% {
  if (length(taus) > 1) {
    currentA <- array(0, dim = c(p, length(taus),2, 1))

  }else{
    return(NA)
  }
  setTxtProgressBar(pb,i)
  if (bmethod == "jack") {   
    mboot <- 2 * ceiling(sqrt(n))
    }else{ mboot <- n}
  
  if (bmethod == "jack") {
    sample <- sample(1:n, mboot);
    
  }
  else if (bmethod == "xy-pair") {
    sample <- table(sample(1:n, mboot, replace = TRUE));
    
  }
  else if (bmethod == "Bose") {
    sample <- rexp(n);
  }
  else stop("invalid bmethod for boot.crq")
    currentA[, ,1,1] <- boot.crq.fixed.sample(x = x1, y = y1, c = c1, taus = taus, method = method, ctype = ctype, bmethod = bmethod, sample = sample, mboot = mboot)$A;
    currentA[, ,2,1] <- boot.crq.fixed.sample(x = x2, y = y2, c = c2, taus = taus, method = method, ctype = ctype, bmethod = bmethod, sample = sample, mboot = mboot)$A;
  return(currentA)  
}

stopCluster(cl)
A1 = res[,,1,]
A2 = res[,,2,]
A = A1 - A2
indices.NA = apply(A,3,function(x){any(is.na(x))})
A = A[,,!indices.NA]


write(paste("rows of data = ",nrow(data) ,"\n",
            "Boostrap iteration R =  ",R,"\n",
            "QR call :",toString(rq1$call), sep = ""), file = file.output,append = TRUE)

#save A 
saveRDS(A, file = paste(current.comparison.string,"_A.rds",sep = "")) #save
#A <- readRDS(paste(current.comparison.string,"_A.rds",sep = "")) #load
#graphical solutions
m05 = matrix(nrow = nrow(A), ncol = ncol(A))
m95 = matrix(nrow = nrow(A), ncol = ncol(A))
mean = matrix(nrow = nrow(A), ncol = ncol(A))
for(i in 1:nrow(A)){
  m05[i,] = apply((A[i,,]), 1,quantile,probs = 0.05)
  m95[i,] = apply((A[i,,]), 1,quantile,probs = 0.95)
  mean[i,] = apply((A[i,,]), 1,mean)
}

names = colnames(x1)
for(j in 1:nrow(t)){
  
  data.current <- data.frame(cbind(taus,
                                   t[j,],
                                   m05[j,],
                                   m95[j,],
                                   mean[j,]
  ))
  
  colnames(data.current) <- c("taus","beta","q05","q95","mean")
  
  p <- ggplot(data.current, aes(x =taus)) +
    geom_line(aes(y=beta), color="red") +
    geom_line(aes(y=q05), color = "black") +
    geom_line(aes(y=q95), color = "black") +
    geom_line(aes(y=mean), linetype = "dashed", color = "black") +
    ggtitle(names[j])
  #print(p)
  #ggsave(p, paste(current.comparison.string,"bandwidth_",names[j],".png"))
  saveRDS(p,paste(current.comparison.string,"bandwidth_",names[j],".rds",sep = ""))
}

# plot all simulated process
# for(i in 1:nrow(t)){
#   #png(file=paste(current.comparison.string,"lines_comparison_",names[i],".png"),width=1200, height=700)
#   matplot(taus, A[i,,],type = "l", main=names[i])
#   lines(taus,t[i,], col="red",lwd=3)
#   #dev.off()
# }
# compute the integral passing the "values" that a function assumes in "points", 
# the points and a possible weight function


# compute the integral passing the "values" that a function assumes in "points", 
# the points and a possible weight function
integral= function(values, points, weight = NULL,...){
  if(length(points) == length(values)){
    values = head(values,-1)
  }
  if(length(points)-1 != length(values) ){
    print("dimension issue")
    return(NULL)
  }
  if(is.null(weight)){
    valutated.increments = tail(points,-1) - head(points,-1)
  }else{
    points = sapply(points,weight,...)
    valutated.increments = tail(points,-1) - head(points,-1)
  }
  return(sum(values * valutated.increments))
}


pvalue = function(T2,T2STAR){
  tryCatch({
    sum( T2STAR> T2 ) / length(T2STAR)
    
  }, warning = function(warning_condition) {
    return(NA)
  }, error = function(error_condition) {
    return(NA)
  }, finally={
    
  })
}

T2 = integral(values = apply( sqrt(n)* (t) ,2,pracma::Norm, p=2),
              points = tail(seq(nu,tauU,step),-1)) #T2 statistic


#print(paste("sup_tau || [ sqrt(n) * ( (hat(beta1) - hat(beta2)) )(tau) ] || :",round(T2,3)))


# compute the T2 statistic for the bootstrap betas (namely T2STAR)
N = pracma::size(A,3) # number of available iterations
T2STAR = array(dim = N) #array of integrals
for(i in 1:N){
  T2STAR[i] = integral(values = apply( sqrt(n) * ( as.matrix(A[,,i])  - t) , 2,pracma::Norm, p=2),
                       points = tail(seq(nu,tauU,step),-1))
}


## compute pvalue of T2 wrt T2STAR  
pv = pvalue(T2,T2STAR)

#print(paste("H0: beta1 = beta2. p-value = ", pv))

write(paste("H0: beta1 = beta2. p-value = ", pv), file = file.output,append = TRUE)
#plot the distrubution of the norm among processes
# p = ggplot() + aes(T2STAR) + geom_histogram() +
#   geom_vline(xintercept = T2, linetype="dotted", color = "red", size=1.5) +
#   ggtitle(paste("norm distrubution pvalue = ",pv))
# print(p)
# ggsave(p,filename = paste(current.comparison.string,"histogram_of_integral_of_norms.png",sep = ""))
# 
# setwd("..")

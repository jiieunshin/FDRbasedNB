library(zoo)
library(dplyr)
library(changepoint)
library(doParallel)

rm(list=ls())
gc()

source("C:\\Users\\uos\\OneDrive - UOS\\바탕 화면\\논문\\R code\\functions.R")


N.tr = 1000
p = 5000

N.te = 1000
N.Rep = 100

#i=1
#j=1

# tmp_list = list()
# sel.var = rep(0, p)

detectCores()
n.clu = 5
cl <- makeCluster(n.clu)
registerDoParallel(cl)

# j=1;i=1
tmp_list = list()
sel.var = rep(0, p)

result = foreach(i = 1:N.Rep, .packages = c('changepoint', 'dplyr')) %dopar%{

    print(i)
    print(Sys.time())
    
    set.seed(i)
    seed = i
    tr = sim.data(N.tr, p)
    te = sim.data(N.te, p)
    
    fit = nbayes(tr$X, tr$Y)
    err.no = predict.nbayes(te$X, te$Y, fit$prob, fit$prior)$miss
    
    # 시간 재기
    time = system.time({ try({sel = cpt.method(tr$X, tr$Y)}, silent = TRUE) })[3]
    
    if(is.na(sel$var[1]) == TRUE){
      err.sel = NA
      n.sel = NA
      tp = NA
      fp = NA
      precision = NA
      recall = NA
      f1.score = NA
      fdr.hat = NA
    } else{
      
      
      pred = predict.nbayes(te$X[,sel$var], te$Y, sel$prob, sel$prior)
      err.sel = pred$miss
      sel.var[sel$var] = sel.var[sel$var] + 1
      
      n.sel = length(sel$var)
      tp = sum(sel$var <= 10)
      fp = sum(sel$var > 10)
      
      result.tab = table(te$Y, pred$Yhat)
      prec.class = diag(result.tab)/colSums(result.tab)
      precision =  mean(prec.class)
      
      recal.class = diag(result.tab)/rowSums(result.tab)
      recall = mean(recal.class)
      f1.score = 2 * (precision * recall) / (precision + recall)
      
      fdr.hat = 1 - precision

    tmp_list[[i]] = c(seed, err.no, err.sel, n.sel, tp, fp, precision, recall, f1.score, fdr.hat, time)

  }
}
stopCluster(cl)

# result
tmp_dat = result %>% do.call(rbind, .) %>% data.frame(., stringsAsFactors = F)

names(tmp_dat) = c('seed', 'err.no' , 'err.sel', 'n.sel', 'tp', 'fp', 'precision', 'recall', 'f1.score', 'fdr.hat', 'time')
apply(tmp_dat, 2, function(x) mean(x, na.rm = T)) %>% round(., 4)
apply(tmp_dat, 2, function(x) sd(x, na.rm = T)/sqrt(N.Rep)) %>% round(., 4)

# path = "C:\\Users\\uos\\Desktop\\논문\\R code\\cpt_result\\"
# 
# save(tmp_dat, file = paste0(path, "cpt_", N.tr, "_", p, "_result.Rdata"))
# save(sel.var, file = paste0(path, "cpt_", N.tr, "_", p, "_var.Rdata"))

Sys.time()


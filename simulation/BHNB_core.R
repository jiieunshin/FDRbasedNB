# library(zoo)
library(dplyr)
library(doParallel)
library(DescTools)
library(Metrics)

rm(list=ls())
gc()
source("C:\\Users\\uos\\OneDrive - UOS\\바탕 화면\\논문\\R code\\functions.R")

# 

N_tr = 100
p = 100
n_sig = 50

N_te = 1000
N_Rep = 100

n_clu = 5
cl <- makeCluster(n_clu)
registerDoParallel(cl)

# j=1;i=1
out_dat = list()
sel_var = rep(0, p)
sel_var_list = list() 
result = foreach(j = 1:N_Rep, .packages = c('zoo')) %dopar%{
  # for(j in 1:N_Rep){
    cat('iter :', j, '\n')
    print(Sys.time())
    
    set.seed(j)
    seed = j
    tr = sim_data(N_tr, p)
    te = sim_data(N_te, p)
    
    fit = nbayes(tr$X, tr$Y)
    err_no = predict_nbayes(te$X, te$Y, fit$prob, fit$prior)$miss
    
    # 시간 재기
    time = system.time({ try({sel = bhfdr(tr$X, tr$Y)}, silent=TRUE) })[3]
    opt_alpha = sel$opt_alpha
    
    pred = predict_nbayes(te$X[,sel$var], te$Y, sel$prob, sel$prior)
    err_sel = pred$miss
    sel_var[sel$var] = sel_var[sel$var] + 1
    
    n_sel = length(sel$var)
    tp = sum(sel$var <= n_sig)
    fp = sum(sel$var > n_sig)
    
    
    true_num_var = rep(c(1, 0), c(n_sig, p - n_sig))
    est_num_var = rep(0, p)
    est_num_var[sel$var] = 1
    
    result_tab = table(true_num_var, est_num_var)

    # result.tab = table(te$Y, pred$Yhat)
    prec_class = diag(result_tab)/colSums(result_tab)
    precision =  mean(prec_class)

    recal_class = diag(result_tab)/rowSums(result_tab)
    recall = mean(recal_class)
    f1_score = 2 * (precision * recall) / (precision + recall)

    out_dat[[j]] = c(seed, err_no, err_sel, n_sel, tp, fp, precision, recall, f1_score, opt_alpha, time)
    sel_var_list[[j]] = sel_var
  }

stopCluster(cl)
Sys.time()

# result

tmp_dat = result %>% do.call(rbind, .) %>% data.frame(., stringsAsFactors = F)

names(tmp_dat) = c('seed', 'err_no' , 'err_sel', 'n_sel', 'tp', 'fp', 'precision', 'recall',
                   'f1_score', 'opt_alpha', 'time')
apply(tmp_dat, 2, function(x) mean(x, na.rm = T)) %>% round(., 4)
apply(tmp_dat, 2, function(x) sd(x, na.rm = T)/sqrt(N_Rep)) %>% round(., 4)

# path = "C:\\Users\\uos\\Desktop\\논문\\R code\\cpt_result\\"
# 
# save(tmp_dat, file = paste0(path, "cpt_", N.tr, "_", p, "_result.Rdata"))
# save(sel.var, file = paste0(path, "cpt_", N.tr, "_", p, "_var.Rdata"))




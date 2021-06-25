# install library 
# library(dplyr)
# library(doParallel)
# library(DescTools)
# library(Metrics)
# library(doSNOW)
# library(Rmpi)

# setting the initial value

N_tr = 100
p = 10000
n_sig = 50

N_te = 1000
N_Rep = 100

out_dat = list()
sel_var = rep(0, p)

n_clu = 5
cl = makeSOCKcluster(n_clu)
registerDoSNOW(cl)

pb = txtProgressBar(max = N_Rep, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

# multiple processing
result = foreach(j = 1:N_Rep, .packages = c('zoo'), .options.snow = opts) %dopar%{
  
  set.seed(j)
  seed = j
  tr = sim_data(N_tr, p)
  te = sim_data(N_te, p)

  fit = nbayes(tr$X, tr$Y)
  err_no = predict_nbayes(te$X, te$Y, fit$prob, fit$prior)$miss
  
  # save the running time
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
  if(dim(result_tab)[2] == 1){
     result_tab = cbind(0, result_tab)}

  # result.tab = table(te$Y, pred$Yhat)
  prec_class = diag(result_tab)/colSums(result_tab)
  prec_class[is.na(prec_class)] = 0
  precision =  mean(prec_class)
  
  recal_class = diag(result_tab)/rowSums(result_tab)
  recall = mean(recal_class)
  f1_score = 2 * (precision * recall) / (precision + recall)
  
  out_dat[[j]] = c(seed, err_no, err_sel, n_sel, tp, fp, precision, recall, f1_score, opt_alpha, time)
  #return(list(out_dat = out_dat, sel_var_list = sel_var))  
  
}

stopCluster(cl)
# Sys.time()

# summary the result
tmp_dat = result %>% do.call(rbind, .) %>% data.frame(., stringsAsFactors = F)
names(tmp_dat) = c('seed', 'err_no' , 'err_sel', 'n_sel', 'tp', 'fp', 'precision', 'recall',
                   'f1_score', 'opt_alpha', 'time')
apply(tmp_dat, 2, function(x) mean(x, na.rm = T)) %>% round(., 4)
apply(tmp_dat, 2, function(x) sd(x, na.rm = T)/sqrt(N_Rep)) %>% round(., 4)

#변수생성

sim_data = function(N, p)
{
  Y = rbinom(N, size = 1, prob = 0.5)
  N_1 = sum(Y == 1)
  N_0 = N - N_1

  X = matrix(0, N, p)
  r = 10
  for (i in 1:r) {
    X[Y == 0, i] = sample(1:3, size = N_0, replace = TRUE, prob = c(0.3, 0.3, 0.4))
    X[Y == 1, i] = sample(1:3, size = N_1, replace = TRUE, prob = c(0.5, 0.4, 0.2))

    X[Y == 0, r + i] = sample(1:3, size = N_0, replace = TRUE, prob = c(0.4, 0.2, 0.4))
    X[Y == 1, r + i] = sample(1:3, size = N_1, replace = TRUE, prob = c(0.3, 0.4, 0.3))

    X[Y == 0, 2 * r + i] = sample(1:3, size = N_0, replace = TRUE, prob = c(0.2, 0.5, 0.3))
    X[Y == 1, 2 * r + i] = sample(1:3, size = N_1, replace = TRUE, prob = c(0.1, 0.5, 0.4))

    X[Y == 0, 3 * r + i] = sample(1:3, size = N_0, replace = TRUE, prob = c(0.2, 0.3, 0.5))
    X[Y == 1, 3 * r + i] = sample(1:3, size = N_1, replace = TRUE, prob = c(0.2, 0.5, 0.3))

    X[Y == 0, 4 * r + i] = sample(1:3, size = N_0, replace = TRUE, prob = c(0.3, 0.3, 0.4))
    X[Y == 1, 4 * r + i] = sample(1:3, size = N_1, replace = TRUE, prob = c(0.3, 0.4, 0.3))
  }
  for (i in 51:p)
    X[, i] = sample(1:3, size = N, replace = TRUE)
  Y = Y + 1
  return(list(X = X, Y = Y))
}



######################################################################################################

nbayes = function(X, Y)
{ 
  X = as.matrix(X)
  N = length(Y)
  p = ncol(X)
  N_class = length(unique(Y))
  Y_class = sort(unique(c(Y)))
  X_class = sort(unique(c(X)))
  
  if(length(p) == 0){
    result = list(ind = 0)  
  } else{
    cnt = prob = list()
    df = chi = p_val = rep(0, p)
    prior = as.numeric(table(Y) / N)
    #i=1
    for (i in 1:p)
    {
      cnt[[i]] = table(Y, X[, i])                              # p번째 변수와 Y의 table
      
      for(q in 1:length(X_class))
      {
        if(q == 1){
          if(colnames(cnt[[i]])[q] != q || is.na(colnames(cnt[[i]])[q]))
          { cname = colnames(cnt[[i]])
          cnt[[i]] = cbind(0, cnt[[i]])
          colnames(cnt[[i]]) = sort(c(cname, q))
          }
        }
        
        if(q != 1 & q != length(X_class)){
          if(colnames(cnt[[i]])[q] != q ||is.na(colnames(cnt[[i]])[q]))
          { cname = colnames(cnt[[i]])
          cnt[[i]] = cbind(cnt[[i]][, c(1:(q - 1))], 0, cnt[[i]][, -c(1:(q - 1))])
          colnames(cnt[[i]]) = sort(c(cname, q))
          }    
        }
        
        if(q == length(X_class)){
          if(colnames(cnt[[i]])[q] != q ||is.na(colnames(cnt[[i]])[q]))
          { cname = colnames(cnt[[i]])
          cnt[[i]] = cbind(cnt[[i]][, c(1:(q - 1))], 0)
          colnames(cnt[[i]]) = sort(c(cname, q))
          }
        }
      }
      
    }
    
    cnt = lapply(1:p, function(x) cnt[[x]] + 1)     # correct zero probability
  
    ##
    for(i in 1:p) 
    {
      df[i] = (N_class - 1) * (length(unique(X[, i])) - 1)     # 자유도
      chi[i] = ifelse(any(dim(cnt[[i]]) %in% 1) , 0, chisq.test(cnt[[i]])$statistic)
      p_val[i] = 1 - pchisq(chi[i], df[i])                       # p-value 단측
      prob[[i]] = cnt[[i]] / apply(cnt[[i]], 1, sum)             # 비율 table
    }

    ##
    ind = sort.int(chi, decreasing = TRUE, index.return = TRUE)$ix
    sort_chi = sort.int(chi, decreasing = TRUE, index.return = TRUE)$x
    
    result = list(prior = prior, prob = prob, df = df,
                  index = ind, chi = sort_chi, pval = p_val[ind])  #정렬해서 출력
  }
  return(result)
}


######################################################################################################

cptnb = function(X, Y)
{
  X = as.matrix(X)
  N = nrow(X)
  p = ncol(X)
  
  fit = nbayes(X, Y)
  cfit = cpt.var(fit$pval, method = "PELT")@cpts
  if(length(cfit) == 0){return(list(var = NA))
  }else{
    
    fit_index = fit$index
    cv_err = c()
    for (r in 1:length(cfit)) 
    { 
      a = cfit[r]
      var_index = fit_index[1:a]
      
      fold = ifelse(min(table(Y)) < 10, min(table(Y))/2, 10)
      
      cv_err[r] = mean(cv_naive(X[, var_index], Y, n_fold = fold)$test_error)
    }

    sel_cpt = cfit[which.min(cv_err)[1]]
    sel_var = fit$index[1:sel_cpt]
    
    sel_fit = nbayes(X[, sel_var], Y)
    
    return(list(prior = sel_fit$prior, prob = sel_fit$prob, var = sel_var, chi_stat = sel_fit$pval))
    
  }
}


############################################################################################

cv_naive = function(X, Y, n_fold = 10)
{
  X = as.matrix(X)
  N = nrow(X)
  p = ncol(X)
  
  if (p == 0){return(list(test_error = 0))
  } else{
    
    fold = data_split(Y, n_fold = n_fold)
    cv = c()

    for (i in 1:n_fold)
    {
      fit = nbayes(X[fold != i,], Y[fold != i]) #train:test = 9:1
      cv[i] = predict_nbayes(X[fold == i,], Y[fold == i], fit$prob, fit$prior)$miss
      cv_se = sd(cv)/sqrt(n_fold)
    }

    return(list(fold = fold, test_error = cv, test_error_se = cv_se))  #test error는 평균
  }
}

######################################################################################################

data_split = function(Y, n_fold)
{
  
  k = length(unique(Y))
  N = length(Y)
  if (min(Y) == 0)
    y = Y + 1
  else 
    y = Y
    class_size = table(y)
    ran = rep(0, N) 
  if ( (min(class_size) < n_fold) & (n_fold != N) )
  {
    warning('decrease your fold size!\n')
    return(NULL)
  }
  if ( min(class_size) >= n_fold )
  {
    for (j in 1:k)
      ran[y == j] = ceiling(sample(class_size[j])/(class_size[j] + 1) * n_fold) 
  } else if (n_fold == N)
    ran = 1:N
  
  return(ran)
}


######################################################################################################

predict_nbayes = function(X, Y, prob, prior)
{ 
  X = as.matrix(X)
  N = nrow(X)
  p = ncol(X)
  N_class = length(unique(Y))
  Y_class = sort(unique(c(Y)))
  
  est = matrix(log(prior), N, N_class, byrow=TRUE)
  for (i in 1:N) 
  {
    for (j in 1:p)    
    {
      est[i, ] = est[i, ] + as.numeric(log(prob[[j]][, X[i, j]]))
    }
  }
  
  Yhat = apply(est, 1, which.max) - 1

  sub_Yhat = Yhat
  for(i in 1:length(sort(unique(Yhat)))){
    which_uni = which(Yhat == sort(unique(Yhat))[i])
    sub_Yhat[which_uni] = Y_class[i]
  }
  
  miss = mean(Y != sub_Yhat)
  
  return(list(Yhat = Yhat, miss = miss))
}

# table(Y)


######################################################################################################

fdr_fit = function(X, Y, Prop = 0.5, alpha)
{
  result = list()
  X = as.matrix(X)
  N = nrow(X)
  p = ncol(X)
  fit = nbayes(X, Y)
  
  #FDR 
  sort_var = sort(fit$pval)               #p-value 오름차순 정렬
  ind = c(1:p)
  adj_alpha = ind * alpha / N #조정된 alpha
  #  plot(fit$pval, cex=.5, pch=16, xlab="p", ylab="p-value", main="p=1000, n=100")  
  #  lines(adj.alpha[1:100])
  
  #FDR result
  thr = which(!(sort_var <= adj_alpha))[1] - 1
  thr = ifelse(is.na(thr), p, thr)
  thr = ifelse(thr != 0, thr, 0)
  # 유의한 변수 갯수
  
  if(thr > 0){ 
    sel_var = fit$index[1:thr] 
    sel_fit = nbayes(X[, sel_var], Y)
    result = list(prior = sel_fit$prior, prob = sel_fit$prob, alpha = alpha,
                  selected_val = thr, var = sel_var, chi = sel_fit$pval)
  }else{
    result = list(selected_val = 0, var = 0)}
  
  return(result)
}



############################################################################################

bhfdr = function(X, Y)
{

  # 알파 선택
  alpha_pilot = exp(seq(log(1e-4), log(1), length.out = 20))
  
  the_num_alpha = 0
  BH_list = list()   # 변수선택한 결과 저장하는 리스트
  
  for(i in 1:length(alpha_pilot)){

    BH_list[[i]] = fdr_fit(X, Y, alpha = alpha_pilot[i])
    the_num_alpha[i] = ifelse(BH_list[[i]]$selected_val > 1, 1, 0)
  }
  #########################################

  pilot_id = which(the_num_alpha == 1)
  pilot_id = pilot_id[!is.na(pilot_id)]
  
  #모형선택
  var_list = list()   # 선택한 변수 저장하는 리스트
  cv_err = c(0, length(pilot_id))
  for (i in 1:length(pilot_id)){
    var_list[[i]] = BH_list[[pilot_id[i]]]$var
    cv_err[i] = mean(cv_naive(X[, var_list[[i]]], Y, n_fold = 10)$test_error)
  }
  
  # optimal variable
  id = min(which(cv_err == min(cv_err)))
  sel_var = var_list[[id]]
  
  #1 standard error rule
  err_se = cv_naive(X[, sel_var], Y, n_fold = 10)$test_error_se 
  st1_err = cv_err[id] + err_se # minimum vd err + standard error of optimal model
  
  id = min(which(cv_err <= st1_err))    # 최종 선택
  
  opt_alpha = alpha_pilot[pilot_id][id]
  
  # feature selection
  sel_fit = nbayes(X[, sel_var], Y) 
  return(list(prior = sel_fit$prior, prob = sel_fit$prob, var = sel_var, p_val = sel_fit$pval, 
              opt_alpha = opt_alpha))
}



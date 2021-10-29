####Analysis of Covariation#########################
####Michelle Nixon##################################
####April 12, 2021##################################

var.explained.Y <- function(posterior, otu.closed){
  posterior = to_proportions(posterior)
  posterior = to_alr(posterior, ncategories(posterior))
  Y.pred = predict(posterior, response = "LambdaX", from_scratch = FALSE)
  Eta.inv = array(0, dim = c(dim(otu.closed), 2000))
  for(i in 1:dim(Y.pred)[3]){
    Eta.inv[,,i] = alrInv_array(Y.pred[,,i], d = nrow(Y.pred) + 1, 1)
  }
  Ypred = array(0, dim=dim(Eta.inv))
  size = as.data.frame(colSums(posterior$Y))
  for (i in 1:dim(Ypred)[3]){
    for (j in 1:dim(Ypred)[2]){
      Ypred[,j,i] <- rmultinom(1, size[j,] ,Eta.inv[,j,i])
    }
  }
  exp.mat = matrix(NA, nrow = dim(Ypred)[3], ncol = 1)
  for(i in 1:dim(Y.pred)[3]){
    Y.pred.closed = miniclo_array(Ypred[,,i], parts = 1)
    var_fit = sum(apply(Y.pred.closed, MARGIN = 1, FUN = "var"))
    var_res = sum(apply(Y.pred.closed - otu.closed, MARGIN = 1, FUN = "var"))
    exp.mat[i,1] = var_fit/ (var_fit + var_res)
    
  }
  return(exp.mat)
}#end of function


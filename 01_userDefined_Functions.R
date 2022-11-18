####Analysis of Covariation#########################
####Michelle Nixon##################################
####April 12, 2021##################################
var.explained.eta <- function(posterior){
  eta.data <- posterior$Eta
  eta.pred <- predict(posterior, response = "LambdaX")
  exp.mat <- matrix(NA, nrow=dim(eta.pred)[3],ncol = 1)
  
  for(i in 1:dim(eta.pred)[3]){
    var_fit = sum(apply(eta.pred[,,i], MARGIN = 1, FUN = "var"))
    var_res = sum(apply(eta.pred[,,i] - eta.data[,,i], MARGIN = 1, FUN = "var"))
    var_total = sum(apply(eta.data[,,i], MARGIN = 1, FUN = "var"))
    exp.mat[i,1] = var_fit/ (var_fit + var_res)
  }
  return(exp.mat)
}


summary.mcmc_hsar <- function(object, ...)
{
  x<-object
  cat("\nCall:\n")
  print(x$call)
  
  cat("Type:", ' hsar ', "\n")
  
  cat("\nCoefficients:\n")
  print( put_labels_to_coefficients(x$Mbetas, x$labels) )
  
  cat("\n Spatial Coefficients:\n")
  print( cbind( rho= x$Mrho, lambda=x$Mlambda) )
  
  cat("\n Diagnostics \n")
  cat("Deviance information criterion (DIC):", x$DIC, "\n")
  cat("Effective number of parameters (pd):", x$pD, "\n")
  cat("Log likelihood:", x$Log_Likelihood, "\n")
  cat("Pseudo R squared:", x$R_Squared, "\n")
  
  cat("\n Impacts:\n")
  df <- as.data.frame( cbind( t(x$impact_direct), t(x$impact_indirect), t(x$impact_total) ) )
  names(df)<-c("direct","indirect","total")
  row.names(df)<- x$labels
  
  print( df )
  
  invisible(x)
}

summary.mcmc_sar <- function(object, ...)
{
  x<-object
  cat("\nCall:\n")
  print(x$call)
  cat("Type:", ' sar ', "\n")
  
  cat("\nCoefficients:\n")
  print( put_labels_to_coefficients(x$Mbetas, x$labels) )
  
  rho<-x$Mrho
  names(rho)<-'rho'
  cat("\n Spatial Coefficients:\n")
  print( rho )
  
  cat("\n Diagnostics \n")
  cat("Deviance information criterion (DIC):", x$DIC, "\n")
  cat("Effective number of parameters (pd):", x$pD, "\n")
  cat("Log likelihood:", x$Log_Likelihood, "\n")
  cat("Pseudo R squared:", x$R_Squared, "\n")
  
  cat("\n Impacts:\n")
  df <- as.data.frame( cbind( t(x$impact_direct), t(x$impact_indirect), t(x$impact_total) ) )
  names(df)<-c("direct","indirect","total")
  row.names(df)<- x$labels
  
  print( df )
  
  invisible(x)
}

summary.mcmc_hsar_rho_0 <- function(object, ...)
{
  x <- object 
  cat("\nCall:\n")
  print(x$call)
  cat("Type:", ' hsar with rho = 0 ', "\n")
  
  cat("\nCoefficients:\n")
  print( put_labels_to_coefficients(x$Mbetas, x$labels) )
  
  lambda<-x$Mlambda
  names(lambda)<-'lambda'
  cat("\n Spatial Coefficients:\n")
  print( lambda )
  
  cat("\n Diagnostics \n")
  cat("Deviance information criterion (DIC):", x$DIC, "\n")
  cat("Effective number of parameters (pd):", x$pD, "\n")
  cat("Log likelihood:", x$Log_Likelihood, "\n")
  cat("Pseudo R squared:", x$R_Squared, "\n")
  invisible(x)
}

summary.mcmc_hsar_lambda_0 <- function(object, ...)
{
  x <- object
  cat("\nCall:\n")
  print(x$call)
  cat("Type:", ' hsar with lambda = 0 ', "\n")
  
  cat("\nCoefficients:\n")
  print( put_labels_to_coefficients(x$Mbetas, x$labels) )
  
  rho<-x$Mrho
  names(rho)<-'rho'
  cat("\n Spatial Coefficients:\n")
  print( rho )
  
  cat("\n Diagnostics \n")
  cat("Deviance information criterion (DIC):", x$DIC, "\n")
  cat("Effective number of parameters (pd):", x$pD, "\n")
  cat("Log likelihood:", x$Log_Likelihood, "\n")
  cat("Pseudo R squared:", x$R_Squared, "\n")
  invisible(x)
}


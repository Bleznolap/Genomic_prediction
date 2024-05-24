
#' @param N The size of individuals in each population 
#' @param ... An argument that contains the genetic variance and gwa object [SNP effects and their variance]. 
#' Note that the two objects are arranged in that order.
#' @returns a list of the following objects:
#' \itemize{
#'     \item{\code{results}} {A list of objects including Z score and p-value for each population}
#'     \item{\code{Z_Meta}} {A vector of combined z scores}
#'     \item{\code{p_Meta}} {A vector of combined p-value}
#'     \item{\code{beta_meta}} {A vector of combined effect size estimates}
#'     \item{\code{SE_meta}} {A vector of combined standard error}
#'     \item{\code{ghat}} {A vector of combined SNP effects}
#'     \item{\code{var_ghat}} {A vector of combined SNP effects variances}
#'     \item{\code{table}} {A table of standard error, z score, p-value, beta, and SNP effect}
#'     \item{\code{g_var_Meta}} {A vector of combined genetic variances}
#' }

#' @details
#' This function performs a meta-analysis


meta_gp <- function(N, ...){
  num <- N
  objs <- list(...)
  n_combinations <- length(objs) / 2
  
  results <- list()
  all_beta <- NULL
  all_wgt <- NULL
  all_var <- NULL
  all_num <- NULL
  
  for (i in 1:n_combinations) {
    g_var <- objs[[2 * i - 1]]
    gwa <- objs[[2 * i]]
    num_y <- num[[i]]-1
    
    varB<- (g_var)^2/gwa$var_ghat
    wgt <- 1/varB
    wgt[is.na(wgt)] <- 0
    SE <- 1/sqrt(wgt)
    beta <- (gwa$ghat*varB)/g_var
    beta[is.na(beta)] <- 0
    Z_score<- gwa$ghat/sqrt(gwa$var_ghat)
    Z_score[is.na(Z_score)] <- 0
    p_value <- 2 * (1 - pnorm(abs(Z_score)))
    
    if (is.null(all_beta)) {
      all_beta <- beta * wgt
    } else {
      all_beta <- all_beta + (beta * wgt)
    }
    
    if (is.null(all_wgt)) {
      all_wgt <- wgt
    } else {
      all_wgt <- all_wgt + wgt
    }
    
    if (is.null(all_var)) {
      all_var <- num_y*g_var
    } else {
      all_var <- all_var + (num_y*g_var)
    }
    
    if (is.null(all_num)) {
      all_num <- num_y
    } else {
      all_num <- all_num + num_y
    }
    
    results[[i]] <- list(varB=varB,wgt=wgt, SE=SE, beta=beta,
                         Z_score=Z_score, p_value=p_value, g_var=g_var)
  }
  
  var_g <- all_var/all_num
  beta_Meta <- all_beta/all_wgt
  SE_Meta <- sqrt(1/all_wgt)
  ghat_Meta <- (var_g)*(beta_Meta/(SE_Meta^2)) 
  var_ghat_Meta <- var_g^2/(SE_Meta^2)
  Z_Meta <- beta_Meta/SE_Meta
  p_Meta <- 2 * (1 - pnorm(abs(Z_Meta)))
  table_M <- {combined_table <- cbind(head(SE_Meta), head(Z_Meta),
                                      head(p_Meta), head(beta_Meta), head(ghat_Meta))
  
  kable(combined_table, escape = FALSE, caption = "Meta analysis table (First 6 rows)", 
        col.names = c("SE_meta", "Z_meta",
                      "p_meta","beta_meta", "ghat_meta")) %>%
    kable_styling(latex_options = "scale down")}
  
  return(list(results=results, Z_Meta=Z_Meta, p_Meta=p_Meta, beta_Meta=beta_Meta, SE_Meta=SE_Meta, 
              ghat= ghat_Meta, var_ghat=var_ghat_Meta,table=table_M, g_var_Meta=var_g))
}

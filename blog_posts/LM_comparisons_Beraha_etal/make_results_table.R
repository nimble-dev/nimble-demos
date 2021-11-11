load("results_C_BFG.RData")
load("results_C_better.RData")
load("results_conj_bin_BFG.RData")
load("results_conj_bin_better.RData")
load("results_WI.RData")
load("results_WI_better.RData")
load("results_NI.RData")
load("results_NI_better.RData")
load("results_L.RData")
load("results_L_better.RData")


summarize_one_result <- function(result) {
  cols <- colnames(result$Psamples[[1]])
  beta_cols <- grepl("beta", cols)
  meanESS <- mean(coda::effectiveSize(result$Psamples[[1]][,beta_cols]))
  Ns <- nrow(result$Psamples[[1]]) # Number of saved samples
  Nit <- result$niter # Different cases have different ratio of niter / nburnin
  if(is.null(Nit)) Nit <- 10000
  c(ESS_per_Ns = meanESS/Ns, Nit_per_time = Nit/sum(result$t2[1:2]), meanESS_per_time = meanESS/sum(result$t2[1:2]))
}


##

res_conj <- cbind(
  do.call('rbind', lapply(results_C_BFG, summarize_one_result)),
  do.call('rbind', lapply(results_C_better, summarize_one_result))
)

res_conj_bin <- cbind(
  do.call('rbind', lapply(results_conj_bin_BFG, summarize_one_result)),
  do.call('rbind', lapply(results_conj_bin_better, summarize_one_result))
)

res_WI <- cbind(
  do.call('rbind', lapply(results_WI, summarize_one_result)),
  do.call('rbind', lapply(results_WI_better, summarize_one_result))
)

res_NI <- cbind(
  do.call('rbind', lapply(results_NI, summarize_one_result)),
  do.call('rbind', lapply(results_NI_better, summarize_one_result))
)

res_L <- cbind(
  do.call('rbind', lapply(results_L, summarize_one_result)),
  do.call('rbind', lapply(results_L_better, summarize_one_result))
)

all_results <- 
  do.call('rbind',
          list(res_conj,
               res_conj_bin,
               res_WI,
               res_NI,
               res_L))


all_results <- cbind(all_results, all_results[,6]/all_results[,3])

colnames(all_results) <- c(rep( c("ESS/Ns", "Nit/t", "ESS/t"), 2), "Better by")


## -----------------------------------------------------------------------------
library(kableExtra)
kbl(all_results, digits = 2) %>% 
  kable_classic() %>%
  add_header_above(c(" " = 1, "BFG" = 3, "Better code" = 3, "Improvement" = 1)) %>%
  pack_rows("LM-C", 1, 6) %>%
  pack_rows("LM-C Bin", 7, 11) %>%
  pack_rows("LM-WI", 12, 15)  %>%
  pack_rows("LM-NI", 16, 20)  %>%
  pack_rows("LM-Lasso", 21, 28)

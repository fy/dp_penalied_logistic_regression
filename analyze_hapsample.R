###############################################################################
## read in data
###############################################################################

case_path = "./data/case_genotypes_Nov09_interaction1_MAF025_1"
control_path = "./data/anticase_genotypes_Nov09_interaction1_MAF025_1"


case_data = read.table(case_path, header=FALSE, sep='\t', as.is=TRUE)
case_data = data.frame(case_data)
rownames(case_data) = case_data[, 2]
control_data = read.table(control_path, header=FALSE, sep='\t', as.is=TRUE)
control_data = data.frame(control_data)


## transpose the data then combine cases and controls into a single dataframe
case_nn = ncol(case_data)-4
snp_nn = nrow(case_data)

control_data_t = data.frame(t(control_data[,  -(1:4)]))
names(control_data_t) = control_data[, 2]
case_data_t = data.frame(t(case_data[,  -(1:4)]))
names(case_data_t) = case_data[, 2]

all_data = rbind(control_data_t,
                 case_data_t[, match(names(case_data_t),
                                     names(control_data_t))])
all_data = cbind(status=as.factor(rep(c("control", "case"),
                                      c(case_nn, case_nn))), all_data)

###############################################################################
## calculate the chi-square statistics and p-values for each snp
###############################################################################

chisq_test_results = list()
for (ii in 1:(ncol(all_data)-1)) {
  if (length(unique(all_data[, 1+ii])) == 1) {
    chisq_test_results[[ii]] = NULL
  } else {
    chisq_test_results[[ii]] = chisq.test(all_data[, 1+ii],
                                          all_data[, 1],
                                          correct=FALSE)
  }
}
chisq_stat = sapply(chisq_test_results, function(ee) {
  if (is.null(ee)) return(NA)
  if (ee$parameter[['df']] != 2) return(NA)
  return(ee$statistic)
})
chisq_stat = as.vector(chisq_stat)
chisq_pval = sapply(chisq_test_results, function(ee) {
  if (is.null(ee)) return(NA)
  if (ee$parameter[['df']] != 2) return(NA)
  return(ee[['p.value']])
})
chisq_pval = as.vector(chisq_pval)


###############################################################################
## sample training and testing data
###############################################################################

if (TRUE) {
  set.seed(101)
  training_indiv = sample(1:case_nn, case_nn / 2)
}
testing_indiv = c(1:case_nn)[-training_indiv]

training_data = all_data[c(training_indiv, case_nn + training_indiv), ]
testing_data = all_data[c(testing_indiv, case_nn + testing_indiv), ]

###############################################################################
## get valid snps and valid chisquare statistics
###############################################################################
valid_snps = names(all_data)[-1][!is.na(chisq_stat)]
valid_snp_chisq = chisq_stat[!is.na(chisq_stat)]



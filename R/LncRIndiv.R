LncRIndiv <-
function (lncRNA, normal_expr, cancer_expr, freq = 0.95) {
    Lgene = dim(normal_expr)[1]
	
    normal_expr_rank = NULL
    for (i in 1:dim(normal_expr)[2]) {
        tmp_expr = sort(normal_expr[, i])
        tmp_rank = match(normal_expr[, i], tmp_expr)
        normal_expr_rank = cbind(normal_expr_rank, tmp_rank)
    }
    cancer_expr_rank = NULL
    for (j in 1:dim(cancer_expr)[2]) {
        tmp_expr = sort(cancer_expr[, j])
        tmp_rank = match(cancer_expr[, j], tmp_expr)
        cancer_expr_rank = cbind(cancer_expr_rank, tmp_rank)
    }
	
    GenePair_sig = matrix(1, nrow = Lgene, ncol = Lgene)
    la = lfactorial(1:(ncol(normal_expr)+ncol(cancer_expr)))
    for (k in 1:(Lgene - 1)) {
		N_tmp = matrix(t(normal_expr[k, ] - t(normal_expr[-c(1:k), ])), nrow = Lgene - k, byrow = F)
        C_tmp = matrix(t(cancer_expr[k, ] - t(cancer_expr[-c(1:k), ])), nrow = Lgene - k, byrow = F)
        tmp_reverse = cbind(rowSums(N_tmp > 0), rowSums(N_tmp < 0), rowSums(C_tmp > 0), rowSums(C_tmp  < 0))
        GenePair_sig_tmp = vector(length=dim(tmp_reverse)[1])
          
        GenePair_sig_tmp = .C("FisherExactTest",as.integer(t(tmp_reverse)),as.integer(length(tmp_reverse)),as.numeric(c(0,la)),p.value=as.numeric(rep(0,nrow(tmp_reverse))))$p.value
        
        GenePair_sig[k, (k + 1):Lgene] = GenePair_sig_tmp
        GenePair_sig[(k + 1):Lgene, k] = GenePair_sig_tmp
    }
    rm(N_tmp, C_tmp, tmp_reverse)
	
    p = NULL
    for (m in 1:(Lgene - 1)) {
        p = c(p, GenePair_sig[m, (m + 1):Lgene])
    }
    fdr_result = p.adjust(p, "BH")
	
    FDR_TMP = matrix(0, nrow = Lgene, ncol = Lgene)
    for (n in 1:(Lgene - 1)) {
        FDR_TMP[n, (n + 1):Lgene] = fdr_result[1:(Lgene - n)]
        fdr_result = fdr_result[-c(1:(Lgene - n))]
    }
    FDR_result = FDR_TMP + t(FDR_TMP)
    for (s in 1:Lgene) {
        FDR_result[s, s] = 1
    }
	
    pair_info = NULL
    for (k in 1:Lgene) {
        Nnorm = normal_expr[k, ]
        Ncase = cancer_expr[k, ]
        direction = sign(mean(Ncase) - mean(Nnorm))
        N_tmp = matrix(rep(Nnorm, Lgene), ncol = dim(normal_expr)[2], by = T) - normal_expr
        if (direction == -1) {
            control_p = rowSums(N_tmp > 0)
            NR = control_p/dim(normal_expr)[2]
            Nloc_up = which(NR > freq)
            tmp_FDR = FDR_result[, k]
            tmp_FDR = tmp_FDR[Nloc_up]
            tmp_P = GenePair_sig[, k]
            tmp_P = tmp_P[Nloc_up]
            tmp_lnc = lncRNA[Nloc_up]
            if (length(Nloc_up) == 1) {
                tmp_norm_expr = matrix(normal_expr[Nloc_up, ], nrow = 1)
                tmp_case_expr = matrix(cancer_expr[Nloc_up, ], nrow = 1)
                tmp_norm_rank = matrix(normal_expr_rank[Nloc_up, ], nrow = 1)
                tmp_case_rank = matrix(cancer_expr_rank[Nloc_up, ], nrow = 1)
                tmp_direction = sign(rowMeans(tmp_case_expr) - rowMeans(tmp_norm_expr))
                index = which(tmp_direction == -1 & tmp_FDR < 0.1)
                if (length(index) > 0) {
                  tmp_FDR = tmp_FDR[index]
                  tmp_lnc = tmp_lnc[index]
                  tmp_P = tmp_P[index]
                  if (length(index) == 1) {
                    tmp_norm_rank = matrix(tmp_norm_rank[index, ], nrow = 1)
                    tmp_case_rank = matrix(tmp_case_rank[index, ], nrow = 1)
                  } else if (length(index) > 1) {
                    tmp_norm_rank = tmp_norm_rank[index, ]
                    tmp_case_rank = tmp_case_rank[index, ]
                  }
                  tmp_expr_rank = cbind(tmp_norm_rank, tmp_case_rank)
                  tmp_pair_info = matrix(NA, 3, 5)
                  if (length(tmp_lnc) > 0 & length(tmp_lnc) <= 3) {
                    tmp_pair_info[1:length(tmp_lnc), 1] = rep(lncRNA[k], length(tmp_lnc))
                    tmp_pair_info[1:length(tmp_lnc), 2] = tmp_lnc
                    tmp_pair_info[1:length(tmp_lnc), 3] = rep("down", length(tmp_lnc))
                    tmp_pair_info[1:length(tmp_lnc), 4] = tmp_P
                    tmp_pair_info[1:length(tmp_lnc), 5] = tmp_FDR
                    pair_info = rbind(pair_info, tmp_pair_info[1:length(tmp_lnc), ])
                  } else if (length(index) > 3) {
                    coef_var = apply(tmp_expr_rank, 1, var)
                    increase_order = sort(coef_var)
                    increase_order = match(increase_order, coef_var)
                    index = increase_order[1:3]
                    tmp_FDR = tmp_FDR[index]
                    tmp_lnc = tmp_lnc[index]
                    tmp_P = tmp_P[index]
                    tmp_pair_info[, 1] = rep(lncRNA[k], 3)
                    tmp_pair_info[, 2] = tmp_lnc
                    tmp_pair_info[, 3] = rep("down", 3)
                    tmp_pair_info[, 4] = tmp_P
                    tmp_pair_info[, 5] = tmp_FDR
                    pair_info = rbind(pair_info, tmp_pair_info)
                  }
                }
            } else if (length(Nloc_up) > 1) {
                tmp_norm_expr = normal_expr[Nloc_up, ]
                tmp_case_expr = cancer_expr[Nloc_up, ]
                tmp_norm_rank = normal_expr_rank[Nloc_up, ]
                tmp_case_rank = cancer_expr_rank[Nloc_up, ]
                tmp_direction = sign(rowMeans(tmp_case_expr) - rowMeans(tmp_norm_expr))
				
                index = which(tmp_direction == -1 & tmp_FDR < 0.1)
                if (length(index) > 0) {
                  tmp_FDR = tmp_FDR[index]
                  tmp_lnc = tmp_lnc[index]
                  tmp_P = tmp_P[index]
                  if (length(index) == 1) {
                    tmp_norm_rank = matrix(tmp_norm_rank[index, ], nrow = 1)
                    tmp_case_rank = matrix(tmp_case_rank[index, ], nrow = 1)
                  } else if (length(index) > 1) {
                    tmp_norm_rank = tmp_norm_rank[index, ]
                    tmp_case_rank = tmp_case_rank[index, ]
                  }
                  tmp_expr_rank = cbind(tmp_norm_rank, tmp_case_rank)
				  
                  tmp_pair_info = matrix(NA, 3, 5)
                  if (length(tmp_lnc) > 0 & length(tmp_lnc) <= 3) {
                    tmp_pair_info[1:length(tmp_lnc), 1] = rep(lncRNA[k], length(tmp_lnc))
                    tmp_pair_info[1:length(tmp_lnc), 2] = tmp_lnc
                    tmp_pair_info[1:length(tmp_lnc), 3] = rep("down", length(tmp_lnc))
                    tmp_pair_info[1:length(tmp_lnc), 4] = tmp_P
                    tmp_pair_info[1:length(tmp_lnc), 5] = tmp_FDR
                    pair_info = rbind(pair_info, tmp_pair_info[1:length(tmp_lnc), ])
                  } else if (length(index) > 3) {
                    coef_var = apply(tmp_expr_rank, 1, var)
                    increase_order = sort(coef_var)
                    increase_order = match(increase_order, coef_var)
                    index = increase_order[1:3]
                    tmp_FDR = tmp_FDR[index]
                    tmp_lnc = tmp_lnc[index]
                    tmp_P = tmp_P[index]
                    tmp_pair_info[, 1] = rep(lncRNA[k], 3)
                    tmp_pair_info[, 2] = tmp_lnc
                    tmp_pair_info[, 3] = rep("down", 3)
                    tmp_pair_info[, 4] = tmp_P
                    tmp_pair_info[, 5] = tmp_FDR
                    pair_info = rbind(pair_info, tmp_pair_info)
                  }
                }
            }
        } else if (direction == 1) {
            control_pp = rowSums(N_tmp < 0)
            NR = control_pp/dim(normal_expr)[2]
            Nloc_down = which(NR > freq)
            tmp_FDR = FDR_result[, k]
            tmp_FDR = tmp_FDR[Nloc_down]
            tmp_P = GenePair_sig[, k]
            tmp_P = tmp_P[Nloc_down]
            tmp_lnc = lncRNA[Nloc_down]
            if (length(Nloc_down) == 1) {
                tmp_norm_expr = matrix(normal_expr[Nloc_down, ], nrow = 1)
                tmp_case_expr = matrix(cancer_expr[Nloc_down, ], nrow = 1)
                tmp_norm_rank = matrix(normal_expr_rank[Nloc_down, ], nrow = 1)
                tmp_case_rank = matrix(cancer_expr_rank[Nloc_down, ], nrow = 1)
                tmp_direction = sign(rowMeans(tmp_case_expr) - rowMeans(tmp_norm_expr))
                index = which(tmp_direction == 1 & tmp_FDR < 0.1)
                if (length(index) > 0) {
                  tmp_FDR = tmp_FDR[index]
                  tmp_lnc = tmp_lnc[index]
                  tmp_P = tmp_P[index]
                  if (length(index) == 1) {
                    tmp_norm_rank = matrix(tmp_norm_rank[index, ], nrow = 1)
                    tmp_case_rank = matrix(tmp_case_rank[index, ], nrow = 1)
                  } else if (length(index) > 1) {
                    tmp_norm_rank = tmp_norm_rank[index, ]
                    tmp_case_rank = tmp_case_rank[index, ]
                  }
                  tmp_expr_rank = cbind(tmp_norm_rank, tmp_case_rank)
                  tmp_pair_info = matrix(NA, 3, 5)
                  if (length(tmp_lnc) > 0 & length(tmp_lnc) <= 3) {
                    tmp_pair_info[1:length(tmp_lnc), 1] = rep(lncRNA[k], length(tmp_lnc))
                    tmp_pair_info[1:length(tmp_lnc), 2] = tmp_lnc
                    tmp_pair_info[1:length(tmp_lnc), 3] = rep("up", length(tmp_lnc))
                    tmp_pair_info[1:length(tmp_lnc), 4] = tmp_P
                    tmp_pair_info[1:length(tmp_lnc), 5] = tmp_FDR
                    pair_info = rbind(pair_info, tmp_pair_info[1:length(tmp_lnc), ])
                  } else if (length(index) > 3) {
                    coef_var = apply(tmp_expr_rank, 1, var)
                    increase_order = sort(coef_var)
                    increase_order = match(increase_order, coef_var)
                    index = increase_order[1:3]
                    tmp_FDR = tmp_FDR[index]
                    tmp_lnc = tmp_lnc[index]
                    tmp_P = tmp_P[index]
                    tmp_pair_info[, 1] = rep(lncRNA[k], 3)
                    tmp_pair_info[, 2] = tmp_lnc
                    tmp_pair_info[, 3] = rep("up", 3)
                    tmp_pair_info[, 4] = tmp_P
                    tmp_pair_info[, 5] = tmp_FDR
                    pair_info = rbind(pair_info, tmp_pair_info)
                  }
                }
            } else if (length(Nloc_down) > 1) {
                tmp_norm_expr = normal_expr[Nloc_down, ]
                tmp_case_expr = cancer_expr[Nloc_down, ]
                tmp_norm_rank = normal_expr_rank[Nloc_down, ]
                tmp_case_rank = cancer_expr_rank[Nloc_down, ]
                tmp_direction = sign(rowMeans(tmp_case_expr) - rowMeans(tmp_norm_expr))
                index = which(tmp_direction == 1 & tmp_FDR < 0.1)
                if (length(index) > 0) {
                  tmp_FDR = tmp_FDR[index]
                  tmp_lnc = tmp_lnc[index]
                  tmp_P = tmp_P[index]
                  if (length(index) == 1) {
                    tmp_norm_rank = matrix(tmp_norm_rank[index, ], nrow = 1)
                    tmp_case_rank = matrix(tmp_case_rank[index, ], nrow = 1)
                  } else if (length(index) > 1) {
                    tmp_norm_rank = tmp_norm_rank[index, ]
                    tmp_case_rank = tmp_case_rank[index, ]
                  }
                  tmp_expr_rank = cbind(tmp_norm_rank, tmp_case_rank)
                  tmp_pair_info = matrix(NA, 3, 5)
                  if (length(tmp_lnc) > 0 & length(tmp_lnc) <= 3) {
                    tmp_pair_info[1:length(tmp_lnc), 1] = rep(lncRNA[k], length(tmp_lnc))
                    tmp_pair_info[1:length(tmp_lnc), 2] = tmp_lnc
                    tmp_pair_info[1:length(tmp_lnc), 3] = rep("up", length(tmp_lnc))
                    tmp_pair_info[1:length(tmp_lnc), 4] = tmp_P
                    tmp_pair_info[1:length(tmp_lnc), 5] = tmp_FDR
                    pair_info = rbind(pair_info, tmp_pair_info[1:length(tmp_lnc), ])
                  } else if (length(index) > 3) {
                    coef_var = apply(tmp_expr_rank, 1, var)
                    increase_order = sort(coef_var)
                    increase_order = match(increase_order, coef_var)
                    index = increase_order[1:3]
                    tmp_FDR = tmp_FDR[index]
                    tmp_lnc = tmp_lnc[index]
                    tmp_P = tmp_P[index]
                    tmp_pair_info[, 1] = rep(lncRNA[k], 3)
                    tmp_pair_info[, 2] = tmp_lnc
                    tmp_pair_info[, 3] = rep("up", 3)
                    tmp_pair_info[, 4] = tmp_P
                    tmp_pair_info[, 5] = tmp_FDR
                    pair_info = rbind(pair_info, tmp_pair_info)
                  }
                }
            }
        }
    }
	if(length(pair_info)==0){
		return(NULL)
	}else{
		colnames(pair_info) = c("DE lncRNAs", "paired lncRNAs", "DE direction", "P-value", "FDR")
		return(pair_info)
	}
}

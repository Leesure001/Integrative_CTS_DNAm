# 1. Install Requierd Packages
## 1.1 For Cran packages
install_cran_pkgs <- function(pkgs){
  not_installed <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(not_installed)){
    install.packages(not_installed, dependencies = TRUE) # you may get a warning for ewas here, just ignore that
    if("ewastools" %in% pkgs & "ewastools" %in% not_installed){
      devtools::install_github("hhhh5/ewastools@master")
    }
  }else{
    message("Cran packages you need are installed")
  }
}

## 1.2 For bioconductor packages
install_bioconductor_pkgs <- function(pkgs){
  bioc_not_installed <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(bioc_not_installed)){
    BiocManager::install(c(bioc_not_installed))
  }else{
    message("Bioconductor packages you need are installed")
  }
}

## 1.3 Run test to check if all packages are installed
check_installed <- function(pkgs){
  not_installed <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(not_installed)){
    message("Package(s) not installed: ", not_installed)
  }else{
    message("All good ...")
  }
}
# If not, figure out the problem and make sure all packages are istalled
# Ideas are from https://github.com/PGC-PTSD-EWAS/EPIC_QC/tree/main

# 2. Modified Version For HIRE (initialized using inputting cell type proportions)
HIRE_V2 <- function(Ometh, X, num_celltype, props_true, tol = 10^(-5), num_iter=1000, alpha=0.01){
  library("HIREewas")
  library("quadprog")
	#generate samples from a Dirichlet distribution
	rDirichlet <- function(alpha_vec){
		num <- length(alpha_vec)
		temp <- rgamma(num, shape = alpha_vec, rate = 1)
		return(temp / sum(temp))
	}

	#initialize P_t

	CorDescent <- function(MethMatr, num_celltype, props_true, tol = 0.01, showIter = FALSE){
		err0 <- 0
		err1 <- 1000
		m <- nrow(MethMatr) #number of CpG sites
		n <- ncol(MethMatr) #number of samples
		K <- num_celltype
		if(m < n){
			stop("The CpG site number must be larger than the sample number!")
		}
		#initialize P_matr
		P_matr_t <- vapply(seq_len(n), function(i){ rDirichlet(rep(2,K)) }, FUN.VALUE=rep(-1,K))
		while(abs(err1 - err0) >= tol){
			err0 <- err1
			#update U_matr
			Dmat <- 2*P_matr_t%*%t(P_matr_t)
			Amat <- cbind(diag(rep(1,K)), diag(rep(-1,K)))
			bvec <- c(rep(0,K), rep(-1,K))
			U_matr_t <- t( vapply(seq_len(m), function(j){
							dvec <- 2*P_matr_t %*% as.numeric(MethMatr[j,])

							solu <- solve.QP(Dmat, dvec, Amat, bvec)
							solu$solution
					}, FUN.VALUE=rep(-1,K)) )
		
			#update P_matr
			Dmat <- 2*t(U_matr_t) %*% U_matr_t
			Amat <- cbind(matrix(1, K, K), diag(rep(1,K)))
			bvec <- c(rep(1, K), rep(0, K))
			P_matr_t <- vapply(seq_len(n), function(i){
						dvec <- 2 * t(U_matr_t) %*% as.numeric(MethMatr[ ,i])
						solu <- solve.QP(Dmat, dvec, Amat, bvec, meq = K)
						solu$solution 
					}, FUN.VALUE=rep(-1,K))
					
			#calculate errors
			err1 <- sum((MethMatr - U_matr_t %*% P_matr_t)^2)
			if(showIter == TRUE){
				message("  ", err1, "\n")
			}
		}
        P_matr_t <- props_true
		# true_prop matrix: cell types for each row, sample names for each column
	
		return(list(U=U_matr_t, P=P_matr_t))	
	}

	Initialize <- function(Ometh, num_celltype, props_true){
		K <- num_celltype
		sdrow <- apply(Ometh, 1, sd)
		ind <- order(sdrow, decreasing = TRUE)
		m <- nrow(Ometh)
		if(m <= 1000){
			num_cpg_for_init <- m
		}else{
			num_cpg_for_init <- max(3*ncol(Ometh), floor(m/10))
			if(num_cpg_for_init > m){
				num_cpg_for_init <- m
			}
		}

		Ometh_part <- Ometh[ind[seq_len(num_cpg_for_init)],] #select CpG sites with the most num_cpg_for_init variant methylation levels 

		#result <- CorDescent(Ometh_part, num_celltype=K, props_true, tol = 0.1, showIter = FALSE)
		#P_initial <- result$P
    P_initial <- props_true

		mu_initial <- vapply(seq_len(m), function(j){
				if(K > 2){
					fit <- lm(Ometh[j,]~as.matrix(t(P_initial[-1, ])))
				}else{
					fit <- lm(Ometh[j,]~as.numeric(P_initial[-1, ]))
				}
				tmp <- as.numeric(summary(fit)$coeff[ ,1])
				tmp[-1] <- tmp[1] + tmp[-1]
				tmp
			}, FUN.VALUE = rep(-1,K) )
		return(list(P_initial, mu_initial))
	}

	EmEwasRcallC <- function(Ometh, X, P_t, mu_t, beta_t, sig_sqTiss_t, sig_sqErr_t, tol, num_iter){
		args <- list("P_init"=as.numeric(P_t), "mu_init"=as.numeric(mu_t), "beta_init"=as.numeric(beta_t),
				"beta_init_dim"=as.integer(dim(beta_t)), "Ometh_r"=as.numeric(Ometh),
				"Ometh_r_dim"=as.integer(dim(Ometh)), "X_r"=as.numeric(X), "X_r_dim"=as.integer(dim(X)),
				"sig_sqTiss_init"=as.numeric(sig_sqTiss_t), "sig_sqErr_init"=as.numeric(sig_sqErr_t),
				"tol_r" = as.numeric(tol), "num_iter" = as.integer(num_iter))
		ret_list <- .Call("EmEwasRcallC", args)		
		return(ret_list)
	}


	m <- nrow(Ometh) #CpG site number
	n <- ncol(Ometh) #sample number
	p <- nrow(X)
	K <- num_celltype
	
	P_t <- matrix(-1, K, n)
	mu_t <- matrix(-1, m, K)
	beta_t <- array(-1, dim=c(m, K, p))
	sig_sqTiss_t <- matrix(-1, m, K)
	sig_sqErr_t <- rep(-1, m)

	init <- Initialize(Ometh, K, props_true)
	P_t <- init[[1]]
	mu_t <- t(init[[2]])
	message("  Initialization Done.\n")
	message("  Implementing EM algorithm... \n")
	ret_list <- EmEwasRcallC(Ometh, X, P_t, mu_t, beta_t, sig_sqTiss_t, sig_sqErr_t, tol, num_iter)
	message("  Done! \n")
	
	message("  Calculating p-values...\n")
	tmp <- NULL
	for(ell in seq_len(p)){
		tmp <- cbind(tmp, X[ell, ]*t(ret_list$P_t))
	}
	x_matr <- cbind( tmp, t(ret_list$P_t)[, seq(2,K)])
	x_matr <- as.matrix(x_matr)

	pvalues <- t(vapply(seq_len(m), function(j){
					y_vec <- Ometh[j,]
					fit <- lm(y_vec~x_matr)
					summary(fit)$coef[seq(2, (1+p*K)),4]
				}, FUN.VALUE = rep(-1,p*K)))
	message("  Done!\n")
	ret_list[[7]] <- pvalues
	names(ret_list)[7] <- "pvalues"
	names(ret_list)[6] <- "pBIC"
	d <- (K-1)*n + m*(1+2*K+p*K) #number of parameters
	d0 <- sum(ret_list$pvalues > alpha/(p*m*K)) #number of zero parameters
	ret_list[[6]] <- ret_list[[6]] - log(n)*d + log(n)*(d-d0) 
	
	return(ret_list)
}

# 3. Calculate counts
count_non_zero <- function(v){
    count <- sum(v != 0)
    return(count)
}

count_one <- function(v){
    count <- sum(v == 1)
    return(count)
}

count_neg_one <- function(v){
    count <- sum(v == -1)
    return(count)
}

# 4. Get DMCT list
get_dmct <- function(coe, adjPThresh = 0.05){
    library(matrixStats)
    dmct <- matrix(rep(0, length(coe) * nrow(coe[[1]])), ncol = length(coe))
    idx <- which(sapply(coe, "[[", "adjP") < adjPThresh)
    dmct[idx] <- sign(sapply(coe, "[[", "Estimate")[idx])
    dmc <- ifelse(rowAlls(dmct == 0), 0, 1)
    colnames(dmct) <- names(coe)
    rownames(dmct) <- rownames(coe[[1]])
    dmct_total <- matrix(apply(dmct, 2, count_non_zero), nrow = 1)
    dmct_hyper <- matrix(apply(dmct, 2, count_one), nrow = 1)
    dmct_hypo <- matrix(apply(dmct, 2, count_neg_one), nrow = 1)
    dmct_count <- rbind(dmct_total, dmct_hyper, dmct_hypo)
    colnames(dmct_count) <- names(coe)
    rownames(dmct_count) <- c("dmct_total", "dmct_hyper", "dmct_hypo")
    # all(dmct_hypo + dmct_hyper == dmct_total)
    dmct <- cbind(dmc, dmct)
    colnames(dmct) <- c("DMC", names(coe))
    rownames(dmct) <- rownames(coe[[1]])
    return(list(dmct = dmct, dmct_count = dmct_count))
}

# 5. Extract differentially methylated CpGs from DMCT list (The cell types in this function need to be modified case by case.)
extract_cpgs <- function (dmct_list){
    dmct <- dmct_list$dmct
    dmct_count <- dmct_list$dmct_count
    cell_type <- c("Epi", "Fib", "CD4T", "Mono", "B")
    library(assertthat)
    assert_that(all(colnames(dmct_count) == cell_type), msg = "Can not extract results! You need to modify cell types in `extract_cpgs` function based on your case!")
    dmc <- rownames(dmct)[which(dmct[, 1] == 1)]

    epi_all <- rownames(dmct)[which(dmct[, 2] != 0)]
    epi_hyper <- rownames(dmct)[which(dmct[, 2] == 1)]
    epi_hypo <- rownames(dmct)[which(dmct[, 2] == -1)]
    epi <- list(epi_all = epi_all, epi_hyper = epi_hyper, epi_hypo = epi_hypo)

    fib_all <- rownames(dmct)[which(dmct[, 3] != 0)]
    fib_hyper <- rownames(dmct)[which(dmct[, 3] == 1)]
    fib_hypo <- rownames(dmct)[which(dmct[, 3] == -1)]
    fib <- list(fib_all = fib_all, fib_hyper = fib_hyper, fib_hypo = fib_hypo)

    cd4t_all <- rownames(dmct)[which(dmct[, 4] != 0)]
    cd4t_hyper <- rownames(dmct)[which(dmct[, 4] == 1)]
    cd4t_hypo <- rownames(dmct)[which(dmct[, 4] == -1)]
    cd4t <- list(cd4t_all = cd4t_all, cd4t_hyper = cd4t_hyper, cd4t_hypo = cd4t_hypo)

    mono_all <- rownames(dmct)[which(dmct[, 5] != 0)]
    mono_hyper <- rownames(dmct)[which(dmct[, 5] == 1)]
    mono_hypo <- rownames(dmct)[which(dmct[, 5] == -1)]
    mono <- list(mono_all = mono_all, mono_hyper = mono_hyper, mono_hypo = mono_hypo)

    b_all <- rownames(dmct)[which(dmct[, 6] != 0)]
    b_hyper <- rownames(dmct)[which(dmct[, 6] == 1)]
    b_hypo <- rownames(dmct)[which(dmct[, 6] == -1)]
    b <- list(b_all = b_all, b_hyper = b_hyper, b_hypo = b_hypo)

    all_result <- list(dmct = dmct, dmct_count = dmct_count, dmc = dmc, epi = epi, fib = fib, cd4t = cd4t, mono = mono, b = b)
    return(all_result)
}

# 6. Generate coe list for TCA (The cell types in this function need to be modified case by case.)
generate_tca_coe <-function(tca_output, adjPMethod = "fdr"){
    cell_type <- c("Epi", "Fib", "CD4T", "Mono", "B")
    dmct_tca_estimate <- tca_output$gammas_hat
    dmct_tca_pvalue <- tca_output$gammas_hat_pvals
    dmct_tca_estimate <- dmct_tca_estimate[, seq(1, dim(dmct_tca_estimate)[2], dim(dmct_tca_estimate)[2] / length(cell_type))]
    dmct_tca_pvalue <- dmct_tca_pvalue[, seq(1, dim(dmct_tca_pvalue)[2], dim(dmct_tca_pvalue)[2] / length(cell_type))]
    tca_cell_type <- sub("\\.phe$", "", colnames(dmct_tca_estimate))
    library(assertthat)
    assert_that(all(tca_cell_type == cell_type), msg = "Can not extract results! You need to modify cell types in `generate_tca_coe` function based on your case!")

    for (i in seq_along(cell_type)){
        assign(cell_type[i], cbind(dmct_tca_estimate[, i], dmct_tca_pvalue[, i]))
    }
    for (i in seq_along(cell_type)){
        tmp <- get(cell_type[i])
        tmp <- cbind(tmp, p.adjust(tmp[, 2], method = adjPMethod))
        assign(cell_type[i], tmp)
    }
    col.names <- c("Estimate", "p", "adjP")
    for (i in seq_along(cell_type)){
        tmp <- get(cell_type[i])
        colnames(tmp) <- col.names
        tmp <- data.frame(tmp)
        assign(cell_type[i], tmp)
    }
    return(list(Epi = Epi, Fib = Fib, CD4T = CD4T, Mono = Mono, B = B))
}

# 7. Generate coe list for TOAST (The cell types in this function need to be modified case by case.)
generate_toast_coe <-function(toast_output, adjPMethod = "fdr"){
    load('sample_450k.RData') # Load your own data to adjust CpG order from TOAST output.
    cell_type <- c("Epi", "Fib", "CD4T", "Mono", "B")
    library(assertthat)
    assert_that(all(names(toast_output)[-length(names(toast_output))] == cell_type), msg = "Can not extract results! You need to modify cell types in `generate_toast_coe` function based on your case!")
    beta <- c()
    for (i in seq_along(cell_type)){
      beta <- cbind(beta, toast_output[[i]][rownames(sample_450k), 1])
    }
    rownames(beta) <- rownames(sample_450k)
    colnames(beta) <- cell_type

    p.value <- c()
    for (i in seq_along(cell_type)){
      p.value <- cbind(p.value, toast_output[[i]][rownames(sample_450k), 6])
    }
    rownames(p.value) <- rownames(sample_450k)
    colnames(p.value) <- cell_type
    
    cell_type <- colnames(beta)
    for (i in seq_along(cell_type)){
        assign(cell_type[i], cbind(beta[, cell_type[i]], p.value[, cell_type[i]]))
    }
    for (i in seq_along(cell_type)){
        tmp <- get(cell_type[i])
        tmp <- cbind(tmp, p.adjust(tmp[, 2], method = adjPMethod))
        assign(cell_type[i], tmp)
    }
    col.names <- c("Estimate", "p", "adjP")
    for (i in seq_along(cell_type)){
        tmp <- get(cell_type[i])
        colnames(tmp) <- col.names
        tmp <- data.frame(tmp)
        assign(cell_type[i], tmp)
    }
    return(list(Epi = Epi, Fib = Fib, CD4T = CD4T, Mono = Mono, B = B))
}

# 8. Generate coe list for CeDAR (The cell types in this function need to be modified case by case.)
generate_cedar_coe <- function(cedar_output, adjPMethod = "fdr"){
    cell_type <- c("Epi", "Fib", "CD4T", "Mono", "B")
    library(assertthat)
    assert_that(all(names(cedar_output$toast_res) == cell_type), msg = "Can not extract results! You need to modify cell types in `generate_cedar_coe` function based on your case!")
    beta <- c()
    row_names <- rownames(cedar_output$toast_res[[1]])
    for (i in seq_along(cell_type)){
      beta <- cbind(beta, cedar_output$toast_res[[i]][row_names, 1])
    }
    rownames(beta) <- row_names
    colnames(beta) <- cell_type

    find_pp_cpgs <- function(cell.pp){
        temp1 <- sort(cell.pp)
        temp2 <- cumsum(temp1)
        temp3 <- which(temp2 < 0.05) # Change here if you want to change FDR control for CeDAR
        #temp3 <- which(temp2 < 0.01)
        pp_cpgs <- names(temp3)
        return (pp_cpgs)
    }
    temp.pp <- 1 - as.matrix(cedar_output$tree_res$full$pp)
    p.value <- matrix(1, nrow = nrow(temp.pp), ncol = ncol(temp.pp))
    rownames(p.value) <- rownames(temp.pp)
    colnames(p.value) <- colnames(temp.pp)
    epi_cpgs <- find_pp_cpgs(temp.pp[, 1])
    fib_cpgs <- find_pp_cpgs(temp.pp[, 2])
    cd4t_cpgs <- find_pp_cpgs(temp.pp[, 3])
    mono_cpgs <- find_pp_cpgs(temp.pp[, 4])
    b_cpgs <- find_pp_cpgs(temp.pp[, 5])
    p.value[epi_cpgs, 1] =  p.value[fib_cpgs, 2] = p.value[cd4t_cpgs, 3] = p.value[mono_cpgs, 4] = p.value[b_cpgs, 5] <- 0.01 # Change this number less than FDR value above.
    #p.value[epi_cpgs, 1] =  p.value[fib_cpgs, 2] = p.value[cd4t_cpgs, 3] = p.value[mono_cpgs, 4] = p.value[b_cpgs, 5] <- 0.001

    cell_type <- colnames(beta)
    for (i in seq_along(cell_type)){
        assign(cell_type[i], cbind(beta[, cell_type[i]], p.value[, cell_type[i]]))# the second column p is already FDR p values
    }
    for (i in seq_along(cell_type)){
        tmp <- get(cell_type[i])
        tmp <- cbind(tmp, temp.pp[, i])
        assign(cell_type[i], tmp)
    }
    col.names <- c("Estimate", "adjP", "p") # Since the second column is already adjP, so we set it column as adjP, the third column we set 1 - pp.
    for (i in seq_along(cell_type)){
        tmp <- get(cell_type[i])
        colnames(tmp) <- col.names
        tmp <- data.frame(tmp)
        assign(cell_type[i], tmp)
    }
    return(list(Epi = Epi, Fib = Fib, CD4T = CD4T, Mono = Mono, B = B))
}

# 9. Generate coe list for HIRE V2 (The cell types in this function need to be modified case by case.)
generate_hire_v2_coe <- function(hire_output, adjPMethod = "fdr"){
    load('sample_450k.RData') # Load your own data to adjust CpG order from HIRE output.
    cell_type <- c("Epi", "Fib", "CD4T", "Mono", "B")
    library(assertthat)
    assert_that(ncol(hire_output$beta_t) == length(cell_type), msg = "Can not extract results! You need to modify cell types in `generate_hire_v2_coe` function based on your case!")

    beta_t_hire <- hire_output$beta_t
    beta_t_hire_disease <- beta_t_hire[, , 1]
    rownames(beta_t_hire_disease) <- rownames(sample_450k)
    colnames(beta_t_hire_disease) <- cell_type

    p_value_hire <- hire_output$pvalues[, 1:length(cell_type)]
    rownames(p_value_hire) <- rownames(sample_450k)
    colnames(p_value_hire) <- cell_type

    for (i in seq_along(cell_type)){
        assign(cell_type[i], cbind(beta_t_hire_disease[, cell_type[i]], p_value_hire[, cell_type[i]]))
    }
    for (i in seq_along(cell_type)){
        tmp <- get(cell_type[i])
        tmp <- cbind(tmp, p.adjust(tmp[, 2], method = adjPMethod))
        assign(cell_type[i], tmp)
    }
    col.names <- c("Estimate", "p", "adjP")
    for (i in seq_along(cell_type)){
        tmp <- get(cell_type[i])
        colnames(tmp) <- col.names
        tmp <- data.frame(tmp)
        assign(cell_type[i], tmp)
    }
    return(list(Epi = Epi, Fib = Fib, CD4T = CD4T, Mono = Mono, B = B))
}

# 10. Integrated analysis with averaging p values - V1 (The cell types in this function need to be modified case by case.)
generate_overall_coe_v1 <-function(celldmc_coe, tca_coe, toast_coe, cedar_coe, hire_v2_coe, select = 5, adjPMethod = "fdr"){
    cell_type <- c("Epi", "Fib", "CD4T", "Mono", "B")
    library(assertthat)
    assert_that(all(names(celldmc_coe) == cell_type), msg = "Can not extract results! You need to modify cell types in `generate_overall_coe_v1` function based on your case!")

    est_coeff <- c()
    row_names <- rownames(celldmc_coe[[1]])
    for (i in seq_along(cell_type)){
      est_coeff <- cbind(est_coeff, celldmc_coe[[i]][row_names, 1])
    }
    rownames(est_coeff) <- row_names
    colnames(est_coeff) <- cell_type

    if (select == 5){
        p.vals.epi <- cbind(celldmc_coe$Epi$p, tca_coe$Epi$p, toast_coe$Epi$p, cedar_coe$Epi$p, hire_v2_coe$Epi$p)
        p.vals.fib <- cbind(celldmc_coe$Fib$p, tca_coe$Fib$p, toast_coe$Fib$p, cedar_coe$Fib$p, hire_v2_coe$Fib$p)
        p.vals.cd4t <- cbind(celldmc_coe$CD4T$p, tca_coe$CD4T$p, toast_coe$CD4T$p, cedar_coe$CD4T$p, hire_v2_coe$CD4T$p)
        p.vals.mono <- cbind(celldmc_coe$Mono$p, tca_coe$Mono$p, toast_coe$Mono$p, cedar_coe$Mono$p, hire_v2_coe$Mono$p)
        p.vals.b <- cbind(celldmc_coe$B$p, tca_coe$B$p, toast_coe$B$p, cedar_coe$B$p, hire_v2_coe$B$p)
    }

    if (select == 4){
        p.vals.epi <- cbind(celldmc_coe$Epi$p, tca_coe$Epi$p, toast_coe$Epi$p, cedar_coe$Epi$p)
        p.vals.fib <- cbind(celldmc_coe$Fib$p, tca_coe$Fib$p, toast_coe$Fib$p, cedar_coe$Fib$p)
        p.vals.cd4t <- cbind(celldmc_coe$CD4T$p, tca_coe$CD4T$p, toast_coe$CD4T$p, cedar_coe$CD4T$p)
        p.vals.mono <- cbind(celldmc_coe$Mono$p, tca_coe$Mono$p, toast_coe$Mono$p, cedar_coe$Mono$p)
        p.vals.b <- cbind(celldmc_coe$B$p, tca_coe$B$p, toast_coe$B$p, cedar_coe$B$p)
    }

    if (select == 31){
        p.vals.epi <- cbind(celldmc_coe$Epi$p, tca_coe$Epi$p, toast_coe$Epi$p)
        p.vals.fib <- cbind(celldmc_coe$Fib$p, tca_coe$Fib$p, toast_coe$Fib$p)
        p.vals.cd4t <- cbind(celldmc_coe$CD4T$p, tca_coe$CD4T$p, toast_coe$CD4T$p)
        p.vals.mono <- cbind(celldmc_coe$Mono$p, tca_coe$Mono$p, toast_coe$Mono$p)
        p.vals.b <- cbind(celldmc_coe$B$p, tca_coe$B$p, toast_coe$B$p)
    }

    if (select == 32){
        p.vals.epi <- cbind(celldmc_coe$Epi$p, tca_coe$Epi$p, cedar_coe$Epi$p)
        p.vals.fib <- cbind(celldmc_coe$Fib$p, tca_coe$Fib$p, cedar_coe$Fib$p)
        p.vals.cd4t <- cbind(celldmc_coe$CD4T$p, tca_coe$CD4T$p, cedar_coe$CD4T$p)
        p.vals.mono <- cbind(celldmc_coe$Mono$p, tca_coe$Mono$p, cedar_coe$Mono$p)
        p.vals.b <- cbind(celldmc_coe$B$p, tca_coe$B$p, cedar_coe$B$p)
    }

    if (select == 2){
        p.vals.epi <- cbind(celldmc_coe$Epi$p, tca_coe$Epi$p)
        p.vals.fib <- cbind(celldmc_coe$Fib$p, tca_coe$Fib$p)
        p.vals.cd4t <- cbind(celldmc_coe$CD4T$p, tca_coe$CD4T$p)
        p.vals.mono <- cbind(celldmc_coe$Mono$p, tca_coe$Mono$p)
        p.vals.b <- cbind(celldmc_coe$B$p, tca_coe$B$p)
    }

    q.vals.epi <- qnorm(p.vals.epi)
    q.vals.fib <- qnorm(p.vals.fib)
    q.vals.cd4t <- qnorm(p.vals.cd4t)
    q.vals.mono <- qnorm(p.vals.mono)
    q.vals.b <- qnorm(p.vals.b)


    meanq.vals.epi <- apply(q.vals.epi, 1, mean)
    meanq.vals.fib <- apply(q.vals.fib, 1, mean)
    meanq.vals.cd4t <- apply(q.vals.cd4t, 1, mean)
    meanq.vals.mono <- apply(q.vals.mono, 1, mean)
    meanq.vals.b <- apply(q.vals.b, 1, mean)


    p.vals <- cbind(pnorm(meanq.vals.epi), pnorm(meanq.vals.fib), pnorm(meanq.vals.cd4t), pnorm(meanq.vals.mono), pnorm(meanq.vals.b))
    rownames(p.vals) <- row_names
    colnames(p.vals) <- cell_type

    for (i in seq_along(cell_type)){
        assign(cell_type[i], cbind(est_coeff[, cell_type[i]], p.vals[, cell_type[i]]))
    }
    for (i in seq_along(cell_type)){
        tmp <- get(cell_type[i])
        tmp <- cbind(tmp, p.adjust(tmp[, 2], method = adjPMethod))
        assign(cell_type[i], tmp)
    }
    col.names <- c("Estimate", "p", "adjP")
    for (i in seq_along(cell_type)){
        tmp <- get(cell_type[i])
        colnames(tmp) <- col.names
        tmp <- data.frame(tmp)
        assign(cell_type[i], tmp)
    }
    return(list(Epi = Epi, Fib = Fib, CD4T = CD4T, Mono = Mono, B = B))
}

# 11. Integrated analysis with minimizing p values - V2 (The cell types in this function need to be modified case by case.)
generate_overall_coe_v2 <-function(celldmc_coe, tca_coe, toast_coe, cedar_coe, hire_v2_coe, select = 5, adjPMethod = "fdr"){
    cell_type <- c("Epi", "Fib", "CD4T", "Mono", "B")
    library(assertthat)
    assert_that(all(names(celldmc_coe) == cell_type), msg = "Can not extract results! You need to modify cell types in `generate_overall_coe_v2` function based on your case!")

    est_coeff <- c()
    row_names <- rownames(celldmc_coe[[1]])
    for (i in seq_along(cell_type)){
      est_coeff <- cbind(est_coeff, celldmc_coe[[i]][row_names, 1])
    }
    rownames(est_coeff) <- row_names
    colnames(est_coeff) <- cell_type

    if (select == 5){
        p.vals.epi <- cbind(celldmc_coe$Epi$p, tca_coe$Epi$p, toast_coe$Epi$p, cedar_coe$Epi$p, hire_v2_coe$Epi$p)
        p.vals.fib <- cbind(celldmc_coe$Fib$p, tca_coe$Fib$p, toast_coe$Fib$p, cedar_coe$Fib$p, hire_v2_coe$Fib$p)
        p.vals.cd4t <- cbind(celldmc_coe$CD4T$p, tca_coe$CD4T$p, toast_coe$CD4T$p, cedar_coe$CD4T$p, hire_v2_coe$CD4T$p)
        p.vals.mono <- cbind(celldmc_coe$Mono$p, tca_coe$Mono$p, toast_coe$Mono$p, cedar_coe$Mono$p, hire_v2_coe$Mono$p)
        p.vals.b <- cbind(celldmc_coe$B$p, tca_coe$B$p, toast_coe$B$p, cedar_coe$B$p, hire_v2_coe$B$p)
    }

    if (select == 4){
        p.vals.epi <- cbind(celldmc_coe$Epi$p, tca_coe$Epi$p, toast_coe$Epi$p, cedar_coe$Epi$p)
        p.vals.fib <- cbind(celldmc_coe$Fib$p, tca_coe$Fib$p, toast_coe$Fib$p, cedar_coe$Fib$p)
        p.vals.cd4t <- cbind(celldmc_coe$CD4T$p, tca_coe$CD4T$p, toast_coe$CD4T$p, cedar_coe$CD4T$p)
        p.vals.mono <- cbind(celldmc_coe$Mono$p, tca_coe$Mono$p, toast_coe$Mono$p, cedar_coe$Mono$p)
        p.vals.b <- cbind(celldmc_coe$B$p, tca_coe$B$p, toast_coe$B$p, cedar_coe$B$p)
    }

    if (select == 31){
        p.vals.epi <- cbind(celldmc_coe$Epi$p, tca_coe$Epi$p, toast_coe$Epi$p)
        p.vals.fib <- cbind(celldmc_coe$Fib$p, tca_coe$Fib$p, toast_coe$Fib$p)
        p.vals.cd4t <- cbind(celldmc_coe$CD4T$p, tca_coe$CD4T$p, toast_coe$CD4T$p)
        p.vals.mono <- cbind(celldmc_coe$Mono$p, tca_coe$Mono$p, toast_coe$Mono$p)
        p.vals.b <- cbind(celldmc_coe$B$p, tca_coe$B$p, toast_coe$B$p)
    }

    if (select == 32){
        p.vals.epi <- cbind(celldmc_coe$Epi$p, tca_coe$Epi$p, cedar_coe$Epi$p)
        p.vals.fib <- cbind(celldmc_coe$Fib$p, tca_coe$Fib$p, cedar_coe$Fib$p)
        p.vals.cd4t <- cbind(celldmc_coe$CD4T$p, tca_coe$CD4T$p, cedar_coe$CD4T$p)
        p.vals.mono <- cbind(celldmc_coe$Mono$p, tca_coe$Mono$p, cedar_coe$Mono$p)
        p.vals.b <- cbind(celldmc_coe$B$p, tca_coe$B$p, cedar_coe$B$p)
    }

    if (select == 2){
        p.vals.epi <- cbind(celldmc_coe$Epi$p, tca_coe$Epi$p)
        p.vals.fib <- cbind(celldmc_coe$Fib$p, tca_coe$Fib$p)
        p.vals.cd4t <- cbind(celldmc_coe$CD4T$p, tca_coe$CD4T$p)
        p.vals.mono <- cbind(celldmc_coe$Mono$p, tca_coe$Mono$p)
        p.vals.b <- cbind(celldmc_coe$B$p, tca_coe$B$p)
    }

    p.vals <- cbind(apply(p.vals.epi, 1, function(x){pbeta(min(x), 1, length(x))}), apply(p.vals.fib, 1, function(x){pbeta(min(x), 1, length(x))}), apply(p.vals.cd4t, 1, function(x){pbeta(min(x), 1, length(x))}), apply(p.vals.mono, 1, function(x){pbeta(min(x), 1, length(x))}), apply(p.vals.b, 1, function(x){pbeta(min(x), 1, length(x))}))
    rownames(p.vals) <- row_names
    colnames(p.vals) <- cell_type

    for (i in seq_along(cell_type)){
        assign(cell_type[i], cbind(est_coeff[, cell_type[i]], p.vals[, cell_type[i]]))
    }
    for (i in seq_along(cell_type)){
        tmp <- get(cell_type[i])
        tmp <- cbind(tmp, p.adjust(tmp[, 2], method = adjPMethod))
        assign(cell_type[i], tmp)
    }
    col.names <- c("Estimate", "p", "adjP")
    for (i in seq_along(cell_type)){
        tmp <- get(cell_type[i])
        colnames(tmp) <- col.names
        tmp <- data.frame(tmp)
        assign(cell_type[i], tmp)
    }
    return(list(Epi = Epi, Fib = Fib, CD4T = CD4T, Mono = Mono, B = B))
}

# 12. Generate UpSet Plot (The cell types in this function need to be modified case by case.)
getList <- function(celldmc_coe, tca_coe, hire_coe, toast_coe, cedar_coe, FDR = 0.05){
    cell_type <- c("Epi", "Fib", "CD4T", "Mono", "B")
    library(assertthat)
    assert_that(all(names(celldmc_coe) == cell_type), msg = "Can not generate UpSet plot! You need to modify cell types in `getList` function based on your case!")
    cpgnames <- rownames(celldmc_coe[[1]])
    list1 <- list()
    for(i in seq_along(cell_type)){
        CellDMC = TCA = HIRE = TOAST = CeDAR <- rep(0, dim(celldmc_coe[[i]])[1])
        CellDMC[which(celldmc_coe[[i]]$adjP <= FDR)] <- 1
        TCA[which(tca_coe[[i]]$adjP <= FDR)] <- 1
        HIRE[which(hire_coe[[i]]$adjP <= FDR)] <- 1
        TOAST[which(toast_coe[[i]]$adjP <= FDR)] <- 1
        CeDAR[which(cedar_coe[[i]]$adjP <= FDR)] <- 1
        list1[[i]] <- list(CellDMC = cpgnames[CellDMC == 1], TCA = cpgnames[TCA == 1], HIRE = cpgnames[HIRE == 1], TOAST = cpgnames[TOAST == 1], CeDAR = cpgnames[CeDAR == 1])
    }
    names(list1) <- cell_type
    return(list1)
}
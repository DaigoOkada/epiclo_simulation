#snalysis (2023.5.1)
#source("/Users/dokada/Dropbox/analysis/2022.5/methy_ra_matome0517.R") #Run 5.23
set.seed(1000)
library(igraph)
library(glmnet)
library(parallel)
out_path <- "/Users/dokada/Desktop/work/methy_ra_matome0517/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}

Methy_data_simu2 <- function(n_markers=100, n_cells=40, ages=20:80, n_tissues=c(100, 100, 100, 100), n_sel_cell=10, n_cores=9, 
                        network_size=10, NSS=0.8, power=1, 
                        subset_prop=list(c(0.1, 0.4, 0.5), c(0.1, 0.4, 0.5), c(0.1, 0.4, 0.5), c(0.1, 0.4, 0.5)),
                        tissue_sens=list(c(1,1,1), c(1,1,1),c(1,1,1), c(1,1,1)),
                        r_com=0.1,
                        dam_prob1=5 * 10^(-3), dam_prob2=1 * 10^(-3), r_sensitive=0.5,
                        alpha=10, noise_sd_dyn=5,
                        max_cpg=500, cpg_window=100, init_cell_noise=1, measure_noise=10,
                        br_st=0.8, br_freq=0.1,
                        myseed=1000, out_file="/Users/dokada/Desktop/work/methy_rs_hard_ntpara/test.rds"){

    #sample setting
    set.seed(myseed)
    n_samples <- sum(n_tissues)
    n_tis_types <- length(n_tissues)
    n_cs_types <- unique(sapply(subset_prop,length))
    com_edge_size <- network_size * NSS
    spe_edge_size <- network_size * (1 - NSS)

    #Set damage sensevity of eacg genes(modify 2023/2/2)
    n_sen1 <- as.integer(round(n_markers * r_sensitive))
    n_sen2 <- n_markers - n_sen1
    fixed_prob_vec <- sample(c(rep(dam_prob1, n_sen1), rep(dam_prob2, n_sen2)))
    n_difsen_idx <- sample(n_markers, as.integer(round(n_markers * (1 - r_com))))
    prob_list <- list()
    for(i in 1: n_tis_types){
        tmp <- list()
        tmp[[1]] <- fixed_prob_vec * tissue_sens[[i]][1]
        for(j in 2:n_cs_types){
            prob_tmp <- fixed_prob_vec
            if(length(n_difsen_idx)>1){ #n_difsen_idxが整数1つのスカラーの時sample関数の挙動が変わる
                prob_tmp[n_difsen_idx] <- sample(fixed_prob_vec[n_difsen_idx])
                prob_tmp <- prob_tmp * tissue_sens[[i]][j]
            }
            tmp[[j]] <- prob_tmp
        }
        prob_list[[i]] <- tmp
    }

    #Epigenetic acc for samples
    dam_speed <- rnorm(n_samples, mean=1, sd=0.05)
    dam_speed[dam_speed < 0.75] <- 0.75
    dam_speed[dam_speed > 1.25] <- 1.25


    #Network of subsets
    #common network
    if(com_edge_size != 0){
        g_com <- sample_pa(n = n_markers, power = power, m=com_edge_size, directed = F)
        ad <- as.matrix(as_adjacency_matrix(g_com))
        adj_com <- t(ad) + ad
    }else{
        adj_com <- matrix(0, n_markers, n_markers)
    }
    off_mat_com <- matrix(rbinom(n_markers^2,size=1,prob=0.5), n_markers, n_markers)
    mgn_com <- matrix(rnorm(n_markers^2), n_markers, n_markers)
    causal_mat_com <- adj_com * off_mat_com * mgn_com

    #common subset specific network
    if(spe_edge_size != 0){
        g_comsubspe <- sample_pa(n = n_markers, power = power, m=spe_edge_size, directed = F)
        ad <- as.matrix(as_adjacency_matrix(g_comsubspe))
        adj_com <- t(ad) + ad
    }else{
        adj_com <- matrix(0, n_markers, n_markers)
    }
    off_mat_com <- matrix(rbinom(n_markers^2,size=1,prob=0.5), n_markers, n_markers)
    mgn_com <- matrix(rnorm(n_markers^2), n_markers, n_markers)
    causal_mat_g_comsubspe <- adj_com * off_mat_com * mgn_com

    #common subset network
    cmat_com_all <- causal_mat_com + causal_mat_g_comsubspe

    #subset -specific network
    cmat_cs_list <- list()
    for(i in 1: n_tis_types){
        tmp <- list()
        tmp[[1]] <- cmat_com_all
        for(j in 2:n_cs_types){
            if(spe_edge_size != 0){
                g <- sample_pa(n = n_markers, power = power, m=spe_edge_size, directed = F)
                ad <- as.matrix(as_adjacency_matrix(g))
                adj <- t(ad) + ad
            }else{
                adj <- matrix(0, n_markers, n_markers)
            }
            mgn <- matrix(rnorm(n_markers^2), n_markers, n_markers)
            off_mat <- matrix(rbinom(n_markers^2,size=1,prob=0.5), n_markers, n_markers)
            causal_mat <- adj * off_mat * mgn
            cmat_cs <- causal_mat_com + causal_mat
            tmp[[j]] <- cmat_cs
        }
        cmat_cs_list[[i]] <- tmp
    }
    #hist(rowSums(abs( cmat_cs_list[[1]][[2]])))

    #Set initial methylation values of each subset
    init_min <- as.integer(round(max_cpg/2)) - cpg_window
    init_max <- as.integer(round(max_cpg/2)) + cpg_window
    fixed_init_vec <- sample(init_min:init_max, size=n_markers, replace=T)
    init_list <- list()
    for(i in 1: n_tis_types){
        tmp <- list()
        tmp[[1]] <- fixed_init_vec
        for(j in 2:n_cs_types){
            init_tmp <- fixed_init_vec
            if(length(n_difsen_idx)>1){ #n_difsen_idxが整数1つのスカラーの時sample関数の挙動が変わる
                init_tmp[n_difsen_idx] <- sample(fixed_init_vec[n_difsen_idx])
            }
            tmp[[j]] <- init_tmp
        }
        init_list[[i]] <- tmp
    }
    #deb1 <- init_tmp[n_difsen_idx]
    #deb2 <- sample(fixed_init_vec[n_difsen_idx])
    #cat("n_difsen_idx:", n_difsen_idx, "\n")
    #cat("init_tmp:", init_tmp, "\n")
    #cat("fixed_init_vec:",fixed_init_vec, "\n")
    #cat("deb1:", deb1, "\n")
    #cat("deb2:",deb2, "\n")
    #init_list[[1]][[1]]==init_list[[2]][[1]]

    #Simulation
    sc_all_dam <- NULL #add single cel analysis
    sc_all_mat <- NULL #add single cel analysis
    tissue_vec <- rep(1:n_tis_types, times=n_tissues)
    n_samples <- length(tissue_vec)
    n_ages <- length(ages)
    int_dir <- paste0(out_path, "int_dir/")
    dir.create(int_dir)
    wrapper_simu <- function(k){
        set.seed(1000*myseed + k)
        out_file1 <- paste0(int_dir, "sample_bulk", k, ".rds")
        out_file2 <- paste0(int_dir, "sample_sc_dam", k, ".rds")
        out_file3 <- paste0(int_dir, "sample_sc_mat", k, ".rds")
        out_file4 <- paste0(int_dir, "dam_bulk", k, ".rds")
        tmp_tissue <- tissue_vec[k]
        tmp_acc <- dam_speed[k]
        celltypes <- sample(x=c(1, 2, 3), size=n_cells, replace=T, prob=subset_prop[[tmp_tissue]])

        #one human's analysis
        dam_sc_rec <- matrix(0, n_cells, n_ages) #add single cel analysis
        dat_onehuman <- array(rep(NA, n_cells*n_markers*n_ages), dim = c(n_cells, n_markers, n_ages))
        for(i in 1:n_cells){
            tmp_celltype <- celltypes[i]
            init_cpg_status <- init_list[[tmp_tissue]][[tmp_celltype]] + as.integer(round(rnorm(n_markers, 0, init_cell_noise)))
            dam_prob <- prob_list[[tmp_tissue]][[tmp_celltype]] * tmp_acc
            cmat_init <- cmat_cs_list[[tmp_tissue]][[tmp_celltype]] 
            init_cpg_status[init_cpg_status < 0] <- 0
            init_cpg_status[init_cpg_status > max_cpg] <- max_cpg
            dat_onehuman[i,,1] <- init_cpg_status
            dam_sc_rec[i,] <- 0 #add single cel analysis

            cmat_tmp <- cmat_init
            for(j in 1:(n_ages - 1)){    
                prev_val <- dat_onehuman[i,,j]

                #random de-methylation
                dam_cpgs_vec_tmp <- sapply(1:n_markers, function(g1){rbinom(n=1, size=prev_val[g1], prob=dam_prob[g1])}) 
                next_val <- prev_val - dam_cpgs_vec_tmp
                next_val[next_val < 0] <- 0
                next_val[next_val > max_cpg] <- max_cpg
                dam_cpgs_vec <- prev_val - next_val
                dam_sc_rec[i,j+1] <- dam_sc_rec[i,j] + sum(dam_cpgs_vec) #add single cel analysis
                if(any(dam_cpgs_vec < 0)) stop(paste("error on step", k, i))

                #causal network drive
                dam_cpgs_mat <- matrix(dam_cpgs_vec, nrow=1, ncol=length(dam_cpgs_vec))
                next_val2 <- next_val - as.integer(round(alpha * (dam_cpgs_mat %*% cmat_tmp) + rnorm(n_markers, 0, noise_sd_dyn)))
                next_val2[next_val2 < 0] <- 0
                next_val2[next_val2 > max_cpg] <- max_cpg
                dat_onehuman[i,,(j+1)] <- next_val2

                #causal network decline
                broken_edge <- matrix(sample(c(1, br_st), size=n_markers^2, replace=T, prob=c(1-br_freq, br_freq)), n_markers, n_markers)
                cmat_tmp <- cmat_tmp * broken_edge

                #cat(j,"\n")

            }
        }

        #dat_onehuman <- log1p(dat_onehuman)
        #sc_all[[k]] <- dat_onehuman
        bulk <- matrix(NA, n_ages, n_markers)
        for(j in 1:n_ages){
            bulk[j,] <- colMeans(dat_onehuman[,,j,drop=F]) + as.integer(round(rnorm(n_markers, 0, measure_noise)))
        }
        saveRDS(bulk, out_file1)

        #bulk damage record
        dam_bulk <- colSums(dam_sc_rec)
        saveRDS(dam_bulk, out_file4)
        
        #Single cell record
        n_sel_cell <- min(10, n_cells)
        tmp_idx <- sample(1:n_ages, 1)
        selected_cell <- sample(1:n_cells, n_sel_cell)
        tmp_dam <- cbind(dam_sc_rec[selected_cell,tmp_idx], rep(ages[tmp_idx], n_sel_cell), rep(tmp_tissue, n_sel_cell), celltypes[selected_cell], rep(k, n_sel_cell))
        tmp_mat <- dat_onehuman[selected_cell,,tmp_idx]
        saveRDS(tmp_dam, out_file2)
        saveRDS(tmp_mat, out_file3)
        return(0)
    }

    #Run
    check <- mclapply(1:n_samples, wrapper_simu, mc.cores=n_cores)
    if(!all(check==0)) stop(paste("error : wrapper_simu return non-zero value", myseed))

    #Get bulk result
    bulk_all <- matrix(NA, n_samples, n_markers)
    ages_all <- rep(NA, n_samples)
    sc_all_dam <- NULL
    sc_all_mat <- NULL
    bulk_dam_all <- rep(NA, n_samples)
    for(k in 1:n_samples){

        #path
        out_file1 <- paste0(int_dir, "sample_bulk", k, ".rds")
        out_file2 <- paste0(int_dir, "sample_sc_dam", k, ".rds")
        out_file3 <- paste0(int_dir, "sample_sc_mat", k, ".rds")
        out_file4 <- paste0(int_dir, "dam_bulk", k, ".rds")

        #bulk
        bulk <- readRDS(out_file1)
        bulk_dam <- readRDS(out_file4)
        tmp_idx <- sample(1:n_ages, 1)
        bulk_all[k, ] <- bulk[tmp_idx, ]/max_cpg
        ages_all[k] <- ages[tmp_idx]
        bulk_dam_all[k] <- bulk_dam[tmp_idx]

        #single cell
        sc_dam <- readRDS(out_file2)
        sc_mat <- readRDS(out_file3)
        sc_all_dam <- rbind(sc_all_dam, sc_dam)
        sc_all_mat <- rbind(sc_all_mat, sc_mat)
    }
    colnames(sc_all_dam) <- c("cell_damage", "donor_age", "donor_tissue", "cell_subset", "sample_row_idx")
    sc_all_mat <- sc_all_mat/max_cpg

    #All record
    dataset <- list(bulk_all, ages_all, bulk_dam_all, dam_speed, tissue_vec, causal_mat_com, cmat_cs_list, prob_list, init_list, sc_all_dam, sc_all_mat, subset_prop)
    names(dataset) <- c("bulk_all", "ages_all", "bulk_dam_all", "dam_speed", "tissue_vec", "causal_mat_com", "cmat_cs_list", "prob_list", "init_list", "sc_all_dam", "sc_all_mat", "subset_prop")
    saveRDS(dataset, out_file)
    unlink(int_dir, recursive=TRUE)
    return(out_file)
}

Epiclo_est2 <- function(out_file_list, out_dir){
    #Epigenetic clock
    types <- c("mixture", "separate")
    for(t in 1:2){
        out_res <- NULL
        for(i in 1:length(out_file_list)){
            set.seed(i)
            out_file <- out_file_list[[i]]
            seed_id <- strsplit(out_file,"/")[[1]][8]
            dataset <- readRDS(out_file)
            bulk_all <- dataset$bulk_all
            ages_all <- dataset$ages_all
            dam_speed <- dataset$dam_speed
            tissue_vec <- dataset$tissue_vec
            cmat_cs_list <- dataset$cmat_cs_list
            prob_list <- dataset$prob_list
            bulk_dam_all <- dataset$bulk_dam_all

            #Split data
            n_samples <- nrow(bulk_all)
            if(types[t]=="mixture"){
                n_train <- as.integer(round(n_samples*0.75))
                train_idx <- sample(n_samples, n_train)
                test_idx <- setdiff(1:n_samples, train_idx)
            }else{
                train_idx <- which(tissue_vec!=4)
                test_idx <- which(tissue_vec==4)
            }

            #Calculation
            X_train <- bulk_all[train_idx,]
            y_train <- ages_all[train_idx]
            X_test <- bulk_all[test_idx,]
            y_test <- ages_all[test_idx]
            m <- cv.glmnet(x = X_train, y = y_train, family = "gaussian", alpha = 0.5)
            best.lambda <- m$lambda.min
            en.model <- glmnet(x = X_train, y = y_train, family = "gaussian", lambda = best.lambda, alpha = 0.5)
            beta <- as.numeric(en.model$beta)
            est_y_train <- as.vector(predict(en.model, newx = X_train, s = best.lambda))
            est_y_test <- as.vector(predict(en.model, newx = X_test, s = best.lambda))
            epigenetic_acc_test <-  as.vector(est_y_test) - y_test
            mycor_test <- cor(est_y_test, y_test)
            mycor_train <- cor(est_y_train, y_train)
            dam_cor_ch <- cor(bulk_dam_all[test_idx], y_test)
            dam_cor_est <- cor(bulk_dam_all[test_idx], est_y_test)
            err_test <-  median(abs(as.vector(est_y_test) - y_test))
            err_train <-  median(abs(as.vector(est_y_train) - y_train))
            dfs_test <- dam_speed[test_idx]
            mycor_acc <- cor(dfs_test, epigenetic_acc_test)
            med_acc <- median(epigenetic_acc_test)


            #single cell damage check
            sc_all_dam <- dataset$sc_all_dam
            sc_all_mat <- dataset$sc_all_mat
            sc_age <- as.vector(predict(en.model, newx = sc_all_mat, s = best.lambda))
            sc_res <- cbind(sc_age, sc_all_dam)
            sc_res_test <- sc_res[as.integer(sc_res[,"sample_row_idx"]) %in% test_idx, ]
            mycor_ds <- cor(sc_res_test[,"donor_age"], sc_res_test[,"sc_age"])
            mycor_dc <- cor(sc_res_test[,"donor_age"], sc_res_test[,"cell_damage"])
            mycor_sc <- cor(sc_res_test[,"sc_age"], sc_res_test[,"cell_damage"])

            #Record
            tmp <- c(mycor_train, mycor_test, err_train, err_test, med_acc, mycor_acc, mycor_ds, mycor_dc, mycor_sc, dam_cor_ch, dam_cor_est)
            out_res <- rbind(out_res, tmp)
        }

        #record
        rownames(out_res) <- sapply(out_file_list, function(chr){strsplit(chr,"/")[[1]][8]})
        colnames(out_res) <- c("mycor_train", "mycor_test", "err_train", "err_test", "medacc_test", "mycor_acc", "mycor_ds", "mycor_dc", "mycor_sc", "dam_cor_ch", "dam_cor_est")
        #out_res_list[[j]] <- out_res
        write.csv(out_res, file=paste0(out_dir, types[t], ".out_res.csv"))
    }
    return(0)
}

#Set n_rep
n_rep <- 100

#Mixture
myseed1 <- 2000
set.seed(myseed1)
out_dir <- paste0(out_path, "standard/")
if(file.exists(out_dir)){
    unlink(out_dir, recursive=TRUE)
    dir.create(out_dir)
}else{
    dir.create(out_dir)
}
seeds_can <- 1:n_rep
out_file_list <- list()
#options(warn = 2)
for(i in 1:length(seeds_can)){
    myseed <- myseed1 + seeds_can[i]*2
    tmp_NSS <- runif(1)
    tmp_rcom <- runif(1)
    com_subset_r <- runif(1)
    subset_prop1 <- lapply(1:4, function(x){z<-runif(2);return(c(com_subset_r, z*(1 - com_subset_r)/sum(z)))})
    out_file <- paste0(out_dir, "para", ".ds", i, ".rds")
    out_file_list[[i]] <- Methy_data_simu2(n_markers=100, myseed=myseed, NSS=tmp_NSS, r_com=tmp_rcom, out_file=out_file, subset_prop=subset_prop1)
    cat(i, "\n")
}
#options(warn = 0)
res <- Epiclo_est2(out_file_list, out_dir)

#Additiomnal hardest (a=10. default)
myseed1 <- 3000
set.seed(myseed1)
n_rep <- 100
out_dir <- paste0(out_path, "hardest_a10/")
if(file.exists(out_dir)){
    unlink(out_dir, recursive=TRUE)
    dir.create(out_dir)
}else{
    dir.create(out_dir)
}
seeds_can <- 1:n_rep
out_file_list <- list()
for(i in 1:length(seeds_can)){
    myseed <- myseed1 + seeds_can[i]*2
    com_subset_r <- 0
    subset_prop1 <- lapply(1:4, function(x){z<-runif(2);return(c(com_subset_r, z*(1 - com_subset_r)/sum(z)))})
    out_file <- paste0(out_dir, "para", ".ds", i, ".rds")
    out_file_list[[i]] <- Methy_data_simu2(n_markers=100, myseed=myseed, NSS=0, r_com=0, out_file=out_file, subset_prop=subset_prop1, alpha=10)
    cat(i, "\n")
}
res <- Epiclo_est2(out_file_list, out_dir)


#Additiomnal hardest (a=25)
myseed1 <- 4000
set.seed(myseed1)
n_rep <- 100
out_dir <- paste0(out_path, "hardest_a25/")
if(file.exists(out_dir)){
    unlink(out_dir, recursive=TRUE)
    dir.create(out_dir)
}else{
    dir.create(out_dir)
}
seeds_can <- 1:n_rep
out_file_list <- list()
for(i in 1:length(seeds_can)){
    myseed <- myseed1 + seeds_can[i]*2
    com_subset_r <- 0
    subset_prop1 <- lapply(1:4, function(x){z<-runif(2);return(c(com_subset_r, z*(1 - com_subset_r)/sum(z)))})
    out_file <- paste0(out_dir, "para", ".ds", i, ".rds")
    out_file_list[[i]] <- Methy_data_simu2(n_markers=100, myseed=myseed, NSS=0, r_com=0, out_file=out_file, subset_prop=subset_prop1, alpha=25)
    cat(i, "\n")
}
res <- Epiclo_est2(out_file_list, out_dir)

#Additiomnal 20 (a=1)
myseed1 <- 5000
set.seed(myseed1)
n_rep <- 100
out_dir <- paste0(out_path, "hardest_a1/")
if(file.exists(out_dir)){
    unlink(out_dir, recursive=TRUE)
    dir.create(out_dir)
}else{
    dir.create(out_dir)
}
seeds_can <- 1:n_rep
out_file_list <- list()
for(i in 1:length(seeds_can)){
    myseed <- myseed1 + seeds_can[i]*2
    com_subset_r <- 0
    subset_prop1 <- lapply(1:4, function(x){z<-runif(2);return(c(com_subset_r, z*(1 - com_subset_r)/sum(z)))})
    out_file <- paste0(out_dir, "para", ".ds", i, ".rds")
    out_file_list[[i]] <- Methy_data_simu2(n_markers=100, myseed=myseed, NSS=0, r_com=0, out_file=out_file, subset_prop=subset_prop1, alpha=1)
    cat(i, "\n")
}
res <- Epiclo_est2(out_file_list, out_dir)

#Separate comsub 1
myseed1 <- 6000
set.seed(myseed1)
n_rep <- 100
out_dir <- paste0(out_path, "separate_comfix/")
if(file.exists(out_dir)){
    unlink(out_dir, recursive=TRUE)
    dir.create(out_dir)
}else{
    dir.create(out_dir)
}
com_subset_r_can <- runif(n_rep)
seeds_can <- 1:n_rep
out_file_list <- list()
for(i in 1:length(seeds_can)){
    myseed <- myseed1 + seeds_can[i]*2
    com_subset_r <- com_subset_r_can[i]
    subset_prop1 <- lapply(1:4, function(x){z<-runif(2);return(c(com_subset_r, z*(1 - com_subset_r)/sum(z)))})
    out_file <- paste0(out_dir, "para", ".ds", i, ".rds")
    out_file_list[[i]] <- Methy_data_simu2(n_markers=100, myseed=myseed, NSS=0, r_com=0, out_file=out_file, subset_prop=subset_prop1)
    cat(i, "\n")
}
res <- Epiclo_est2(out_file_list, out_dir)
write.csv(com_subset_r_can, file=paste0(out_dir, "com_subset_r_can.csv"))

#Separate comsub 2
myseed1 <- 7000
set.seed(myseed1)
n_rep <- 100
out_dir <- paste0(out_path, "separate_comdev/")
if(file.exists(out_dir)){
    unlink(out_dir, recursive=TRUE)
    dir.create(out_dir)
}else{
    dir.create(out_dir)
}
dev_test_scores1 <- dev_test_scores2 <- rep(NA, n_rep)
seeds_can <- 1:n_rep
out_file_list <- list()
for(i in 1:length(seeds_can)){
    myseed <- myseed1 + seeds_can[i]*2
    com_subset_r_all <- runif(4)
    subset_prop1 <- lapply(1:4, function(x){z<-runif(2);return(c(com_subset_r_all[x], z*(1 - com_subset_r_all[x])/sum(z)))})
    out_file <- paste0(out_dir, "para", ".ds", i, ".rds")
    out_file_list[[i]] <- Methy_data_simu2(n_markers=100, myseed=myseed, NSS=0, r_com=0, out_file=out_file,  subset_prop=subset_prop1)
    dev_test_scores1[i] <- mean(sapply(1:3,function(v1){sqrt(sum((subset_prop1[[4]][1] - subset_prop1[[v1]][1])^2))}))
    dev_test_scores2[i] <- subset_prop1[[4]][1]
    cat(i, "\n")
}
mat <- cbind(dev_test_scores1, dev_test_scores2)
colnames(mat) <- c("dev_test_scores1", "dev_test_scores2")
res <- Epiclo_est2(out_file_list, out_dir)
write.csv(mat, file=paste0(out_dir, "mat.csv"))

#共通サブセットがゼロの時にも時計が確立するか？
myseed1 <- 8000
set.seed(myseed1)
out_dir <- paste0(out_path, "separate_comzero/")
if(file.exists(out_dir)){
    unlink(out_dir, recursive=TRUE)
    dir.create(out_dir)
}else{
    dir.create(out_dir)
}
r_com_can <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
nss_can <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
paras_pair <- expand.grid(r_com_can, nss_can)
paras_can <- rbind(paras_pair, paras_pair, paras_pair, paras_pair, paras_pair)
colnames(paras_can) <- c("r_com", "nss")
n_rep <- nrow(paras_can)
seeds_can <- 1:n_rep
out_file_list <- list()
for(i in 1:length(seeds_can)){
    myseed <- myseed1 + seeds_can[i]*2
    com_subset_r <- 0
    subset_prop1 <- lapply(1:4, function(x){z<-runif(2);return(c(com_subset_r, z*(1 - com_subset_r)/sum(z)))})
    out_file <- paste0(out_dir, "para", ".ds", i, ".rds")
    out_file_list[[i]] <- Methy_data_simu2(n_markers=100, myseed=myseed, NSS=paras_can[i,"nss"], r_com=paras_can[i,"r_com"], out_file=out_file, subset_prop=subset_prop1)
}
res <- Epiclo_est2(out_file_list, out_dir)
write.csv(paras_can, file=paste0(out_dir, "paras_can.csv"))


#NSSのみ変えた場合にどうなるか？
myseed1 <- 9000
set.seed(myseed1)
n_rep <- 100
out_dir <- paste0(out_path, "separate_nss/")
if(file.exists(out_dir)){
    unlink(out_dir, recursive=TRUE)
    dir.create(out_dir)
}else{
    dir.create(out_dir)
}
nss_can <- runif(n_rep)
seeds_can <- 1:n_rep
out_file_list <- list()
for(i in 1:length(seeds_can)){
    myseed <- myseed1 + seeds_can[i]*2
    com_subset_r <- 0.25
    subset_prop1 <- lapply(1:4, function(x){z<-runif(2);return(c(com_subset_r, z*(1 - com_subset_r)/sum(z)))})
    out_file <- paste0(out_dir, "para", ".ds", i, ".rds")
    out_file_list[[i]] <- Methy_data_simu2(n_markers=100, myseed=myseed, NSS=nss_can[i], r_com=0.5, out_file=out_file, subset_prop=subset_prop1)
    cat(i, "\n")
}
res <- Epiclo_est2(out_file_list, out_dir)
write.csv(nss_can, file=paste0(out_dir, "nss_can.csv"))


#R_comのみ変えた場合にどうなるか？
myseed1 <- 10000
set.seed(myseed1)
n_rep <- 100
out_dir <- paste0(out_path, "separate_rcom/")
if(file.exists(out_dir)){
    unlink(out_dir, recursive=TRUE)
    dir.create(out_dir)
}else{
    dir.create(out_dir)
}
r_com_can <- runif(n_rep)
seeds_can <- 1:n_rep
out_file_list <- list()
for(i in 1:length(seeds_can)){
    myseed <- myseed1 + seeds_can[i]*2
    com_subset_r <- 0.25
    subset_prop1 <- lapply(1:4, function(x){z<-runif(2);return(c(com_subset_r, z*(1 - com_subset_r)/sum(z)))})
    out_file <- paste0(out_dir, "para", ".ds", i, ".rds")
    out_file_list[[i]] <- Methy_data_simu2(n_markers=100, myseed=myseed, NSS=0.5, r_com=r_com_can[i], out_file=out_file, subset_prop=subset_prop1)
    cat(i, "\n")
}
res <- Epiclo_est2(out_file_list, out_dir)
write.csv(r_com_can, file=paste0(out_dir, "r_com_can.csv"))

#epigenetic時計が不成立の場合
myseed1 <- 11000
set.seed(myseed1)
n_rep <- 100
out_dir <- paste0(out_path, "failcase/")
if(file.exists(out_dir)){
    unlink(out_dir, recursive=TRUE)
    dir.create(out_dir)
}else{
    dir.create(out_dir)
}
seeds_can <- 1:n_rep
out_file_list <- list()
for(i in 1:length(seeds_can)){
    myseed <- myseed1 + seeds_can[i]*2
    com_subset_r <- 0
    subset_prop1 <- lapply(1:4, function(x){z<-runif(2);return(c(com_subset_r, z*(1 - com_subset_r)/sum(z)))})
    out_file <- paste0(out_dir, "para", ".ds", i, ".rds")
    out_file_list[[i]] <- Methy_data_simu2(n_markers=100, myseed=myseed, NSS=0.2, r_com=0.2, out_file=out_file, subset_prop=subset_prop1)
    cat(i, "\n")
}
res <- Epiclo_est2(out_file_list, out_dir)




#Resul visualization (2923.5.8)
#labelの文字サイズだけ変更
#source("/Users/dokada/Dropbox/analysis/2022.5/methy_viz0517_2.R") #Run.2023.5.25
library(igraph)
library(glmnet)
library(parallel)
library(gplots)
out_path <- "/Users/dokada/Desktop/work/methy_viz0517_2/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}

#Example of epigenetic clock if 1st example
set.seed(1)
dataset <- readRDS("/Users/dokada/Desktop/work/methy_ra_matome0517/standard/para.ds1.rds")
bulk_all <- dataset$bulk_all
ages_all <- dataset$ages_all
dam_speed <- dataset$dam_speed
tissue_vec <- dataset$tissue_vec
cmat_cs_list <- dataset$cmat_cs_list
prob_list <- dataset$prob_list
bulk_dam_all <- dataset$bulk_dam_all

#Split data
n_samples <- nrow(bulk_all)
n_train <- as.integer(round(n_samples*0.75))
train_idx <- sample(n_samples, n_train)
test_idx <- setdiff(1:n_samples, train_idx)

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
err_test <-  median(abs(as.vector(est_y_test) - y_test))
err_train <-  median(abs(as.vector(est_y_train) - y_train))
dfs_test <- dam_speed[test_idx]
mycor_acc <- cor(dfs_test, epigenetic_acc_test)
med_acc <- median(epigenetic_acc_test)
dam_cor <- cor(bulk_dam_all[test_idx], y_test)

#standard files
out_dir <- paste0(out_path, "example_standard/")
dir.create(out_dir)

#Train
png(paste0(out_dir, "train_example.png"), width=960, height=960)
par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
xlab = "Chronological Age"
ylab = "Epigenetic Age"
main = paste0("Training, ", "cor=", round(mycor_train,2), ",", "err=", round(err_train,2))
plot(y_train, est_y_train, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, pch=16)
abline(0, 1, col="red", lwd=3)
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=4, adj=0)
dev.off()

#Test
png(paste0(out_dir, "test_example.png"), width=960, height=960)
par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
xlab = "Chronological Age"
ylab = "Epigenetic Age"
main = paste0("Test, ", "cor=", round(mycor_test,2), ",", "err=", round(err_test,2))
plot(y_test, est_y_test, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, pch=16)
abline(0, 1, col="red", lwd=3)
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=4, adj=0)
dev.off()

#Epigenetic acc vs dam speed 
png(paste0(out_dir, "acc_example.png"), width=960, height=960)
par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
xlab1 = "Epigenetic Acc"
ylab1 = "Damage speed"
cr <- round(cor(epigenetic_acc_test, dfs_test), 2)
main = paste0("Epigenetic Acc vs Damage speed, ", "cor = ", cr)
plot(epigenetic_acc_test, dfs_test, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, pch=16)
mtext(xlab1, side=1, line=6, cex=4)
mtext(ylab1, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=3, adj=0)
dev.off()


#Epigenetic acc vs bulk dam
png(paste0(out_dir, "bulkdam_example.png"), width=960, height=960)
par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
xlab = "Epigenetic Age"
ylab = "Accumlated Damage"
dam_test <-  bulk_dam_all[test_idx]
cr <- round(cor(est_y_test,dam_test), 2)
main = paste0("Epigenetic Age vs Accumlated Damage, ", "cor = ", cr)
plot(est_y_test, dam_test, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, pch=16)
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=3, adj=0)
dev.off()

#single cell damage check
sc_all_dam <- dataset$sc_all_dam
sc_all_mat <- dataset$sc_all_mat
sc_age <- as.vector(predict(en.model, newx = sc_all_mat, s = best.lambda))
sc_res <- cbind(sc_age, sc_all_dam)
sc_res_test <- sc_res[as.integer(sc_res[,"sample_row_idx"]) %in% test_idx, ]
mycor_ds <- cor(sc_res_test[,"donor_age"], sc_res_test[,"sc_age"])
mycor_dc <- cor(sc_res_test[,"donor_age"], sc_res_test[,"cell_damage"])
mycor_sc <- cor(sc_res_test[,"sc_age"], sc_res_test[,"cell_damage"])

#Single cells
png(paste0(out_dir, "sc_example1.png"), width=960, height=960)
par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
xlab = "Donor Chronological Age"
ylab = "Cellular damage"
xx <- sc_res_test[,"donor_age"]
yy <- sc_res_test[,"cell_damage"]
cr <- round(cor(xx, yy), 2)
main = paste0("Donor Age vs Cellular damage, ", "cor = ", cr)
plot(xx, yy, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, pch=16)
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=3, adj=0)
dev.off()

png(paste0(out_dir, "sc_example2.png"), width=960, height=960)
par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
xlab = "Donor Chronological Age"
ylab = "SC Epigenetic Age"
xx <- sc_res_test[,"donor_age"]
yy <- sc_res_test[,"sc_age"]
cr <- round(cor(xx, yy), 2)
main = paste0("Donor Age vs SC Epigenetic Age, ", "cor = ", cr)
plot(xx, yy, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, pch=16)
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=3, adj=0)
dev.off()

png(paste0(out_dir, "sc_example3.png"), width=960, height=960)
par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
xlab = "SC Epigenetic Age"
ylab = "Cellular damage"
xx <- sc_res_test[,"sc_age"]
yy <- sc_res_test[,"cell_damage"]
cr <- round(cor(xx, yy), 2)
main = paste0("SC Epigenetic Age vs Cellular damage, ", "cor = ", cr)
plot(xx, yy, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, pch=16)
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=3, adj=0)
dev.off()

#Guisograms of degree
cnt <- 1
for(i in 1:4){
    for(j in 1:3){
        if(i==1 && j==1){
            commat <- cmat_cs_list[[1]][[1]]
            Connecitivy <- rowSums(abs(commat))
            png(paste0(out_dir, "Common.png"), width=960, height=960)
            par(mar = c(9, 9, 6, 4)) ##bottom, left, top, right
            xlab = "Connecivity"
            ylab = "Frequency"
            main = "Common"
            hist(Connecitivy, cex.axis=4, cex.lab=4, cex.main=4,xlab="", ylab="", main="")
            mtext(xlab, side=1, line=6, cex=4)
            mtext(ylab, side=2, line=6, cex=4)
            mtext(main, side=3, line=3, cex=4)
            dev.off()
        }else if(j != 1){
            spemat <- cmat_cs_list[[i]][[j]]
            Connecitivy <- rowSums(abs(spemat))
            png(paste0(out_dir, "Sp", cnt, ".png"), width=960, height=960)
            par(mar = c(9, 9, 6, 4)) ##bottom, left, top, right
            xlab = "Connecivity"
            ylab = "Frequency"
            main = paste0("Specific", cnt)
            hist(Connecitivy, cex.axis=4, cex.lab=4, cex.main=4,xlab="", ylab="", main="")
            mtext(xlab, side=1, line=6, cex=4)
            mtext(ylab, side=2, line=6, cex=4)
            mtext(main, side=3, line=3, cex=4)
            dev.off()
            cnt <- cnt + 1
        }
    }
}
system(paste0("cd ", out_dir, ";","montage Common.png Sp1.png Sp2.png Sp3.png Sp4.png Sp5.png Sp6.png Sp7.png Sp8.png -tile 3x3 -geometry +2+2 ex_degree.png"))


#micture_standard_all (All 100 replicates)
dirs <- c("standard", "hardest_a1", "hardest_a10", "hardest_a25")
for(i in 1:length(dirs)){

    cans <- c("mixture", "separate")
    dir.create(paste0(out_path, dirs[i]))
    for(k in 1:2){
        dir <- dirs[i]
        out_dir <- paste0(out_path, dir, "/", cans[k], "/")
        dir.create(out_dir, recursive=T)
        path <- paste0("/Users/dokada/Desktop/work/methy_ra_matome0517/", dir, "/", cans[k], ".out_res.csv")
        res <- read.csv(path, header=T, row.names=1)

        main <- "Cor"
        png(paste0(out_dir, main, ".png"), width=960, height=960)
        par(mar = c(9, 12, 9, 4)) ##bottom, left, top, right
        lab_sub <- c("Train", "Test")
        tab <- res[,c("mycor_train", "mycor_test")]
        colnames(tab) <- lab_sub
        boxplot(tab, cex.axis=4, cex.lab=4, main=main, cex.main=4)
        mtext(main, side=2, line=6, cex=4)
        dev.off()

        main <- "Error"
        png(paste0(out_dir, main, ".png"), width=960, height=960)
        par(mar = c(9, 12, 9, 4)) ##bottom, left, top, right
        lab_sub <- c("Train", "Test")
        tab <- res[,c("err_train", "err_test")]
        colnames(tab) <- lab_sub
        boxplot(tab, cex.axis=4, cex.lab=4, main=main, cex.main=4)
        mtext(main, side=2, line=6, cex=4)
        dev.off()

        features <- c("mycor_acc",  "mycor_ds",  "mycor_dc",  "mycor_sc", "dam_cor_est")
        f_names <-  c("Epigenetic Acc vs Damage speed",  "Donor Age vs SC Epigenetic Age",  "Donor Age vs Cellular damage", "SC Epigenetic Age vs Cellular damage", "Epigenetic Age vs Accumlated Damage")
        y_labs <-  c("Cor",  "Cor",  "Cor",  "Cor", "Cor")
        for(j in 1:length(features)){
            png(paste0(out_dir, features[j], ".png"), width=960, height=960)
            par(mar = c(9, 12, 9, 4)) ##bottom, left, top, right
            tab <- res[,features[j]]
            boxplot(tab, cex.axis=4, cex.lab=4, main=f_names[j], cex.main=3, xlab="")
            mtext(y_labs[j], side=2, line=6, cex=4)
            dev.off()
        }
    }

    #Example histograms
    dataset <- readRDS(paste0("/Users/dokada/Desktop/work/methy_ra_matome0517/", dirs[i], "/para.ds1.rds"))
    bulk_all <- dataset$bulk_all
    ages_all <- dataset$ages_all
    cors <- rep(NA, ncol(bulk_all))
    for(s in 1:ncol(bulk_all)){
        cors[s] <- cor(bulk_all[,s], ages_all)
    }
    out_dir <- paste0(out_path, dir, "/")
    png(paste0(out_dir, dir, ".cor_hist.png"), width=960, height=960)
    par(mar = c(9, 9, 6, 4)) ##bottom, left, top, right
    xlab = "Correlation coef"
    ylab = "Frequency"
    main = "Histogram of Cor"
    hist(cors, cex.axis=4, cex.lab=4, cex.main=4,xlab="", ylab="", main="")
    mtext(xlab, side=1, line=6, cex=4)
    mtext(ylab, side=2, line=6, cex=4)
    mtext(main, side=3, line=3, cex=4)
    dev.off()
}

#共通サブセットが重要
#common subset1
#range(res[,"dam_cor_est"]- res[,"mycor_test"]) #similar but dufferent
res <- read.csv("/Users/dokada/Desktop/work/methy_ra_matome0517/separate_comfix/separate.out_res.csv", header=T, row.names=1)
para <- read.csv("/Users/dokada/Desktop/work/methy_ra_matome0517/separate_comfix/com_subset_r_can.csv", header=T, row.names=1)
term <- c("mycor_test", "err_test", "mycor_acc","dam_cor_est")
mains <- c("Performance(Cor)", "Performance(Error)", "Acc", "Damage")
ylabs <- c("Cor", "Error", "Cor", "Cor")
out_dir <- paste0(out_path, "comsub1/")
dir.create(out_dir)
for(i in 1:length(mains)){
    xx <- para[,1]
    yy <- res[,term[i]]
    png(paste0(out_dir, mains[i], "comsubset.png"), width=960, height=960)
    par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
    xlab = "Common_Subset_Proportion"
    ylab = mains[i]
    cr <- round(cor(xx, yy, method="spearman"), 2)
    main = paste0(mains[i], ", ", "cor = ", cr)
    plot(xx, yy, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, pch=16)
    mtext(xlab, side=1, line=6, cex=4)
    mtext(ylab, side=2, line=6, cex=4)
    mtext(main, side=3, line=3, cex=4, adj=0)
    dev.off()
}




#Common subset2
res <- read.csv("/Users/dokada/Desktop/work/methy_ra_matome0517/separate_comdev/separate.out_res.csv", header=T, row.names=1)
mat <- read.csv("/Users/dokada/Desktop/work/methy_ra_matome0517/separate_comdev/mat.csv", header=T, row.names=1)
term2 <- c("dev_test_scores1", "dev_test_scores2")
term <- c("mycor_test", "err_test", "mycor_acc","dam_cor_est")
mains <- c("Performance(Cor)", "Performance(Error)", "Acc", "Damage")
mains2 <- c("Deviation of Common_Subset_Proportion","Common_Subset_Proportion in Test tissue")
ylabs <- c("Cor", "Error", "Cor", "Cor")
out_dir <- paste0(out_path, "comsub2/")
dir.create(out_dir)
for(j in 1:2){
    for(i in 1:length(mains)){
        xx <- mat[,term2[j]]
        yy <- res[,term[i]]
        png(paste0(out_dir, mains2[j], ".", mains[i], "comsubset.png"), width=960, height=960)
        par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
        xlab = mains2[j]
        ylab = mains[i]
        cr <- round(cor(xx, yy, method="spearman"), 2)
        main = paste0(mains[i], ", ", "cor = ", cr)
        plot(xx, yy, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, pch=16)
        mtext(xlab, side=1, line=6, cex=3)
        mtext(ylab, side=2, line=6, cex=4)
        mtext(main, side=3, line=3, cex=4, adj=0)
        dev.off()
    }
}

#共通サブセットがゼロの時にも時計が確立するか？
res <- read.csv("/Users/dokada/Desktop/work/methy_ra_matome0517/separate_comzero/separate.out_res.csv", header=T, row.names=1)
para <- read.csv("/Users/dokada/Desktop/work/methy_ra_matome0517/separate_comzero/paras_can.csv", header=T, row.names=1)
term <- c("err_test", "mycor_test", "dam_cor_est")
mains <- c("Error", "Cor", "Cor")
r_com_can <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
nss_can <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
out_dir <- paste0(out_path, "comzero", "/")
dir.create(out_dir)
for(k in 1:length(term)){
    per_mat <- matrix(NA, length(r_com_can), length(nss_can))
    rownames(per_mat) <- r_com_can
    colnames(per_mat) <- nss_can
    for(i in 1:length(r_com_can)){
        for(j in 1:length(nss_can)){
            idxs <- which((para[,"r_com"]==r_com_can[i]) & (para[,"nss"]==nss_can[j]))
            per_mat[i,j] <- median(res[idxs,term[k]])
        }
    }
    png(paste0(out_dir, term[k], ".png"), width=480, height=480)
    heatmap.2(per_mat, dendrogram="none", Rowv=FALSE, Colv=FALSE, trace="none",density.info="none", xlab="NSS",ylab="R_COM", main=mains[k]) #xlabがcolumn
    dev.off()
}

#NSSのみ変えた場合にどうなるか？
res <- read.csv("/Users/dokada/Desktop/work/methy_ra_matome0517/separate_nss/separate.out_res.csv", header=T, row.names=1)
para <- read.csv("/Users/dokada/Desktop/work/methy_ra_matome0517/separate_nss/nss_can.csv", header=T, row.names=1)
term <- c("mycor_test", "err_test", "mycor_acc","dam_cor_est")
mains <- c("Performance(Cor)", "Performance(Error)", "Acc", "Damage")
ylabs <- c("Cor", "Error", "Cor", "Cor")
out_dir <- paste0(out_path, "nss/")
dir.create(out_dir)
for(i in 1:length(mains)){
    xx <- para[,1]
    yy <- res[,term[i]]
    png(paste0(out_dir, mains[i], ".NSS.png"), width=960, height=960)
    par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
    xlab = "Network_Similarity_Score"
    ylab = mains[i]
    cr <- round(cor(xx, yy, method="spearman"), 2)
    main = paste0(mains[i], ", ", "cor = ", cr)
    plot(xx, yy, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, pch=16)
    mtext(xlab, side=1, line=6, cex=4)
    mtext(ylab, side=2, line=6, cex=4)
    mtext(main, side=3, line=3, cex=4, adj=0)
    dev.off()
}

#R_comのみ変えた場合にどうなるか？
res <- read.csv("/Users/dokada/Desktop/work/methy_ra_matome0517/separate_rcom/separate.out_res.csv", header=T, row.names=1)
para <- read.csv("/Users/dokada/Desktop/work/methy_ra_matome0517/separate_rcom/r_com_can.csv", header=T, row.names=1)
term <- c("mycor_test", "err_test", "mycor_acc","dam_cor_est")
mains <- c("Performance(Cor)", "Performance(Error)", "Acc", "Damage")
ylabs <- c("Cor", "Error", "Cor", "Cor")
out_dir <- paste0(out_path, "rcom/")
dir.create(out_dir)
for(i in 1:length(mains)){
    xx <- para[,1]
    yy <- res[,term[i]]
    png(paste0(out_dir, mains[i], "r_com.png"), width=960, height=960)
    par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
    xlab = "Common_Marker_Proportion"
    ylab = mains[i]
    cr <- round(cor(xx, yy, method="spearman"), 2)
    main = paste0(mains[i], ", ", "cor = ", cr)
    plot(xx, yy, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, pch=16)
    mtext(xlab, side=1, line=6, cex=4)
    mtext(ylab, side=2, line=6, cex=4)
    mtext(main, side=3, line=3, cex=4, adj=0)
    dev.off()
}

#epigenetic時計が不成立の場合
out_dir <- paste0(out_path, "failbias/")
dir.create(out_dir)
res <- read.csv("/Users/dokada/Desktop/work/methy_ra_matome0517/failcase/separate.out_res.csv", header=T, row.names=1)
xx <- res[,"err_test"]
yy <- res[,"medacc_test"]
png(paste0(out_dir, "fail_bias.png"), width=960, height=960)
par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
xlab = "Error"
ylab = "Med Acc"
main=""
plot(xx, yy, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, pch=16)
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=2, adj=0)
dev.off()

#Histograms of cor
png(paste0(out_dir, ".perform.cor_hist.png"), width=960, height=960)
par(mar = c(9, 9, 6, 4)) ##bottom, left, top, right
xlab = "Cor"
ylab = "Frequency"
main = ""
hist(res[,"mycor_test"], cex.axis=4, cex.lab=4, cex.main=4,xlab="", ylab="", main="")
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=4)
dev.off()

#下振れと上振れ
term <- c("Up", "Down")
for(i in 1:2){
    if(term[i]=="Up"){
        fn <- rownames(res)[which.max(res[,"medacc_test"])]
    }else{
        fn <- rownames(res)[which.min(res[,"medacc_test"])]
    }
    
    dataset <- readRDS(paste0("/Users/dokada/Desktop/work/methy_ra_matome0517/failcase/", fn))
    bulk_all <- dataset$bulk_all
    ages_all <- dataset$ages_all
    dam_speed <- dataset$dam_speed
    tissue_vec <- dataset$tissue_vec
    cmat_cs_list <- dataset$cmat_cs_list
    prob_list <- dataset$prob_list
    bulk_dam_all <- dataset$bulk_dam_all


    #Split data
    train_idx <- which(tissue_vec!=4)
    test_idx <- which(tissue_vec==4)

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
    err_test <-  median(abs(as.vector(est_y_test) - y_test))
    err_train <-  median(abs(as.vector(est_y_train) - y_train))
    dfs_test <- dam_speed[test_idx]
    mycor_acc <- cor(dfs_test, epigenetic_acc_test)
    med_acc <- median(epigenetic_acc_test)
    dam_cor <- cor(bulk_dam_all[test_idx], y_test)

    #Train
    png(paste0(out_dir, term[i], ".ex.train_example.png"), width=960, height=960)
    par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
    xlab = "Chronological Age"
    ylab = "Epigenetic Age"
    main = paste0("Training, ", "cor=", round(mycor_train,2), ",", "err=", round(err_train,2))
    plot(y_train, est_y_train, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, pch=16)
    abline(0, 1, col="red", lwd=3)
    mtext(xlab, side=1, line=6, cex=4)
    mtext(ylab, side=2, line=6, cex=4)
    mtext(main, side=3, line=3, cex=4, adj=0)
    dev.off()

    #Test
    png(paste0(out_dir,  term[i], ".ex.test_example.png"), width=960, height=960)
    par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
    xlab = "Chronological Age"
    ylab = "Epigenetic Age"
    main = paste0("Test, ", "cor=", round(mycor_test,2), ",", "err=", round(err_test,2))
    plot(y_test, est_y_test, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, pch=16)
    abline(0, 1, col="red", lwd=3)
    mtext(xlab, side=1, line=6, cex=4)
    mtext(ylab, side=2, line=6, cex=4)
    mtext(main, side=3, line=3, cex=4, adj=0)
    dev.off()
}


#Subsetdrift
res <- read.csv("/Users/dokada/Desktop/work/methy_ra_drift0518/drift_com2/separate.out_res.csv", header=T, row.names=1)
para <- read.csv("/Users/dokada/Desktop/work/methy_ra_drift0518/drift_com2/w_can.csv", header=T, row.names=1)
term <- c("mycor_test", "err_test", "mycor_acc","dam_cor_est", "mycor_sc")
mains <- c("Performance(Cor)", "Performance(Error)", "Acc", "Damage", "SC_Age")
ylabs <- c("Cor", "Error", "Cor", "Cor", "Cor")
out_dir <- paste0(out_path, "drift_com/")
dir.create(out_dir)
pvals <- rep(NA, length(mains))
names(pvals) <- mains
for(i in 1:length(mains)){
    xx <- para[,1]
    yy1 <- res[1:50,term[i]]
    yy2 <- res[51:100,term[i]]
    zz <- list(yy1,yy2)
    png(paste0(out_dir, mains[i], ".drift_com.png"), width=960, height=960)
    par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
    names(zz) <- c("Normal", "Drift")
    mm <- paste0(mains[i])
    boxplot(zz, cex.axis=4, cex.lab=4, main=mm, cex.main=4)
    dev.off()
    pvals[i] <- t.test(yy1, yy2)$p.value
}
write.csv(pvals, file=paste0(out_dir, "pvaldrift.csv"))
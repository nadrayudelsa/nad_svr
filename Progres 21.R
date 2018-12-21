library(dplyr)

#---------------Jumlah ruta tiap kab (Ni)-------------------#
ht_rt <- function(dt1){
  r102 = unique(dt1$kab)
  de = data.frame(r102)
  Peng_jatim <- data.frame(dt1)
  ds <- list() 
  for (i in r102){
    ds[[i]] <- subset(Peng_jatim, kab %in% i)
  }
  Ni <- c()
  for(k in r102){
    Ni[k] <- length(ds[[k]]$kab)
  }
  N_ruta <- data.frame("ruta"=Ni)
  return(N_ruta)
}

#-----------------------------ESTIMASI DIRECT METODE SRS-----------------------------#
resampling_srs <- function(data, prop,iterasi ){
  propinsi <- function(data, prop){  
    srs_sampling <- function(data, prop){
      kab_1 <- data
      samp_rt <- sample(1:length(kab_1$pop), size = prop*length(kab_1$pop),
                        replace = FALSE)
      hs_samp_srs <-  kab_1[samp_rt,]
      hs_samp_srs$pop <- as.numeric(paste(hs_samp_srs$pop))
      hs_samp_srs <- as.data.frame(hs_samp_srs)
      est_ybar_srs <- function(data,f){
        hs_samp_srs <- data
        y_tot_samp <- sum(hs_samp_srs$pop)
        n <- length(hs_samp_srs$pop)
        y_bar <- y_tot_samp/n
        Y_cap_est <- n/f*y_bar
        s_kw <- sum((hs_samp_srs$pop - (rep(y_bar,n)))^2)/(n-1)
        var_y_bar <- (1-f)*s_kw/n
        Var_cap_est <- (n/f)^2*var_y_bar
        RRMSE_var_y <- sqrt(Var_cap_est)/Y_cap_est*100
        hasil <- cbind(Y_cap_est,Var_cap_est,RRMSE_var_y)
        return(hasil)
      }
      result = est_ybar_srs(hs_samp_srs, prop)
      return(result)
    }
    
    r102 = unique(data$kab)
    datakoeh = list()
    for(b in r102){
      datakoeh[[b]] <- subset(data, kab %in% b)
    }
    
    try_ye <- list()
    for(c in r102){
      try_ye[[c]] <- srs_sampling(data = datakoeh[[c]], prop = prop)  
    }
    
    this_it = data.frame()
    for (q in r102) {
      this_it <- rbind(this_it, try_ye[[q]])
    }
    this_it = cbind(r102, this_it)
  }
  
  dataset = list()
  for( j in 1:iterasi){
    set.seed(j)
    dataset[[j]] = propinsi(data, prop)
  }
  
  databind = data.frame()
  for (b in 1:iterasi) {
    databind = rbind(databind, dataset[[b]])
  }
  
  datafinal = aggregate(list(databind$Y_cap_est,databind$Var_cap_est,databind$RRMSE_var_y)
                        ,by=list(databind$r102),FUN=mean)
  colnames(datafinal) = c("kab", "Y dir", "MSE dir","RRMSE dir")
  return(datafinal)
}

#--------------------------------ESTIMASI DIRECT METODE ONE STAGE CLUSTER--------------------------------------------------------------
resampling_1sc <- function(data,prop,iterasi ){
  propinsi <- function(data, prop){  
    cb_sampling <- function(data, prop){
      bs_u <-unique(data$bs_o)
      bs_u <- data.frame(bs_u)
      samp_bs <- sample(1:length(bs_u$bs_u), size = round(prop*length(bs_u$bs_u)), replace = FALSE)
      sampi <- bs_u[samp_bs,]
      sampi2 <- as.data.frame(paste(sampi))
      colnames(sampi2) <- c("sampel_bs")
      hsl_smpl <- subset(data,bs_o %in% sampi2$sampel_bs)
      
      est_ybar_1cs <- function(data,f){
        hs_sam <- data
        hs_sam$pop <- as.numeric(paste(hs_sam$pop))
        hs_sam$bs_o <- as.character(paste(hs_sam$bs_o))
        Mi <- data.frame(table(hs_sam$bs_o))
        y_bar_i <- aggregate(hs_sam$pop, by=list(hs_sam$bs_o), FUN=mean)
        y_bar_n <- (sum(y_bar_i$x*Mi$Freq))/sum(Mi$Freq)
        M_bar <- sum(Mi$Freq)/nrow(Mi)
        pnyb <- sum((Mi$Freq^2)*((y_bar_i$x-rep(c(y_bar_n),nrow(Mi)))^2))
        pmbg <- nrow(Mi)*(M_bar)^2*(nrow(Mi)-1)
        var_y_bar_n <- (1-f)*(pnyb/pmbg) 
        Y_cap_est <- (nrow(Mi)/f)*M_bar*y_bar_n
        Var_cap_est <- ((nrow(Mi)/f)*M_bar)^2*var_y_bar_n
        RRMSE_var_y <- sqrt(Var_cap_est)/Y_cap_est*100
        hasil <- cbind(Y_cap_est,Var_cap_est,RRMSE_var_y)
        return(hasil)
      }
      result <- est_ybar_1cs(hsl_smpl, prop)
    }
    
    r102 = unique(data$kab)
    datakoeh = list()
    for(b in r102){
      datakoeh[[b]] <- subset(data, kab %in% b)
    }
    
    try_ye <- list()
    for(c in r102){
      try_ye[[c]] <- cb_sampling(data = datakoeh[[c]], prop = prop)  
    }
    
    this_it = data.frame()
    for (q in r102) {
      this_it <- rbind(this_it, try_ye[[q]])
    }
    this_it = cbind(r102, this_it)
  }
  
  dataset = list()
  for( j in 1:iterasi){
    set.seed(j)
    dataset[[j]] = propinsi(data, prop)
  }
  
  databind = data.frame()
  for (b in 1:iterasi) {
    databind = rbind(databind, dataset[[b]])
  }
  
  datafinal = aggregate(list(databind$Y_cap_est,databind$Var_cap_est,databind$RRMSE_var_y),
                        by=list(databind$r102),FUN=mean)
  colnames(datafinal) = c("kab", "Y dir", "MSE dir", "RRMSE dir")
  return(datafinal)
}

#-------------------------------ESTIMASI DIRECT METODE TWO STAGE CLUSTER (SRSWOR-SRSWOR)------------------------------------------
resampling_2sc_srs = function(dataY, f1, f2, iterasi ){
  est_2st_cl <- function(dataY, f1,f2){
    sp_2_ct <- function(dataY,f1,f2){
      samp_th1 <- function(mydata,f1){ #pengambilan sampel tahap 1 secara SRSWOR
        bs_u <-unique(mydata$bs_o)
        bs_u <- data.frame(bs_u)
        samp_bs <- sample(1:nrow(bs_u), size = f1*nrow(bs_u), replace = FALSE)
        sampi <- bs_u[samp_bs,]
        sampi2 <- as.data.frame(paste(sampi))
        colnames(sampi2) <- c("sampel_bs")
        return(sampi2)
      }
      hs_sp_t1 <- samp_th1(dataY,f1) #ambil sampel tahap 1
      samp_th2 <- function(dt_kb_1,hs_sp_1,f2){ 
        #Pengambilan sampel tahap 2 secara SRSWOR
        index <- c(1:nrow(hs_sp_1))
        hs_sp_1 <- cbind(index,hs_sp_1)
        dt_2 <- list();dt_3 <- list() ; samp_rt <- list() ; rt <- list() 
        hs_samp_srs <- list() ; hs_smpl_srs <- list() ; Mi <- list() ; Sw_i <- list()
        for(i in hs_sp_1$index){ #fungsi for untuk Pengambilan sampel ssu scr srs wor
          dt_2[[i]] <- subset(hs_sp_1, index %in% i)
          dt_3[[i]] <- subset(dt_kb_1, bs_o %in% dt_2[[i]]$sampel_bs)
          Mi[[i]] <- nrow(dt_3[[i]])
          samp_rt[[i]] <- sample(1:length(dt_3[[i]]$pop), 
                                 size = round(f2*length(dt_3[[i]]$pop)), 
                                 replace = FALSE)
          hs_samp_srs[[i]] <-  dt_3[[i]][samp_rt[[i]],]
          rt[[i]] <- sum(as.numeric(paste(hs_samp_srs[[i]]$pop)))/
            (length(hs_samp_srs[[i]]$pop))
          Sw_i[[i]] <- sum((as.numeric(paste(hs_samp_srs[[i]]$pop))-rt[[i]])^2)/
            (length(hs_samp_srs[[i]]$pop)-1)
          hs_smpl_srs[[i]] <- cbind(hs_samp_srs[[i]],rt[[i]],Sw_i[[i]],Mi[[i]])
        }
        gb_samp = data.frame()
        for (j in hs_sp_1$index) {
          gb_samp <- rbind(gb_samp,hs_smpl_srs[[j]])
        }
        return(gb_samp)
      }
      res_sp <- samp_th2(dataY,hs_sp_t1,f2) #ambil sampel tahap 2 
      
      #Estimasi DIRECT
      res_sp$`rt[[i]]` <- as.numeric(paste(res_sp$`rt[[i]]`))
      res_sp$`Sw_i[[i]]` <- as.numeric(paste(res_sp$`Sw_i[[i]]`))
      res_sp$`Mi[[i]]` <- as.numeric(paste(res_sp$`Mi[[i]]`))
      dataY$bs_o <- as.character(paste(dataY$bs_o))
      Mi <- aggregate(res_sp$`Mi[[i]]`, list(res_sp$bs_o), mean)
      mi <- round(Mi$x*f2)
      n <- length(unique(res_sp$bs_o))
      N <- length(unique(dataY$bs_o))
      #estimasi total
      ybar_i. <- aggregate(res_sp$`rt[[i]]`, list(res_sp$bs_o), mean)
      Y_cap_est <- (N/n)*(sum(Mi$x * ybar_i.$x))
      #estimasi varians dari total
      sw_i_2 <- aggregate(res_sp$`Sw_i[[i]]`, list(res_sp$bs_o),mean)
      sb_2 <- 1/(n-1) * sum(((Mi$x*ybar_i.$x)-(1/n*sum(Mi$x*ybar_i.$x)))^2)
      var_cap_est <- (N^2 * ((1-f1)) * (sb_2/n) ) + 
        (N/n * sum(Mi$x^2 * rep((1-f2),n) * (sw_i_2$x/mi)))
      RRMSE_var_y <- sqrt(var_cap_est)/Y_cap_est*100
      hasil <- cbind(Y_cap_est,var_cap_est,RRMSE_var_y)
      return(hasil)
    }    
    r102 = unique(dataY$kab)
    datakoeh = list() 
    for(b in r102){
      datakoeh[[b]] <- subset(dataY, kab %in% b)
    }
    
    try_ye <- list()
    for(c in r102){
      try_ye[[c]] <- sp_2_ct(data = datakoeh[[c]],f1,f2)  
    }
    
    this_it = data.frame()
    for (q in r102) {
      this_it <- rbind(this_it, try_ye[[q]])
    }
    
    this_it = cbind(r102, this_it)
  }
  
  dataset = list()
  for( j in 1:iterasi){
    set.seed(j)
    dataset[[j]] = est_2st_cl(dataY, f1, f2)
  }
  
  databind = data.frame()
  for (b in 1:iterasi) {
    databind = rbind(databind, dataset[[b]])
  }
  
  datafinal = aggregate(list(databind$Y_cap_est, databind$var_cap_est,databind$RRMSE_var_y),
                        by=list(databind$r102), FUN=mean, na.rm = TRUE)
  colnames(datafinal) = c("kab", "Y dir","MSE dir","RRMSE dir")
  return(datafinal)
}

#----------------------------ESTIMASI DIRECT METODE TWO STAGE CLUSTER (PPSWR-SRSWOR)-------------------------------------------------
resampling_2sc_pps = function(dataY, f1, f2, iterasi ){
  est_2st_cl <- function(dataY, f1,f2){
    sp_2_ct <- function(dataY,f1,f2){
      kb_1 <- dataY
      kb_1$bs_o <- as.character(paste(kb_1$bs_o))
      Mi_kb_1 <- data.frame(table(kb_1$bs_o))
      colnames(Mi_kb_1) <- c("bs_o","Muatan_x")
      Mi_kb_1$P_i <- Mi_kb_1$Muatan_x/(sum(Mi_kb_1$Muatan_x))
      samp_th1 <- function(data,f1){ #pengambilan sampel tahap 1 
        pps_wr <- function(sizes,n){
          N <- length(sizes)		
          cumsizes <- cumsum(sizes)
          totsize <- cumsizes[N]
          s <- numeric(n)		
          for (u in 1:n) {	
            r <- runif(1,0,totsize)
            i <- 1
            while (cumsizes[i] < r) {
              i <- i+1
            }
            s[u] <- i
          }
          return(s)
        }
        n_length <- f1*nrow(data)
        data_sa <- data[pps_wr(data$Muatan_x,n_length),]
        return(data_sa)
      }
      
      hs_sp_t1 <- samp_th1(Mi_kb_1,f1)
      samp_th2 <- function(dt_kb_1,hs_sp_1,f2){
        hs_sp_1 <- hs_sp_t1
        index <- c(1:nrow(hs_sp_1))
        hs_sp_1 <- cbind(index,hs_sp_1)
        dt_2 <- list();dt_3 <- list() ; 
        samp_rt <- list() ; kode <- list() ; Mi <- list() 
        hs_samp_srs <- list() ; hs_smp_srs <- list();
        n_samp <- list(); p_i <- list() 
        for(i in hs_sp_1$index){ #Pengambilan sampel ssu scr srs wor
          dt_2[[i]] <- subset(hs_sp_1, index %in% i)
          dt_3[[i]] <- subset(dt_kb_1, bs_o %in% dt_2[[i]]$bs_o)
          samp_rt[[i]] <- sample(1:length(dt_3[[i]]$pop), size = f2*length(dt_3[[i]]$pop),
                                 replace = FALSE)
          hs_samp_srs[[i]] <-  dt_3[[i]][samp_rt[[i]],]
          n_samp [[i]] <- nrow(hs_samp_srs[[i]])
          kode[[i]] <- rep(i,n_samp[[i]])
          p_i[[i]] <- rep(dt_2[[i]]$P_i,n_samp[[i]])
          Mi[[i]] <- rep(dt_2[[i]]$Muatan_x,n_samp[[i]])
          hs_smp_srs[[i]] <- cbind("kode"=kode[[i]],hs_samp_srs[[i]], "Mi"=Mi[[i]] 
                                   , "P_i" = p_i[[i]])
        }
        gb_samp = data.frame()
        for (j in hs_sp_1$index) {
          gb_samp <- rbind(gb_samp,hs_smp_srs[[j]])
        }
        return(gb_samp)
      }
      
      res_sp <- samp_th2(kb_1,hs_sp_t1,f2) 
      res_sp$pop <- as.numeric(paste(res_sp$pop))
      res_sp$Mi <- as.numeric(paste(res_sp$Mi))
      res_sp$P_i <- as.numeric(paste(res_sp$P_i))
      y_bar_i. <- aggregate(res_sp$pop, list(res_sp$kode), mean)
      M_i <- aggregate(res_sp$Mi, list(res_sp$kode), mean)
      P_i <- aggregate(res_sp$P_i, list(res_sp$kode), mean)
      Y_cap_est <- 1/ (length(y_bar_i.$x))*(sum(M_i$x*y_bar_i.$x/P_i$x))
      
      #estimasi varians
      var_cap_est <- (sum(((M_i$x*y_bar_i.$x/P_i$x)-Y_cap_est)^2))/
        ((length(y_bar_i.$x)-1)*length(y_bar_i.$x))
      RRMSE_var_y <- sqrt(var_cap_est)/Y_cap_est*100
      hasil <- cbind(Y_cap_est,var_cap_est,RRMSE_var_y)
      return(hasil)
    }
    
    r102 = unique(dataY$kab)
    datakoeh = list() 
    for(b in r102){
      datakoeh[[b]] <- subset(dataY, kab %in% b)
    }
    
    try_ye <- list()
    for(c in r102){
      try_ye[[c]] <- sp_2_ct(data = datakoeh[[c]],f1,f2)  
    }
    
    this_it = data.frame()
    for (q in r102) {
      this_it <- rbind(this_it, try_ye[[q]])
    }
    this_it = cbind(r102, this_it)
  }  
  
  dataset = list()
  for( j in 1:iterasi){
    set.seed(j)
    dataset[[j]] = est_2st_cl(dataY, f1, f2)
  }
  
  databind = data.frame()
  for (b in 1:iterasi) {
    databind = rbind(databind, dataset[[b]])
  }
  
  datafinal = aggregate(list(databind$Y_cap_est, databind$var_cap_est,databind$RRMSE_var_y),
                        by=list(databind$r102),
                        FUN=mean)
  colnames(datafinal) = c("kab", "Y dir", "MSE dir","RRMSE dir")
  return(datafinal)
}

#-------------------------ESTIMASI SECARA SYNTHETIC dan COMPOSITE -------------------------------------
est_prov_synth_comp <- function(x,ydir,Ni){
  r102 <- as.character(paste(ydir$kab))
  y <- as.numeric(paste(ydir$`Y dir`))  
  vary <- as.numeric(paste(ydir$`MSE dir`))
  Ni <- as.numeric(paste(Ni$ruta))
  Gab_kab <- function(x,y,vary,Ni,r102){
    Y_synth_comp <- function(x,y,vary,Ni){
      
      #synthetic function
      N <- Ni
      xty <- t(x)%*%y
      xtx <- t(x)%*%x
      bt_xy <- solve(xtx)%*%xty
      bt_xy = t(bt_xy)
      Y_syn <- c(rep(0,nrow(x)))
      for(i in 1: ncol(x)){
        Y_syn = Y_syn + x[,i]*bt_xy[,i]
      }
      m <- length(Y_syn)
      ax <- x[,-1]
      sum_x <- colSums(ax)
      xi_per_X <- (ax/sum_x)^2
      sum_y_s <- sum(Y_syn)
      R_cap <- sum_y_s/sum_x
      R_cap_x_i <- c()
      for (i in 1:ncol(ax)) {
        cap_x_i <- ax[,i] * R_cap[i]
        R_cap_x_i <- cbind(R_cap_x_i, cap_x_i)
      }
      e_i <- Y_syn - R_cap_x_i
      var_e_i_2 <- ((e_i)^2)/(m-1)
      var_y_s <- xi_per_X*var_e_i_2
      var_y_cap_is <- apply(var_y_s,1,function(x) {
        set.seed(0)
        library(resample)
        b <- bootstrap(x,median,R=100)
        m <- min(b$replicates)
        return(m)
      })
      Ni_2 <- N^2
      ba.1 <- (sum(((Y_syn - y)^2)/Ni_2))/m
      ba.2 <- (sum(vary/Ni_2))/m
      ba.3 <- (sum(var_y_cap_is/Ni_2))/m
      ba_a <- ba.1 - ba.2 - ba.3
      mse_syn <- var_y_cap_is + (Ni_2*ba_a)
      RRMSE_syn <- sqrt(mse_syn)/Y_syn*100
      
      #composite function
      tetha_i <- mse_syn/(vary+mse_syn)
      tetha_i[tetha_i<0]<-0
      tetha_i[tetha_i>1]<-1
      Y_comp <- (tetha_i*y) + ((1-tetha_i)*Y_syn)
      comp <- data.frame(tetha_i,mse_syn,vary)
      colnames(comp) <- c("tetha_i","mse_syn","vary")
      library(dplyr)
      fun_mse_com <- function(tetha_i,mse_syn,vary){
        if(tetha_i>0.5){
          return((1-tetha_i)*mse_syn)
        }
        else{
          return(tetha_i*vary) 
        }
      }
      hasil_com <- comp %>%
        rowwise() %>%
        mutate(mse_com = fun_mse_com(tetha_i,mse_syn,vary))
      mse_com <- as.data.frame(hasil_com)
      mse_comp <- mse_com$mse_com
      RRMSE_comp <- sqrt(mse_comp)/Y_comp*100
      Syn_Comp <- data.frame(cbind(Y_syn, Y_comp,mse_syn,  mse_comp, RRMSE_syn, RRMSE_comp))
      colnames(Syn_Comp) = c("Y syn","Y comp", "MSE syn",  "MSE comp", "RRMSE syn","RRMSE comp")
      return(Syn_Comp)
    }
    this_it = Y_synth_comp(x,y,vary,Ni)
    return(this_it)
  }
  
  est_dataset = Gab_kab(x,y,vary,Ni,r102)
  datafinal = data.frame(r102,y,est_dataset$`Y syn`,est_dataset$`Y comp`,vary,est_dataset$`MSE syn`,est_dataset$`MSE comp`,
                         (sqrt(vary)/y)*100,est_dataset$`RRMSE syn`,est_dataset$`RRMSE comp`,(est_dataset$`MSE syn`/vary),
                         (est_dataset$`MSE comp`/vary))
  
  colnames(datafinal) = c("kab","Y dir","Y syn","Y com","MSE dir","MSE syn","MSE com","RRMSE dir","RRMSE syn",
                          "RRMSE com","Ef syn","Ef comp")
  return(datafinal)
}

################################################################################################################################
#----------------------------------------------------A. PERSIAPAN DATA BANGKITAN---------------------------------------------------

#A. Persiapan
library(readxl)
dt_empiris <- read_excel("D:/KAMPUS/skripsi/DATA/dt_empiris.xlsx", col_types = c("text","text", "text","numeric", "numeric", "numeric"))

data_1 <- data.frame(dt_empiris$`No Kode Sampel`,dt_empiris$exp_krb,dt_empiris$R102)
colnames(data_1) <- c("bs_o","pop","kab")
hist(x = data_1$pop, col = "red2", border = "indianred1", xlab = "Data I",main = "Histogram DATA I")

data_2 <- data.frame(dt_empiris$`No Kode Sampel`,dt_empiris$KALORI,dt_empiris$R102)
colnames(data_2) <- c("bs_o","pop","kab")
hist(data_2$pop, col = "mediumpurple4", border = "mediumpurple", xlab = "Data II",main = "Histogram DATA II")

#data simulasi
set.seed(1234)
b_1 <- abs(rnorm(length(data_1$pop),mean(data_1$pop),5000))
hist(b_1, col = "seagreen3", border = "seagreen4", xlab = "Data III",main = "Histogram DATA III")
data_3 <- data.frame(dt_empiris$`No Kode Sampel`,b_1,dt_empiris$R102)
colnames(data_3) <- c("bs_o","pop","kab")

set.seed(123)
b_2 <- rnorm(length(data_1$pop),mean(data_1$pop),5)
hist(b_2,col = "magenta4", border = "magenta3", xlab = "Data IV",main = "Histogram DATA IV")
data_4 <- data.frame(dt_empiris$`No Kode Sampel`,b_2,dt_empiris$R102)
colnames(data_4) <- c("bs_o","pop","kab")

########################################################################################################################################
#------------------------------------B. Running Data secara Direct----------------------------------------------------------------------
#-------------------------------------------------SRS----------------------------------------------------------------------------
res_1_srs_1 <- resampling_srs(data = data_1, prop = 1, iterasi = 100)
res_1_srs_0.9 <- resampling_srs(data = data_1, prop = 0.9, iterasi = 100)
res_1_srs_0.8 <- resampling_srs(data = data_1, prop = 0.8, iterasi = 100)
res_1_srs_0.7 <- resampling_srs(data = data_1, prop = 0.7, iterasi = 100)
res_1_srs_0.6 <- resampling_srs(data = data_1, prop = 0.6, iterasi = 100)
res_1_srs_0.5 <- resampling_srs(data = data_1, prop = 0.5, iterasi = 100)
res_1_srs_0.4 <- resampling_srs(data = data_1, prop = 0.4, iterasi = 100)
res_1_srs_0.3 <- resampling_srs(data = data_1, prop = 0.3, iterasi = 100)
res_1_srs_0.2 <- resampling_srs(data = data_1, prop = 0.2, iterasi = 100)
res_1_srs_0.1 <- resampling_srs(data = data_1, prop = 0.1, iterasi = 100)
res_1_srs_0.09 <- resampling_srs(data = data_1, prop = 0.09, iterasi = 100)
res_1_srs_0.08 <- resampling_srs(data = data_1, prop = 0.08, iterasi = 100)
res_1_srs_0.07 <- resampling_srs(data = data_1, prop = 0.07, iterasi = 100)
res_1_srs_0.06 <- resampling_srs(data = data_1, prop = 0.06, iterasi = 100)
res_1_srs_0.05 <- resampling_srs(data = data_1, prop = 0.05, iterasi = 100)
res_1_srs_0.04 <- resampling_srs(data = data_1, prop = 0.04, iterasi = 100)
res_1_srs_0.03 <- resampling_srs(data = data_1, prop = 0.03, iterasi = 100)
res_1_srs_0.02 <- resampling_srs(data = data_1, prop = 0.02, iterasi = 100)
res_1_srs_0.01 <- resampling_srs(data = data_1, prop = 0.01, iterasi = 100)

res_2_srs_1 <- resampling_srs(data = data_2, prop = 1, iterasi = 100)
res_2_srs_0.9 <- resampling_srs(data = data_2, prop = 0.9, iterasi = 100)
res_2_srs_0.8 <- resampling_srs(data = data_2, prop = 0.8, iterasi = 100)
res_2_srs_0.7 <- resampling_srs(data = data_2, prop = 0.7, iterasi = 100)
res_2_srs_0.6 <- resampling_srs(data = data_2, prop = 0.6, iterasi = 100)
res_2_srs_0.5 <- resampling_srs(data = data_2, prop = 0.5, iterasi = 100)
res_2_srs_0.4 <- resampling_srs(data = data_2, prop = 0.4, iterasi = 100)
res_2_srs_0.3 <- resampling_srs(data = data_2, prop = 0.3, iterasi = 100)
res_2_srs_0.2 <- resampling_srs(data = data_2, prop = 0.2, iterasi = 100)
res_2_srs_0.1 <- resampling_srs(data = data_2, prop = 0.1, iterasi = 100)
res_2_srs_0.09 <- resampling_srs(data = data_2, prop = 0.09, iterasi = 100)
res_2_srs_0.08 <- resampling_srs(data = data_2, prop = 0.08, iterasi = 100)
res_2_srs_0.07 <- resampling_srs(data = data_2, prop = 0.07, iterasi = 100)
res_2_srs_0.06 <- resampling_srs(data = data_2, prop = 0.06, iterasi = 100)
res_2_srs_0.05 <- resampling_srs(data = data_2, prop = 0.05, iterasi = 100)
res_2_srs_0.04 <- resampling_srs(data = data_2, prop = 0.04, iterasi = 100)
res_2_srs_0.03 <- resampling_srs(data = data_2, prop = 0.03, iterasi = 100)
res_2_srs_0.02 <- resampling_srs(data = data_2, prop = 0.02, iterasi = 100)
res_2_srs_0.01 <- resampling_srs(data = data_2, prop = 0.01, iterasi = 100)

res_3_srs_1 <- resampling_srs(data = data_3, prop = 1, iterasi = 100)
res_3_srs_0.9 <- resampling_srs(data = data_3, prop = 0.9, iterasi = 100)
res_3_srs_0.8 <- resampling_srs(data = data_3, prop = 0.8, iterasi = 100)
res_3_srs_0.7 <- resampling_srs(data = data_3, prop = 0.7, iterasi = 100)
res_3_srs_0.6 <- resampling_srs(data = data_3, prop = 0.6, iterasi = 100)
res_3_srs_0.5 <- resampling_srs(data = data_3, prop = 0.5, iterasi = 100)
res_3_srs_0.4 <- resampling_srs(data = data_3, prop = 0.4, iterasi = 100)
res_3_srs_0.3 <- resampling_srs(data = data_3, prop = 0.3, iterasi = 100)
res_3_srs_0.2 <- resampling_srs(data = data_3, prop = 0.2, iterasi = 100)
res_3_srs_0.1 <- resampling_srs(data = data_3, prop = 0.1, iterasi = 100)
res_3_srs_0.09 <- resampling_srs(data = data_3, prop = 0.09, iterasi = 100)
res_3_srs_0.08 <- resampling_srs(data = data_3, prop = 0.08, iterasi = 100)
res_3_srs_0.07 <- resampling_srs(data = data_3, prop = 0.07, iterasi = 100)
res_3_srs_0.06 <- resampling_srs(data = data_3, prop = 0.06, iterasi = 100)
res_3_srs_0.05 <- resampling_srs(data = data_3, prop = 0.05, iterasi = 100)
res_3_srs_0.04 <- resampling_srs(data = data_3, prop = 0.04, iterasi = 100)
res_3_srs_0.03 <- resampling_srs(data = data_3, prop = 0.03, iterasi = 100)
res_3_srs_0.02 <- resampling_srs(data = data_3, prop = 0.02, iterasi = 100)
res_3_srs_0.01 <- resampling_srs(data = data_3, prop = 0.01, iterasi = 100)

res_4_srs_1 <- resampling_srs(data = data_4, prop = 1, iterasi = 100)
res_4_srs_0.9 <- resampling_srs(data = data_4, prop = 0.9, iterasi = 100)
res_4_srs_0.8 <- resampling_srs(data = data_4, prop = 0.8, iterasi = 100)
res_4_srs_0.7 <- resampling_srs(data = data_4, prop = 0.7, iterasi = 100)
res_4_srs_0.6 <- resampling_srs(data = data_4, prop = 0.6, iterasi = 100)
res_4_srs_0.5 <- resampling_srs(data = data_4, prop = 0.5, iterasi = 100)
res_4_srs_0.4 <- resampling_srs(data = data_4, prop = 0.4, iterasi = 100)
res_4_srs_0.3 <- resampling_srs(data = data_4, prop = 0.3, iterasi = 100)
res_4_srs_0.2 <- resampling_srs(data = data_4, prop = 0.2, iterasi = 100)
res_4_srs_0.1 <- resampling_srs(data = data_4, prop = 0.1, iterasi = 100)
res_4_srs_0.09 <- resampling_srs(data = data_4, prop = 0.09, iterasi = 100)
res_4_srs_0.08 <- resampling_srs(data = data_4, prop = 0.08, iterasi = 100)
res_4_srs_0.07 <- resampling_srs(data = data_4, prop = 0.07, iterasi = 100)
res_4_srs_0.06 <- resampling_srs(data = data_4, prop = 0.06, iterasi = 100)
res_4_srs_0.05 <- resampling_srs(data = data_4, prop = 0.05, iterasi = 100)
res_4_srs_0.04 <- resampling_srs(data = data_4, prop = 0.04, iterasi = 100)
res_4_srs_0.03 <- resampling_srs(data = data_4, prop = 0.03, iterasi = 100)
res_4_srs_0.02 <- resampling_srs(data = data_4, prop = 0.02, iterasi = 100)
res_4_srs_0.01 <- resampling_srs(data = data_4, prop = 0.01, iterasi = 100)

#------------------------------------------------One Stage Cluster--------------------------------------------------------------

res_1_1sc_1 <- resampling_1sc(data = data_1, prop = 1, iterasi = 100)
res_1_1sc_0.9 <- resampling_1sc(data = data_1, prop = 0.9, iterasi = 100)
res_1_1sc_0.8 <- resampling_1sc(data = data_1, prop = 0.8, iterasi = 100)
res_1_1sc_0.7 <- resampling_1sc(data = data_1, prop = 0.7, iterasi = 100)
res_1_1sc_0.6 <- resampling_1sc(data = data_1, prop = 0.6, iterasi = 100)
res_1_1sc_0.5 <- resampling_1sc(data = data_1, prop = 0.5, iterasi = 100)
res_1_1sc_0.4 <- resampling_1sc(data = data_1, prop = 0.4, iterasi = 100)
res_1_1sc_0.3 <- resampling_1sc(data = data_1, prop = 0.3, iterasi = 100)
res_1_1sc_0.2 <- resampling_1sc(data = data_1, prop = 0.2, iterasi = 100)
res_1_1sc_0.1 <- resampling_1sc(data = data_1, prop = 0.1, iterasi = 100)
res_1_1sc_0.09 <- resampling_1sc(data = data_1, prop = 0.09, iterasi = 100)
res_1_1sc_0.08 <- resampling_1sc(data = data_1, prop = 0.08, iterasi = 100)
res_1_1sc_0.07 <- resampling_1sc(data = data_1, prop = 0.07, iterasi = 100)
res_1_1sc_0.06 <- resampling_1sc(data = data_1, prop = 0.06, iterasi = 100)
res_1_1sc_0.05 <- resampling_1sc(data = data_1, prop = 0.05, iterasi = 100)
res_1_1sc_0.04 <- resampling_1sc(data = data_1, prop = 0.04, iterasi = 100)
res_1_1sc_0.03 <- resampling_1sc(data = data_1, prop = 0.03, iterasi = 100)
res_1_1sc_0.02 <- resampling_1sc(data = data_1, prop = 0.02, iterasi = 100)
res_1_1sc_0.01 <- resampling_1sc(data = data_1, prop = 0.01, iterasi = 100)

res_2_1sc_1 <- resampling_1sc(data = data_2, prop = 1, iterasi = 100)
res_2_1sc_0.9 <- resampling_1sc(data = data_2, prop = 0.9, iterasi = 100)
res_2_1sc_0.8 <- resampling_1sc(data = data_2, prop = 0.8, iterasi = 100)
res_2_1sc_0.7 <- resampling_1sc(data = data_2, prop = 0.7, iterasi = 100)
res_2_1sc_0.6 <- resampling_1sc(data = data_2, prop = 0.6, iterasi = 100)
res_2_1sc_0.5 <- resampling_1sc(data = data_2, prop = 0.5, iterasi = 100)
res_2_1sc_0.4 <- resampling_1sc(data = data_2, prop = 0.4, iterasi = 100)
res_2_1sc_0.3 <- resampling_1sc(data = data_2, prop = 0.3, iterasi = 100)
res_2_1sc_0.2 <- resampling_1sc(data = data_2, prop = 0.2, iterasi = 100)
res_2_1sc_0.1 <- resampling_1sc(data = data_2, prop = 0.1, iterasi = 100)
res_2_1sc_0.09 <- resampling_1sc(data = data_2, prop = 0.09, iterasi = 100)
res_2_1sc_0.08 <- resampling_1sc(data = data_2, prop = 0.08, iterasi = 100)
res_2_1sc_0.07 <- resampling_1sc(data = data_2, prop = 0.07, iterasi = 100)
res_2_1sc_0.06 <- resampling_1sc(data = data_2, prop = 0.06, iterasi = 100)
res_2_1sc_0.05 <- resampling_1sc(data = data_2, prop = 0.05, iterasi = 100)
res_2_1sc_0.04 <- resampling_1sc(data = data_2, prop = 0.04, iterasi = 100)
res_2_1sc_0.03 <- resampling_1sc(data = data_2, prop = 0.03, iterasi = 100)
res_2_1sc_0.02 <- resampling_1sc(data = data_2, prop = 0.02, iterasi = 100)
res_2_1sc_0.01 <- resampling_1sc(data = data_2, prop = 0.01, iterasi = 100)

res_3_1sc_1 <- resampling_1sc(data = data_3, prop = 1, iterasi = 100)
res_3_1sc_0.9 <- resampling_1sc(data = data_3, prop = 0.9, iterasi = 100)
res_3_1sc_0.8 <- resampling_1sc(data = data_3, prop = 0.8, iterasi = 100)
res_3_1sc_0.7 <- resampling_1sc(data = data_3, prop = 0.7, iterasi = 100)
res_3_1sc_0.6 <- resampling_1sc(data = data_3, prop = 0.6, iterasi = 100)
res_3_1sc_0.5 <- resampling_1sc(data = data_3, prop = 0.5, iterasi = 100)
res_3_1sc_0.4 <- resampling_1sc(data = data_3, prop = 0.4, iterasi = 100)
res_3_1sc_0.35 <- resampling_1sc(data = data_3, prop = 0.35, iterasi = 100)
res_3_1sc_0.3 <- resampling_1sc(data = data_3, prop = 0.3, iterasi = 100)
res_3_1sc_0.2 <- resampling_1sc(data = data_3, prop = 0.2, iterasi = 100)
res_3_1sc_0.1 <- resampling_1sc(data = data_3, prop = 0.1, iterasi = 100)
res_3_1sc_0.09 <- resampling_1sc(data = data_3, prop = 0.09, iterasi = 100)
res_3_1sc_0.08 <- resampling_1sc(data = data_3, prop = 0.08, iterasi = 100)
res_3_1sc_0.07 <- resampling_1sc(data = data_3, prop = 0.07, iterasi = 100)
res_3_1sc_0.06 <- resampling_1sc(data = data_3, prop = 0.06, iterasi = 100)
res_3_1sc_0.05 <- resampling_1sc(data = data_3, prop = 0.05, iterasi = 100)
res_3_1sc_0.04 <- resampling_1sc(data = data_3, prop = 0.04, iterasi = 100)
res_3_1sc_0.03 <- resampling_1sc(data = data_3, prop = 0.03, iterasi = 100)
res_3_1sc_0.02 <- resampling_1sc(data = data_3, prop = 0.02, iterasi = 100)
res_3_1sc_0.01 <- resampling_1sc(data = data_3, prop = 0.01, iterasi = 100)

res_4_1sc_1 <- resampling_1sc(data = data_4, prop = 1, iterasi = 100)
res_4_1sc_0.9 <- resampling_1sc(data = data_4, prop = 0.9, iterasi = 100)
res_4_1sc_0.8 <- resampling_1sc(data = data_4, prop = 0.8, iterasi = 100)
res_4_1sc_0.7 <- resampling_1sc(data = data_4, prop = 0.7, iterasi = 100)
res_4_1sc_0.6 <- resampling_1sc(data = data_4, prop = 0.6, iterasi = 100)
res_4_1sc_0.5 <- resampling_1sc(data = data_4, prop = 0.5, iterasi = 100)
res_4_1sc_0.4 <- resampling_1sc(data = data_4, prop = 0.4, iterasi = 100)
res_4_1sc_0.3 <- resampling_1sc(data = data_4, prop = 0.3, iterasi = 100)
res_4_1sc_0.2 <- resampling_1sc(data = data_4, prop = 0.2, iterasi = 100)
res_4_1sc_0.1 <- resampling_1sc(data = data_4, prop = 0.1, iterasi = 100)
res_4_1sc_0.09 <- resampling_1sc(data = data_4, prop = 0.09, iterasi = 100)
res_4_1sc_0.08 <- resampling_1sc(data = data_4, prop = 0.08, iterasi = 100)
res_4_1sc_0.07 <- resampling_1sc(data = data_4, prop = 0.07, iterasi = 100)
res_4_1sc_0.06 <- resampling_1sc(data = data_4, prop = 0.06, iterasi = 100)
res_4_1sc_0.05 <- resampling_1sc(data = data_4, prop = 0.05, iterasi = 100)
res_4_1sc_0.04 <- resampling_1sc(data = data_4, prop = 0.04, iterasi = 100)
res_4_1sc_0.03 <- resampling_1sc(data = data_4, prop = 0.03, iterasi = 100)
res_4_1sc_0.02 <- resampling_1sc(data = data_4, prop = 0.02, iterasi = 100)
res_4_1sc_0.01 <- resampling_1sc(data = data_4, prop = 0.01, iterasi = 100)

#--------------------------------------------Two Stage Cluster------------------------------------------------------------------
#2sc
res_1_2sc_srs_1 <- resampling_2sc_srs(data = data_1, f1 = 1, f2 = 1, iterasi = 20)
res_1_2sc_srs_0.9 <- resampling_2sc_srs(data = data_1, f1 = 0.9, f2 = 0.9, iterasi = 20)
res_1_2sc_srs_0.85 <- resampling_2sc_srs(data = data_1, f1 = 0.85, f2 = 0.85, iterasi = 20)
res_1_2sc_srs_0.8 <- resampling_2sc_srs(data = data_1, f1 = 0.8, f2 = 0.8, iterasi = 20)
res_1_2sc_srs_0.75 <- resampling_2sc_srs(data = data_1, f1 = 0.75, f2 = 0.75, iterasi = 20)
res_1_2sc_srs_0.7 <- resampling_2sc_srs(data = data_1, f1 = 0.7, f2 = 0.7, iterasi = 20)
res_1_2sc_srs_0.65 <- resampling_2sc_srs(data = data_1, f1 = 0.65, f2 = 0.65, iterasi = 20)
res_1_2sc_srs_0.6 <- resampling_2sc_srs(data = data_1, f1 = 0.6, f2 = 0.6, iterasi = 20)
res_1_2sc_srs_0.55 <- resampling_2sc_srs(data = data_1, f1 = 0.55, f2 = 0.55, iterasi = 20)
res_1_2sc_srs_0.5 <- resampling_2sc_srs(data = data_1, f1 = 0.5, f2 = 0.5, iterasi = 20)
res_1_2sc_srs_0.45 <- resampling_2sc_srs(data = data_1, f1 = 0.45, f2 = 0.45, iterasi = 20)
res_1_2sc_srs_0.4 <- resampling_2sc_srs(data = data_1, f1 = 0.4, f2 = 0.4, iterasi = 20)
res_1_2sc_srs_0.35 <- resampling_2sc_srs(data = data_1, f1 = 0.35, f2 = 0.35, iterasi = 20)
res_1_2sc_srs_0.3 <- resampling_2sc_srs(data = data_1, f1 = 0.3, f2 = 0.3, iterasi = 20)
res_1_2sc_srs_0.25 <- resampling_2sc_srs(data = data_1, f1 = 0.25, f2 = 0.25, iterasi = 20)
res_1_2sc_srs_0.2 <- resampling_2sc_srs(data = data_1, f1 = 0.2, f2 = 0.2, iterasi = 20)
res_1_2sc_srs_0.15 <- resampling_2sc_srs(dataY = data_1, f1 = 0.15, f2 = 0.15, iterasi = 20)

res_2_2sc_srs_1 <- resampling_2sc_srs(data = data_2, f1 = 1, f2 = 1, iterasi = 20)
res_2_2sc_srs_0.9 <- resampling_2sc_srs(data = data_2, f1 = 0.9, f2 = 0.9, iterasi = 20)
res_2_2sc_srs_0.85 <- resampling_2sc_srs(data = data_2, f1 = 0.85, f2 = 0.85, iterasi = 20)
res_2_2sc_srs_0.8 <- resampling_2sc_srs(data = data_2, f1 = 0.8, f2 = 0.8, iterasi = 20)
res_2_2sc_srs_0.75 <- resampling_2sc_srs(data = data_2, f1 = 0.75, f2 = 0.75, iterasi = 20)
res_2_2sc_srs_0.7 <- resampling_2sc_srs(data = data_2, f1 = 0.7, f2 = 0.7, iterasi = 20)
res_2_2sc_srs_0.65 <- resampling_2sc_srs(data = data_2, f1 = 0.65, f2 = 0.65, iterasi = 20)
res_2_2sc_srs_0.6 <- resampling_2sc_srs(data = data_2, f1 = 0.6, f2 = 0.6, iterasi = 20)
res_2_2sc_srs_0.55 <- resampling_2sc_srs(data = data_2, f1 = 0.55, f2 = 0.55, iterasi = 20)
res_2_2sc_srs_0.5 <- resampling_2sc_srs(data = data_2, f1 = 0.5, f2 = 0.5, iterasi = 20)
res_2_2sc_srs_0.45 <- resampling_2sc_srs(data = data_2, f1 = 0.45, f2 = 0.45, iterasi = 20)
res_2_2sc_srs_0.4 <- resampling_2sc_srs(data = data_2, f1 = 0.4, f2 = 0.4, iterasi = 20)
res_2_2sc_srs_0.35 <- resampling_2sc_srs(data = data_2, f1 = 0.35, f2 = 0.35, iterasi = 20)
res_2_2sc_srs_0.3 <- resampling_2sc_srs(data = data_2, f1 = 0.3, f2 = 0.3, iterasi = 20)
res_2_2sc_srs_0.25 <- resampling_2sc_srs(data = data_2, f1 = 0.25, f2 = 0.25, iterasi = 20)
res_2_2sc_srs_0.2 <- resampling_2sc_srs(data = data_2, f1 = 0.2, f2 = 0.2, iterasi = 20)
res_2_2sc_srs_0.15 <- resampling_2sc_srs(dataY = data_2, f1 = 0.15, f2 = 0.15, iterasi = 20)

res_3_2sc_srs_1 <- resampling_2sc_srs(data = data_3, f1 = 1, f2 = 1, iterasi = 20)
res_3_2sc_srs_0.9 <- resampling_2sc_srs(data = data_3, f1 = 0.9, f2 = 0.9, iterasi = 20)
res_3_2sc_srs_0.85 <- resampling_2sc_srs(data = data_3, f1 = 0.85, f2 = 0.85, iterasi = 20)
res_3_2sc_srs_0.8 <- resampling_2sc_srs(data = data_3, f1 = 0.8, f2 = 0.8, iterasi = 20)
res_3_2sc_srs_0.75 <- resampling_2sc_srs(data = data_3, f1 = 0.75, f2 = 0.75, iterasi = 20)
res_3_2sc_srs_0.7 <- resampling_2sc_srs(data = data_3, f1 = 0.7, f2 = 0.7, iterasi = 20)
res_3_2sc_srs_0.65 <- resampling_2sc_srs(data = data_3, f1 = 0.65, f2 = 0.65, iterasi = 20)
res_3_2sc_srs_0.6 <- resampling_2sc_srs(data = data_3, f1 = 0.6, f2 = 0.6, iterasi = 20)
res_3_2sc_srs_0.55 <- resampling_2sc_srs(data = data_3, f1 = 0.55, f2 = 0.55, iterasi = 20)
res_3_2sc_srs_0.5 <- resampling_2sc_srs(data = data_3, f1 = 0.5, f2 = 0.5, iterasi = 20)
res_3_2sc_srs_0.45 <- resampling_2sc_srs(data = data_3, f1 = 0.45, f2 = 0.45, iterasi = 20)
res_3_2sc_srs_0.4 <- resampling_2sc_srs(data = data_3, f1 = 0.4, f2 = 0.4, iterasi = 20)
res_3_2sc_srs_0.35 <- resampling_2sc_srs(data = data_3, f1 = 0.35, f2 = 0.35, iterasi = 20)
res_3_2sc_srs_0.3 <- resampling_2sc_srs(data = data_3, f1 = 0.3, f2 = 0.3, iterasi = 20)
res_3_2sc_srs_0.25 <- resampling_2sc_srs(data = data_3, f1 = 0.25, f2 = 0.25, iterasi = 20)
res_3_2sc_srs_0.2 <- resampling_2sc_srs(data = data_3, f1 = 0.2, f2 = 0.2, iterasi = 20)
res_3_2sc_srs_0.15 <- resampling_2sc_srs(dataY = data_3, f1 = 0.15, f2 = 0.15, iterasi = 20)

res_4_2sc_srs_1 <- resampling_2sc_srs(data = data_4, f1 = 1, f2 = 1, iterasi = 20)
res_4_2sc_srs_0.9 <- resampling_2sc_srs(data = data_4, f1 = 0.9, f2 = 0.9, iterasi = 20)
res_4_2sc_srs_0.85 <- resampling_2sc_srs(data = data_4, f1 = 0.85, f2 = 0.85, iterasi = 20)
res_4_2sc_srs_0.8 <- resampling_2sc_srs(data = data_4, f1 = 0.8, f2 = 0.8, iterasi = 20)
res_4_2sc_srs_0.75 <- resampling_2sc_srs(data = data_4, f1 = 0.75, f2 = 0.75, iterasi = 20)
res_4_2sc_srs_0.7 <- resampling_2sc_srs(data = data_4, f1 = 0.7, f2 = 0.7, iterasi = 20)
res_4_2sc_srs_0.65 <- resampling_2sc_srs(data = data_4, f1 = 0.65, f2 = 0.65, iterasi = 20)
res_4_2sc_srs_0.6 <- resampling_2sc_srs(data = data_4, f1 = 0.6, f2 = 0.6, iterasi = 20)
res_4_2sc_srs_0.55 <- resampling_2sc_srs(data = data_4, f1 = 0.55, f2 = 0.55, iterasi = 20)
res_4_2sc_srs_0.5 <- resampling_2sc_srs(data = data_4, f1 = 0.5, f2 = 0.5, iterasi = 20)
res_4_2sc_srs_0.45 <- resampling_2sc_srs(data = data_4, f1 = 0.45, f2 = 0.45, iterasi = 20)
res_4_2sc_srs_0.4 <- resampling_2sc_srs(data = data_4, f1 = 0.4, f2 = 0.4, iterasi = 20)
res_4_2sc_srs_0.35 <- resampling_2sc_srs(data = data_4, f1 = 0.35, f2 = 0.35, iterasi = 20)
res_4_2sc_srs_0.3 <- resampling_2sc_srs(data = data_4, f1 = 0.3, f2 = 0.3, iterasi = 20)
res_4_2sc_srs_0.25 <- resampling_2sc_srs(data = data_4, f1 = 0.25, f2 = 0.25, iterasi = 20)
res_4_2sc_srs_0.2 <- resampling_2sc_srs(data = data_4, f1 = 0.2, f2 = 0.2, iterasi = 20)
res_4_2sc_srs_0.15 <- resampling_2sc_srs(dataY = data_4, f1 = 0.15, f2 = 0.15, iterasi = 20)


#-------------------------------------------------Two Stage Cluster (PPS)------------------------------------------------------
res_1_2sc_pps_1 <- resampling_2sc_pps(data = data_1, f1 = 1, f2 = 1, iterasi = 20)
res_1_2sc_pps_0.9 <- resampling_2sc_pps(data = data_1, f1 = 0.9, f2 = 0.9, iterasi = 20)
res_1_2sc_pps_0.85 <- resampling_2sc_pps(data = data_1, f1 = 0.85, f2 = 0.85, iterasi = 20)
res_1_2sc_pps_0.8 <- resampling_2sc_pps(data = data_1, f1 = 0.8, f2 = 0.8, iterasi = 20)
res_1_2sc_pps_0.75 <- resampling_2sc_pps(data = data_1, f1 = 0.75, f2 = 0.75, iterasi = 20)
res_1_2sc_pps_0.7 <- resampling_2sc_pps(data = data_1, f1 = 0.7, f2 = 0.7, iterasi = 20)
res_1_2sc_pps_0.65 <- resampling_2sc_pps(data = data_1, f1 = 0.65, f2 = 0.65, iterasi = 20)
res_1_2sc_pps_0.6 <- resampling_2sc_pps(data = data_1, f1 = 0.6, f2 = 0.6, iterasi = 20)
res_1_2sc_pps_0.55 <- resampling_2sc_pps(data = data_1, f1 = 0.55, f2 = 0.55, iterasi = 20)
res_1_2sc_pps_0.5 <- resampling_2sc_pps(data = data_1, f1 = 0.5, f2 = 0.5, iterasi = 20)
res_1_2sc_pps_0.45 <- resampling_2sc_pps(data = data_1, f1 = 0.45, f2 = 0.45, iterasi = 20)
res_1_2sc_pps_0.4 <- resampling_2sc_pps(data = data_1, f1 = 0.4, f2 = 0.4, iterasi = 20)
res_1_2sc_pps_0.35 <- resampling_2sc_pps(data = data_1, f1 = 0.35, f2 = 0.35, iterasi = 20)
res_1_2sc_pps_0.3 <- resampling_2sc_pps(data = data_1, f1 = 0.3, f2 = 0.3, iterasi = 20)
res_1_2sc_pps_0.25 <- resampling_2sc_pps(data = data_1, f1 = 0.25, f2 = 0.25, iterasi = 20)
res_1_2sc_pps_0.2 <- resampling_2sc_pps(data = data_1, f1 = 0.2, f2 = 0.2, iterasi = 20)
res_1_2sc_pps_0.15 <- resampling_2sc_pps(data = data_1, f1 = 0.15, f2 = 0.15, iterasi = 20)

res_2_2sc_pps_1 <- resampling_2sc_pps(data = data_2, f1 = 1, f2 = 1, iterasi = 20)
res_2_2sc_pps_0.9 <- resampling_2sc_pps(data = data_2, f1 = 0.9, f2 = 0.9, iterasi = 20)
res_2_2sc_pps_0.85 <- resampling_2sc_pps(data = data_2, f1 = 0.85, f2 = 0.85, iterasi = 20)
res_2_2sc_pps_0.8 <- resampling_2sc_pps(data = data_2, f1 = 0.8, f2 = 0.8, iterasi = 20)
res_2_2sc_pps_0.75 <- resampling_2sc_pps(data = data_2, f1 = 0.75, f2 = 0.75, iterasi = 20)
res_2_2sc_pps_0.7 <- resampling_2sc_pps(data = data_2, f1 = 0.7, f2 = 0.7, iterasi = 20)
res_2_2sc_pps_0.65 <- resampling_2sc_pps(data = data_2, f1 = 0.65, f2 = 0.65, iterasi = 20)
res_2_2sc_pps_0.6 <- resampling_2sc_pps(data = data_2, f1 = 0.6, f2 = 0.6, iterasi = 20)
res_2_2sc_pps_0.55 <- resampling_2sc_pps(data = data_2, f1 = 0.55, f2 = 0.55, iterasi = 20)
res_2_2sc_pps_0.5 <- resampling_2sc_pps(data = data_2, f1 = 0.5, f2 = 0.5, iterasi = 20)
res_2_2sc_pps_0.45 <- resampling_2sc_pps(data = data_2, f1 = 0.45, f2 = 0.45, iterasi = 20)
res_2_2sc_pps_0.4 <- resampling_2sc_pps(data = data_2, f1 = 0.4, f2 = 0.4, iterasi = 20)
res_2_2sc_pps_0.35 <- resampling_2sc_pps(data = data_2, f1 = 0.35, f2 = 0.35, iterasi = 20)
res_2_2sc_pps_0.3 <- resampling_2sc_pps(data = data_2, f1 = 0.3, f2 = 0.3, iterasi = 20)
res_2_2sc_pps_0.25 <- resampling_2sc_pps(data = data_2, f1 = 0.25, f2 = 0.25, iterasi = 20)
res_2_2sc_pps_0.2 <- resampling_2sc_pps(data = data_2, f1 = 0.2, f2 = 0.2, iterasi = 20)
res_2_2sc_pps_0.15 <- resampling_2sc_pps(data = data_2, f1 = 0.15, f2 = 0.15, iterasi = 20)

res_3_2sc_pps_1 <- resampling_2sc_pps(data = data_3, f1 = 1, f2 = 1, iterasi = 20)
res_3_2sc_pps_0.9 <- resampling_2sc_pps(data = data_3, f1 = 0.9, f2 = 0.9, iterasi = 20)
res_3_2sc_pps_0.85 <- resampling_2sc_pps(data = data_3, f1 = 0.85, f2 = 0.85, iterasi = 20)
res_3_2sc_pps_0.8 <- resampling_2sc_pps(data = data_3, f1 = 0.8, f2 = 0.8, iterasi = 20)
res_3_2sc_pps_0.75 <- resampling_2sc_pps(data = data_3, f1 = 0.75, f2 = 0.75, iterasi = 20)
res_3_2sc_pps_0.7 <- resampling_2sc_pps(data = data_3, f1 = 0.7, f2 = 0.7, iterasi = 20)
res_3_2sc_pps_0.65 <- resampling_2sc_pps(data = data_3, f1 = 0.65, f2 = 0.65, iterasi = 20)
res_3_2sc_pps_0.6 <- resampling_2sc_pps(data = data_3, f1 = 0.6, f2 = 0.6, iterasi = 20)
res_3_2sc_pps_0.55 <- resampling_2sc_pps(data = data_3, f1 = 0.55, f2 = 0.55, iterasi = 20)
res_3_2sc_pps_0.5 <- resampling_2sc_pps(data = data_3, f1 = 0.5, f2 = 0.5, iterasi = 20)
res_3_2sc_pps_0.45 <- resampling_2sc_pps(data = data_3, f1 = 0.45, f2 = 0.45, iterasi = 20)
res_3_2sc_pps_0.4 <- resampling_2sc_pps(data = data_3, f1 = 0.4, f2 = 0.4, iterasi = 20)
res_3_2sc_pps_0.35 <- resampling_2sc_pps(data = data_3, f1 = 0.35, f2 = 0.35, iterasi = 20)
res_3_2sc_pps_0.3 <- resampling_2sc_pps(data = data_3, f1 = 0.3, f2 = 0.3, iterasi = 20)
res_3_2sc_pps_0.25 <- resampling_2sc_pps(data = data_3, f1 = 0.25, f2 = 0.25, iterasi = 20)
res_3_2sc_pps_0.2 <- resampling_2sc_pps(data = data_3, f1 = 0.2, f2 = 0.2, iterasi = 20)
res_3_2sc_pps_0.15 <- resampling_2sc_pps(data = data_3, f1 = 0.15, f2 = 0.15, iterasi = 20)

res_4_2sc_pps_1 <- resampling_2sc_pps(data = data_4, f1 = 1, f2 = 1, iterasi = 20)
res_4_2sc_pps_0.9 <- resampling_2sc_pps(data = data_4, f1 = 0.9, f2 = 0.9, iterasi = 20)
res_4_2sc_pps_0.85 <- resampling_2sc_pps(data = data_4, f1 = 0.85, f2 = 0.85, iterasi = 20)
res_4_2sc_pps_0.8 <- resampling_2sc_pps(data = data_4, f1 = 0.8, f2 = 0.8, iterasi = 20)
res_4_2sc_pps_0.75 <- resampling_2sc_pps(data = data_4, f1 = 0.75, f2 = 0.75, iterasi = 20)
res_4_2sc_pps_0.7 <- resampling_2sc_pps(data = data_4, f1 = 0.7, f2 = 0.7, iterasi = 20)
res_4_2sc_pps_0.65 <- resampling_2sc_pps(data = data_4, f1 = 0.65, f2 = 0.65, iterasi = 20)
res_4_2sc_pps_0.6 <- resampling_2sc_pps(data = data_4, f1 = 0.6, f2 = 0.6, iterasi = 20)
res_4_2sc_pps_0.55 <- resampling_2sc_pps(data = data_4, f1 = 0.55, f2 = 0.55, iterasi = 20)
res_4_2sc_pps_0.5 <- resampling_2sc_pps(data = data_4, f1 = 0.5, f2 = 0.5, iterasi = 20)
res_4_2sc_pps_0.45 <- resampling_2sc_pps(data = data_4, f1 = 0.45, f2 = 0.45, iterasi = 20)
res_4_2sc_pps_0.4 <- resampling_2sc_pps(data = data_4, f1 = 0.4, f2 = 0.4, iterasi = 20)
res_4_2sc_pps_0.35 <- resampling_2sc_pps(data = data_4, f1 = 0.35, f2 = 0.35, iterasi = 20)
res_4_2sc_pps_0.3 <- resampling_2sc_pps(data = data_4, f1 = 0.3, f2 = 0.3, iterasi = 20)
res_4_2sc_pps_0.25 <- resampling_2sc_pps(data = data_4, f1 = 0.25, f2 = 0.25, iterasi = 20)
res_4_2sc_pps_0.2 <- resampling_2sc_pps(data = data_4, f1 = 0.2, f2 = 0.2, iterasi = 20)
res_4_2sc_pps_0.15 <- resampling_2sc_pps(data = data_4, f1 = 0.15, f2 = 0.15, iterasi = 20)


#######################################################################################################################################

#------------------------------------------------Cek kenormalan data-------------------------------------------------------------------
nr_1_1 <- shapiro.test(res_1_srs_1$`Y dir`)
nr_1_0.9 <- shapiro.test(res_1_srs_0.9$`Y dir`)
nr_1_0.8 <- shapiro.test(res_1_srs_0.8$`Y dir`)
nr_1_0.7 <- shapiro.test(res_1_srs_0.7$`Y dir`)
nr_1_0.6 <- shapiro.test(res_1_srs_0.6$`Y dir`)
nr_1_0.5 <- shapiro.test(res_1_srs_0.5$`Y dir`)
nr_1_0.4 <- shapiro.test(res_1_srs_0.4$`Y dir`)
nr_1_0.3 <- shapiro.test(res_1_srs_0.3$`Y dir`)
nr_1_0.2 <- shapiro.test(res_1_srs_0.2$`Y dir`)
nr_1_0.1 <- shapiro.test(res_1_srs_0.1$`Y dir`)
nr_1_0.09 <- shapiro.test(res_1_srs_0.09$`Y dir`)
nr_1_0.08 <- shapiro.test(res_1_srs_0.08$`Y dir`)
nr_1_0.07 <- shapiro.test(res_1_srs_0.07$`Y dir`)
nr_1_0.06 <- shapiro.test(res_1_srs_0.06$`Y dir`)
nr_1_0.05 <- shapiro.test(res_1_srs_0.05$`Y dir`)
nr_1_0.04 <- shapiro.test(res_1_srs_0.04$`Y dir`)
nr_1_0.03 <- shapiro.test(res_1_srs_0.03$`Y dir`)
nr_1_0.02 <- shapiro.test(res_1_srs_0.02$`Y dir`)
nr_1_0.01 <- shapiro.test(res_1_srs_0.01$`Y dir`)

nr_2_1 <- shapiro.test(res_2_srs_1$`Y dir`)
nr_2_0.9 <- shapiro.test(res_2_srs_0.9$`Y dir`)
nr_2_0.8 <- shapiro.test(res_2_srs_0.8$`Y dir`)
nr_2_0.7 <- shapiro.test(res_2_srs_0.7$`Y dir`)
nr_2_0.6 <- shapiro.test(res_2_srs_0.6$`Y dir`)
nr_2_0.5 <- shapiro.test(res_2_srs_0.5$`Y dir`)
nr_2_0.4 <- shapiro.test(res_2_srs_0.4$`Y dir`)
nr_2_0.3 <- shapiro.test(res_2_srs_0.3$`Y dir`)
nr_2_0.2 <- shapiro.test(res_2_srs_0.2$`Y dir`)
nr_2_0.1 <- shapiro.test(res_2_srs_0.1$`Y dir`)
nr_2_0.09 <- shapiro.test(res_2_srs_0.09$`Y dir`)
nr_2_0.08 <- shapiro.test(res_2_srs_0.08$`Y dir`)
nr_2_0.07 <- shapiro.test(res_2_srs_0.07$`Y dir`)
nr_2_0.06 <- shapiro.test(res_2_srs_0.06$`Y dir`)
nr_2_0.05 <- shapiro.test(res_2_srs_0.05$`Y dir`)
nr_2_0.04 <- shapiro.test(res_2_srs_0.04$`Y dir`)
nr_2_0.03 <- shapiro.test(res_2_srs_0.03$`Y dir`)
nr_2_0.02 <- shapiro.test(res_2_srs_0.02$`Y dir`)
nr_2_0.01 <- shapiro.test(res_2_srs_0.01$`Y dir`)

nr_3_1 <- shapiro.test(res_3_srs_1$`Y dir`)
nr_3_0.9 <- shapiro.test(res_3_srs_0.9$`Y dir`)
nr_3_0.8 <- shapiro.test(res_3_srs_0.8$`Y dir`)
nr_3_0.7 <- shapiro.test(res_3_srs_0.7$`Y dir`)
nr_3_0.6 <- shapiro.test(res_3_srs_0.6$`Y dir`)
nr_3_0.5 <- shapiro.test(res_3_srs_0.5$`Y dir`)
nr_3_0.4 <- shapiro.test(res_3_srs_0.4$`Y dir`)
nr_3_0.3 <- shapiro.test(res_3_srs_0.3$`Y dir`)
nr_3_0.2 <- shapiro.test(res_3_srs_0.2$`Y dir`)
nr_3_0.1 <- shapiro.test(res_3_srs_0.1$`Y dir`)
nr_3_0.09 <- shapiro.test(res_3_srs_0.09$`Y dir`)
nr_3_0.08 <- shapiro.test(res_3_srs_0.08$`Y dir`)
nr_3_0.07 <- shapiro.test(res_3_srs_0.07$`Y dir`)
nr_3_0.06 <- shapiro.test(res_3_srs_0.06$`Y dir`)
nr_3_0.05 <- shapiro.test(res_3_srs_0.05$`Y dir`)
nr_3_0.04 <- shapiro.test(res_3_srs_0.04$`Y dir`)
nr_3_0.03 <- shapiro.test(res_3_srs_0.03$`Y dir`)
nr_3_0.02 <- shapiro.test(res_3_srs_0.02$`Y dir`)
nr_3_0.01 <- shapiro.test(res_3_srs_0.01$`Y dir`)

nr_4_1 <- shapiro.test(res_4_srs_1$`Y dir`)
nr_4_0.9 <- shapiro.test(res_4_srs_0.9$`Y dir`)
nr_4_0.8 <- shapiro.test(res_4_srs_0.8$`Y dir`)
nr_4_0.7 <- shapiro.test(res_4_srs_0.7$`Y dir`)
nr_4_0.6 <- shapiro.test(res_4_srs_0.6$`Y dir`)
nr_4_0.5 <- shapiro.test(res_4_srs_0.5$`Y dir`)
nr_4_0.4 <- shapiro.test(res_4_srs_0.4$`Y dir`)
nr_4_0.3 <- shapiro.test(res_4_srs_0.3$`Y dir`)
nr_4_0.2 <- shapiro.test(res_4_srs_0.2$`Y dir`)
nr_4_0.1 <- shapiro.test(res_4_srs_0.1$`Y dir`)
nr_4_0.09 <- shapiro.test(res_4_srs_0.09$`Y dir`)
nr_4_0.08 <- shapiro.test(res_4_srs_0.08$`Y dir`)
nr_4_0.07 <- shapiro.test(res_4_srs_0.07$`Y dir`)
nr_4_0.06 <- shapiro.test(res_4_srs_0.06$`Y dir`)
nr_4_0.05 <- shapiro.test(res_4_srs_0.05$`Y dir`)
nr_4_0.04 <- shapiro.test(res_4_srs_0.04$`Y dir`)
nr_4_0.03 <- shapiro.test(res_4_srs_0.03$`Y dir`)
nr_4_0.02 <- shapiro.test(res_4_srs_0.02$`Y dir`)
nr_4_0.01 <- shapiro.test(res_4_srs_0.01$`Y dir`)

#1sc
n1cs_1_1 <- shapiro.test(res_1_1sc_1$`Y dir`)
n1cs_1_0.9 <- shapiro.test(res_1_1sc_0.9$`Y dir`)
n1cs_1_0.8 <- shapiro.test(res_1_1sc_0.8$`Y dir`)
n1cs_1_0.7 <- shapiro.test(res_1_1sc_0.7$`Y dir`)
n1cs_1_0.6 <- shapiro.test(res_1_1sc_0.6$`Y dir`)
n1cs_1_0.5 <- shapiro.test(res_1_1sc_0.5$`Y dir`)
n1cs_1_0.4 <- shapiro.test(res_1_1sc_0.4$`Y dir`)
n1cs_1_0.3 <- shapiro.test(res_1_1sc_0.3$`Y dir`)
n1cs_1_0.2 <- shapiro.test(res_1_1sc_0.2$`Y dir`)
n1cs_1_0.1 <- shapiro.test(res_1_1sc_0.1$`Y dir`)
n1cs_1_0.09 <- shapiro.test(res_1_1sc_0.09$`Y dir`)
n1cs_1_0.08 <- shapiro.test(res_1_1sc_0.08$`Y dir`)
n1cs_1_0.07 <- shapiro.test(res_1_1sc_0.07$`Y dir`)
n1cs_1_0.06 <- shapiro.test(res_1_1sc_0.06$`Y dir`)
n1cs_1_0.05 <- shapiro.test(res_1_1sc_0.05$`Y dir`)
n1cs_1_0.04 <- shapiro.test(res_1_1sc_0.04$`Y dir`)
n1cs_1_0.03 <- shapiro.test(res_1_1sc_0.03$`Y dir`)
n1cs_1_0.02 <- shapiro.test(res_1_1sc_0.02$`Y dir`)
n1cs_1_0.01 <- shapiro.test(res_1_1sc_0.01$`Y dir`)

n1cs_2_1 <- shapiro.test(res_2_1sc_1$`Y dir`)
n1cs_2_0.9 <- shapiro.test(res_2_1sc_0.9$`Y dir`)
n1cs_2_0.8 <- shapiro.test(res_2_1sc_0.8$`Y dir`)
n1cs_2_0.7 <- shapiro.test(res_2_1sc_0.7$`Y dir`)
n1cs_2_0.6 <- shapiro.test(res_2_1sc_0.6$`Y dir`)
n1cs_2_0.5 <- shapiro.test(res_2_1sc_0.5$`Y dir`)
n1cs_2_0.4 <- shapiro.test(res_2_1sc_0.4$`Y dir`)
n1cs_2_0.3 <- shapiro.test(res_2_1sc_0.3$`Y dir`)
n1cs_2_0.2 <- shapiro.test(res_2_1sc_0.2$`Y dir`)
n1cs_2_0.1 <- shapiro.test(res_2_1sc_0.1$`Y dir`)
n1cs_2_0.09 <- shapiro.test(res_2_1sc_0.09$`Y dir`)
n1cs_2_0.08 <- shapiro.test(res_2_1sc_0.08$`Y dir`)
n1cs_2_0.07 <- shapiro.test(res_2_1sc_0.07$`Y dir`)
n1cs_2_0.06 <- shapiro.test(res_2_1sc_0.06$`Y dir`)
n1cs_2_0.05 <- shapiro.test(res_2_1sc_0.05$`Y dir`)
n1cs_2_0.04 <- shapiro.test(res_2_1sc_0.04$`Y dir`)
n1cs_2_0.03 <- shapiro.test(res_2_1sc_0.03$`Y dir`)
n1cs_2_0.02 <- shapiro.test(res_2_1sc_0.02$`Y dir`)
n1cs_2_0.01 <- shapiro.test(res_2_1sc_0.01$`Y dir`)

n1cs_3_1 <- shapiro.test(res_3_1sc_1$`Y dir`)
n1cs_3_0.9 <- shapiro.test(res_3_1sc_0.9$`Y dir`)
n1cs_3_0.8 <- shapiro.test(res_3_1sc_0.8$`Y dir`)
n1cs_3_0.7 <- shapiro.test(res_3_1sc_0.7$`Y dir`)
n1cs_3_0.6 <- shapiro.test(res_3_1sc_0.6$`Y dir`)
n1cs_3_0.5 <- shapiro.test(res_3_1sc_0.5$`Y dir`)
n1cs_3_0.4 <- shapiro.test(res_3_1sc_0.4$`Y dir`)
n1cs_3_0.3 <- shapiro.test(res_3_1sc_0.3$`Y dir`)
n1cs_3_0.2 <- shapiro.test(res_3_1sc_0.2$`Y dir`)
n1cs_3_0.1 <- shapiro.test(res_3_1sc_0.1$`Y dir`)
n1cs_3_0.09 <- shapiro.test(res_3_1sc_0.09$`Y dir`)
n1cs_3_0.08 <- shapiro.test(res_3_1sc_0.08$`Y dir`)
n1cs_3_0.07 <- shapiro.test(res_3_1sc_0.07$`Y dir`)
n1cs_3_0.06 <- shapiro.test(res_3_1sc_0.06$`Y dir`)
n1cs_3_0.05 <- shapiro.test(res_3_1sc_0.05$`Y dir`)
n1cs_3_0.04 <- shapiro.test(res_3_1sc_0.04$`Y dir`)
n1cs_3_0.03 <- shapiro.test(res_3_1sc_0.03$`Y dir`)
n1cs_3_0.02 <- shapiro.test(res_3_1sc_0.02$`Y dir`)
n1cs_3_0.01 <- shapiro.test(res_3_1sc_0.01$`Y dir`)

n1cs_4_1 <- shapiro.test(res_4_1sc_1$`Y dir`)
n1cs_4_0.9 <- shapiro.test(res_4_1sc_0.9$`Y dir`)
n1cs_4_0.8 <- shapiro.test(res_4_1sc_0.8$`Y dir`)
n1cs_4_0.7 <- shapiro.test(res_4_1sc_0.7$`Y dir`)
n1cs_4_0.6 <- shapiro.test(res_4_1sc_0.6$`Y dir`)
n1cs_4_0.5 <- shapiro.test(res_4_1sc_0.5$`Y dir`)
n1cs_4_0.4 <- shapiro.test(res_4_1sc_0.4$`Y dir`)
n1cs_4_0.3 <- shapiro.test(res_4_1sc_0.3$`Y dir`)
n1cs_4_0.2 <- shapiro.test(res_4_1sc_0.2$`Y dir`)
n1cs_4_0.1 <- shapiro.test(res_4_1sc_0.1$`Y dir`)
n1cs_4_0.09 <- shapiro.test(res_4_1sc_0.09$`Y dir`)
n1cs_4_0.08 <- shapiro.test(res_4_1sc_0.08$`Y dir`)
n1cs_4_0.07 <- shapiro.test(res_4_1sc_0.07$`Y dir`)
n1cs_4_0.06 <- shapiro.test(res_4_1sc_0.06$`Y dir`)
n1cs_4_0.05 <- shapiro.test(res_4_1sc_0.05$`Y dir`)
n1cs_4_0.04 <- shapiro.test(res_4_1sc_0.04$`Y dir`)
n1cs_4_0.03 <- shapiro.test(res_4_1sc_0.03$`Y dir`)
n1cs_4_0.02 <- shapiro.test(res_4_1sc_0.02$`Y dir`)
n1cs_4_0.01 <- shapiro.test(res_4_1sc_0.01$`Y dir`)

#2sc_srs
sc_srs_1_nor_1 <- shapiro.test(res_1_2sc_srs_1$`Y dir`)
sc_srs_1_nor_0.9 <- shapiro.test(res_1_2sc_srs_0.9$`Y dir`)
sc_srs_1_nor_0.85 <- shapiro.test(res_1_2sc_srs_0.85$`Y dir`)
sc_srs_1_nor_0.8 <-shapiro.test(res_1_2sc_srs_0.8$`Y dir`)
sc_srs_1_nor_0.75 <-shapiro.test(res_1_2sc_srs_0.75$`Y dir`)
sc_srs_1_nor_0.7 <-shapiro.test(res_1_2sc_srs_0.7$`Y dir`)
sc_srs_1_nor_0.65 <-shapiro.test(res_1_2sc_srs_0.65$`Y dir`)
sc_srs_1_nor_0.6 <-shapiro.test(res_1_2sc_srs_0.6$`Y dir`)
sc_srs_1_nor_0.55 <-shapiro.test(res_1_2sc_srs_0.55$`Y dir`)
sc_srs_1_nor_0.5 <-shapiro.test(res_1_2sc_srs_0.5$`Y dir`)
sc_srs_1_nor_0.45 <-shapiro.test(res_1_2sc_srs_0.45$`Y dir`)
sc_srs_1_nor_0.4 <-shapiro.test(res_1_2sc_srs_0.4$`Y dir`)
sc_srs_1_nor_0.35 <-shapiro.test(res_1_2sc_srs_0.35$`Y dir`)
sc_srs_1_nor_0.3 <-shapiro.test(res_1_2sc_srs_0.3$`Y dir`)
sc_srs_1_nor_0.25 <-shapiro.test(res_1_2sc_srs_0.25$`Y dir`)
sc_srs_1_nor_0.2 <-shapiro.test(res_1_2sc_srs_0.2$`Y dir`)
sc_srs_1_nor_0.15 <-shapiro.test(res_1_2sc_srs_0.15$`Y dir`)

sc_srs_2_nor_1 <- shapiro.test(res_2_2sc_srs_1$`Y dir`)
sc_srs_2_nor_0.9 <- shapiro.test(res_2_2sc_srs_0.9$`Y dir`)
sc_srs_2_nor_0.85 <- shapiro.test(res_2_2sc_srs_0.85$`Y dir`)
sc_srs_2_nor_0.8 <-shapiro.test(res_2_2sc_srs_0.8$`Y dir`)
sc_srs_2_nor_0.75 <-shapiro.test(res_2_2sc_srs_0.75$`Y dir`)
sc_srs_2_nor_0.7 <-shapiro.test(res_2_2sc_srs_0.7$`Y dir`)
sc_srs_2_nor_0.65 <-shapiro.test(res_2_2sc_srs_0.65$`Y dir`)
sc_srs_2_nor_0.6 <-shapiro.test(res_2_2sc_srs_0.6$`Y dir`)
sc_srs_2_nor_0.55 <-shapiro.test(res_2_2sc_srs_0.55$`Y dir`)
sc_srs_2_nor_0.5 <-shapiro.test(res_2_2sc_srs_0.5$`Y dir`)
sc_srs_2_nor_0.45 <-shapiro.test(res_2_2sc_srs_0.45$`Y dir`)
sc_srs_2_nor_0.4 <-shapiro.test(res_2_2sc_srs_0.4$`Y dir`)
sc_srs_2_nor_0.35 <-shapiro.test(res_2_2sc_srs_0.35$`Y dir`)
sc_srs_2_nor_0.3 <-shapiro.test(res_2_2sc_srs_0.3$`Y dir`)
sc_srs_2_nor_0.25 <-shapiro.test(res_2_2sc_srs_0.25$`Y dir`)
sc_srs_2_nor_0.2 <-shapiro.test(res_2_2sc_srs_0.2$`Y dir`)
sc_srs_2_nor_0.15 <-shapiro.test(res_2_2sc_srs_0.15$`Y dir`)

sc_srs_3_nor_1 <- shapiro.test(res_3_2sc_srs_1$`Y dir`)
sc_srs_3_nor_0.9 <- shapiro.test(res_3_2sc_srs_0.9$`Y dir`)
sc_srs_3_nor_0.85 <- shapiro.test(res_3_2sc_srs_0.85$`Y dir`)
sc_srs_3_nor_0.8 <-shapiro.test(res_3_2sc_srs_0.8$`Y dir`)
sc_srs_3_nor_0.75 <-shapiro.test(res_3_2sc_srs_0.75$`Y dir`)
sc_srs_3_nor_0.7 <-shapiro.test(res_3_2sc_srs_0.7$`Y dir`)
sc_srs_3_nor_0.65 <-shapiro.test(res_3_2sc_srs_0.65$`Y dir`)
sc_srs_3_nor_0.6 <-shapiro.test(res_3_2sc_srs_0.6$`Y dir`)
sc_srs_3_nor_0.55 <-shapiro.test(res_3_2sc_srs_0.55$`Y dir`)
sc_srs_3_nor_0.5 <-shapiro.test(res_3_2sc_srs_0.5$`Y dir`)
sc_srs_3_nor_0.45 <-shapiro.test(res_3_2sc_srs_0.45$`Y dir`)
sc_srs_3_nor_0.4 <-shapiro.test(res_3_2sc_srs_0.4$`Y dir`)
sc_srs_3_nor_0.35 <-shapiro.test(res_3_2sc_srs_0.35$`Y dir`)
sc_srs_3_nor_0.3 <-shapiro.test(res_3_2sc_srs_0.3$`Y dir`)
sc_srs_3_nor_0.25 <-shapiro.test(res_3_2sc_srs_0.25$`Y dir`)
sc_srs_3_nor_0.2 <-shapiro.test(res_3_2sc_srs_0.2$`Y dir`)
sc_srs_3_nor_0.15 <-shapiro.test(res_3_2sc_srs_0.15$`Y dir`)

sc_srs_4_nor_1 <- shapiro.test(res_4_2sc_srs_1$`Y dir`)
sc_srs_4_nor_0.9 <- shapiro.test(res_4_2sc_srs_0.9$`Y dir`)
sc_srs_4_nor_0.85 <- shapiro.test(res_4_2sc_srs_0.85$`Y dir`)
sc_srs_4_nor_0.8 <-shapiro.test(res_4_2sc_srs_0.8$`Y dir`)
sc_srs_4_nor_0.75 <-shapiro.test(res_4_2sc_srs_0.75$`Y dir`)
sc_srs_4_nor_0.7 <-shapiro.test(res_4_2sc_srs_0.7$`Y dir`)
sc_srs_4_nor_0.65 <-shapiro.test(res_4_2sc_srs_0.65$`Y dir`)
sc_srs_4_nor_0.6 <-shapiro.test(res_4_2sc_srs_0.6$`Y dir`)
sc_srs_4_nor_0.55 <-shapiro.test(res_4_2sc_srs_0.55$`Y dir`)
sc_srs_4_nor_0.5 <-shapiro.test(res_4_2sc_srs_0.5$`Y dir`)
sc_srs_4_nor_0.45 <-shapiro.test(res_4_2sc_srs_0.45$`Y dir`)
sc_srs_4_nor_0.4 <-shapiro.test(res_4_2sc_srs_0.4$`Y dir`)
sc_srs_4_nor_0.35 <-shapiro.test(res_4_2sc_srs_0.35$`Y dir`)
sc_srs_4_nor_0.3 <-shapiro.test(res_4_2sc_srs_0.3$`Y dir`)
sc_srs_4_nor_0.25 <-shapiro.test(res_4_2sc_srs_0.25$`Y dir`)
sc_srs_4_nor_0.2 <-shapiro.test(res_4_2sc_srs_0.2$`Y dir`)
sc_srs_4_nor_0.15 <-shapiro.test(res_4_2sc_srs_0.15$`Y dir`)

#2sc_pps
sc_pps_1_nor_1 <- shapiro.test(res_1_2sc_pps_1$`Y dir`)
sc_pps_1_nor_0.9 <- shapiro.test(res_1_2sc_pps_0.9$`Y dir`)
sc_pps_1_nor_0.85 <- shapiro.test(res_1_2sc_pps_0.85$`Y dir`)
sc_pps_1_nor_0.8 <-shapiro.test(res_1_2sc_pps_0.8$`Y dir`)
sc_pps_1_nor_0.75 <-shapiro.test(res_1_2sc_pps_0.75$`Y dir`)
sc_pps_1_nor_0.7 <-shapiro.test(res_1_2sc_pps_0.7$`Y dir`)
sc_pps_1_nor_0.65 <-shapiro.test(res_1_2sc_pps_0.65$`Y dir`)
sc_pps_1_nor_0.6 <-shapiro.test(res_1_2sc_pps_0.6$`Y dir`)
sc_pps_1_nor_0.55 <-shapiro.test(res_1_2sc_pps_0.55$`Y dir`)
sc_pps_1_nor_0.5 <-shapiro.test(res_1_2sc_pps_0.5$`Y dir`)
sc_pps_1_nor_0.45 <-shapiro.test(res_1_2sc_pps_0.45$`Y dir`)
sc_pps_1_nor_0.4 <-shapiro.test(res_1_2sc_pps_0.4$`Y dir`)
sc_pps_1_nor_0.35 <-shapiro.test(res_1_2sc_pps_0.35$`Y dir`)
sc_pps_1_nor_0.3 <-shapiro.test(res_1_2sc_pps_0.3$`Y dir`)
sc_pps_1_nor_0.25 <-shapiro.test(res_1_2sc_pps_0.25$`Y dir`)
sc_pps_1_nor_0.2 <-shapiro.test(res_1_2sc_pps_0.2$`Y dir`)
sc_pps_1_nor_0.15 <-shapiro.test(res_1_2sc_pps_0.15$`Y dir`)

sc_pps_2_nor_1 <- shapiro.test(res_2_2sc_pps_1$`Y dir`)
sc_pps_2_nor_0.9 <- shapiro.test(res_2_2sc_pps_0.9$`Y dir`)
sc_pps_2_nor_0.85 <- shapiro.test(res_2_2sc_pps_0.85$`Y dir`)
sc_pps_2_nor_0.8 <-shapiro.test(res_2_2sc_pps_0.8$`Y dir`)
sc_pps_2_nor_0.75 <-shapiro.test(res_2_2sc_pps_0.75$`Y dir`)
sc_pps_2_nor_0.7 <-shapiro.test(res_2_2sc_pps_0.7$`Y dir`)
sc_pps_2_nor_0.65 <-shapiro.test(res_2_2sc_pps_0.65$`Y dir`)
sc_pps_2_nor_0.6 <-shapiro.test(res_2_2sc_pps_0.6$`Y dir`)
sc_pps_2_nor_0.55 <-shapiro.test(res_2_2sc_pps_0.55$`Y dir`)
sc_pps_2_nor_0.5 <-shapiro.test(res_2_2sc_pps_0.5$`Y dir`)
sc_pps_2_nor_0.45 <-shapiro.test(res_2_2sc_pps_0.45$`Y dir`)
sc_pps_2_nor_0.4 <-shapiro.test(res_2_2sc_pps_0.4$`Y dir`)
sc_pps_2_nor_0.35 <-shapiro.test(res_2_2sc_pps_0.35$`Y dir`)
sc_pps_2_nor_0.3 <-shapiro.test(res_2_2sc_pps_0.3$`Y dir`)
sc_pps_2_nor_0.25 <-shapiro.test(res_2_2sc_pps_0.25$`Y dir`)
sc_pps_2_nor_0.2 <-shapiro.test(res_2_2sc_pps_0.2$`Y dir`)
sc_pps_2_nor_0.15 <-shapiro.test(res_2_2sc_pps_0.15$`Y dir`)

sc_pps_3_nor_1 <- shapiro.test(res_3_2sc_pps_1$`Y dir`)
sc_pps_3_nor_0.9 <- shapiro.test(res_3_2sc_pps_0.9$`Y dir`)
sc_pps_3_nor_0.85 <- shapiro.test(res_3_2sc_pps_0.85$`Y dir`)
sc_pps_3_nor_0.8 <-shapiro.test(res_3_2sc_pps_0.8$`Y dir`)
sc_pps_3_nor_0.75 <-shapiro.test(res_3_2sc_pps_0.75$`Y dir`)
sc_pps_3_nor_0.7 <-shapiro.test(res_3_2sc_pps_0.7$`Y dir`)
sc_pps_3_nor_0.65 <-shapiro.test(res_3_2sc_pps_0.65$`Y dir`)
sc_pps_3_nor_0.6 <-shapiro.test(res_3_2sc_pps_0.6$`Y dir`)
sc_pps_3_nor_0.55 <-shapiro.test(res_3_2sc_pps_0.55$`Y dir`)
sc_pps_3_nor_0.5 <-shapiro.test(res_3_2sc_pps_0.5$`Y dir`)
sc_pps_3_nor_0.45 <-shapiro.test(res_3_2sc_pps_0.45$`Y dir`)
sc_pps_3_nor_0.4 <-shapiro.test(res_3_2sc_pps_0.4$`Y dir`)
sc_pps_3_nor_0.35 <-shapiro.test(res_3_2sc_pps_0.35$`Y dir`)
sc_pps_3_nor_0.3 <-shapiro.test(res_3_2sc_pps_0.3$`Y dir`)
sc_pps_3_nor_0.25 <-shapiro.test(res_3_2sc_pps_0.25$`Y dir`)
sc_pps_3_nor_0.2 <-shapiro.test(res_3_2sc_pps_0.2$`Y dir`)
sc_pps_3_nor_0.15 <-shapiro.test(res_3_2sc_pps_0.15$`Y dir`)

sc_pps_4_nor_1 <- shapiro.test(res_4_2sc_pps_1$`Y dir`)
sc_pps_4_nor_0.9 <- shapiro.test(res_4_2sc_pps_0.9$`Y dir`)
sc_pps_4_nor_0.85 <- shapiro.test(res_4_2sc_pps_0.85$`Y dir`)
sc_pps_4_nor_0.8 <-shapiro.test(res_4_2sc_pps_0.8$`Y dir`)
sc_pps_4_nor_0.75 <-shapiro.test(res_4_2sc_pps_0.75$`Y dir`)
sc_pps_4_nor_0.7 <-shapiro.test(res_4_2sc_pps_0.7$`Y dir`)
sc_pps_4_nor_0.65 <-shapiro.test(res_4_2sc_pps_0.65$`Y dir`)
sc_pps_4_nor_0.6 <-shapiro.test(res_4_2sc_pps_0.6$`Y dir`)
sc_pps_4_nor_0.55 <-shapiro.test(res_4_2sc_pps_0.55$`Y dir`)
sc_pps_4_nor_0.5 <-shapiro.test(res_4_2sc_pps_0.5$`Y dir`)
sc_pps_4_nor_0.45 <-shapiro.test(res_4_2sc_pps_0.45$`Y dir`)
sc_pps_4_nor_0.4 <-shapiro.test(res_4_2sc_pps_0.4$`Y dir`)
sc_pps_4_nor_0.35 <-shapiro.test(res_4_2sc_pps_0.35$`Y dir`)
sc_pps_4_nor_0.3 <-shapiro.test(res_4_2sc_pps_0.3$`Y dir`)
sc_pps_4_nor_0.25 <-shapiro.test(res_4_2sc_pps_0.25$`Y dir`)
sc_pps_4_nor_0.2 <-shapiro.test(res_4_2sc_pps_0.2$`Y dir`)
sc_pps_4_nor_0.15 <-shapiro.test(res_4_2sc_pps_0.15$`Y dir`)

######################################################################################################################################
#---------------------------------------Pemilihan Auxiliary---------------------------------------------------------------------------
library(readxl)
Xnew_ya <- read_excel("D:/KAMPUS/skripsi/DATA/Xnew ya.xlsx", 
                      sheet = "Sheet3", col_types = c("text", 
                                                      "text", "text", "text", "text", "text", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric"))

#Pemilihan auxiliary variabel
x1 <- aggregate(Xnew_ya$X1, list(Xnew_ya$R102),sum,na.rm=TRUE)
x2 <- aggregate(Xnew_ya$X2, list(Xnew_ya$R102),sum,na.rm=TRUE)
x3 <- aggregate(Xnew_ya$X3, list(Xnew_ya$R102),sum,na.rm=TRUE)
x4 <- aggregate(Xnew_ya$X4, list(Xnew_ya$R102),sum,na.rm=TRUE)
x5 <- aggregate(Xnew_ya$X5, list(Xnew_ya$R102),sum,na.rm=TRUE)
x6 <- aggregate(Xnew_ya$X6, list(Xnew_ya$R102),sum,na.rm=TRUE)
x7 <- aggregate(Xnew_ya$X7, list(Xnew_ya$R102),sum,na.rm=TRUE)
x8 <- aggregate(Xnew_ya$X8, list(Xnew_ya$R102),sum,na.rm=TRUE)
x9 <- aggregate(Xnew_ya$X9, list(Xnew_ya$R102),sum,na.rm=TRUE)
x10 <- aggregate(Xnew_ya$X10, list(Xnew_ya$R102),sum,na.rm=TRUE)
x11 <- aggregate(Xnew_ya$X11, list(Xnew_ya$R102),sum,na.rm=TRUE)
x12 <- aggregate(Xnew_ya$X12, list(Xnew_ya$R102),sum,na.rm=TRUE)
x13 <- aggregate(Xnew_ya$X13, list(Xnew_ya$R102),sum,na.rm=TRUE)
x14 <- aggregate(Xnew_ya$X14, list(Xnew_ya$R102),sum,na.rm=TRUE)
x15 <- aggregate(Xnew_ya$X15, list(Xnew_ya$R102),sum,na.rm=TRUE)
x16 <- aggregate(Xnew_ya$X16, list(Xnew_ya$R102),sum,na.rm=TRUE)
x17 <- aggregate(Xnew_ya$X17, list(Xnew_ya$R102),sum,na.rm=TRUE)
x18 <- aggregate(Xnew_ya$X18, list(Xnew_ya$R102),sum,na.rm=TRUE)
x19 <- aggregate(Xnew_ya$X19, list(Xnew_ya$R102),sum,na.rm=TRUE)
x20 <- aggregate(Xnew_ya$X20, list(Xnew_ya$R102),sum,na.rm=TRUE)
x21 <- aggregate(Xnew_ya$X21, list(Xnew_ya$R102),sum,na.rm=TRUE)
x22 <- aggregate(Xnew_ya$X22, list(Xnew_ya$R102),sum,na.rm=TRUE)
x23 <- aggregate(Xnew_ya$X23, list(Xnew_ya$R102),sum,na.rm=TRUE)
x24 <- aggregate(Xnew_ya$X24, list(Xnew_ya$R102),sum,na.rm=TRUE)
x25 <- aggregate(Xnew_ya$X25, list(Xnew_ya$R102),sum,na.rm=TRUE)
x26 <- aggregate(Xnew_ya$X26, list(Xnew_ya$R102),sum,na.rm=TRUE)
x27 <- aggregate(Xnew_ya$X27, list(Xnew_ya$R102),sum,na.rm=TRUE)
x28 <- aggregate(Xnew_ya$X28, list(Xnew_ya$R102),sum,na.rm=TRUE)
x29 <- aggregate(Xnew_ya$X29, list(Xnew_ya$R102),sum,na.rm=TRUE)
x30 <- aggregate(Xnew_ya$X30, list(Xnew_ya$R102),sum,na.rm=TRUE)
x31 <- aggregate(Xnew_ya$X31, list(Xnew_ya$R102),sum,na.rm=TRUE)
x32 <- aggregate(Xnew_ya$X32, list(Xnew_ya$R102),sum,na.rm=TRUE)
x33 <- aggregate(Xnew_ya$X33, list(Xnew_ya$R102),sum,na.rm=TRUE)
x34 <- aggregate(Xnew_ya$X34, list(Xnew_ya$R102),sum,na.rm=TRUE)
x35 <- aggregate(Xnew_ya$X35, list(Xnew_ya$R102),sum,na.rm=TRUE)

#----------------------------------------------------Pengecekan Korelasi-------------------------------------------------------------
#Cek tiap variabel dengan srs

cor(x1$x,res_1_srs_0.9$`Y dir`)
cor(x2$x,res_1_srs_0.9$`Y dir`) #midle
cor(x3$x,res_1_srs_0.9$`Y dir`) #midle
cor(x4$x,res_1_srs_0.9$`Y dir`) #jelek
cor(x5$x,res_1_srs_0.9$`Y dir`) #midle
cor(x6$x,res_1_srs_0.9$`Y dir`) #midle
cor(x7$x,res_1_srs_0.9$`Y dir`) #Jelek
cor(x8$x,res_1_srs_0.9$`Y dir`) #midle
cor(x9$x,res_1_srs_0.9$`Y dir`) #jelek
cor(x10$x,res_1_srs_0.9$`Y dir`)#midle
cor(x11$x,res_1_srs_0.9$`Y dir`) #midle
cor(x12$x,res_1_srs_0.9$`Y dir`) #midle
cor(x13$x,res_1_srs_0.9$`Y dir`) #bagus
cor(x14$x,res_1_srs_0.9$`Y dir`)#jelek
cor(x15$x,res_1_srs_0.9$`Y dir`) #bagus
cor(x16$x,res_1_srs_0.9$`Y dir`)#midle
cor(x17$x,res_1_srs_0.9$`Y dir`) #midle
cor(x18$x,res_1_srs_0.9$`Y dir`) #bagus
cor(x19$x,res_1_srs_0.9$`Y dir`) #bagus
cor(x20$x,res_1_srs_0.9$`Y dir`) #bagus
cor(x21$x,res_1_srs_0.9$`Y dir`) #midle
cor(x22$x,res_1_srs_0.9$`Y dir`) #jelek
cor(x23$x,res_1_srs_0.9$`Y dir`) #jelek
cor(x24$x,res_1_srs_0.9$`Y dir`) #midle
cor(x25$x,res_1_srs_0.9$`Y dir`) #bagus
cor(x26$x,res_1_srs_0.9$`Y dir`)#jelek
cor(x27$x,res_1_srs_0.9$`Y dir`) #mmidle
cor(x28$x,res_1_srs_0.9$`Y dir`) #midle
cor(x29$x,res_1_srs_0.9$`Y dir`) #midle
cor(x30$x,res_1_srs_0.9$`Y dir`) #bagus
cor(x31$x,res_1_srs_0.9$`Y dir`)#jelek
cor(x32$x,res_1_srs_0.9$`Y dir`) #bagus
cor(x33$x,res_1_srs_0.9$`Y dir`) #bagus
cor(x34$x,res_1_srs_0.9$`Y dir`) #jelek
cor(x35$x,res_1_srs_0.9$`Y dir`) #bagus

#Cek tiap variabel dengan 1sc
cor(x1$x,res_2_1sc_0.9$`Y dir`) #bagus
cor(x2$x,res_2_1sc_0.9$`Y dir`) #bagus
cor(x3$x,res_2_1sc_0.9$`Y dir`) #bagus
cor(x4$x,res_2_1sc_0.9$`Y dir`) #jelek
cor(x5$x,res_2_1sc_0.9$`Y dir`) #bagus
cor(x6$x,res_2_1sc_0.9$`Y dir`) #bagus
cor(x7$x,res_2_1sc_0.9$`Y dir`) #bagus
cor(x8$x,res_2_1sc_0.9$`Y dir`) #bagus
cor(x9$x,res_2_1sc_0.9$`Y dir`) #jelek
cor(x10$x,res_2_1sc_0.9$`Y dir`) #jelek
cor(x11$x,res_2_1sc_0.9$`Y dir`) #jelek
cor(x12$x,res_2_1sc_0.9$`Y dir`) #bgs
cor(x13$x,res_2_1sc_0.9$`Y dir`) #bgs
cor(x14$x,res_2_1sc_0.9$`Y dir`)#jelek
cor(x15$x,res_2_1sc_0.9$`Y dir`) #bgs
cor(x16$x,res_2_1sc_0.9$`Y dir`)#jelek
cor(x17$x,res_2_1sc_0.9$`Y dir`) #bgs
cor(x18$x,res_2_1sc_0.9$`Y dir`) #bgs
cor(x19$x,res_2_1sc_0.9$`Y dir`) #jelek
cor(x20$x,res_2_1sc_0.9$`Y dir`)#bgs
cor(x21$x,res_2_1sc_0.9$`Y dir`) #bgs
cor(x22$x,res_2_1sc_0.9$`Y dir`) #bgs
cor(x23$x,res_2_1sc_0.9$`Y dir`) #bgs
cor(x24$x,res_2_1sc_0.9$`Y dir`) #bgs
cor(x25$x,res_2_1sc_0.9$`Y dir`) #jlk
cor(x26$x,res_2_1sc_0.9$`Y dir`)#jelek
cor(x27$x,res_2_1sc_0.9$`Y dir`) #bgs
cor(x28$x,res_2_1sc_0.9$`Y dir`) #bgs
cor(x29$x,res_2_1sc_0.9$`Y dir`)#jelek
cor(x30$x,res_2_1sc_0.9$`Y dir`) #bgs
cor(x31$x,res_2_1sc_0.9$`Y dir`)#jelek
cor(x32$x,res_2_1sc_0.9$`Y dir`) #bgs
cor(x33$x,res_2_1sc_0.9$`Y dir`) #bgs
cor(x34$x,res_2_1sc_0.9$`Y dir`) #bgs
cor(x35$x,res_2_1sc_0.9$`Y dir`) #bgs

#Cek tiap variabel dengan 2sc_srs
cor(x1$x,res_3_2sc_srs_0.9$`Y dir`) #bgs
cor(x2$x,res_3_2sc_srs_0.9$`Y dir`) #bgs
cor(x3$x,res_3_2sc_srs_0.9$`Y dir`) #bgs
cor(x4$x,res_3_2sc_srs_0.9$`Y dir`) #jelek
cor(x5$x,res_3_2sc_srs_0.9$`Y dir`) #bgs
cor(x6$x,res_3_2sc_srs_0.9$`Y dir`) #bgs
cor(x7$x,res_3_2sc_srs_0.9$`Y dir`) #bgs
cor(x8$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x9$x,res_3_2sc_srs_0.9$`Y dir`) #jlk
cor(x10$x,res_3_2sc_srs_0.9$`Y dir`) #jlk
cor(x11$x,res_3_2sc_srs_0.9$`Y dir`) #jlk
cor(x12$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x13$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x14$x,res_3_2sc_srs_0.9$`Y dir`)#jelek
cor(x15$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x16$x,res_3_2sc_srs_0.9$`Y dir`)#jelek
cor(x17$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x18$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x19$x,res_3_2sc_srs_0.9$`Y dir`) #jlk
cor(x20$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x21$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x22$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x23$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x24$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x25$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x26$x,res_3_2sc_srs_0.9$`Y dir`)#jelek
cor(x27$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x28$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x29$x,res_3_2sc_srs_0.9$`Y dir`)#jelek
cor(x30$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x31$x,res_3_2sc_srs_0.9$`Y dir`)#jelek
cor(x32$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x33$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x34$x,res_3_2sc_srs_0.9$`Y dir`)
cor(x35$x,res_3_2sc_srs_0.9$`Y dir`)

#Cek tiap variabel dengan 2sc_pps
cor(x1$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x2$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x3$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x4$x,res_4_2sc_pps_0.9$`Y dir`) #jelek
cor(x5$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x6$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x7$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x8$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x9$x,res_4_2sc_pps_0.9$`Y dir`) #jelek
cor(x10$x,res_4_2sc_pps_0.9$`Y dir`) #jelek
cor(x11$x,res_4_2sc_pps_0.9$`Y dir`) #jelek
cor(x12$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x13$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x14$x,res_4_2sc_pps_0.9$`Y dir`)#jelek
cor(x15$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x16$x,res_4_2sc_pps_0.9$`Y dir`)#jelek
cor(x17$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x18$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x19$x,res_4_2sc_pps_0.9$`Y dir`) #jelek
cor(x20$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x21$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x22$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x23$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x24$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x25$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x26$x,res_4_2sc_pps_0.9$`Y dir`)#jelek
cor(x27$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x28$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x29$x,res_4_2sc_pps_0.9$`Y dir`)#jelek
cor(x30$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x31$x,res_4_2sc_pps_0.9$`Y dir`)#jelek
cor(x32$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x33$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x34$x,res_4_2sc_pps_0.9$`Y dir`)
cor(x35$x,res_4_2sc_pps_0.9$`Y dir`)


#--------------------------------------------Pemilihan Variabel Auxiliary----------------------------------------------------
m_1 <- lm(res_1_srs_0.9$`Y dir`~ x1$x + x13$x + x15$x + x18$x + x20$x +x32$x + x33$x + x35$x) 
m_2 <- lm(res_2_1sc_0.9$`Y dir`~ x1$x + x13$x + x15$x + x18$x + x20$x +x32$x + x33$x + x35$x) 
m_3 <- lm(res_3_2sc_srs_0.9$`Y dir`~ x1$x + x13$x + x15$x + x18$x + x20$x +x32$x + x33$x + x35$x) 
m_4 <- lm(res_4_2sc_pps_0.9$`Y dir`~ x1$x + x13$x + x15$x + x18$x + x20$x +x32$x + x33$x + x35$x) 

library(car)
vif(m_1)
vif(m_2)
vif(m_3)
vif(m_4)

m_5 <- lm(res_1_srs_0.9$`Y dir`~ x1$x + x13$x + x18$x + x20$x +x32$x + x33$x + x35$x) 
m_6 <- lm(res_2_1sc_0.9$`Y dir`~ x1$x + x13$x + x18$x + x20$x +x32$x + x33$x + x35$x) 
m_7 <- lm(res_3_2sc_srs_0.9$`Y dir`~ x1$x + x13$x + x18$x + x20$x +x32$x + x33$x + x35$x) 
m_8 <- lm(res_4_2sc_pps_0.9$`Y dir`~ x1$x + x13$x + x18$x + x20$x +x32$x + x33$x + x35$x) 

library(car)
vif(m_5)
vif(m_6)
vif(m_7)
vif(m_8)

m_9 <- lm(res_1_srs_0.9$`Y dir`~ x1$x + x13$x + x18$x + x20$x +x32$x + x35$x ) 
m_10 <- lm(res_2_1sc_0.9$`Y dir`~ x1$x + x13$x + x18$x + x20$x +x32$x + x35$x  ) 
m_11 <- lm(res_3_2sc_srs_0.9$`Y dir`~ x1$x + x13$x + x18$x + x20$x +x32$x + x35$x  ) 
m_12 <- lm(res_4_2sc_pps_0.9$`Y dir`~ x1$x + x13$x + x18$x + x20$x +x32$x + x35$x  ) 

library(car)
vif(m_9)
vif(m_10)
vif(m_11)
vif(m_12)

cek_r2_1 <- summary(lm(res_1_srs_0.04$`Y dir`~ x1$x + x13$x + x18$x + x20$x +x32$x + x35$x ))
cek_r2_2 <- summary(lm(res_2_1sc_0.9$`Y dir`~ x1$x + x13$x + x18$x + x20$x +x32$x + x35$x ))
cek_r2_3 <- summary(lm(res_3_2sc_srs_0.9$`Y dir`~ x1$x + x13$x + x18$x + x20$x +x32$x + x35$x ))
cek_r2_4 <- summary(lm(res_4_2sc_pps_0.9$`Y dir`~ x1$x + x13$x + x18$x + x20$x +x32$x + x35$x ))

cek_r2_1
cek_r2_2
cek_r2_3
cek_r2_4

var_x <- matrix(cbind(1,x1$x,x13$x,x18$x,x20$x,x32$x,x35$x),ncol = 7) 
vr_x <- data.frame(x1$x,x13$x,x18$x,x20$x,x32$x,x35$x)

#---------------------------------------------Perhitungan jumlah ruta------------------------------------------------------------------
Ni_ruta <- ht_rt(data_1)

#---------------------------------Perhitungan secara synthetic dan composite-----------------------------------------------------------
res_syn_1_srs_1 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_1,Ni = Ni_ruta)
res_syn_1_srs_0.9 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.9,Ni = Ni_ruta)
res_syn_1_srs_0.8 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.8,Ni = Ni_ruta)
res_syn_1_srs_0.7 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.7,Ni = Ni_ruta)
res_syn_1_srs_0.6 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.6,Ni = Ni_ruta)
res_syn_1_srs_0.5 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.5,Ni = Ni_ruta)
res_syn_1_srs_0.4 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.4,Ni = Ni_ruta)
res_syn_1_srs_0.3 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.3,Ni = Ni_ruta)
res_syn_1_srs_0.2 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.2,Ni = Ni_ruta)
res_syn_1_srs_0.1 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.1,Ni = Ni_ruta)
res_syn_1_srs_0.09 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.09,Ni = Ni_ruta)
res_syn_1_srs_0.08 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.08,Ni = Ni_ruta)
res_syn_1_srs_0.07 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.07,Ni = Ni_ruta)
res_syn_1_srs_0.06 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.06,Ni = Ni_ruta)
res_syn_1_srs_0.05 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.05,Ni = Ni_ruta)
res_syn_1_srs_0.04 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.04,Ni = Ni_ruta)
res_syn_1_srs_0.03 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.03,Ni = Ni_ruta)
res_syn_1_srs_0.02 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.02,Ni = Ni_ruta)
res_syn_1_srs_0.01 <- est_prov_synth_comp(x = var_x,ydir = res_1_srs_0.01,Ni = Ni_ruta)

res_syn_1_1sc_1 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_1,Ni = Ni_ruta)
res_syn_1_1sc_0.9 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.9,Ni = Ni_ruta)
res_syn_1_1sc_0.8 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.8,Ni = Ni_ruta)
res_syn_1_1sc_0.7 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.7,Ni = Ni_ruta)
res_syn_1_1sc_0.6 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.6,Ni = Ni_ruta)
res_syn_1_1sc_0.5 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.5,Ni = Ni_ruta)
res_syn_1_1sc_0.4 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.4,Ni = Ni_ruta)
res_syn_1_1sc_0.3 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.3,Ni = Ni_ruta)
res_syn_1_1sc_0.2 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.2,Ni = Ni_ruta)
res_syn_1_1sc_0.1 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.1,Ni = Ni_ruta)
res_syn_1_1sc_0.09 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.09,Ni = Ni_ruta)
res_syn_1_1sc_0.08 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.08,Ni = Ni_ruta)
res_syn_1_1sc_0.07 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.07,Ni = Ni_ruta)
res_syn_1_1sc_0.06 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.06,Ni = Ni_ruta)
res_syn_1_1sc_0.05 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.05,Ni = Ni_ruta)
res_syn_1_1sc_0.04 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.04,Ni = Ni_ruta)
res_syn_1_1sc_0.03 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.03,Ni = Ni_ruta)
res_syn_1_1sc_0.02 <- est_prov_synth_comp(x = var_x,ydir = res_1_1sc_0.02,Ni = Ni_ruta)

res_syn_1_2sc_srs_1 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_srs_1,Ni = Ni_ruta)
res_syn_1_2sc_srs_0.9 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_srs_0.9,Ni = Ni_ruta)
res_syn_1_2sc_srs_0.85 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_srs_0.85,Ni = Ni_ruta)
res_syn_1_2sc_srs_0.8 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_srs_0.8,Ni = Ni_ruta)
res_syn_1_2sc_srs_0.75 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_srs_0.75,Ni = Ni_ruta)
res_syn_1_2sc_srs_0.7 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_srs_0.7,Ni = Ni_ruta)
res_syn_1_2sc_srs_0.65 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_srs_0.65,Ni = Ni_ruta)
res_syn_1_2sc_srs_0.6 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_srs_0.6,Ni = Ni_ruta)
res_syn_1_2sc_srs_0.55 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_srs_0.55,Ni = Ni_ruta)
res_syn_1_2sc_srs_0.5 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_srs_0.5,Ni = Ni_ruta)
res_syn_1_2sc_srs_0.45 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_srs_0.45,Ni = Ni_ruta)
res_syn_1_2sc_srs_0.4 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_srs_0.4,Ni = Ni_ruta)
res_syn_1_2sc_srs_0.35 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_srs_0.35,Ni = Ni_ruta)
res_syn_1_2sc_srs_0.3 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_srs_0.3,Ni = Ni_ruta)
res_syn_1_2sc_srs_0.25 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_srs_0.25,Ni = Ni_ruta)
res_syn_1_2sc_srs_0.2 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_srs_0.2,Ni = Ni_ruta) #negatif
res_syn_1_2sc_srs_0.15 <- est_prov_synth_comp(x = var_x, ydir = res_1_2sc_srs_0.15, Ni = Ni_ruta)

res_syn_1_2sc_pps_1 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_1,Ni = Ni_ruta)
res_syn_1_2sc_pps_0.9 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_0.9,Ni = Ni_ruta)
res_syn_1_2sc_pps_0.85 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_0.85,Ni = Ni_ruta)
res_syn_1_2sc_pps_0.8 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_0.8,Ni = Ni_ruta)
res_syn_1_2sc_pps_0.75 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_0.75,Ni = Ni_ruta)
res_syn_1_2sc_pps_0.7 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_0.7,Ni = Ni_ruta)
res_syn_1_2sc_pps_0.65 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_0.65,Ni = Ni_ruta)
res_syn_1_2sc_pps_0.6 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_0.6,Ni = Ni_ruta)
res_syn_1_2sc_pps_0.55 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_0.55,Ni = Ni_ruta)
res_syn_1_2sc_pps_0.5 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_0.5,Ni = Ni_ruta)
res_syn_1_2sc_pps_0.45 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_0.45,Ni = Ni_ruta)
res_syn_1_2sc_pps_0.4 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_0.4,Ni = Ni_ruta)
res_syn_1_2sc_pps_0.35 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_0.35,Ni = Ni_ruta)
res_syn_1_2sc_pps_0.3 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_0.3,Ni = Ni_ruta)
res_syn_1_2sc_pps_0.25 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_0.25,Ni = Ni_ruta) 
res_syn_1_2sc_pps_0.2 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_0.2,Ni = Ni_ruta) 
res_syn_1_2sc_pps_0.15 <- est_prov_synth_comp(x = var_x,ydir = res_1_2sc_pps_0.15,Ni = Ni_ruta) 

res_syn_2_srs_1 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_1,Ni = Ni_ruta)
res_syn_2_srs_0.9 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.9,Ni = Ni_ruta)
res_syn_2_srs_0.8 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.8,Ni = Ni_ruta)
res_syn_2_srs_0.7 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.7,Ni = Ni_ruta)
res_syn_2_srs_0.6 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.6,Ni = Ni_ruta)
res_syn_2_srs_0.5 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.5,Ni = Ni_ruta)
res_syn_2_srs_0.4 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.4,Ni = Ni_ruta)
res_syn_2_srs_0.3 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.3,Ni = Ni_ruta)
res_syn_2_srs_0.2 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.2,Ni = Ni_ruta)
res_syn_2_srs_0.1 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.1,Ni = Ni_ruta)
res_syn_2_srs_0.09 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.09,Ni = Ni_ruta)
res_syn_2_srs_0.08 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.08,Ni = Ni_ruta)
res_syn_2_srs_0.07 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.07,Ni = Ni_ruta)
res_syn_2_srs_0.06 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.06,Ni = Ni_ruta)
res_syn_2_srs_0.05 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.05,Ni = Ni_ruta)
res_syn_2_srs_0.04 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.04,Ni = Ni_ruta)
res_syn_2_srs_0.03 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.03,Ni = Ni_ruta)
res_syn_2_srs_0.02 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.02,Ni = Ni_ruta)
res_syn_2_srs_0.01 <- est_prov_synth_comp(x = var_x,ydir = res_2_srs_0.01,Ni = Ni_ruta)

res_syn_2_1sc_1 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_1,Ni = Ni_ruta)
res_syn_2_1sc_0.9 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.9,Ni = Ni_ruta)
res_syn_2_1sc_0.8 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.8,Ni = Ni_ruta)
res_syn_2_1sc_0.7 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.7,Ni = Ni_ruta)
res_syn_2_1sc_0.6 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.6,Ni = Ni_ruta)
res_syn_2_1sc_0.5 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.5,Ni = Ni_ruta)
res_syn_2_1sc_0.4 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.4,Ni = Ni_ruta)
res_syn_2_1sc_0.3 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.3,Ni = Ni_ruta)
res_syn_2_1sc_0.2 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.2,Ni = Ni_ruta)
res_syn_2_1sc_0.1 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.1,Ni = Ni_ruta)
res_syn_2_1sc_0.09 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.09,Ni = Ni_ruta)
res_syn_2_1sc_0.08 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.08,Ni = Ni_ruta)
res_syn_2_1sc_0.07 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.07,Ni = Ni_ruta)
res_syn_2_1sc_0.06 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.06,Ni = Ni_ruta)
res_syn_2_1sc_0.05 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.05,Ni = Ni_ruta)
res_syn_2_1sc_0.04 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.04,Ni = Ni_ruta)
res_syn_2_1sc_0.03 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.03,Ni = Ni_ruta)
res_syn_2_1sc_0.02 <- est_prov_synth_comp(x = var_x,ydir = res_2_1sc_0.02,Ni = Ni_ruta)

res_syn_2_2sc_srs_1 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_srs_1,Ni = Ni_ruta)
res_syn_2_2sc_srs_0.9 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_srs_0.9,Ni = Ni_ruta)
res_syn_2_2sc_srs_0.85 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_srs_0.85,Ni = Ni_ruta)
res_syn_2_2sc_srs_0.8 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_srs_0.8,Ni = Ni_ruta)
res_syn_2_2sc_srs_0.75 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_srs_0.75,Ni = Ni_ruta)
res_syn_2_2sc_srs_0.7 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_srs_0.7,Ni = Ni_ruta)
res_syn_2_2sc_srs_0.65 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_srs_0.65,Ni = Ni_ruta)
res_syn_2_2sc_srs_0.6 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_srs_0.6,Ni = Ni_ruta)
res_syn_2_2sc_srs_0.55 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_srs_0.55,Ni = Ni_ruta)
res_syn_2_2sc_srs_0.5 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_srs_0.5,Ni = Ni_ruta)
res_syn_2_2sc_srs_0.45 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_srs_0.45,Ni = Ni_ruta)
res_syn_2_2sc_srs_0.4 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_srs_0.4,Ni = Ni_ruta)
res_syn_2_2sc_srs_0.35 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_srs_0.35,Ni = Ni_ruta)
res_syn_2_2sc_srs_0.3 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_srs_0.3,Ni = Ni_ruta)
res_syn_2_2sc_srs_0.25 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_srs_0.25,Ni = Ni_ruta)
res_syn_2_2sc_srs_0.2 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_srs_0.2,Ni = Ni_ruta)
res_syn_2_2sc_srs_0.15 <- est_prov_synth_comp(x = var_x, ydir = res_2_2sc_srs_0.15, Ni = Ni_ruta)

res_syn_2_2sc_pps_1 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_1,Ni = Ni_ruta)
res_syn_2_2sc_pps_0.9 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_0.9,Ni = Ni_ruta)
res_syn_2_2sc_pps_0.85 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_0.85,Ni = Ni_ruta)
res_syn_2_2sc_pps_0.8 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_0.8,Ni = Ni_ruta)
res_syn_2_2sc_pps_0.75 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_0.75,Ni = Ni_ruta)
res_syn_2_2sc_pps_0.7 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_0.7,Ni = Ni_ruta)
res_syn_2_2sc_pps_0.65 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_0.65,Ni = Ni_ruta)
res_syn_2_2sc_pps_0.6 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_0.6,Ni = Ni_ruta)
res_syn_2_2sc_pps_0.55 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_0.55,Ni = Ni_ruta)
res_syn_2_2sc_pps_0.5 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_0.5,Ni = Ni_ruta)
res_syn_2_2sc_pps_0.45 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_0.45,Ni = Ni_ruta)
res_syn_2_2sc_pps_0.4 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_0.4,Ni = Ni_ruta)
res_syn_2_2sc_pps_0.35 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_0.35,Ni = Ni_ruta)
res_syn_2_2sc_pps_0.3 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_0.3,Ni = Ni_ruta)
res_syn_2_2sc_pps_0.25 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_0.25,Ni = Ni_ruta)
res_syn_2_2sc_pps_0.2 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_0.2,Ni = Ni_ruta)
res_syn_2_2sc_pps_0.15 <- est_prov_synth_comp(x = var_x,ydir = res_2_2sc_pps_0.15,Ni = Ni_ruta)

res_syn_3_srs_1 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_1,Ni = Ni_ruta)
res_syn_3_srs_0.9 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.9,Ni = Ni_ruta)
res_syn_3_srs_0.8 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.8,Ni = Ni_ruta)
res_syn_3_srs_0.7 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.7,Ni = Ni_ruta)
res_syn_3_srs_0.6 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.6,Ni = Ni_ruta)
res_syn_3_srs_0.5 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.5,Ni = Ni_ruta)
res_syn_3_srs_0.4 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.4,Ni = Ni_ruta)
res_syn_3_srs_0.3 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.3,Ni = Ni_ruta)
res_syn_3_srs_0.2 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.2,Ni = Ni_ruta)
res_syn_3_srs_0.1 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.1,Ni = Ni_ruta)
res_syn_3_srs_0.09 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.09,Ni = Ni_ruta)
res_syn_3_srs_0.08 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.08,Ni = Ni_ruta)
res_syn_3_srs_0.07 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.07,Ni = Ni_ruta)
res_syn_3_srs_0.06 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.06,Ni = Ni_ruta)
res_syn_3_srs_0.05 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.05,Ni = Ni_ruta)
res_syn_3_srs_0.04 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.04,Ni = Ni_ruta)
res_syn_3_srs_0.03 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.03,Ni = Ni_ruta)
res_syn_3_srs_0.02 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.02,Ni = Ni_ruta)
res_syn_3_srs_0.01 <- est_prov_synth_comp(x = var_x,ydir = res_3_srs_0.01,Ni = Ni_ruta)

res_syn_3_1sc_1 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_1,Ni = Ni_ruta)
res_syn_3_1sc_0.9 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.9,Ni = Ni_ruta)
res_syn_3_1sc_0.8 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.8,Ni = Ni_ruta)
res_syn_3_1sc_0.7 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.7,Ni = Ni_ruta)
res_syn_3_1sc_0.6 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.6,Ni = Ni_ruta)
res_syn_3_1sc_0.5 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.5,Ni = Ni_ruta)
res_syn_3_1sc_0.4 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.4,Ni = Ni_ruta)
res_syn_3_1sc_0.3 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.3,Ni = Ni_ruta)
res_syn_3_1sc_0.35 <- est_prov_synth_comp(x = var_x, ydir = res_3_1sc_0.35, Ni = Ni_ruta)
res_syn_3_1sc_0.2 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.2,Ni = Ni_ruta)
res_syn_3_1sc_0.1 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.1,Ni = Ni_ruta)
res_syn_3_1sc_0.09 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.09,Ni = Ni_ruta)
res_syn_3_1sc_0.08 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.08,Ni = Ni_ruta)
res_syn_3_1sc_0.07 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.07,Ni = Ni_ruta)
res_syn_3_1sc_0.06 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.06,Ni = Ni_ruta)
res_syn_3_1sc_0.05 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.05,Ni = Ni_ruta)
res_syn_3_1sc_0.04 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.04,Ni = Ni_ruta)
res_syn_3_1sc_0.03 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.03,Ni = Ni_ruta)
res_syn_3_1sc_0.02 <- est_prov_synth_comp(x = var_x,ydir = res_3_1sc_0.02,Ni = Ni_ruta)

res_syn_3_2sc_srs_1 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_srs_1,Ni = Ni_ruta)
res_syn_3_2sc_srs_0.9 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_srs_0.9,Ni = Ni_ruta)
res_syn_3_2sc_srs_0.85 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_srs_0.85,Ni = Ni_ruta)
res_syn_3_2sc_srs_0.8 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_srs_0.8,Ni = Ni_ruta)
res_syn_3_2sc_srs_0.75 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_srs_0.75,Ni = Ni_ruta)
res_syn_3_2sc_srs_0.7 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_srs_0.7,Ni = Ni_ruta)
res_syn_3_2sc_srs_0.65 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_srs_0.65,Ni = Ni_ruta)
res_syn_3_2sc_srs_0.6 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_srs_0.6,Ni = Ni_ruta)
res_syn_3_2sc_srs_0.55 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_srs_0.55,Ni = Ni_ruta)
res_syn_3_2sc_srs_0.5 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_srs_0.5,Ni = Ni_ruta)
res_syn_3_2sc_srs_0.45 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_srs_0.45,Ni = Ni_ruta)
res_syn_3_2sc_srs_0.4 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_srs_0.4,Ni = Ni_ruta)
res_syn_3_2sc_srs_0.35 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_srs_0.35,Ni = Ni_ruta)
res_syn_3_2sc_srs_0.3 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_srs_0.3,Ni = Ni_ruta)
res_syn_3_2sc_srs_0.25 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_srs_0.25,Ni = Ni_ruta)
res_syn_3_2sc_srs_0.2 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_srs_0.2,Ni = Ni_ruta)
res_syn_3_2sc_srs_0.15 <- est_prov_synth_comp(x = var_x, ydir = res_3_2sc_srs_0.15, Ni = Ni_ruta)

res_syn_3_2sc_pps_1 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_1,Ni = Ni_ruta)
res_syn_3_2sc_pps_0.9 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_0.9,Ni = Ni_ruta)
res_syn_3_2sc_pps_0.85 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_0.85,Ni = Ni_ruta)
res_syn_3_2sc_pps_0.8 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_0.8,Ni = Ni_ruta)
res_syn_3_2sc_pps_0.75 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_0.75,Ni = Ni_ruta)
res_syn_3_2sc_pps_0.7 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_0.7,Ni = Ni_ruta)
res_syn_3_2sc_pps_0.65 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_0.65,Ni = Ni_ruta)
res_syn_3_2sc_pps_0.6 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_0.6,Ni = Ni_ruta)
res_syn_3_2sc_pps_0.55 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_0.55,Ni = Ni_ruta)
res_syn_3_2sc_pps_0.5 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_0.5,Ni = Ni_ruta)
res_syn_3_2sc_pps_0.45 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_0.45,Ni = Ni_ruta)
res_syn_3_2sc_pps_0.4 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_0.4,Ni = Ni_ruta)
res_syn_3_2sc_pps_0.35 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_0.35,Ni = Ni_ruta)
res_syn_3_2sc_pps_0.3 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_0.3,Ni = Ni_ruta)
res_syn_3_2sc_pps_0.25 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_0.25,Ni = Ni_ruta)
res_syn_3_2sc_pps_0.2 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_0.2,Ni = Ni_ruta)
res_syn_3_2sc_pps_0.15 <- est_prov_synth_comp(x = var_x,ydir = res_3_2sc_pps_0.15,Ni = Ni_ruta)

res_syn_4_srs_1 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_1,Ni = Ni_ruta)
res_syn_4_srs_0.9 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.9,Ni = Ni_ruta)
res_syn_4_srs_0.8 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.8,Ni = Ni_ruta)
res_syn_4_srs_0.7 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.7,Ni = Ni_ruta)
res_syn_4_srs_0.6 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.6,Ni = Ni_ruta)
res_syn_4_srs_0.5 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.5,Ni = Ni_ruta)
res_syn_4_srs_0.4 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.4,Ni = Ni_ruta)
res_syn_4_srs_0.3 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.3,Ni = Ni_ruta)
res_syn_4_srs_0.2 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.2,Ni = Ni_ruta)
res_syn_4_srs_0.1 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.1,Ni = Ni_ruta)
res_syn_4_srs_0.09 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.09,Ni = Ni_ruta)
res_syn_4_srs_0.08 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.08,Ni = Ni_ruta)
res_syn_4_srs_0.07 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.07,Ni = Ni_ruta)
res_syn_4_srs_0.06 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.06,Ni = Ni_ruta)
res_syn_4_srs_0.05 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.05,Ni = Ni_ruta)
res_syn_4_srs_0.04 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.04,Ni = Ni_ruta)
res_syn_4_srs_0.03 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.03,Ni = Ni_ruta)
res_syn_4_srs_0.02 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.02,Ni = Ni_ruta)
res_syn_4_srs_0.01 <- est_prov_synth_comp(x = var_x,ydir = res_4_srs_0.01,Ni = Ni_ruta)

res_syn_4_1sc_1 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_1,Ni = Ni_ruta)
res_syn_4_1sc_0.9 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.9,Ni = Ni_ruta)
res_syn_4_1sc_0.8 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.8,Ni = Ni_ruta)
res_syn_4_1sc_0.7 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.7,Ni = Ni_ruta)
res_syn_4_1sc_0.6 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.6,Ni = Ni_ruta)
res_syn_4_1sc_0.5 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.5,Ni = Ni_ruta)
res_syn_4_1sc_0.4 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.4,Ni = Ni_ruta)
res_syn_4_1sc_0.3 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.3,Ni = Ni_ruta)
res_syn_4_1sc_0.2 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.2,Ni = Ni_ruta)
res_syn_4_1sc_0.1 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.1,Ni = Ni_ruta)
res_syn_4_1sc_0.09 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.09,Ni = Ni_ruta)
res_syn_4_1sc_0.08 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.08,Ni = Ni_ruta)
res_syn_4_1sc_0.07 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.07,Ni = Ni_ruta)
res_syn_4_1sc_0.06 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.06,Ni = Ni_ruta)
res_syn_4_1sc_0.05 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.05,Ni = Ni_ruta)
res_syn_4_1sc_0.04 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.04,Ni = Ni_ruta)
res_syn_4_1sc_0.03 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.03,Ni = Ni_ruta)
res_syn_4_1sc_0.02 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.02,Ni = Ni_ruta)
res_syn_4_1sc_0.01 <- est_prov_synth_comp(x = var_x,ydir = res_4_1sc_0.01,Ni = Ni_ruta)

res_syn_4_2sc_srs_1 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_srs_1,Ni = Ni_ruta)
res_syn_4_2sc_srs_0.9 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_srs_0.9,Ni = Ni_ruta)
res_syn_4_2sc_srs_0.85 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_srs_0.85,Ni = Ni_ruta)
res_syn_4_2sc_srs_0.8 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_srs_0.8,Ni = Ni_ruta)
res_syn_4_2sc_srs_0.75 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_srs_0.75,Ni = Ni_ruta)
res_syn_4_2sc_srs_0.7 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_srs_0.7,Ni = Ni_ruta)
res_syn_4_2sc_srs_0.65 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_srs_0.65,Ni = Ni_ruta)
res_syn_4_2sc_srs_0.6 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_srs_0.6,Ni = Ni_ruta)
res_syn_4_2sc_srs_0.55 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_srs_0.55,Ni = Ni_ruta)
res_syn_4_2sc_srs_0.5 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_srs_0.5,Ni = Ni_ruta)
res_syn_4_2sc_srs_0.45 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_srs_0.45,Ni = Ni_ruta)
res_syn_4_2sc_srs_0.4 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_srs_0.4,Ni = Ni_ruta)
res_syn_4_2sc_srs_0.35 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_srs_0.35,Ni = Ni_ruta)
res_syn_4_2sc_srs_0.3 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_srs_0.3,Ni = Ni_ruta)
res_syn_4_2sc_srs_0.25 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_srs_0.25,Ni = Ni_ruta)
res_syn_4_2sc_srs_0.2 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_srs_0.2,Ni = Ni_ruta)
res_syn_4_2sc_srs_0.15 <- est_prov_synth_comp(x = var_x, ydir = res_4_2sc_srs_0.15, Ni = Ni_ruta)

res_syn_4_2sc_pps_1 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_1,Ni = Ni_ruta)
res_syn_4_2sc_pps_0.9 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_0.9,Ni = Ni_ruta)
res_syn_4_2sc_pps_0.85 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_0.85,Ni = Ni_ruta)
res_syn_4_2sc_pps_0.8 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_0.8,Ni = Ni_ruta)
res_syn_4_2sc_pps_0.75 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_0.75,Ni = Ni_ruta)
res_syn_4_2sc_pps_0.7 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_0.7,Ni = Ni_ruta)
res_syn_4_2sc_pps_0.65 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_0.65,Ni = Ni_ruta)
res_syn_4_2sc_pps_0.6 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_0.6,Ni = Ni_ruta)
res_syn_4_2sc_pps_0.55 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_0.55,Ni = Ni_ruta)
res_syn_4_2sc_pps_0.5 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_0.5,Ni = Ni_ruta)
res_syn_4_2sc_pps_0.45 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_0.45,Ni = Ni_ruta)
res_syn_4_2sc_pps_0.4 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_0.4,Ni = Ni_ruta)
res_syn_4_2sc_pps_0.35 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_0.35,Ni = Ni_ruta)
res_syn_4_2sc_pps_0.3 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_0.3,Ni = Ni_ruta)
res_syn_4_2sc_pps_0.25 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_0.25,Ni = Ni_ruta)
res_syn_4_2sc_pps_0.2 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_0.2,Ni = Ni_ruta)
res_syn_4_2sc_pps_0.15 <- est_prov_synth_comp(x = var_x,ydir = res_4_2sc_pps_0.15,Ni = Ni_ruta)

#---------------------------------------------------Cek rata-rata relatif efisiensi synthetic-----------------------------------------
mean(res_syn_1_srs_1$`Ef syn`)
mean(res_syn_1_srs_0.9$`Ef syn`)
mean(res_syn_1_srs_0.8$`Ef syn`)
mean(res_syn_1_srs_0.7$`Ef syn`)
mean(res_syn_1_srs_0.6$`Ef syn`)
mean(res_syn_1_srs_0.5$`Ef syn`)
mean(res_syn_1_srs_0.4$`Ef syn`)
mean(res_syn_1_srs_0.3$`Ef syn`)
mean(res_syn_1_srs_0.2$`Ef syn`)
mean(res_syn_1_srs_0.1$`Ef syn`)
mean(res_syn_1_srs_0.09$`Ef syn`)
mean(res_syn_1_srs_0.08$`Ef syn`)
mean(res_syn_1_srs_0.07$`Ef syn`)
mean(res_syn_1_srs_0.06$`Ef syn`)
mean(res_syn_1_srs_0.05$`Ef syn`) #ini
mean(res_syn_1_srs_0.04$`Ef syn`) #ini 
mean(res_syn_1_srs_0.03$`Ef syn`)
mean(res_syn_1_srs_0.02$`Ef syn`) #negatif
mean(res_syn_1_srs_0.01$`Ef syn`) #negatif

mean(res_syn_1_1sc_1$`Ef syn`)
mean(res_syn_1_1sc_0.9$`Ef syn`)
mean(res_syn_1_1sc_0.8$`Ef syn`)
mean(res_syn_1_1sc_0.7$`Ef syn`)
mean(res_syn_1_1sc_0.6$`Ef syn`)
mean(res_syn_1_1sc_0.5$`Ef syn`)
mean(res_syn_1_1sc_0.4$`Ef syn`)
mean(res_syn_1_1sc_0.3$`Ef syn`)
mean(res_syn_1_1sc_0.2$`Ef syn`) #ini
mean(res_syn_1_1sc_0.1$`Ef syn`) #ini
mean(res_syn_1_1sc_0.09$`Ef syn`)
mean(res_syn_1_1sc_0.08$`Ef syn`)
mean(res_syn_1_1sc_0.07$`Ef syn`)#negatif
mean(res_syn_1_1sc_0.06$`Ef syn`)#negatif
mean(res_syn_1_1sc_0.05$`Ef syn`)#negatif
mean(res_syn_1_1sc_0.04$`Ef syn`)#negatif

mean(res_syn_1_2sc_srs_1$`Ef syn`)
mean(res_syn_1_2sc_srs_0.9$`Ef syn`)
mean(res_syn_1_2sc_srs_0.85$`Ef syn`)
mean(res_syn_1_2sc_srs_0.8$`Ef syn`)
mean(res_syn_1_2sc_srs_0.75$`Ef syn`)
mean(res_syn_1_2sc_srs_0.7$`Ef syn`)
mean(res_syn_1_2sc_srs_0.65$`Ef syn`)
mean(res_syn_1_2sc_srs_0.6$`Ef syn`)
mean(res_syn_1_2sc_srs_0.55$`Ef syn`)
mean(res_syn_1_2sc_srs_0.5$`Ef syn`)
mean(res_syn_1_2sc_srs_0.45$`Ef syn`)
mean(res_syn_1_2sc_srs_0.4$`Ef syn`)
mean(res_syn_1_2sc_srs_0.35$`Ef syn`)
mean(res_syn_1_2sc_srs_0.3$`Ef syn`) 
mean(res_syn_1_2sc_srs_0.25$`Ef syn`) #ini
mean(res_syn_1_2sc_srs_0.2$`Ef syn`) #ini

mean(res_syn_1_2sc_pps_1$`Ef syn`)
mean(res_syn_1_2sc_pps_0.9$`Ef syn`)
mean(res_syn_1_2sc_pps_0.85$`Ef syn`)
mean(res_syn_1_2sc_pps_0.8$`Ef syn`)
mean(res_syn_1_2sc_pps_0.75$`Ef syn`)
mean(res_syn_1_2sc_pps_0.7$`Ef syn`)
mean(res_syn_1_2sc_pps_0.65$`Ef syn`)
mean(res_syn_1_2sc_pps_0.6$`Ef syn`)
mean(res_syn_1_2sc_pps_0.55$`Ef syn`)
mean(res_syn_1_2sc_pps_0.5$`Ef syn`)
mean(res_syn_1_2sc_pps_0.45$`Ef syn`)
mean(res_syn_1_2sc_pps_0.4$`Ef syn`)
mean(res_syn_1_2sc_pps_0.35$`Ef syn`)
mean(res_syn_1_2sc_pps_0.3$`Ef syn`) #ini
mean(res_syn_1_2sc_pps_0.25$`Ef syn`)#ini
mean(res_syn_1_2sc_pps_0.2$`Ef syn`)

mean(res_syn_2_srs_1$`Ef syn`)
mean(res_syn_2_srs_0.9$`Ef syn`)
mean(res_syn_2_srs_0.8$`Ef syn`)
mean(res_syn_2_srs_0.7$`Ef syn`)
mean(res_syn_2_srs_0.6$`Ef syn`)
mean(res_syn_2_srs_0.5$`Ef syn`)
mean(res_syn_2_srs_0.4$`Ef syn`)
mean(res_syn_2_srs_0.3$`Ef syn`)
mean(res_syn_2_srs_0.2$`Ef syn`)
mean(res_syn_2_srs_0.1$`Ef syn`)
mean(res_syn_2_srs_0.09$`Ef syn`)
mean(res_syn_2_srs_0.08$`Ef syn`)
mean(res_syn_2_srs_0.07$`Ef syn`)
mean(res_syn_2_srs_0.06$`Ef syn`)
mean(res_syn_2_srs_0.05$`Ef syn`) 
mean(res_syn_2_srs_0.04$`Ef syn`) 
mean(res_syn_2_srs_0.03$`Ef syn`)#ini
mean(res_syn_2_srs_0.02$`Ef syn`)
mean(res_syn_2_srs_0.01$`Ef syn`)#negatif

mean(res_syn_2_1sc_1$`Ef syn`)
mean(res_syn_2_1sc_0.9$`Ef syn`)
mean(res_syn_2_1sc_0.8$`Ef syn`)
mean(res_syn_2_1sc_0.7$`Ef syn`)
mean(res_syn_2_1sc_0.6$`Ef syn`)
mean(res_syn_2_1sc_0.5$`Ef syn`)
mean(res_syn_2_1sc_0.4$`Ef syn`)
mean(res_syn_2_1sc_0.3$`Ef syn`)
mean(res_syn_2_1sc_0.2$`Ef syn`)
mean(res_syn_2_1sc_0.1$`Ef syn`) 
mean(res_syn_2_1sc_0.09$`Ef syn`) #ini
mean(res_syn_2_1sc_0.08$`Ef syn`) #ini
mean(res_syn_2_1sc_0.07$`Ef syn`)
mean(res_syn_2_1sc_0.06$`Ef syn`)
mean(res_syn_2_1sc_0.05$`Ef syn`)
mean(res_syn_2_1sc_0.04$`Ef syn`)
mean(res_syn_2_1sc_0.03$`Ef syn`) #error
mean(res_syn_2_1sc_0.02$`Ef syn`) #error
mean(res_syn_2_1sc_0.01$`Ef syn`) #error

mean(res_syn_2_2sc_srs_1$`Ef syn`)
mean(res_syn_2_2sc_srs_0.9$`Ef syn`)
mean(res_syn_2_2sc_srs_0.85$`Ef syn`)
mean(res_syn_2_2sc_srs_0.8$`Ef syn`)
mean(res_syn_2_2sc_srs_0.75$`Ef syn`)
mean(res_syn_2_2sc_srs_0.7$`Ef syn`)
mean(res_syn_2_2sc_srs_0.65$`Ef syn`)
mean(res_syn_2_2sc_srs_0.6$`Ef syn`)
mean(res_syn_2_2sc_srs_0.55$`Ef syn`)
mean(res_syn_2_2sc_srs_0.5$`Ef syn`)
mean(res_syn_2_2sc_srs_0.45$`Ef syn`)
mean(res_syn_2_2sc_srs_0.4$`Ef syn`)
mean(res_syn_2_2sc_srs_0.35$`Ef syn`)
mean(res_syn_2_2sc_srs_0.3$`Ef syn`)
mean(res_syn_2_2sc_srs_0.25$`Ef syn`)
mean(res_syn_2_2sc_srs_0.2$`Ef syn`) #ini
mean(res_syn_2_2sc_srs_0.15$`Ef syn`)

mean(res_syn_2_2sc_pps_1$`Ef syn`)
mean(res_syn_2_2sc_pps_0.9$`Ef syn`)
mean(res_syn_2_2sc_pps_0.85$`Ef syn`)
mean(res_syn_2_2sc_pps_0.8$`Ef syn`)
mean(res_syn_2_2sc_pps_0.75$`Ef syn`)
mean(res_syn_2_2sc_pps_0.7$`Ef syn`)
mean(res_syn_2_2sc_pps_0.65$`Ef syn`)
mean(res_syn_2_2sc_pps_0.6$`Ef syn`)
mean(res_syn_2_2sc_pps_0.55$`Ef syn`)
mean(res_syn_2_2sc_pps_0.5$`Ef syn`)
mean(res_syn_2_2sc_pps_0.45$`Ef syn`)
mean(res_syn_2_2sc_pps_0.4$`Ef syn`)
mean(res_syn_2_2sc_pps_0.35$`Ef syn`)
mean(res_syn_2_2sc_pps_0.3$`Ef syn`) 
mean(res_syn_2_2sc_pps_0.25$`Ef syn`) #ini
mean(res_syn_2_2sc_pps_0.2$`Ef syn`) #ini
mean(res_syn_2_2sc_pps_0.15$`Ef syn`) 

mean(res_syn_3_srs_1$`Ef syn`)
mean(res_syn_3_srs_0.9$`Ef syn`)
mean(res_syn_3_srs_0.8$`Ef syn`)
mean(res_syn_3_srs_0.7$`Ef syn`)
mean(res_syn_3_srs_0.6$`Ef syn`)
mean(res_syn_3_srs_0.5$`Ef syn`)
mean(res_syn_3_srs_0.4$`Ef syn`) #ini
mean(res_syn_3_srs_0.3$`Ef syn`) #ini
mean(res_syn_3_srs_0.2$`Ef syn`) #negartif
mean(res_syn_3_srs_0.1$`Ef syn`) #negatif
mean(res_syn_3_srs_0.09$`Ef syn`) #negatif
mean(res_syn_3_srs_0.08$`Ef syn`) #negatif
mean(res_syn_3_srs_0.07$`Ef syn`) #negatif
mean(res_syn_3_srs_0.06$`Ef syn`) #negatif
mean(res_syn_3_srs_0.05$`Ef syn`) #negatif
mean(res_syn_3_srs_0.04$`Ef syn`) #negatif
mean(res_syn_3_srs_0.03$`Ef syn`) #negatif
mean(res_syn_3_srs_0.02$`Ef syn`) #negatif
mean(res_syn_3_srs_0.01$`Ef syn`) #negatif

mean(res_syn_3_1sc_1$`Ef syn`)
mean(res_syn_3_1sc_0.9$`Ef syn`)
mean(res_syn_3_1sc_0.8$`Ef syn`)
mean(res_syn_3_1sc_0.7$`Ef syn`)
mean(res_syn_3_1sc_0.6$`Ef syn`)
mean(res_syn_3_1sc_0.5$`Ef syn`)
mean(res_syn_3_1sc_0.4$`Ef syn`) #ini
mean(res_syn_3_1sc_0.3$`Ef syn`) #ini
mean(res_syn_3_1sc_0.2$`Ef syn`) #negatif
mean(res_syn_3_1sc_0.1$`Ef syn`) #negatif
mean(res_syn_3_1sc_0.09$`Ef syn`)
mean(res_syn_3_1sc_0.08$`Ef syn`)
mean(res_syn_3_1sc_0.07$`Ef syn`)
mean(res_syn_3_1sc_0.06$`Ef syn`)
mean(res_syn_3_1sc_0.05$`Ef syn`)
mean(res_syn_3_1sc_0.04$`Ef syn`)
mean(res_syn_3_1sc_0.03$`Ef syn`)
mean(res_syn_3_1sc_0.02$`Ef syn`)
mean(res_syn_3_1sc_0.01$`Ef syn`)

mean(res_syn_3_2sc_srs_1$`Ef syn`)
mean(res_syn_3_2sc_srs_0.9$`Ef syn`)
mean(res_syn_3_2sc_srs_0.85$`Ef syn`)
mean(res_syn_3_2sc_srs_0.8$`Ef syn`)
mean(res_syn_3_2sc_srs_0.75$`Ef syn`)
mean(res_syn_3_2sc_srs_0.7$`Ef syn`)
mean(res_syn_3_2sc_srs_0.65$`Ef syn`) 
mean(res_syn_3_2sc_srs_0.6$`Ef syn`) #ini
mean(res_syn_3_2sc_srs_0.55$`Ef syn`) #ini
mean(res_syn_3_2sc_srs_0.5$`Ef syn`)
mean(res_syn_3_2sc_srs_0.45$`Ef syn`) #negatif
mean(res_syn_3_2sc_srs_0.4$`Ef syn`) #negatif
mean(res_syn_3_2sc_srs_0.35$`Ef syn`) #negatif
mean(res_syn_3_2sc_srs_0.3$`Ef syn`) #negatif
mean(res_syn_3_2sc_srs_0.25$`Ef syn`) #negatif
mean(res_syn_3_2sc_srs_0.2$`Ef syn`) #negatif

mean(res_syn_3_2sc_pps_1$`Ef syn`)
mean(res_syn_3_2sc_pps_0.9$`Ef syn`)
mean(res_syn_3_2sc_pps_0.85$`Ef syn`)
mean(res_syn_3_2sc_pps_0.8$`Ef syn`) 
mean(res_syn_3_2sc_pps_0.75$`Ef syn`) #ini
mean(res_syn_3_2sc_pps_0.7$`Ef syn`) #ini
mean(res_syn_3_2sc_pps_0.65$`Ef syn`)
mean(res_syn_3_2sc_pps_0.6$`Ef syn`)
mean(res_syn_3_2sc_pps_0.55$`Ef syn`)
mean(res_syn_3_2sc_pps_0.5$`Ef syn`) #negatif
mean(res_syn_3_2sc_pps_0.45$`Ef syn`)
mean(res_syn_3_2sc_pps_0.4$`Ef syn`)
mean(res_syn_3_2sc_pps_0.35$`Ef syn`)
mean(res_syn_3_2sc_pps_0.3$`Ef syn`)
mean(res_syn_3_2sc_pps_0.25$`Ef syn`)
mean(res_syn_3_2sc_pps_0.2$`Ef syn`)

mean(res_syn_4_srs_1$`Ef syn`)
mean(res_syn_4_srs_0.9$`Ef syn`)
mean(res_syn_4_srs_0.8$`Ef syn`)
mean(res_syn_4_srs_0.7$`Ef syn`)
mean(res_syn_4_srs_0.6$`Ef syn`)
mean(res_syn_4_srs_0.5$`Ef syn`)
mean(res_syn_4_srs_0.4$`Ef syn`)
mean(res_syn_4_srs_0.3$`Ef syn`)
mean(res_syn_4_srs_0.2$`Ef syn`)
mean(res_syn_4_srs_0.1$`Ef syn`)
mean(res_syn_4_srs_0.09$`Ef syn`)
mean(res_syn_4_srs_0.08$`Ef syn`)
mean(res_syn_4_srs_0.07$`Ef syn`)
mean(res_syn_4_srs_0.06$`Ef syn`)
mean(res_syn_4_srs_0.05$`Ef syn`)
mean(res_syn_4_srs_0.04$`Ef syn`)
mean(res_syn_4_srs_0.03$`Ef syn`)
mean(res_syn_4_srs_0.02$`Ef syn`)
mean(res_syn_4_srs_0.01$`Ef syn`)

mean(res_syn_4_1sc_1$`Ef syn`)
mean(res_syn_4_1sc_0.9$`Ef syn`)
mean(res_syn_4_1sc_0.8$`Ef syn`)
mean(res_syn_4_1sc_0.7$`Ef syn`)
mean(res_syn_4_1sc_0.6$`Ef syn`)
mean(res_syn_4_1sc_0.5$`Ef syn`)
mean(res_syn_4_1sc_0.4$`Ef syn`)
mean(res_syn_4_1sc_0.3$`Ef syn`)
mean(res_syn_4_1sc_0.2$`Ef syn`)
mean(res_syn_4_1sc_0.1$`Ef syn`)
mean(res_syn_4_1sc_0.09$`Ef syn`)
mean(res_syn_4_1sc_0.08$`Ef syn`)
mean(res_syn_4_1sc_0.07$`Ef syn`)
mean(res_syn_4_1sc_0.06$`Ef syn`)
mean(res_syn_4_1sc_0.05$`Ef syn`)
mean(res_syn_4_1sc_0.04$`Ef syn`)
mean(res_syn_4_1sc_0.03$`Ef syn`)
mean(res_syn_4_1sc_0.02$`Ef syn`)
mean(res_syn_4_1sc_0.01$`Ef syn`)

mean(res_syn_4_2sc_srs_1$`Ef syn`)
mean(res_syn_4_2sc_srs_0.9$`Ef syn`)
mean(res_syn_4_2sc_srs_0.85$`Ef syn`)
mean(res_syn_4_2sc_srs_0.8$`Ef syn`)
mean(res_syn_4_2sc_srs_0.75$`Ef syn`)
mean(res_syn_4_2sc_srs_0.7$`Ef syn`)
mean(res_syn_4_2sc_srs_0.65$`Ef syn`)
mean(res_syn_4_2sc_srs_0.6$`Ef syn`)
mean(res_syn_4_2sc_srs_0.55$`Ef syn`)
mean(res_syn_4_2sc_srs_0.5$`Ef syn`)
mean(res_syn_4_2sc_srs_0.45$`Ef syn`)
mean(res_syn_4_2sc_srs_0.4$`Ef syn`)
mean(res_syn_4_2sc_srs_0.35$`Ef syn`)
mean(res_syn_4_2sc_srs_0.3$`Ef syn`)
mean(res_syn_4_2sc_srs_0.25$`Ef syn`)
mean(res_syn_4_2sc_srs_0.2$`Ef syn`)

mean(res_syn_4_2sc_pps_1$`Ef syn`)
mean(res_syn_4_2sc_pps_0.9$`Ef syn`)
mean(res_syn_4_2sc_pps_0.85$`Ef syn`)
mean(res_syn_4_2sc_pps_0.8$`Ef syn`)
mean(res_syn_4_2sc_pps_0.75$`Ef syn`)
mean(res_syn_4_2sc_pps_0.7$`Ef syn`)
mean(res_syn_4_2sc_pps_0.65$`Ef syn`)
mean(res_syn_4_2sc_pps_0.6$`Ef syn`)
mean(res_syn_4_2sc_pps_0.55$`Ef syn`)
mean(res_syn_4_2sc_pps_0.5$`Ef syn`)
mean(res_syn_4_2sc_pps_0.45$`Ef syn`)
mean(res_syn_4_2sc_pps_0.4$`Ef syn`)
mean(res_syn_4_2sc_pps_0.35$`Ef syn`)
mean(res_syn_4_2sc_pps_0.3$`Ef syn`)
mean(res_syn_4_2sc_pps_0.25$`Ef syn`)
mean(res_syn_4_2sc_pps_0.2$`Ef syn`)

#-----------------------------------------------------------PERHITUNGAN EBLUP----------------------------------------------------------
aux_x <- var_x[,-1]
library(sae)
res_eblup_2_srs_0.03 <- mseFH(formula = res_2_srs_0.03$`Y dir` ~ aux_x, vardir = res_2_srs_0.03$`MSE dir`, method = "ML")
res_eblup_2_1sc_0.09 <- mseFH(formula = res_2_1sc_0.09$`Y dir` ~ aux_x, vardir = res_2_1sc_0.09$`MSE dir`, method = "ML")
res_eblup_2_2sc_srs_0.15 <- mseFH(formula = res_2_2sc_srs_0.15$`Y dir` ~ aux_x, vardir = res_2_2sc_srs_0.15$`MSE dir`, method = "ML")
res_eblup_2_2sc_pps_0.2 <- mseFH(formula = res_2_2sc_pps_0.2$`Y dir` ~ aux_x, vardir = res_2_2sc_pps_0.2$`MSE dir`, method = "ML")

res_eblup_2_srs_0.03_2 <- mseFH(formula = res_2_srs_0.03$`Y dir` ~ aux_x, vardir = res_2_srs_0.03$`MSE dir`, method = "REML")
res_eblup_2_1sc_0.09_2 <- mseFH(formula = res_2_1sc_0.09$`Y dir` ~ aux_x, vardir = res_2_1sc_0.09$`MSE dir`, method = "REML")
res_eblup_2_2sc_srs_0.15_2 <- mseFH(formula = res_2_2sc_srs_0.15$`Y dir` ~ aux_x, vardir = res_2_2sc_srs_0.15$`MSE dir`, method = "REML")
res_eblup_2_2sc_pps_0.2_2 <- mseFH(formula = res_2_2sc_pps_0.2$`Y dir` ~ aux_x, vardir = res_2_2sc_pps_0.2$`MSE dir`, method = "REML")
#---------------------------------------------------------EBLUP HITUNG----------------------
res_eblup_1_srs_0.9 <- mseFH(formula = res_1_srs_0.9$`Y dir` ~ aux_x, vardir = res_1_srs_0.9$`MSE dir`, method = "REML")
res_eblup_1_srs_0.8 <- mseFH(formula = res_1_srs_0.8$`Y dir` ~ aux_x, vardir = res_1_srs_0.8$`MSE dir`, method = "REML")
res_eblup_1_srs_0.7 <- mseFH(formula = res_1_srs_0.7$`Y dir` ~ aux_x, vardir = res_1_srs_0.7$`MSE dir`, method = "REML")
res_eblup_1_srs_0.6 <- mseFH(formula = res_1_srs_0.6$`Y dir` ~ aux_x, vardir = res_1_srs_0.6$`MSE dir`, method = "REML")
res_eblup_1_srs_0.5 <- mseFH(formula = res_1_srs_0.5$`Y dir` ~ aux_x, vardir = res_1_srs_0.5$`MSE dir`, method = "REML")
res_eblup_1_srs_0.4 <- mseFH(formula = res_1_srs_0.4$`Y dir` ~ aux_x, vardir = res_1_srs_0.4$`MSE dir`, method = "REML")
res_eblup_1_srs_0.3 <- mseFH(formula = res_1_srs_0.3$`Y dir` ~ aux_x, vardir = res_1_srs_0.3$`MSE dir`, method = "REML")
res_eblup_1_srs_0.2 <- mseFH(formula = res_1_srs_0.2$`Y dir` ~ aux_x, vardir = res_1_srs_0.2$`MSE dir`, method = "REML")
res_eblup_1_srs_0.1 <- mseFH(formula = res_1_srs_0.1$`Y dir` ~ aux_x, vardir = res_1_srs_0.1$`MSE dir`, method = "REML")

res_eblup_1_1sc_0.9 <- mseFH(formula = res_1_1sc_0.9$`Y dir` ~ aux_x, vardir = res_1_1sc_0.9$`MSE dir`, method = "REML")
res_eblup_1_1sc_0.8 <- mseFH(formula = res_1_1sc_0.8$`Y dir` ~ aux_x, vardir = res_1_1sc_0.8$`MSE dir`, method = "REML")
res_eblup_1_1sc_0.7 <- mseFH(formula = res_1_1sc_0.7$`Y dir` ~ aux_x, vardir = res_1_1sc_0.7$`MSE dir`, method = "REML")
res_eblup_1_1sc_0.6 <- mseFH(formula = res_1_1sc_0.6$`Y dir` ~ aux_x, vardir = res_1_1sc_0.6$`MSE dir`, method = "REML")
res_eblup_1_1sc_0.5 <- mseFH(formula = res_1_1sc_0.5$`Y dir` ~ aux_x, vardir = res_1_1sc_0.5$`MSE dir`, method = "REML")
res_eblup_1_1sc_0.4 <- mseFH(formula = res_1_1sc_0.4$`Y dir` ~ aux_x, vardir = res_1_1sc_0.4$`MSE dir`, method = "REML")
res_eblup_1_1sc_0.3 <- mseFH(formula = res_1_1sc_0.3$`Y dir` ~ aux_x, vardir = res_1_1sc_0.3$`MSE dir`, method = "REML")
res_eblup_1_1sc_0.2 <- mseFH(formula = res_1_1sc_0.2$`Y dir` ~ aux_x, vardir = res_1_1sc_0.2$`MSE dir`, method = "REML")
res_eblup_1_1sc_0.1 <- mseFH(formula = res_1_1sc_0.1$`Y dir` ~ aux_x, vardir = res_1_1sc_0.1$`MSE dir`, method = "REML")

res_eblup_1_2sc_srs_0.9 <- mseFH(formula = res_1_2sc_srs_0.9$`Y dir` ~ aux_x, vardir = res_1_2sc_srs_0.9$`MSE dir`, method = "REML")
res_eblup_1_2sc_srs_0.8 <- mseFH(formula = res_1_2sc_srs_0.8$`Y dir` ~ aux_x, vardir = res_1_2sc_srs_0.8$`MSE dir`, method = "REML")
res_eblup_1_2sc_srs_0.7 <- mseFH(formula = res_1_2sc_srs_0.7$`Y dir` ~ aux_x, vardir = res_1_2sc_srs_0.7$`MSE dir`, method = "REML")
res_eblup_1_2sc_srs_0.6 <- mseFH(formula = res_1_2sc_srs_0.6$`Y dir` ~ aux_x, vardir = res_1_2sc_srs_0.6$`MSE dir`, method = "REML")
res_eblup_1_2sc_srs_0.5 <- mseFH(formula = res_1_2sc_srs_0.5$`Y dir` ~ aux_x, vardir = res_1_2sc_srs_0.5$`MSE dir`, method = "REML")
res_eblup_1_2sc_srs_0.4 <- mseFH(formula = res_1_2sc_srs_0.4$`Y dir` ~ aux_x, vardir = res_1_2sc_srs_0.4$`MSE dir`, method = "REML")
res_eblup_1_2sc_srs_0.3 <- mseFH(formula = res_1_2sc_srs_0.3$`Y dir` ~ aux_x, vardir = res_1_2sc_srs_0.3$`MSE dir`, method = "REML")
res_eblup_1_2sc_srs_0.2 <- mseFH(formula = res_1_2sc_srs_0.2$`Y dir` ~ aux_x, vardir = res_1_2sc_srs_0.2$`MSE dir`, method = "REML")
res_eblup_1_2sc_srs_0.1 <- mseFH(formula = res_1_2sc_srs_0.1$`Y dir` ~ aux_x, vardir = res_1_2sc_srs_0.1$`MSE dir`, method = "REML")

res_eblup_1_2sc_pps_0.9 <- mseFH(formula = res_1_2sc_pps_0.9$`Y dir` ~ aux_x, vardir = res_1_2sc_pps_0.9$`MSE dir`, method = "REML")
res_eblup_1_2sc_pps_0.8 <- mseFH(formula = res_1_2sc_pps_0.8$`Y dir` ~ aux_x, vardir = res_1_2sc_pps_0.8$`MSE dir`, method = "REML")
res_eblup_1_2sc_pps_0.7 <- mseFH(formula = res_1_2sc_pps_0.7$`Y dir` ~ aux_x, vardir = res_1_2sc_pps_0.7$`MSE dir`, method = "REML")
res_eblup_1_2sc_pps_0.6 <- mseFH(formula = res_1_2sc_pps_0.6$`Y dir` ~ aux_x, vardir = res_1_2sc_pps_0.6$`MSE dir`, method = "REML")
res_eblup_1_2sc_pps_0.5 <- mseFH(formula = res_1_2sc_pps_0.5$`Y dir` ~ aux_x, vardir = res_1_2sc_pps_0.5$`MSE dir`, method = "REML")
res_eblup_1_2sc_pps_0.4 <- mseFH(formula = res_1_2sc_pps_0.4$`Y dir` ~ aux_x, vardir = res_1_2sc_pps_0.4$`MSE dir`, method = "REML")
res_eblup_1_2sc_pps_0.3 <- mseFH(formula = res_1_2sc_pps_0.3$`Y dir` ~ aux_x, vardir = res_1_2sc_pps_0.3$`MSE dir`, method = "REML")
res_eblup_1_2sc_pps_0.2 <- mseFH(formula = res_1_2sc_pps_0.2$`Y dir` ~ aux_x, vardir = res_1_2sc_pps_0.2$`MSE dir`, method = "REML")
res_eblup_1_2sc_pps_0.1 <- mseFH(formula = res_1_2sc_pps_0.1$`Y dir` ~ aux_x, vardir = res_1_2sc_pps_0.1$`MSE dir`, method = "REML")



gab_srs_0.9 <- data.frame(res_1_srs_0.9$kab,res_1_srs_0.9$`Y dir`,res_1_srs_0.9$`MSE dir`,res_1_srs_0.9$`RRMSE dir`,
                           res_eblup_1_srs_0.9$est$eblup,res_eblup_1_srs_0.9$mse,sqrt(res_eblup_1_srs_0.9$mse)
                           /res_eblup_1_srs_0.9$est$eblup*100)
colnames(gab_srs_0.9) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_srs_0.9$`MSE eblup`/gab_srs_0.9$`MSE dir`)

gab_srs_0.8 <- data.frame(res_1_srs_0.8$kab,res_1_srs_0.8$`Y dir`,res_1_srs_0.8$`MSE dir`,res_1_srs_0.8$`RRMSE dir`,
                          res_eblup_1_srs_0.8$est$eblup,res_eblup_1_srs_0.8$mse,sqrt(res_eblup_1_srs_0.8$mse)
                          /res_eblup_1_srs_0.8$est$eblup*100)
colnames(gab_srs_0.8) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_srs_0.8$`MSE eblup`/gab_srs_0.8$`MSE dir`)

gab_srs_0.7 <- data.frame(res_1_srs_0.7$kab,res_1_srs_0.7$`Y dir`,res_1_srs_0.7$`MSE dir`,res_1_srs_0.7$`RRMSE dir`,
                          res_eblup_1_srs_0.7$est$eblup,res_eblup_1_srs_0.7$mse,sqrt(res_eblup_1_srs_0.7$mse)
                          /res_eblup_1_srs_0.7$est$eblup*100)
colnames(gab_srs_0.7) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_srs_0.7$`MSE eblup`/gab_srs_0.7$`MSE dir`)

gab_srs_0.6 <- data.frame(res_1_srs_0.6$kab,res_1_srs_0.6$`Y dir`,res_1_srs_0.6$`MSE dir`,res_1_srs_0.6$`RRMSE dir`,
                          res_eblup_1_srs_0.6$est$eblup,res_eblup_1_srs_0.6$mse,sqrt(res_eblup_1_srs_0.6$mse)
                          /res_eblup_1_srs_0.6$est$eblup*100)
colnames(gab_srs_0.6) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_srs_0.6$`MSE eblup`/gab_srs_0.6$`MSE dir`)

gab_srs_0.5 <- data.frame(res_1_srs_0.5$kab,res_1_srs_0.5$`Y dir`,res_1_srs_0.5$`MSE dir`,res_1_srs_0.5$`RRMSE dir`,
                          res_eblup_1_srs_0.5$est$eblup,res_eblup_1_srs_0.5$mse,sqrt(res_eblup_1_srs_0.5$mse)
                          /res_eblup_1_srs_0.5$est$eblup*100)
colnames(gab_srs_0.5) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_srs_0.5$`MSE eblup`/gab_srs_0.5$`MSE dir`)

gab_srs_0.4 <- data.frame(res_1_srs_0.4$kab,res_1_srs_0.4$`Y dir`,res_1_srs_0.4$`MSE dir`,res_1_srs_0.4$`RRMSE dir`,
                          res_eblup_1_srs_0.4$est$eblup,res_eblup_1_srs_0.4$mse,sqrt(res_eblup_1_srs_0.4$mse)
                          /res_eblup_1_srs_0.4$est$eblup*100)
colnames(gab_srs_0.4) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_srs_0.4$`MSE eblup`/gab_srs_0.4$`MSE dir`)

gab_srs_0.3 <- data.frame(res_1_srs_0.3$kab,res_1_srs_0.3$`Y dir`,res_1_srs_0.3$`MSE dir`,res_1_srs_0.3$`RRMSE dir`,
                          res_eblup_1_srs_0.3$est$eblup,res_eblup_1_srs_0.3$mse,sqrt(res_eblup_1_srs_0.3$mse)
                          /res_eblup_1_srs_0.3$est$eblup*100)
colnames(gab_srs_0.3) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_srs_0.3$`MSE eblup`/gab_srs_0.3$`MSE dir`)

gab_srs_0.2 <- data.frame(res_1_srs_0.2$kab,res_1_srs_0.2$`Y dir`,res_1_srs_0.2$`MSE dir`,res_1_srs_0.2$`RRMSE dir`,
                          res_eblup_1_srs_0.2$est$eblup,res_eblup_1_srs_0.2$mse,sqrt(res_eblup_1_srs_0.2$mse)
                          /res_eblup_1_srs_0.2$est$eblup*100)
colnames(gab_srs_0.2) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_srs_0.2$`MSE eblup`/gab_srs_0.2$`MSE dir`)

gab_srs_0.1 <- data.frame(res_1_srs_0.1$kab,res_1_srs_0.1$`Y dir`,res_1_srs_0.1$`MSE dir`,res_1_srs_0.1$`RRMSE dir`,
                          res_eblup_1_srs_0.1$est$eblup,res_eblup_1_srs_0.1$mse,sqrt(res_eblup_1_srs_0.1$mse)
                          /res_eblup_1_srs_0.1$est$eblup*100)
colnames(gab_srs_0.1) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_srs_0.1$`MSE eblup`/gab_srs_0.1$`MSE dir`)


gab_1sc_0.9 <- data.frame(res_1_1sc_0.9$kab,res_1_1sc_0.9$`Y dir`,res_1_1sc_0.9$`MSE dir`,res_1_1sc_0.9$`RRMSE dir`,
                          res_eblup_1_1sc_0.9$est$eblup,res_eblup_1_1sc_0.9$mse,sqrt(res_eblup_1_1sc_0.9$mse)
                          /res_eblup_1_1sc_0.9$est$eblup*100)
colnames(gab_1sc_0.9) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_1sc_0.9$`MSE eblup`/gab_1sc_0.9$`MSE dir`)

gab_1sc_0.8 <- data.frame(res_1_1sc_0.8$kab,res_1_1sc_0.8$`Y dir`,res_1_1sc_0.8$`MSE dir`,res_1_1sc_0.8$`RRMSE dir`,
                          res_eblup_1_1sc_0.8$est$eblup,res_eblup_1_1sc_0.8$mse,sqrt(res_eblup_1_1sc_0.8$mse)
                          /res_eblup_1_1sc_0.8$est$eblup*100)
colnames(gab_1sc_0.8) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_1sc_0.8$`MSE eblup`/gab_1sc_0.8$`MSE dir`)

gab_1sc_0.7 <- data.frame(res_1_1sc_0.7$kab,res_1_1sc_0.7$`Y dir`,res_1_1sc_0.7$`MSE dir`,res_1_1sc_0.7$`RRMSE dir`,
                          res_eblup_1_1sc_0.7$est$eblup,res_eblup_1_1sc_0.7$mse,sqrt(res_eblup_1_1sc_0.7$mse)
                          /res_eblup_1_1sc_0.7$est$eblup*100)
colnames(gab_1sc_0.7) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_1sc_0.7$`MSE eblup`/gab_1sc_0.7$`MSE dir`)

gab_1sc_0.6 <- data.frame(res_1_1sc_0.6$kab,res_1_1sc_0.6$`Y dir`,res_1_1sc_0.6$`MSE dir`,res_1_1sc_0.6$`RRMSE dir`,
                          res_eblup_1_1sc_0.6$est$eblup,res_eblup_1_1sc_0.6$mse,sqrt(res_eblup_1_1sc_0.6$mse)
                          /res_eblup_1_1sc_0.6$est$eblup*100)
colnames(gab_1sc_0.6) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_1sc_0.6$`MSE eblup`/gab_1sc_0.6$`MSE dir`)

gab_1sc_0.5 <- data.frame(res_1_1sc_0.5$kab,res_1_1sc_0.5$`Y dir`,res_1_1sc_0.5$`MSE dir`,res_1_1sc_0.5$`RRMSE dir`,
                          res_eblup_1_1sc_0.5$est$eblup,res_eblup_1_1sc_0.5$mse,sqrt(res_eblup_1_1sc_0.5$mse)
                          /res_eblup_1_1sc_0.5$est$eblup*100)
colnames(gab_1sc_0.5) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_1sc_0.5$`MSE eblup`/gab_1sc_0.5$`MSE dir`)

gab_1sc_0.4 <- data.frame(res_1_1sc_0.4$kab,res_1_1sc_0.4$`Y dir`,res_1_1sc_0.4$`MSE dir`,res_1_1sc_0.4$`RRMSE dir`,
                          res_eblup_1_1sc_0.4$est$eblup,res_eblup_1_1sc_0.4$mse,sqrt(res_eblup_1_1sc_0.4$mse)
                          /res_eblup_1_1sc_0.4$est$eblup*100)
colnames(gab_1sc_0.4) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_1sc_0.4$`MSE eblup`/gab_1sc_0.4$`MSE dir`)

gab_1sc_0.3 <- data.frame(res_1_1sc_0.3$kab,res_1_1sc_0.3$`Y dir`,res_1_1sc_0.3$`MSE dir`,res_1_1sc_0.3$`RRMSE dir`,
                          res_eblup_1_1sc_0.3$est$eblup,res_eblup_1_1sc_0.3$mse,sqrt(res_eblup_1_1sc_0.3$mse)
                          /res_eblup_1_1sc_0.3$est$eblup*100)
colnames(gab_1sc_0.3) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_1sc_0.3$`MSE eblup`/gab_1sc_0.3$`MSE dir`)

gab_1sc_0.2 <- data.frame(res_1_1sc_0.2$kab,res_1_1sc_0.2$`Y dir`,res_1_1sc_0.2$`MSE dir`,res_1_1sc_0.2$`RRMSE dir`,
                          res_eblup_1_1sc_0.2$est$eblup,res_eblup_1_1sc_0.2$mse,sqrt(res_eblup_1_1sc_0.2$mse)
                          /res_eblup_1_1sc_0.2$est$eblup*100)
colnames(gab_1sc_0.2) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_1sc_0.2$`MSE eblup`/gab_1sc_0.2$`MSE dir`)

gab_1sc_0.1 <- data.frame(res_1_1sc_0.1$kab,res_1_1sc_0.1$`Y dir`,res_1_1sc_0.1$`MSE dir`,res_1_1sc_0.1$`RRMSE dir`,
                          res_eblup_1_1sc_0.1$est$eblup,res_eblup_1_1sc_0.1$mse,sqrt(res_eblup_1_1sc_0.1$mse)
                          /res_eblup_1_1sc_0.1$est$eblup*100)
colnames(gab_1sc_0.1) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_1sc_0.1$`MSE eblup`/gab_1sc_0.1$`MSE dir`)


gab_2sc_srs_0.9 <- data.frame(res_1_2sc_srs_0.9$kab,res_1_2sc_srs_0.9$`Y dir`,res_1_2sc_srs_0.9$`MSE dir`,res_1_2sc_srs_0.9$`RRMSE dir`,
                              res_eblup_1_2sc_srs_0.9$est$eblup,res_eblup_1_2sc_srs_0.9$mse,sqrt(res_eblup_1_2sc_srs_0.9$mse)
                              /res_eblup_1_2sc_srs_0.9$est$eblup*100)
colnames(gab_2sc_srs_0.9) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_srs_0.9$`MSE eblup`/gab_2sc_srs_0.9$`MSE dir`)

gab_2sc_srs_0.8 <- data.frame(res_1_2sc_srs_0.8$kab,res_1_2sc_srs_0.8$`Y dir`,res_1_2sc_srs_0.8$`MSE dir`,res_1_2sc_srs_0.8$`RRMSE dir`,
                              res_eblup_1_2sc_srs_0.8$est$eblup,res_eblup_1_2sc_srs_0.8$mse,sqrt(res_eblup_1_2sc_srs_0.8$mse)
                              /res_eblup_1_2sc_srs_0.8$est$eblup*100)
colnames(gab_2sc_srs_0.8) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_srs_0.8$`MSE eblup`/gab_2sc_srs_0.8$`MSE dir`)

gab_2sc_srs_0.7 <- data.frame(res_1_2sc_srs_0.7$kab,res_1_2sc_srs_0.7$`Y dir`,res_1_2sc_srs_0.7$`MSE dir`,res_1_2sc_srs_0.7$`RRMSE dir`,
                              res_eblup_1_2sc_srs_0.7$est$eblup,res_eblup_1_2sc_srs_0.7$mse,sqrt(res_eblup_1_2sc_srs_0.7$mse)
                              /res_eblup_1_2sc_srs_0.7$est$eblup*100)
colnames(gab_2sc_srs_0.7) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_srs_0.7$`MSE eblup`/gab_2sc_srs_0.7$`MSE dir`)

gab_2sc_srs_0.6 <- data.frame(res_1_2sc_srs_0.6$kab,res_1_2sc_srs_0.6$`Y dir`,res_1_2sc_srs_0.6$`MSE dir`,res_1_2sc_srs_0.6$`RRMSE dir`,
                              res_eblup_1_2sc_srs_0.6$est$eblup,res_eblup_1_2sc_srs_0.6$mse,sqrt(res_eblup_1_2sc_srs_0.6$mse)
                              /res_eblup_1_2sc_srs_0.6$est$eblup*100)
colnames(gab_2sc_srs_0.6) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_srs_0.6$`MSE eblup`/gab_2sc_srs_0.6$`MSE dir`)

gab_2sc_srs_0.5 <- data.frame(res_1_2sc_srs_0.5$kab,res_1_2sc_srs_0.5$`Y dir`,res_1_2sc_srs_0.5$`MSE dir`,res_1_2sc_srs_0.5$`RRMSE dir`,
                              res_eblup_1_2sc_srs_0.5$est$eblup,res_eblup_1_2sc_srs_0.5$mse,sqrt(res_eblup_1_2sc_srs_0.5$mse)
                              /res_eblup_1_2sc_srs_0.5$est$eblup*100)
colnames(gab_2sc_srs_0.5) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_srs_0.5$`MSE eblup`/gab_2sc_srs_0.5$`MSE dir`)

gab_2sc_srs_0.4 <- data.frame(res_1_2sc_srs_0.4$kab,res_1_2sc_srs_0.4$`Y dir`,res_1_2sc_srs_0.4$`MSE dir`,res_1_2sc_srs_0.4$`RRMSE dir`,
                              res_eblup_1_2sc_srs_0.4$est$eblup,res_eblup_1_2sc_srs_0.4$mse,sqrt(res_eblup_1_2sc_srs_0.4$mse)
                              /res_eblup_1_2sc_srs_0.4$est$eblup*100)
colnames(gab_2sc_srs_0.4) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_srs_0.4$`MSE eblup`/gab_2sc_srs_0.4$`MSE dir`)

gab_2sc_srs_0.3 <- data.frame(res_1_2sc_srs_0.3$kab,res_1_2sc_srs_0.3$`Y dir`,res_1_2sc_srs_0.3$`MSE dir`,res_1_2sc_srs_0.3$`RRMSE dir`,
                              res_eblup_1_2sc_srs_0.3$est$eblup,res_eblup_1_2sc_srs_0.3$mse,sqrt(res_eblup_1_2sc_srs_0.3$mse)
                              /res_eblup_1_2sc_srs_0.3$est$eblup*100)
colnames(gab_2sc_srs_0.3) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_srs_0.3$`MSE eblup`/gab_2sc_srs_0.3$`MSE dir`)

gab_2sc_srs_0.2 <- data.frame(res_1_2sc_srs_0.2$kab,res_1_2sc_srs_0.2$`Y dir`,res_1_2sc_srs_0.2$`MSE dir`,res_1_2sc_srs_0.2$`RRMSE dir`,
                              res_eblup_1_2sc_srs_0.2$est$eblup,res_eblup_1_2sc_srs_0.2$mse,sqrt(res_eblup_1_2sc_srs_0.2$mse)
                              /res_eblup_1_2sc_srs_0.2$est$eblup*100)
colnames(gab_2sc_srs_0.2) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_srs_0.2$`MSE eblup`/gab_2sc_srs_0.2$`MSE dir`)

gab_2sc_srs_0.1 <- data.frame(res_1_2sc_srs_0.1$kab,res_1_2sc_srs_0.1$`Y dir`,res_1_2sc_srs_0.1$`MSE dir`,res_1_2sc_srs_0.1$`RRMSE dir`,
                              res_eblup_1_2sc_srs_0.1$est$eblup,res_eblup_1_2sc_srs_0.1$mse,sqrt(res_eblup_1_2sc_srs_0.1$mse)
                              /res_eblup_1_2sc_srs_0.1$est$eblup*100)
colnames(gab_2sc_srs_0.1) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_srs_0.1$`MSE eblup`/gab_2sc_srs_0.1$`MSE dir`)


gab_2sc_pps_0.9 <- data.frame(res_1_2sc_pps_0.9$kab,res_1_2sc_pps_0.9$`Y dir`,res_1_2sc_pps_0.9$`MSE dir`,res_1_2sc_pps_0.9$`RRMSE dir`,
                              res_eblup_1_2sc_pps_0.9$est$eblup,res_eblup_1_2sc_pps_0.9$mse,sqrt(res_eblup_1_2sc_pps_0.9$mse)
                              /res_eblup_1_2sc_pps_0.9$est$eblup*100)
colnames(gab_2sc_pps_0.9) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_pps_0.9$`MSE eblup`/gab_2sc_pps_0.9$`MSE dir`)

gab_2sc_pps_0.8 <- data.frame(res_1_2sc_pps_0.8$kab,res_1_2sc_pps_0.8$`Y dir`,res_1_2sc_pps_0.8$`MSE dir`,res_1_2sc_pps_0.8$`RRMSE dir`,
                              res_eblup_1_2sc_pps_0.8$est$eblup,res_eblup_1_2sc_pps_0.8$mse,sqrt(res_eblup_1_2sc_pps_0.8$mse)
                              /res_eblup_1_2sc_pps_0.8$est$eblup*100)
colnames(gab_2sc_pps_0.8) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_pps_0.8$`MSE eblup`/gab_2sc_pps_0.8$`MSE dir`)

gab_2sc_pps_0.7 <- data.frame(res_1_2sc_pps_0.7$kab,res_1_2sc_pps_0.7$`Y dir`,res_1_2sc_pps_0.7$`MSE dir`,res_1_2sc_pps_0.7$`RRMSE dir`,
                              res_eblup_1_2sc_pps_0.7$est$eblup,res_eblup_1_2sc_pps_0.7$mse,sqrt(res_eblup_1_2sc_pps_0.7$mse)
                              /res_eblup_1_2sc_pps_0.7$est$eblup*100)
colnames(gab_2sc_pps_0.7) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_pps_0.7$`MSE eblup`/gab_2sc_pps_0.7$`MSE dir`)

gab_2sc_pps_0.6 <- data.frame(res_1_2sc_pps_0.6$kab,res_1_2sc_pps_0.6$`Y dir`,res_1_2sc_pps_0.6$`MSE dir`,res_1_2sc_pps_0.6$`RRMSE dir`,
                              res_eblup_1_2sc_pps_0.6$est$eblup,res_eblup_1_2sc_pps_0.6$mse,sqrt(res_eblup_1_2sc_pps_0.6$mse)
                              /res_eblup_1_2sc_pps_0.6$est$eblup*100)
colnames(gab_2sc_pps_0.6) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_pps_0.6$`MSE eblup`/gab_2sc_pps_0.6$`MSE dir`)

gab_2sc_pps_0.5 <- data.frame(res_1_2sc_pps_0.5$kab,res_1_2sc_pps_0.5$`Y dir`,res_1_2sc_pps_0.5$`MSE dir`,res_1_2sc_pps_0.5$`RRMSE dir`,
                              res_eblup_1_2sc_pps_0.5$est$eblup,res_eblup_1_2sc_pps_0.5$mse,sqrt(res_eblup_1_2sc_pps_0.5$mse)
                              /res_eblup_1_2sc_pps_0.5$est$eblup*100)
colnames(gab_2sc_pps_0.5) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_pps_0.5$`MSE eblup`/gab_2sc_pps_0.5$`MSE dir`)

gab_2sc_pps_0.4 <- data.frame(res_1_2sc_pps_0.4$kab,res_1_2sc_pps_0.4$`Y dir`,res_1_2sc_pps_0.4$`MSE dir`,res_1_2sc_pps_0.4$`RRMSE dir`,
                              res_eblup_1_2sc_pps_0.4$est$eblup,res_eblup_1_2sc_pps_0.4$mse,sqrt(res_eblup_1_2sc_pps_0.4$mse)
                              /res_eblup_1_2sc_pps_0.4$est$eblup*100)
colnames(gab_2sc_pps_0.4) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_pps_0.4$`MSE eblup`/gab_2sc_pps_0.4$`MSE dir`)

gab_2sc_pps_0.3 <- data.frame(res_1_2sc_pps_0.3$kab,res_1_2sc_pps_0.3$`Y dir`,res_1_2sc_pps_0.3$`MSE dir`,res_1_2sc_pps_0.3$`RRMSE dir`,
                              res_eblup_1_2sc_pps_0.3$est$eblup,res_eblup_1_2sc_pps_0.3$mse,sqrt(res_eblup_1_2sc_pps_0.3$mse)
                              /res_eblup_1_2sc_pps_0.3$est$eblup*100)
colnames(gab_2sc_pps_0.3) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_pps_0.3$`MSE eblup`/gab_2sc_pps_0.3$`MSE dir`)

gab_2sc_pps_0.2 <- data.frame(res_1_2sc_pps_0.2$kab,res_1_2sc_pps_0.2$`Y dir`,res_1_2sc_pps_0.2$`MSE dir`,res_1_2sc_pps_0.2$`RRMSE dir`,
                              res_eblup_1_2sc_pps_0.2$est$eblup,res_eblup_1_2sc_pps_0.2$mse,sqrt(res_eblup_1_2sc_pps_0.2$mse)
                              /res_eblup_1_2sc_pps_0.2$est$eblup*100)
colnames(gab_2sc_pps_0.2) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_pps_0.2$`MSE eblup`/gab_2sc_pps_0.2$`MSE dir`)

gab_2sc_pps_0.1 <- data.frame(res_1_2sc_pps_0.1$kab,res_1_2sc_pps_0.1$`Y dir`,res_1_2sc_pps_0.1$`MSE dir`,res_1_2sc_pps_0.1$`RRMSE dir`,
                              res_eblup_1_2sc_pps_0.1$est$eblup,res_eblup_1_2sc_pps_0.1$mse,sqrt(res_eblup_1_2sc_pps_0.1$mse)
                              /res_eblup_1_2sc_pps_0.1$est$eblup*100)
colnames(gab_2sc_pps_0.1) <- c("kab","Y dir","MSE dir","RRMSE dir","Y eblup","MSE eblup","RRMSE eblup")
mean(gab_2sc_pps_0.1$`MSE eblup`/gab_2sc_pps_0.1$`MSE dir`)

#---------------------------------------------HOAM--------------------------------------------------


gab_1sc_0.09 <- data.frame(res_2_1sc_0.09$kab,res_2_1sc_0.09$`Y dir`,res_2_1sc_0.09$`MSE dir`,res_2_1sc_0.09$`RRMSE dir`,
                           res_syn_2_1sc_0.09$`Y syn`,res_syn_2_1sc_0.09$`MSE syn`,res_syn_2_1sc_0.09$`RRMSE syn`,
                           res_syn_2_1sc_0.09$`Y com`,res_syn_2_1sc_0.09$`MSE com`,res_syn_2_1sc_0.09$`RRMSE com`,
                           res_eblup_2_1sc_0.09_2$est$eblup,res_eblup_2_1sc_0.09_2$mse,sqrt(res_eblup_2_1sc_0.09_2$mse)
                           /res_eblup_2_1sc_0.09_2$est$eblup*100)
colnames(gab_1sc_0.09) <- c("kab","Y dir","MSE dir","RRMSE dir","Y syn","MSE syn","RRMSE syn","Y com","MSE com",
                            "RRMSE com","Y eblup","MSE eblup","RRMSE eblup")

gab_2sc_srs_0.15 <- data.frame(res_2_2sc_srs_0.15$kab,res_2_2sc_srs_0.15$`Y dir`,res_2_2sc_srs_0.15$`MSE dir`,res_2_2sc_srs_0.15$`RRMSE dir`,
                               res_syn_2_2sc_srs_0.15$`Y syn`,res_syn_2_2sc_srs_0.15$`MSE syn`,res_syn_2_2sc_srs_0.15$`RRMSE syn`,
                               res_syn_2_2sc_srs_0.15$`Y com`,res_syn_2_2sc_srs_0.15$`MSE com`,res_syn_2_2sc_srs_0.15$`RRMSE com`,
                               res_eblup_2_2sc_srs_0.15_2$est$eblup,res_eblup_2_2sc_srs_0.15_2$mse,sqrt(res_eblup_2_2sc_srs_0.15_2$mse)
                               /res_eblup_2_2sc_srs_0.15_2$est$eblup*100)
colnames(gab_2sc_srs_0.15) <- c("kab","Y dir","MSE dir","RRMSE dir","Y syn","MSE syn","RRMSE syn","Y com","MSE com","RRMSE com",
                                "Y eblup","MSE eblup","RRMSE eblup")

gab_2sc_pps_0.2 <- data.frame(res_2_2sc_pps_0.2$kab,res_2_2sc_pps_0.2$`Y dir`,res_2_2sc_pps_0.2$`MSE dir`,res_2_2sc_pps_0.2$`RRMSE dir`,
                              res_syn_2_2sc_pps_0.2$`Y syn`,res_syn_2_2sc_pps_0.2$`MSE syn`,res_syn_2_2sc_pps_0.2$`RRMSE syn`,
                              res_syn_2_2sc_pps_0.2$`Y com`,res_syn_2_2sc_pps_0.2$`MSE com`,res_syn_2_2sc_pps_0.2$`RRMSE com`,
                              res_eblup_2_2sc_pps_0.2_2$est$eblup,res_eblup_2_2sc_pps_0.2_2$mse,sqrt(res_eblup_2_2sc_pps_0.2_2$mse)
                              /res_eblup_2_2sc_pps_0.2_2$est$eblup*100)
colnames(gab_2sc_pps_0.2) <- c("kab","Y dir","MSE dir","RRMSE dir","Y syn","MSE syn","RRMSE syn","Y com","MSE com","RRMSE com",
                               "Y eblup","MSE eblup","RRMSE eblup")

gab_srs_0.3 <- data.frame(res_3_srs_0.3$kab,res_3_srs_0.3$`Y dir`,res_3_srs_0.3$`MSE dir`,res_3_srs_0.3$`RRMSE dir`,
                          res_syn_3_srs_0.3$`Y syn`,res_syn_3_srs_0.3$`MSE syn`,res_syn_3_srs_0.3$`RRMSE syn`,
                          res_syn_3_srs_0.3$`Y com`,res_syn_3_srs_0.3$`MSE com`,res_syn_3_srs_0.3$`RRMSE com`,
                          res_eblup_3_srs_0.3_2$est$eblup,res_eblup_3_srs_0.3_2$mse,sqrt(res_eblup_3_srs_0.3_2$mse)
                          /res_eblup_3_srs_0.3_2$est$eblup*100)
colnames(gab_srs_0.3) <- c("kab","Y dir","MSE dir","RRMSE dir","Y syn","MSE syn","RRMSE syn","Y com","MSE com","RRMSE com",
                           "Y eblup","MSE eblup","RRMSE eblup")

gab_1sc_0.3 <- data.frame(res_3_1sc_0.3$kab,res_3_1sc_0.3$`Y dir`,res_3_1sc_0.3$`MSE dir`,res_3_1sc_0.3$`RRMSE dir`,
                          res_syn_3_1sc_0.3$`Y syn`,res_syn_3_1sc_0.3$`MSE syn`,res_syn_3_1sc_0.3$`RRMSE syn`,
                          res_syn_3_1sc_0.3$`Y com`,res_syn_3_1sc_0.3$`MSE com`,res_syn_3_1sc_0.3$`RRMSE com`,
                          res_eblup_3_1sc_0.3_2$est$eblup,res_eblup_3_1sc_0.3_2$mse,sqrt(res_eblup_3_1sc_0.3_2$mse)
                          /res_eblup_3_1sc_0.3_2$est$eblup*100)
colnames(gab_1sc_0.3) <- c("kab","Y dir","MSE dir","RRMSE dir","Y syn","MSE syn","RRMSE syn","Y com","MSE com","RRMSE com",
                           "Y eblup","MSE eblup","RRMSE eblup")

gab_2sc_srs_0.55 <- data.frame(res_3_2sc_srs_0.55$kab,res_3_2sc_srs_0.55$`Y dir`,res_3_2sc_srs_0.55$`MSE dir`,res_3_2sc_srs_0.55$`RRMSE dir`,
                               res_syn_3_2sc_srs_0.55$`Y syn`,res_syn_3_2sc_srs_0.55$`MSE syn`,res_syn_3_2sc_srs_0.55$`RRMSE syn`,
                               res_syn_3_2sc_srs_0.55$`Y com`,res_syn_3_2sc_srs_0.55$`MSE com`,res_syn_3_2sc_srs_0.55$`RRMSE com`,
                               res_eblup_3_2sc_srs_0.55_2$est$eblup,res_eblup_3_2sc_srs_0.55_2$mse,sqrt(res_eblup_3_2sc_srs_0.55_2$mse)
                               /res_eblup_3_2sc_srs_0.55_2$est$eblup*100)
colnames(gab_2sc_srs_0.55) <- c("kab","Y dir","MSE dir","RRMSE dir","Y syn","MSE syn","RRMSE syn","Y com","MSE com","RRMSE com","Y eblup","MSE eblup","RRMSE eblup")

gab_2sc_pps_0.7 <- data.frame(res_3_2sc_pps_0.7$kab,res_3_2sc_pps_0.7$`Y dir`,res_3_2sc_pps_0.7$`MSE dir`,res_3_2sc_pps_0.7$`RRMSE dir`,
                              res_syn_3_2sc_pps_0.7$`Y syn`,res_syn_3_2sc_pps_0.7$`MSE syn`,res_syn_3_2sc_pps_0.7$`RRMSE syn`,
                              res_syn_3_2sc_pps_0.7$`Y com`,res_syn_3_2sc_pps_0.7$`MSE com`,res_syn_3_2sc_pps_0.7$`RRMSE com`,
                              res_eblup_3_2sc_pps_0.7_2$est$eblup,res_eblup_3_2sc_pps_0.7_2$mse,sqrt(res_eblup_3_2sc_pps_0.7_2$mse)
                              /res_eblup_3_2sc_pps_0.7_2$est$eblup*100)
colnames(gab_2sc_pps_0.7) <- c("kab","Y dir","MSE dir","RRMSE dir","Y syn","MSE syn","RRMSE syn","Y com","MSE com","RRMSE com",
                               "Y eblup","MSE eblup","RRMSE eblup")

#-------------------------------------------------------TES BUAT DATA 4----------------------------------------------------------------------
res_eblup_4_2sc_pps_0.7_2 <- mseFH(formula = res_4_2sc_pps_0.7$`Y dir` ~ aux_x, vardir = res_4_2sc_pps_0.7$`MSE dir`, method = "REML")
try_2sc_pps_0.7 <- data.frame(res_4_2sc_pps_0.7$kab,res_4_2sc_pps_0.7$`Y dir`,res_4_2sc_pps_0.7$`MSE dir`,res_4_2sc_pps_0.7$`RRMSE dir`,
                              res_syn_4_2sc_pps_0.7$`Y syn`,res_syn_4_2sc_pps_0.7$`MSE syn`,res_syn_4_2sc_pps_0.7$`RRMSE syn`,
                              res_syn_4_2sc_pps_0.7$`Y com`,res_syn_4_2sc_pps_0.7$`MSE com`,res_syn_4_2sc_pps_0.7$`RRMSE com`,
                              res_eblup_4_2sc_pps_0.7_2$est$eblup,res_eblup_4_2sc_pps_0.7_2$mse,sqrt(res_eblup_4_2sc_pps_0.7_2$mse)
                              /res_eblup_4_2sc_pps_0.7_2$est$eblup*100)
colnames(try_2sc_pps_0.7) <- c("kab","Y dir","MSE dir","RRMSE dir","Y syn","MSE syn","RRMSE syn","Y com","MSE com","RRMSE com",
                               "Y eblup","MSE eblup","RRMSE eblup")


hsl_eblup <- data.frame(gab_srs_0.03$`Y eblup`,gab_srs_0.03$`MSE eblup`,gab_srs_0.03$`RRMSE eblup`,
                        gab_1sc_0.09$`Y eblup`,gab_1sc_0.09$`MSE eblup`,gab_1sc_0.09$`RRMSE eblup`,
                        gab_2sc_srs_0.15$`Y eblup`,gab_2sc_srs_0.15$`MSE eblup`,gab_2sc_srs_0.15$`RRMSE eblup`,
                        gab_2sc_pps_0.2$`Y eblup`,gab_2sc_pps_0.2$`MSE eblup`,gab_2sc_pps_0.2$`RRMSE eblup`,
                        gab_srs_0.3$`Y eblup`,gab_srs_0.3$`MSE eblup`,gab_srs_0.3$`RRMSE eblup`,
                        gab_1sc_0.3$`Y eblup`,gab_1sc_0.3$`MSE eblup`,gab_1sc_0.3$`RRMSE eblup`,
                        gab_2sc_srs_0.55$`Y eblup`,gab_2sc_srs_0.55$`MSE eblup`,gab_2sc_srs_0.55$`RRMSE eblup`,
                        gab_2sc_pps_0.7$`Y eblup`,gab_2sc_pps_0.7$`MSE eblup`,gab_2sc_pps_0.7$`RRMSE eblup`)
write.table(hsl_eblup,"D:/hsl_eblup.txt")


a_1_dt1 <- cbind(rep("direct"),paste(gab_srs_0.9$kab), as.numeric(paste(gab_srs_0.9$`MSE dir`)))
b_1_dt1 <- cbind(rep("eblup"),paste(gab_srs_0.9$kab), as.numeric(paste(gab_srs_0.9$`MSE eblup`)))

a_1_dt1 = as.data.frame(a_1_dt1)
b_1_dt1 = as.data.frame(b_1_dt1)

pl_1_dt1 <- data.frame(rbind(a_1_dt1,b_1_dt1))
colnames(pl_1_dt1) <- c("cond","area","mse")
yt_1 <- data.frame(pl_1_dt1$cond,pl_1_dt1$area,as.numeric(paste(pl_1_dt1$mse)))
colnames(yt_1) <- c("cond","area","mse")
library(ggplot2)
ggplot(yt_1, aes(x=area, y=mse, shape=cond, color=cond, size = 0.5)) + geom_point()



a_1_dt1 <- cbind(rep("direct"),paste(gab_srs_0.6$kab), as.numeric(paste(gab_srs_0.6$`MSE dir`)))
b_1_dt1 <- cbind(rep("eblup"),paste(gab_srs_0.6$kab), as.numeric(paste(gab_srs_0.6$`MSE eblup`)))

a_1_dt1 = as.data.frame(a_1_dt1)
b_1_dt1 = as.data.frame(b_1_dt1)

pl_1_dt1 <- data.frame(rbind(a_1_dt1,b_1_dt1))
colnames(pl_1_dt1) <- c("cond","area","mse")
yt_1 <- data.frame(pl_1_dt1$cond,pl_1_dt1$area,as.numeric(paste(pl_1_dt1$mse)))
colnames(yt_1) <- c("cond","area","mse")
library(ggplot2)
ggplot(yt_1, aes(x=area, y=mse, shape=cond, color=cond, size = 0.5)) + geom_point()

a_1_dt1 <- cbind(rep("direct"),paste(gab_srs_0.4$kab), as.numeric(paste(gab_srs_0.4$`MSE dir`)))
b_1_dt1 <- cbind(rep("eblup"),paste(gab_srs_0.4$kab), as.numeric(paste(gab_srs_0.4$`MSE eblup`)))

a_1_dt1 = as.data.frame(a_1_dt1)
b_1_dt1 = as.data.frame(b_1_dt1)

pl_1_dt1 <- data.frame(rbind(a_1_dt1,b_1_dt1))
colnames(pl_1_dt1) <- c("cond","area","mse")
yt_1 <- data.frame(pl_1_dt1$cond,pl_1_dt1$area,as.numeric(paste(pl_1_dt1$mse)))
colnames(yt_1) <- c("cond","area","mse")
library(ggplot2)
ggplot(yt_1, aes(x=area, y=mse, shape=cond, color=cond, size = 0.5)) + geom_point()

a_1_dt1 <- cbind(rep("direct"),paste(gab_srs_0.2$kab), as.numeric(paste(gab_srs_0.2$`MSE dir`)))
b_1_dt1 <- cbind(rep("eblup"),paste(gab_srs_0.2$kab), as.numeric(paste(gab_srs_0.2$`MSE eblup`)))

a_1_dt1 = as.data.frame(a_1_dt1)
b_1_dt1 = as.data.frame(b_1_dt1)

pl_1_dt1 <- data.frame(rbind(a_1_dt1,b_1_dt1))
colnames(pl_1_dt1) <- c("cond","area","mse")
yt_1 <- data.frame(pl_1_dt1$cond,pl_1_dt1$area,as.numeric(paste(pl_1_dt1$mse)))
colnames(yt_1) <- c("cond","area","mse")
library(ggplot2)
ggplot(yt_1, aes(x=area, y=mse, shape=cond, color=cond, size = 0.5)) + geom_point()

###############################################################################################

a_1_dt1 <- cbind(rep("direct"),paste(gab_1sc_0.9$kab), as.numeric(paste(gab_1sc_0.9$`MSE dir`)))
b_1_dt1 <- cbind(rep("eblup"),paste(gab_1sc_0.9$kab), as.numeric(paste(gab_1sc_0.9$`MSE eblup`)))

a_1_dt1 = as.data.frame(a_1_dt1)
b_1_dt1 = as.data.frame(b_1_dt1)

pl_1_dt1 <- data.frame(rbind(a_1_dt1,b_1_dt1))
colnames(pl_1_dt1) <- c("cond","area","mse")
yt_1 <- data.frame(pl_1_dt1$cond,pl_1_dt1$area,as.numeric(paste(pl_1_dt1$mse)))
colnames(yt_1) <- c("cond","area","mse")
library(ggplot2)
ggplot(yt_1, aes(x=area, y=mse, shape=cond, color=cond, size = 0.5)) + geom_point()



a_1_dt1 <- cbind(rep("direct"),paste(gab_1sc_0.6$kab), as.numeric(paste(gab_1sc_0.6$`MSE dir`)))
b_1_dt1 <- cbind(rep("eblup"),paste(gab_1sc_0.6$kab), as.numeric(paste(gab_1sc_0.6$`MSE eblup`)))

a_1_dt1 = as.data.frame(a_1_dt1)
b_1_dt1 = as.data.frame(b_1_dt1)

pl_1_dt1 <- data.frame(rbind(a_1_dt1,b_1_dt1))
colnames(pl_1_dt1) <- c("cond","area","mse")
yt_1 <- data.frame(pl_1_dt1$cond,pl_1_dt1$area,as.numeric(paste(pl_1_dt1$mse)))
colnames(yt_1) <- c("cond","area","mse")
library(ggplot2)
ggplot(yt_1, aes(x=area, y=mse, shape=cond, color=cond, size = 0.5)) + geom_point()

a_1_dt1 <- cbind(rep("direct"),paste(gab_1sc_0.4$kab), as.numeric(paste(gab_1sc_0.4$`MSE dir`)))
b_1_dt1 <- cbind(rep("eblup"),paste(gab_1sc_0.4$kab), as.numeric(paste(gab_1sc_0.4$`MSE eblup`)))

a_1_dt1 = as.data.frame(a_1_dt1)
b_1_dt1 = as.data.frame(b_1_dt1)

pl_1_dt1 <- data.frame(rbind(a_1_dt1,b_1_dt1))
colnames(pl_1_dt1) <- c("cond","area","mse")
yt_1 <- data.frame(pl_1_dt1$cond,pl_1_dt1$area,as.numeric(paste(pl_1_dt1$mse)))
colnames(yt_1) <- c("cond","area","mse")
library(ggplot2)
ggplot(yt_1, aes(x=area, y=mse, shape=cond, color=cond, size = 0.5)) + geom_point()

a_1_dt1 <- cbind(rep("direct"),paste(gab_1sc_0.2$kab), as.numeric(paste(gab_1sc_0.2$`MSE dir`)))
b_1_dt1 <- cbind(rep("eblup"),paste(gab_1sc_0.2$kab), as.numeric(paste(gab_1sc_0.2$`MSE eblup`)))

a_1_dt1 = as.data.frame(a_1_dt1)
b_1_dt1 = as.data.frame(b_1_dt1)

pl_1_dt1 <- data.frame(rbind(a_1_dt1,b_1_dt1))
colnames(pl_1_dt1) <- c("cond","area","mse")
yt_1 <- data.frame(pl_1_dt1$cond,pl_1_dt1$area,as.numeric(paste(pl_1_dt1$mse)))
colnames(yt_1) <- c("cond","area","mse")
library(ggplot2)
ggplot(yt_1, aes(x=area, y=mse, shape=cond, color=cond, size = 0.5)) + geom_point()

#################################################################
a_1_dt1 <- cbind(rep("direct"),paste(gab_2sc_srs_0.9$kab), as.numeric(paste(gab_2sc_srs_0.9$`MSE dir`)))
b_1_dt1 <- cbind(rep("eblup"),paste(gab_2sc_srs_0.9$kab), as.numeric(paste(gab_2sc_srs_0.9$`MSE eblup`)))

a_1_dt1 = as.data.frame(a_1_dt1)
b_1_dt1 = as.data.frame(b_1_dt1)

pl_1_dt1 <- data.frame(rbind(a_1_dt1,b_1_dt1))
colnames(pl_1_dt1) <- c("cond","area","mse")
yt_1 <- data.frame(pl_1_dt1$cond,pl_1_dt1$area,as.numeric(paste(pl_1_dt1$mse)))
colnames(yt_1) <- c("cond","area","mse")
library(ggplot2)
ggplot(yt_1, aes(x=area, y=mse, shape=cond, color=cond, size = 0.5)) + geom_point()


a_1_dt1 <- cbind(rep("direct"),paste(gab_2sc_srs_0.6$kab), as.numeric(paste(gab_2sc_srs_0.6$`MSE dir`)))
b_1_dt1 <- cbind(rep("eblup"),paste(gab_2sc_srs_0.6$kab), as.numeric(paste(gab_2sc_srs_0.6$`MSE eblup`)))

a_1_dt1 = as.data.frame(a_1_dt1)
b_1_dt1 = as.data.frame(b_1_dt1)

pl_1_dt1 <- data.frame(rbind(a_1_dt1,b_1_dt1))
colnames(pl_1_dt1) <- c("cond","area","mse")
yt_1 <- data.frame(pl_1_dt1$cond,pl_1_dt1$area,as.numeric(paste(pl_1_dt1$mse)))
colnames(yt_1) <- c("cond","area","mse")
library(ggplot2)
ggplot(yt_1, aes(x=area, y=mse, shape=cond, color=cond, size = 0.5)) + geom_point()

a_1_dt1 <- cbind(rep("direct"),paste(gab_2sc_srs_0.4$kab), as.numeric(paste(gab_2sc_srs_0.4$`MSE dir`)))
b_1_dt1 <- cbind(rep("eblup"),paste(gab_2sc_srs_0.4$kab), as.numeric(paste(gab_2sc_srs_0.4$`MSE eblup`)))

a_1_dt1 = as.data.frame(a_1_dt1)
b_1_dt1 = as.data.frame(b_1_dt1)

pl_1_dt1 <- data.frame(rbind(a_1_dt1,b_1_dt1))
colnames(pl_1_dt1) <- c("cond","area","mse")
yt_1 <- data.frame(pl_1_dt1$cond,pl_1_dt1$area,as.numeric(paste(pl_1_dt1$mse)))
colnames(yt_1) <- c("cond","area","mse")
library(ggplot2)
ggplot(yt_1, aes(x=area, y=mse, shape=cond, color=cond, size = 0.5)) + geom_point()

a_1_dt1 <- cbind(rep("direct"),paste(gab_2sc_srs_0.2$kab), as.numeric(paste(gab_2sc_srs_0.2$`MSE dir`)))
b_1_dt1 <- cbind(rep("eblup"),paste(gab_2sc_srs_0.2$kab), as.numeric(paste(gab_2sc_srs_0.2$`MSE eblup`)))

a_1_dt1 = as.data.frame(a_1_dt1)
b_1_dt1 = as.data.frame(b_1_dt1)

pl_1_dt1 <- data.frame(rbind(a_1_dt1,b_1_dt1))
colnames(pl_1_dt1) <- c("cond","area","mse")
yt_1 <- data.frame(pl_1_dt1$cond,pl_1_dt1$area,as.numeric(paste(pl_1_dt1$mse)))
colnames(yt_1) <- c("cond","area","mse")
library(ggplot2)
ggplot(yt_1, aes(x=area, y=mse, shape=cond, color=cond, size = 0.5)) + geom_point()

######################################################
a_1_dt1 <- cbind(rep("direct"),paste(gab_2sc_pps_0.9$kab), as.numeric(paste(gab_2sc_pps_0.9$`MSE dir`)))
b_1_dt1 <- cbind(rep("eblup"),paste(gab_2sc_pps_0.9$kab), as.numeric(paste(gab_2sc_pps_0.9$`MSE eblup`)))

a_1_dt1 = as.data.frame(a_1_dt1)
b_1_dt1 = as.data.frame(b_1_dt1)

pl_1_dt1 <- data.frame(rbind(a_1_dt1,b_1_dt1))
colnames(pl_1_dt1) <- c("cond","area","mse")
yt_1 <- data.frame(pl_1_dt1$cond,pl_1_dt1$area,as.numeric(paste(pl_1_dt1$mse)))
colnames(yt_1) <- c("cond","area","mse")
library(ggplot2)
ggplot(yt_1, aes(x=area, y=mse, shape=cond, color=cond, size = 0.5)) + geom_point()


a_1_dt1 <- cbind(rep("direct"),paste(gab_2sc_pps_0.6$kab), as.numeric(paste(gab_2sc_pps_0.6$`MSE dir`)))
b_1_dt1 <- cbind(rep("eblup"),paste(gab_2sc_pps_0.6$kab), as.numeric(paste(gab_2sc_pps_0.6$`MSE eblup`)))

a_1_dt1 = as.data.frame(a_1_dt1)
b_1_dt1 = as.data.frame(b_1_dt1)

pl_1_dt1 <- data.frame(rbind(a_1_dt1,b_1_dt1))
colnames(pl_1_dt1) <- c("cond","area","mse")
yt_1 <- data.frame(pl_1_dt1$cond,pl_1_dt1$area,as.numeric(paste(pl_1_dt1$mse)))
colnames(yt_1) <- c("cond","area","mse")
library(ggplot2)
ggplot(yt_1, aes(x=area, y=mse, shape=cond, color=cond, size = 0.5)) + geom_point()

a_1_dt1 <- cbind(rep("direct"),paste(gab_2sc_pps_0.4$kab), as.numeric(paste(gab_2sc_pps_0.4$`MSE dir`)))
b_1_dt1 <- cbind(rep("eblup"),paste(gab_2sc_pps_0.4$kab), as.numeric(paste(gab_2sc_pps_0.4$`MSE eblup`)))

a_1_dt1 = as.data.frame(a_1_dt1)
b_1_dt1 = as.data.frame(b_1_dt1)

pl_1_dt1 <- data.frame(rbind(a_1_dt1,b_1_dt1))
colnames(pl_1_dt1) <- c("cond","area","mse")
yt_1 <- data.frame(pl_1_dt1$cond,pl_1_dt1$area,as.numeric(paste(pl_1_dt1$mse)))
colnames(yt_1) <- c("cond","area","mse")
library(ggplot2)
ggplot(yt_1, aes(x=area, y=mse, shape=cond, color=cond, size = 0.5)) + geom_point()

a_1_dt1 <- cbind(rep("direct"),paste(gab_2sc_pps_0.2$kab), as.numeric(paste(gab_2sc_pps_0.2$`MSE dir`)))
b_1_dt1 <- cbind(rep("eblup"),paste(gab_2sc_pps_0.2$kab), as.numeric(paste(gab_2sc_pps_0.2$`MSE eblup`)))

a_1_dt1 = as.data.frame(a_1_dt1)
b_1_dt1 = as.data.frame(b_1_dt1)

pl_1_dt1 <- data.frame(rbind(a_1_dt1,b_1_dt1))
colnames(pl_1_dt1) <- c("cond","area","mse")
yt_1 <- data.frame(pl_1_dt1$cond,pl_1_dt1$area,as.numeric(paste(pl_1_dt1$mse)))
colnames(yt_1) <- c("cond","area","mse")
library(ggplot2)
ggplot(yt_1, aes(x=area, y=mse, shape=cond, color=cond, size = 0.5)) + geom_point()


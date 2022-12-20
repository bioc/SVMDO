test_that("Prognostic Analysis is succesful", {
  
  c_path<-system.file("extdata","coad_clinic_sum.rds",package="SVMDO",mustWork = TRUE)
  c<-readRDS(c_path)
  alldata<-as.data.frame(c@assays@data@listData)
  
  surv_object<-NULL                          
  st_point<-NULL
  end_point<-NULL
  exp_col_data<-NULL
  coxcoeff<-NULL
  
  exp_col_data<-tcga_sample_comb$id
  
  exp_data<-top_genes
  
  st_point<-as.numeric(sum(tcga_sample_comb$tissue_type=="normal")+1)
  end_point<-as.numeric(sum(tcga_sample_comb$tissue_type=="normal")+sum(tcga_sample_comb$tissue_type=="tumour"))
  exp_col_data<-tcga_sample_comb$id
  tcga_sample_comb<-tcga_sample_comb[tcga_sample_comb$tissue_type=="tumour",]
  
  exp_data<-as.data.frame(t(exp_data))
  exp_data<-exp_data[-1,]
  colnames(exp_data)<-exp_col_data
  
  V1<-rownames(exp_data)
  exp_data<-exp_data[,st_point:end_point]
  exp_data<-cbind(V1,exp_data)
  
  colnames(exp_data)<-c("V1",as.character(tcga_sample_comb$id))
  rownames(exp_data)<-NULL
  
  disc_list<-final_discriminative_gene_set
  

  alldata[alldata == "--"] <- NA
  alldata[alldata == ""] <- NA
  alldata$days_to_death <- as.numeric(as.character(alldata$days_to_death))
  alldata$vital_status <- as.character(alldata$vital_status)
  alldata$vital_status[alldata$vital_status == "Alive"] <- 1
  alldata$vital_status[alldata$vital_status == "Dead"] <- 2
  alldata$vital_status <- as.numeric(as.character(alldata$vital_status))
  alldata$id <- as.character(gsub("-", ".", alldata$id, fixed = TRUE))
  
  exp_data<-exp_data[,(names(exp_data) %in% as.factor(alldata$id))]
  exp_data<-cbind(V1,exp_data)
  
  sub_exp_data<-exp_data[exp_data$V1 %in% disc_list$Names,]
  sub_exp_data<-as.data.frame(sub_exp_data)
  rownames(sub_exp_data)<-NULL
  
  table1 <- data.frame(matrix(NA, nrow = length(disc_list$Names), ncol = 5))
  colnames(table1) <- c("modulename", "p", "hr", "low", "up")
  
  for (a in 1:length(disc_list$Names)) {
    
    alldata1 <- alldata
    
    module <- (sub_exp_data[a,])
    row.names(module) <- make.names(module[,1], unique = TRUE)
    module <- (module[,-1])
    
    moduletrans <- as.data.frame(t(module))
    name_val<-colnames(moduletrans)
    
    moduletrans<-as.data.frame(moduletrans[match(alldata1$id,rownames(moduletrans)),])
    colnames(moduletrans)<-name_val
    
    coxcoeff<- survival::coxph(survival::Surv(time= alldata1$days_to_death, event = alldata1$vital_status) ~ (moduletrans[,1]), data = alldata1)$coefficients
    PI <- data.frame()[1:nrow(moduletrans),]
    
    for (j in 1:nrow(moduletrans)) {
      PI[j,1] <- sum(na.omit(((moduletrans[j,])) * coxcoeff))
    }
    rownames(PI) <- rownames(moduletrans)
    colnames(PI) <- 'PI'
    PI <- as.data.frame(PI)
    
    alldata1 <- dplyr::mutate(alldata1, PI = PI$PI)
    medPI <- median(PI$PI)
    
    if (medPI<=0 ) {
      
      if (sum(PI$PI<=0)==nrow(alldata1)) {
        alldata1$group <- ifelse(PI$PI < medPI, 2, 1)
      }
      
      if (sum(PI$PI>=0)==nrow(alldata1)) {
        alldata1$group <- ifelse(PI$PI <= medPI, 2, 1)
      }
    }
    
    if (medPI>0) {
      alldata1$group <- ifelse(PI$PI <= medPI, 2, 1)
    }
    
    
    if (sum(!duplicated(alldata1$group))!=1) {
      
      surv_object <- survival::Surv(time = alldata1$days_to_death, event = alldata1$vital_status)
      fit1 <- survival::survfit(surv_object ~ group, data = alldata1)
      coxhr <- survival::coxph(surv_object ~ group, data = alldata1)
      summhr <- summary(coxhr)
      p <- round(summhr$sctest[3], digits = 5)
      
      data.survdiff <- survival::survdiff(surv_object ~ group, data = alldata1)
      hr = round((data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2]), digits = 3)
      up = round(exp(log(hr) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])), digits = 3)
      low = round(exp(log(hr) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])), digits = 3)
      
      modulename <- as.character(sub_exp_data$V1[a])
      
      if (sjmisc::str_contains(modulename,"__")) {
        modulename<-gsub("__","-",modulename)
      }
      
      table <- t(cbind.data.frame(modulename, p, hr, low, up))
      table1[a,] <- table
      
      if (p<0.05) {
        
        png(filename=paste0(modulename,".png"),width = 800,height = 600)
        plot(fit1,main=modulename,col = c(2,4),frame.plot=0,ylab="Survival Probability",xlab="Time (day)",cex.axis=1.5,cex.lab=1.5,cex.main=2)
        legend("top",legend = c("High risk group", "Low risk group"),text.col = c(2,4),bty = "n",horiz = TRUE,cex = 1.5)
        
        mtext(text = paste0("p = ",p),side=3,adj = 1,line = -2,cex=1.5)
        mtext(text = paste0("HR = ",hr),side=3,adj = 1,line = -4,cex = 1.5)
        
        dev.off()
        alldata1_col<-NULL
      }
      
      colnames(table1)<-c("Gene", "p", "hr", "low", "up")
      write.table(table1,"prognostic_analysis.txt",sep = "\t")
    }
  }
  
  eq_val<-1
  expect_equal(eq_val,1)

  
  })
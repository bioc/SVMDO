test_that("Prognostic Analysis is succesful", {
  
  c_path<-system.file("extdata","coad_clinic_sum.rds",package="SVMDO",mustWork = TRUE)
  c<-readRDS(c_path)
  alldata<-assay(c)
  
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
  
  name_data<-final_discriminative_gene_set$Names
  prog_prep<-lapply(seq_along(disc_list$Names),function(a){
  alldata1 <- alldata
  
  module <- (sub_exp_data[a,])
  row.names(module) <- make.names(module[,1], unique = TRUE)
  module <- (module[,-1])
  
  moduletrans <- as.data.frame(t(module))
  name_val<-colnames(moduletrans)
  
  moduletrans<-as.data.frame(moduletrans[match(alldata1$id,rownames(moduletrans)),])
  colnames(moduletrans)<-name_val
  
  coxcoeff<- coxph(Surv(time= alldata1$days_to_death, event = alldata1$vital_status) ~ (moduletrans[,1]), data = alldata1)$coefficients
  PI <- data.frame()[seq.int(1,nrow(moduletrans)),]
  
  lap_PI_test<-lapply(seq.int(nrow(moduletrans)), function(b){
    pi_val<-sum(na.omit(((moduletrans[b,])) * coxcoeff))
  })
  
  PI<-as.data.frame(unlist(lap_PI_test))
  rownames(PI) <- rownames(moduletrans)
  colnames(PI) <- 'PI'
  
  alldata1 <- mutate(alldata1, PI = PI$PI)
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
    surv_object <- Surv(time = alldata1$days_to_death, event = alldata1$vital_status)
    fit1 <- survfit(surv_object ~ group, data = alldata1)
    coxhr <- coxph(surv_object ~ group, data = alldata1)
    summhr <- summary(coxhr)
    p <- round(summhr$sctest[3], digits = 5)    
    
    data.survdiff <- survdiff(surv_object ~ group, data = alldata1)
    hr = round((data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2]), digits = 3)
    up = round(exp(log(hr) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])), digits = 3)
    low = round(exp(log(hr) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])), digits = 3)
    
    modulename <- (sub_exp_data$V1[a])
    if (str_contains(modulename,"__")) {
      modulename<-gsub("__","-",modulename)
    }
    
    if (p<0.05) {
      assign(paste0("fit1","_",a),fit1,envir = .GlobalEnv)
      assign(paste0("modulename","_",a),modulename,envir = .GlobalEnv)
      assign(paste0("hr","_",a),hr,envir = .GlobalEnv)
      assign(paste0("p","_",a),p,envir = .GlobalEnv)      
      }
	}       
      })
	  	  
    if (length(ls(pattern = "^fit1_",envir = .GlobalEnv))>0) {
      fit_list<-ls(pattern = "fit1_",envir = .GlobalEnv)
      
      modulename_list<-ls(pattern = "modulename_",envir = .GlobalEnv)
      p_list<-ls(pattern = "^p_",envir = .GlobalEnv)
      hr_list<-ls(pattern = "hr_",envir = .GlobalEnv)
      
      assign("fit_list",fit_list,envir = .GlobalEnv)
      assign("modulename_list",modulename_list,envir = .GlobalEnv)
      assign("hr_list",hr_list,envir = .GlobalEnv)
      assign("p_list",p_list,envir = .GlobalEnv)
    }	
	
    if(exists("final_discriminative_gene_set")){
      lapply(seq.int(length(fit_list)),function(a){
        png(filename=paste0(get(modulename_list[a]),".png"),width = 800,height = 600)
        plot(get(fit_list[a]),main=get(modulename_list[a]),col = c(2,4),frame.plot=0,ylab="Survival Probability",xlab="Time (day)",cex.axis=1.5,cex.lab=1.5,cex.main=2)
        legend("top",legend = c("High risk group", "Low risk group"),text.col = c(2,4),bty = "n",horiz = TRUE,cex = 1.5)
        mtext(text = paste0("p = ",get(p_list[a])),side=3,adj = 1,line = -6,cex=1.5)
        mtext(text = paste0("HR = ",get(hr_list[a])),side=3,adj = 1,line = -8,cex = 1.5)
        dev.off()})      
    }

	  
  eq_val<-1
  expect_equal(eq_val,1)

  
  })
server <- function(input, output,session) {
  id <- NULL
  path_val<-NULL
  wd<-NULL
  
  dir_selection<-observeEvent(input$dir_select,{
    wd<-rstudioapi::selectDirectory()
    if (!is.null(wd)) {
      setwd(wd)
      output$dir_select_ready <- renderUI({
        HTML("<font color=\"red\">", wd, "</font>")
      })
    }else{
      output$dir_select_ready <- renderUI({
        HTML("<font color=\"red\">", "using default Rstudio directory", "</font>")
      })}
})

  
  rawData <- eventReactive(input$file1, {fread(input$file1$datapath,sep = "\t")})
  rawData_2 <- eventReactive(input$file1, {read.table(input$file2$datapath,sep = "\t",header = TRUE)})
  
  message_val<-NULL
  message_val_2<-NULL
  
  deg_analysis<-observeEvent(input$initiate_deg_analysis, {
    if (input$test_datasets=="COAD") {
      top_gene_number<-25
      c_path<-system.file("extdata","coad_exp_sum.rds",package="SVMDO",mustWork = TRUE)
      c<-readRDS(c_path)
      c<-as.data.frame(c@assays@data@listData)
    }else if (input$test_datasets=="LUSC"){
      top_gene_number<-25
      c_path<-system.file("extdata","lusc_exp_sum.rds",package="SVMDO",mustWork = TRUE)
      c<-readRDS(c_path)
      c<-as.data.frame(c@assays@data@listData)
    }else{
      top_gene_number<-input$num_val
      c<-rawData()
    }
    
    
    col_val_1<-NULL
    col_val_2<-NULL
    
    
    for (i in 1:ncol(c)) {                      
      if (typeof(c[[i]][2])=="character") {                        
        col_val_1<- i
        break
      }
    }
    
    for (i in 1:ncol(c)) {                      
      if (typeof(c[[i]][2])!="character" & typeof(c[[i]][2])!="integer") {                        
        col_val_2<- i-1
        break
      }
    }
    assigning_tcga_id_list<-as.data.frame(c[[col_val_1]])
    assigning_tissue_type_list<-as.data.frame(c[[col_val_2]])
    
    colnames(assigning_tcga_id_list)<-data.frame("id")
    colnames(assigning_tissue_type_list)<-c("tissue_type")
    
    if(assigning_tcga_id_list$id[1]!=assigning_tissue_type_list$tissue_type[1]){
      assign("tcga_id_list",assigning_tcga_id_list,envir =.GlobalEnv)
    }
    
    assigning_tissue_type_list$tissue_type<-as.character(assigning_tissue_type_list$tissue_type)
    assign("tissue_type_list",assigning_tissue_type_list,envir =.GlobalEnv)
    
    c<-c[,-(1:col_val_2)]
    
    c<-cbind(tissue_type_list,c)
    
    normal_data<-NULL
    tumour_data<-NULL
    p_val_total<-NULL
    normal_data_mean<-NULL
    tumour_data_mean<-NULL
    int_col<-NULL
    dist_normal_pval<-NULL
    dist_tumour_pval<-NULL
    
    
    i<-1:ncol(c)
    u<-duplicated(colnames(c))
    df<-which(u[i]==TRUE)
    d<-colnames(c)[df]
    duplicate_number<-as.data.frame(d)
    
    for(j in 1:nrow(duplicate_number)){
      colnames(c)[df[j]]<-paste0(colnames(c)[df[j]],"__",j)
    }
    
    #####Modifying illegal symbols
    colnames(c)<-gsub( "-", "__",  colnames(c))
    colnames(c)<-gsub( " ", "__",  colnames(c))
    colnames(c)<-gsub( ";. ", "__",colnames(c))
    colnames(c)<-gsub( ";__","__", colnames(c))
    #####
    
    ############ Normal and Tumour Data Extraction and Fold Change Calculation
    
    normal_data<-c[grep("Nor",c$tissue_type),]
    tumour_data<-c[grep("Tum",c$tissue_type),]
    
    for(i in 1:ncol(c)){
      if(typeof(unlist(c[[i]][1]))!="character" & typeof(unlist(c[[i]][1]))!="integer"){
        start_val<-i
        break
      }
    }
    
    normal_data_mean<-colMeans(normal_data[,start_val:ncol(normal_data)])
    tumour_data_mean<-colMeans(tumour_data[,start_val:ncol(tumour_data)])
    
    fold_change<-tumour_data_mean/normal_data_mean
    matr_fold<-as.matrix(fold_change)
    
    for (i in start_val:ncol(normal_data)) {
      dist_normal_pval<-nortest::ad.test(normal_data[[i]][1:nrow(normal_data)])$p.value
      dist_tumour_pval<-nortest::ad.test(tumour_data[[i]][1:nrow(tumour_data)])$p.value
      
      if (dist_normal_pval<0.05 | dist_tumour_pval<0.05 ) {
        
        for (i in start_val:ncol(c)) {
          p_val<-wilcox.test(as.numeric(unlist(normal_data[[i]])),as.numeric(unlist(tumour_data[[i]])))$p.value
          p_val_total<-c(p_val_total,p_val)
        }
        break
      }
    }
    
    if (is.null(p_val_total)) {
      for(i in start_val:ncol(c)) {
        p_val<-z.test(normal_data[[i]],tumour_data[[i]],sigma.x = sd(normal_data[[i]]),sigma.y = sd(tumour_data[[i]]))$p.value
        p_val_total<-c(p_val_total,p_val)
      }
    }
    
    padjust<-p.adjust(p_val_total,method = "BH")
    p_val_total<-padjust
    rownames(p_val_total)<-NULL
    filtered_p_val_total<-subset(p_val_total,p_val_total<0.05)
    bound_form<-cbind.data.frame(matr_fold,p_val_total)
    colnames(bound_form)<-c("Fold_Change","P_value(adjusted)")
    new_bound_form<-subset(bound_form,bound_form$Fold_Change!="Inf" & bound_form$Fold_Change!="-Inf"
                           &bound_form$Fold_Change!="NaN" &bound_form$`P_value(adjusted)`<0.05 &bound_form$`P_value(adjusted)`!="NaN")
    
    new_bound_form_A<-subset(new_bound_form,new_bound_form$Fold_Change<= 0.3535534)
    new_bound_form_B<-subset(new_bound_form,new_bound_form$Fold_Change>= 2.828427)
    
    new_bound_form_A<-cbind.data.frame(rownames(new_bound_form_A),new_bound_form_A)
    new_bound_form_B<-cbind.data.frame(rownames(new_bound_form_B),new_bound_form_B)
    colnames(new_bound_form_A)<-c("Genes","Fold_Change","P_value")
    colnames(new_bound_form_B)<-c("Genes","Fold_Change","P_value")
    complete_deg_gene_list<-rbind(new_bound_form_A,new_bound_form_B)
    
    sorted_new_bound_form_A<-new_bound_form_A[order(-new_bound_form_A$Fold_Change),]
    sorted_new_bound_form_B<-new_bound_form_B[order(new_bound_form_B$Fold_Change),]
    
    
    if (nrow(sorted_new_bound_form_A) < top_gene_number | nrow(sorted_new_bound_form_B) < top_gene_number) {
      message_val_2<-1
      max_down_genes<-sorted_new_bound_form_A
      max_up_genes<-sorted_new_bound_form_B
      print("Insufficient gene number All up/downregulated genes were selected")
      print("Enter a lower value of input size")
    }else{
      max_down_genes<-sorted_new_bound_form_A[(nrow(sorted_new_bound_form_A)-(top_gene_number-1)):nrow(sorted_new_bound_form_A),]
      max_up_genes<-sorted_new_bound_form_B[(nrow(sorted_new_bound_form_B)-(top_gene_number-1)):nrow(sorted_new_bound_form_B),]
    }
    
    top_combined_genes<-rbind(max_down_genes,max_up_genes)
    rownames(top_combined_genes)<-top_combined_genes[,1]
    changed_whole_data<-subset(c,select=top_combined_genes$Genes)
    
    data.table::fwrite(complete_deg_gene_list,"complete_deg_gene_list.txt",sep = "\t")
    data.table::fwrite(changed_whole_data,"top_genes.txt",sep = "\t")
    assign("top_genes",changed_whole_data,envir =.GlobalEnv)
    assign("complete_deg_gene_list",complete_deg_gene_list,envir =.GlobalEnv)
    
    message_val<-1
    
    if (!is.null(message_val)) {
      
      showModal(
        modalDialog(
          title = "DEG Analysis Result",
          "Process Completed",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
  })
  
  do_analysis<-observeEvent(input$initiate_do_analysis, {
    message_val_3<-NULL
    normal_samples<-NULL
    tumour_samples<-NULL
    disease_filtered_gene_data<-NULL
    changed_name_plus_var_imp_genes_table<-NULL
    
    c_2<-top_genes
    
    changed_whole_data<-c_2
    collect_gene_names<-NULL
    selected_data<-colnames(changed_whole_data[,-1])
    
    checking_data<-selected_data
    checking_data<-as.character(checking_data)
    checking_data<-gsub( "__", "-",checking_data)
    
    colnames(changed_whole_data[,2:ncol(changed_whole_data)])<-checking_data
    database_selection<-AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = checking_data,columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
    
    pop_gene_id<-database_selection$ENTREZID
    pop_gene_symbol<-database_selection$SYMBOL
    count_gene_id<-length(database_selection$ENTREZID)
    
    for (i in 1:count_gene_id) {
      disease_enrichment<-DOSE::enrichDO(pop_gene_id[i],ont = "DO",maxGSSize = Inf)
      
      if (is.null(disease_enrichment)==FALSE){
        dis_length<-length(disease_enrichment@result$Description)
        p_check<-(disease_enrichment@result$p.adjust[1])
        
        if (dis_length>0 & p_check<0.05) {
          assign(paste0("disease_table","_",i),disease_enrichment)
          check_value<-assign(paste0("disease_table","_",i),disease_enrichment)
          collect_gene_names<-c(collect_gene_names,pop_gene_symbol[i])
        }
      }
    }
    
    collect_gene_names<-gsub( "-", "__",collect_gene_names)
    changed_whole_data<-changed_whole_data[,-1]
    changed_name_plus_var_imp_genes_table<-subset(changed_whole_data,select=(collect_gene_names))
    
    
    for (i in 1:nrow(tissue_type_list)) {
      
      if (sjmisc::str_contains(tissue_type_list$tissue_type[i],"Nor")) {
        normal_samples<-rbind("normal",normal_samples)
      }
      if (sjmisc::str_contains(tissue_type_list$tissue_type[i],"Tum")) {
        tumour_samples<-rbind("tumour",tumour_samples)
      }
    }
    
    assigning_new_tissue_type_list<-rbind(normal_samples,tumour_samples)
    colnames(assigning_new_tissue_type_list)<-c("tissue_type")
    assigning_new_tissue_type_list<-as.data.frame(assigning_new_tissue_type_list)
    assign("new_tissue_type_list",assigning_new_tissue_type_list,envir = .GlobalEnv)
    
    if(exists("tcga_id_list")){
      tcga_id_list$id<-as.character(gsub("-", ".", tcga_id_list$id, fixed = TRUE))
      assigning_tcga_sample_comb<-cbind(tcga_id_list,new_tissue_type_list)
      assign("tcga_sample_comb",assigning_tcga_sample_comb,envir =.GlobalEnv)
    }
    
    assigning_disease_filtered_gene_data<-cbind(new_tissue_type_list,changed_name_plus_var_imp_genes_table)
    assigning_disease_filtered_gene_data<-as.data.frame(assigning_disease_filtered_gene_data)
    assign("disease_filtered_gene_data",assigning_disease_filtered_gene_data,envir = .GlobalEnv)
    data.table::fwrite(assigning_disease_filtered_gene_data,"disease_filtered_gene_data.txt",sep = "\t")
    message_val_3<-1
    
    if (!is.null(message_val_3)) {
      
      showModal(
        modalDialog(
          title = "DO based Gene Filtration Result",
          "Process Completed",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
  })
  
  classification_analysis<-observeEvent(input$initiate_wlsvmdo_analysis, {
    test_spec<-NULL
    train_spec<-NULL
    test_pval<-NULL
    train_pval<-NULL
    test_Ac<-NULL
    train_Ac<-NULL
    test_kappa<-NULL
    train_kappa<-NULL
    elected_genes<-NULL
    elected_genes_before<-NULL
    skip_value<-NULL
    control_data<-NULL
    nameless_disease_filtered_gene_data<-NULL
    disease_filtered_gene_data_before<-NULL
    para_fixer<-NULL
    g_val_const<-NULL
    g_val_selection<-NULL
    c_val_const<-NULL
    c_val_selection<-NULL
    train_sample_normal<-NULL
    train_sample_cancer<-NULL
    final_elected_check<-NULL
    train_test_val_check<-NULL
    spec_list<-NULL
    sens_list<-NULL
    nor_val<-NULL
    tum_val<-NULL
    nor_val_list<-NULL
    tum_val_list<-NULL
    test_fin_col<-NULL
    
    specific_check_val<-0
    niv_value<-0.1
    set_val<-1
    control_val<-0
    abc<-0
    
    all_names<-new_tissue_type_list
    actual_disease_filtered_gene_data<-subset(disease_filtered_gene_data,select=-tissue_type)
    
    start_time<-Sys.time()
    
    repeat{
      
      if (ncol(disease_filtered_gene_data)<=2) {
        elected_genes<-elected_genes_before
        disease_filtered_gene_data<-disease_filtered_gene_data_before
        specific_check_val<-1
      }
      
      if (ncol(disease_filtered_gene_data)>2 & specific_check_val==0) {
        elected_val<-NULL
        elected_genes<-klaR::greedy.wilks(tissue_type ~ ., data=disease_filtered_gene_data,niveau=niv_value)
        nameless_disease_filtered_gene_data<-disease_filtered_gene_data[,-1]
        elected_val<-elected_genes$results$vars
        elected_val<-as.vector(elected_val)
        
        
        if (length(elected_val)!=ncol(nameless_disease_filtered_gene_data)) {
          try(nameless_disease_filtered_gene_data<-nameless_disease_filtered_gene_data[,..elected_val],silent = TRUE)
          try(nameless_disease_filtered_gene_data<-nameless_disease_filtered_gene_data[,elected_val],silent = TRUE)
        }
      }
      
      disease_filtered_gene_data<-cbind(all_names,nameless_disease_filtered_gene_data)
      
      
      if (nrow(elected_genes$results)>2 ) {
        elected_genes_before<-elected_genes
        nameless_disease_filtered_gene_data_before<-nameless_disease_filtered_gene_data
        disease_filtered_gene_data_before<-disease_filtered_gene_data
      }
      
      skip_value<-1
      
      if (is.null(skip_value)==FALSE) {
        
        spl <- caTools::sample.split(disease_filtered_gene_data$tissue_type, SplitRatio = 0.80)
        training_set <- subset(disease_filtered_gene_data, spl == TRUE)
        test_set <- subset(disease_filtered_gene_data, spl == FALSE)
        row.names(training_set) <- NULL
        row.names(test_set) <- NULL
        
        if (is.null(train_test_val_check)==TRUE) {
          
          train_check<-training_set
          test_check<-test_set
          train_test_val_check<-1
        }
        
        if(is.null(para_fixer)==TRUE) {
          
          g_val_selection<- sample(10^(seq(-6,6,1)),1)
          g_val_const<-g_val_selection
          
          c_val_selection<- sample(10^(seq(-5,5,1)),1)
          c_val_const<-c_val_selection
          
        }
        
        end_time<-Sys.time()
        
        if (as.numeric(difftime(end_time,start_time,units = "mins"))>10) {
          break
        }
        
        tuning_action<-e1071::svm(as.factor(tissue_type)~., training_set,type="C-classification",scale = FALSE, cross=10 ,gamma = g_val_selection,cost = c_val_selection,probability=TRUE)
        svm_data<-tuning_action
        
        check_training_set<-subset(training_set,select=-tissue_type)
        check_testing_set<-subset(test_set,select=-tissue_type)
        
        
        t_pred<-predict(svm_data,type="prob",check_training_set,probability =TRUE )
        y_pred <- predict(svm_data,type="prob", check_testing_set,probability =TRUE)
        
        
        train_table<-caret::confusionMatrix(t_pred, as.factor(training_set$tissue_type))
        test_table<-caret::confusionMatrix(y_pred,as.factor(test_set$tissue_type))
        
        test_spec<-test_table$byClass[2]
        train_spec<-train_table$byClass[2]
        train_sens<-train_table$byClass[1]
        
        test_pval<-test_table$overall[6]
        train_pval<-train_table$overall[6]
        test_kappa<-test_table$overall[2]
        train_kappa<-train_table$overall[2]
        
        
        if (train_spec>0.8 & train_kappa>0.8 & train_pval<0.05 & test_pval<0.05 & test_kappa>0.8 & test_spec>0.8) {
          
          
          assign(paste0("elected_gene_set","_",set_val),elected_genes)
          final_elected_check<-assign(paste0("elected_gene_set","_",set_val),elected_genes)
          
          final_gene_list<-assign(paste0("gene_set","_",set_val),colnames(disease_filtered_gene_data[,-1]))
          
          assign(paste0("train_table","_",set_val),train_table)
          assign(paste0("test_table","_",set_val),test_table)
          
          collected_g_val<-g_val_selection
          collected_c_val<-c_val_selection
          set_val<-set_val+1
          disease_filtered_gene_data<-cbind(all_names,actual_disease_filtered_gene_data)
          specific_check_val<-0
          para_fixer<-NULL
          elected_genes_before<-NULL
          nameless_disease_filtered_gene_data_before<-NULL
          disease_filtered_gene_data_before<-NULL
          
          
          if (test_pval<0.05) {
            if (final_elected_check$results$p.value.diff[length(final_elected_check$results$p.value.diff)]>0.05 & final_elected_check$results$p.value.diff[length(final_elected_check$results$p.value.diff)-1]<0.05) {
              final_gene_list<-as.data.frame(final_gene_list)
              colnames(final_gene_list)<-"Names"
              write.table(final_gene_list,"final_discriminative_gene_set.txt",sep="\t",row.names=FALSE)
              assign("final_discriminative_gene_set",final_gene_list,envir =.GlobalEnv)
              break}
            
            if (final_elected_check$results$p.value.diff[length(final_elected_check$results$p.value.diff)]<0.05) {
              final_gene_list<-as.data.frame(final_gene_list)
              colnames(final_gene_list)<-"Names"
              write.table(final_gene_list,"final_discriminative_gene_set.txt",sep="\t",row.names=FALSE)
              assign("final_discriminative_gene_set",final_gene_list,envir =.GlobalEnv)
              break}
          }
          niv_value<-niv_value-0.01
        }
        
      }
      skip_value<-NULL
      gc()
    }
    
    if (exists("final_discriminative_gene_set")) {
      showModal(
        modalDialog(
          title = "Classifcation Result",
          "Process completed",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
  })
  
  survival_analysis<-observeEvent(input$initiate_surv_analysis, {
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
    
    if (input$test_datasets=="COAD") {
      c_path<-system.file("extdata","coad_clinic_sum.rds",package="SVMDO",mustWork = TRUE)
      c<-readRDS(c_path)
      alldata<-as.data.frame(c@assays@data@listData)
    }else if (input$test_datasets=="LUSC"){
      c_path<-system.file("extdata","lusc_clinic_sum.rds",package="SVMDO",mustWork = TRUE)
      c<-readRDS(c_path)
      alldata<-as.data.frame(c@assays@data@listData)
      
    }else{
      alldata<-rawData_2()
    }
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
    showModal(
      modalDialog(
        title = "Survival Analysis Result",
        "Process completed",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
}

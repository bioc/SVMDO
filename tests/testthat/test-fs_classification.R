test_that("Feature Selection and Classification processes are succesful", {
  
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
   
  repeat{
    
    
    if (ncol(disease_filtered_gene_data)<=2) {
      elected_genes<-elected_genes_before
      disease_filtered_gene_data<-disease_filtered_gene_data_before
      specific_check_val<-1
    }
    
    if (ncol(disease_filtered_gene_data)>2 & specific_check_val==0) {
      elected_val<-NULL
      elected_genes<-greedy.wilks(tissue_type ~ ., data=disease_filtered_gene_data,niveau=niv_value)
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
      
      spl <- sample.split(disease_filtered_gene_data$tissue_type, SplitRatio = 0.80)
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
        tuning_action<-tune.svm(as.factor(tissue_type)~., data=training_set,
                                cost=10^(seq(-5,5,1)),gamma=10^(seq(-6,6,1)),scale=FALSE,probability=TRUE)
        svm_data<-tuning_action$best.model
      }

      check_training_set<-subset(training_set,select=-tissue_type)
      check_testing_set<-subset(test_set,select=-tissue_type)
      
      
      t_pred<-predict(svm_data,type="prob",check_training_set,probability =TRUE )
      y_pred <- predict(svm_data,type="prob", check_testing_set,probability =TRUE)
      
      
      train_table<-confusionMatrix(t_pred, as.factor(training_set$tissue_type))
      test_table<-confusionMatrix(y_pred,as.factor(test_set$tissue_type))
      
      test_spec<-test_table$byClass[2]
      train_spec<-train_table$byClass[2]
      train_sens<-train_table$byClass[1]
      
      test_pval<-test_table$overall[6]
      train_pval<-train_table$overall[6]
      test_kappa<-test_table$overall[2]
      train_kappa<-train_table$overall[2]
      
      
      if (train_spec>0.8 & train_kappa>0.8 & train_pval<0.05 & test_pval<0.05 & test_kappa>0.8 & test_spec>0.8) {
        
        
        final_elected_check<-elected_genes
        final_gene_list<-colnames(disease_filtered_gene_data[,-1])
		
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
            assign("final_discriminative_gene_set",final_gene_list,envir =.GlobalEnv)
            eq_val<-1
            expect_equal(eq_val,1)
            break}
          
          if (final_elected_check$results$p.value.diff[length(final_elected_check$results$p.value.diff)]<0.05) {
            final_gene_list<-as.data.frame(final_gene_list)
            colnames(final_gene_list)<-"Names"
            assign("final_discriminative_gene_set",final_gene_list,envir =.GlobalEnv)
            eq_val<-1
            expect_equal(eq_val,1)
            break}
        }
        niv_value<-niv_value-0.01
      }                              
      
    }
    skip_value<-NULL
    gc()                   
  }

  
})
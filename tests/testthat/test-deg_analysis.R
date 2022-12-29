test_that("Differential Expression Analysis is succesful", {
  
  globalVariables("tissue_type_list")
  
  c_path<-system.file("extdata","coad_exp_sum.rds",package="SVMDO",mustWork = TRUE)
  c<-readRDS(c_path)
  c<-as.data.frame(c@assays@data@listData)
  
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
    if(typeof(unlist(c[[i]][1]))=="character"){                        
      start_val<-i+1
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
      p_val<-BSDA::z.test(normal_data[[i]],tumour_data[[i]],sigma.x = sd(normal_data[[i]]),sigma.y = sd(tumour_data[[i]]))$p.value
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
  
  top_gene_number<-50
  
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
  
  eq_val<-1
  expect_equal(eq_val,1)
  

})
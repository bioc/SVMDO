#' @title SVMDO
#' @name deg_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @param rawData expression dataset provided from innerServer_exp_server
#' @param rval Selected radio button information provided from innerServer_rad_server
#' @return Server section of differential gene expression analysis


innerServer_3 <- function(input, output, session,rawData, rval) {
  observeEvent(input$initiate_deg_analysis, {
    
    rval_sel<-rval()
    
    message_val<-NULL
    message_val_2<-NULL
    
    if (rval_sel=="COAD") {
      package_path<-find.package("SVMDO",quiet=TRUE)
      rds_path<-"extdata/coad_exp_sum.rds"
      comb_path<-file.path(package_path,rds_path,fsep = "/")
      c<-readRDS(comb_path)
      c<-assay(c)
      type_num_data<-seq.int(ncol(c))
      type_data<-c[1,type_num_data]
    }else if (rval_sel=="LUSC"){
      package_path<-find.package("SVMDO",quiet=TRUE)
      rds_path<-"extdata/lusc_exp_sum.rds"
      comb_path<-file.path(package_path,rds_path,fsep = "/")
      c<-readRDS(comb_path)
      c<-assay(c)
      type_num_data<-seq.int(ncol(c))
      type_data<-c[1,type_num_data]
    }else{
      c<-rawData()
      type_num_data<-seq.int(ncol(c))
      type_data<-c[1,..type_num_data]
    }
    
    type_apply_1<-lapply(type_data,is.character)
    type_loc_1<-as.numeric(min((which(type_apply_1==TRUE))))
    
    if (!is.character(c[[type_loc_1+1]])) {
      
      type_loc_2<-type_loc_1+1
      col_val_1<-type_loc_1
      col_val_2<-type_loc_2
      
      assigning_tissue_type_list<-as.data.frame(c[[col_val_1]])
      colnames(assigning_tissue_type_list)<-c("tissue_type")
      assign("tissue_type_list",assigning_tissue_type_list,envir =.GlobalEnv)
      
      start_val<-type_loc_1
      num_data<-start_val+1
    }else{
      type_loc_2<-type_loc_1+2
      col_val_1<-type_loc_1
      col_val_2<-type_loc_2
      
      assigning_tcga_id_list<-as.data.frame(c[[col_val_1]])
      assigning_tissue_type_list<-as.data.frame(c[[col_val_1+1]])
      
      colnames(assigning_tcga_id_list)<-data.frame("id")
      colnames(assigning_tissue_type_list)<-c("tissue_type")
      
      assign("tissue_type_list",assigning_tissue_type_list,envir =.GlobalEnv)
      assign("tcga_id_list",assigning_tcga_id_list,envir =.GlobalEnv)
      
      start_val<-type_loc_1+1
      num_data<-start_val+2
      
    }
    
    c<-c[,num_data:ncol(c)]
    c<-cbind(tissue_type_list,c)
    
    normal_data<-NULL
    tumour_data<-NULL
    p_val_total<-NULL
    normal_data_mean<-NULL
    tumour_data_mean<-NULL
    int_col<-NULL
    dist_normal_pval<-NULL
    dist_tumour_pval<-NULL
    
    colnames(c)<-(make.unique(colnames(c)))
    
    #####Modifying illegal symbols
    colnames(c)<-gsub( "-", "__",  colnames(c))
    colnames(c)<-gsub( " ", "__",  colnames(c))
    colnames(c)<-gsub( ";. ", "__",colnames(c))
    colnames(c)<-gsub( ";__","__", colnames(c))
    #####
    
    ############ Normal and Tumour Data Extraction and Fold Change Calculation
    
    normal_data<-c[grep("Nor",c$tissue_type),]
    tumour_data<-c[grep("Tum",c$tissue_type),]
    
    num_data_norm<-(seq.int(start_val,ncol(normal_data)))
    num_data_tum<-(seq.int(start_val,ncol(tumour_data)))
    normal_data_mean<-colMeans(normal_data[,num_data_norm])
    tumour_data_mean<-colMeans(tumour_data[,num_data_tum])
    
    fold_change<-tumour_data_mean/normal_data_mean
    matr_fold<-as.matrix(fold_change)
    
    i<-start_val
    u<-start_val
    s<-start_val
    
    dist_normal_pval<-ad.test(normal_data[[i]][seq.int(ncol(normal_data))])$p.value
    dist_tumour_pval<-ad.test(tumour_data[[i]][seq.int(ncol(tumour_data))])$p.value
    
    repeat{
      if (!is.null(dist_normal_pval) & (dist_normal_pval>=0.05 | dist_tumour_pval>=0.05)) {
        dist_normal_pval<-ad.test(normal_data[[i]][seq.int(ncol(normal_data))])$p.value
        dist_tumour_pval<-ad.test(tumour_data[[i]][seq.int(ncol(tumour_data))])$p.value
        if (i==ncol(normal_data)) {
          p_val<-z.test(normal_data[[s]],tumour_data[[s]],sigma.x = sd(normal_data[[s]]),sigma.y = sd(tumour_data[[s]]))$p.value
          p_val_total<-c(p_val_total,p_val)
          s<-s+1
        }else{
          i<-i+1
        }
      }
      
      if (dist_normal_pval<0.05 | dist_tumour_pval<0.05) {
        p_val<-wilcox.test(as.numeric(unlist(normal_data[[u]])),as.numeric(unlist(tumour_data[[u]])))$p.value
        p_val_total<-c(p_val_total,p_val)
        u<-u+1
      }
      if (u>ncol(normal_data) | s>ncol(normal_data)) {
        break
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
    
    assign("sorted_new_bound_form_A",sorted_new_bound_form_A,envir =.GlobalEnv)
    assign("sorted_new_bound_form_B",sorted_new_bound_form_B,envir =.GlobalEnv)
    assign("complete_deg_gene_list",complete_deg_gene_list,envir =.GlobalEnv)
    assign("total_exp_dataset",c, envir =.GlobalEnv)
    
    if (rval_sel %in% c("COAD","LUSC")) {
      top_gene_number<-50
      if (nrow(sorted_new_bound_form_A) < top_gene_number | nrow(sorted_new_bound_form_B) < top_gene_number) {
        message_val_2<-1
        max_down_genes<-sorted_new_bound_form_A
        max_up_genes<-sorted_new_bound_form_B
        message("Insufficient gene number All up/downregulated genes were selected")
        message("Enter a lower value of input size")
      }else{
        max_down_genes<-sorted_new_bound_form_A[seq.int((nrow(sorted_new_bound_form_A)-(top_gene_number-1)),nrow(sorted_new_bound_form_A)),]
        max_up_genes<-sorted_new_bound_form_B[seq.int((nrow(sorted_new_bound_form_B)-(top_gene_number-1)),nrow(sorted_new_bound_form_B)),]
      }
      top_combined_genes<-rbind(max_down_genes,max_up_genes)
      rownames(top_combined_genes)<-top_combined_genes[,1]
      changed_whole_data<-subset(c,select=top_combined_genes$Genes)    
      assign("top_genes_test",changed_whole_data,envir =.GlobalEnv)
      message_val_2<-1
    }else{
      message_val<-1
    }
    
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
    
    if (!is.null(message_val_2)) {
      showModal(
        modalDialog(
          title = "Test Data-Based DEG Analysis Result",
          "Process Completed",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
  })
}
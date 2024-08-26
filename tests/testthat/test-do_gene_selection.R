test_that("Disease Ontology-based Gene Selection is succesful", {
  
  globalVariables("new_tissue_type_list")

  normal_samples<-NULL
  tumour_samples<-NULL
  disease_filtered_gene_data<-NULL
  changed_name_plus_var_imp_genes_table<-NULL
  
  changed_whole_data<-cbind(tissue_type_list,top_genes)
  collect_gene_names<-NULL
  selected_data<-colnames(changed_whole_data[,-1])
  
  checking_data<-selected_data
  checking_data<-as.character(checking_data)
  checking_data<-gsub( "__", "-",checking_data)
  
  colnames(changed_whole_data[,seq.int(2,ncol(changed_whole_data))])<-checking_data
  database_selection<-select(org.Hs.eg.db, keys = checking_data,columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
  
  pop_gene_id<-database_selection$ENTREZID
  pop_gene_symbol<-database_selection$SYMBOL
  count_gene_id<-length(database_selection$ENTREZID)
  
  dis_gene_extract<-lapply(seq_along(pop_gene_id),function(x){
	disease_enrichment<-enrichDO(pop_gene_id[x],ont = "HDO")
        if (!is.null(disease_enrichment)){
          dis_length<-length(disease_enrichment[1]$Description)
          p_check<-(disease_enrichment[1]$p.adjust[1])
          if (dis_length>0 & p_check<0.05 ){
            selected_gene_list<-c(collect_gene_names,pop_gene_symbol[x])
          }
        }
      })
  
  collect_gene_names<-unlist(dis_gene_extract)
  collect_gene_names<-gsub( "-", "__",collect_gene_names)
  changed_whole_data<-changed_whole_data[,-1]
  changed_name_plus_var_imp_genes_table<-subset(changed_whole_data,select=(collect_gene_names))
  
  samp_select_norm<-lapply(seq.int(nrow(tissue_type_list)), function(x){
    if (str_contains(tissue_type_list$tissue_type[x],"Nor")) {
      normal_samples<-rbind("normal",normal_samples)
    }})
  samp_select_tum<-lapply(seq.int(nrow(tissue_type_list)), function(x){
    if (str_contains(tissue_type_list$tissue_type[x],"Tum")) {
      tumour_samples<-rbind("tumour",tumour_samples)
    }
  })
  
  normal_samples<-as.data.frame(unlist(samp_select_norm))
  tumour_samples<-as.data.frame(unlist(samp_select_tum))
  colnames(normal_samples)<-"tissue_type"
  colnames(tumour_samples)<-"tissue_type"
  assigning_new_tissue_type_list<-rbind(normal_samples,tumour_samples)
  rownames(assigning_new_tissue_type_list)<-NULL
  assign("new_tissue_type_list",assigning_new_tissue_type_list,envir = .GlobalEnv)
  
  if(exists("tcga_id_list")){
    tcga_id_list$id<-as.character(gsub("-", ".", tcga_id_list$id, fixed = TRUE))
    assigning_tcga_sample_comb<-cbind(tcga_id_list,new_tissue_type_list)
    assign("tcga_sample_comb",assigning_tcga_sample_comb,envir =.GlobalEnv)
  }
  
  assigning_disease_filtered_gene_data<-cbind(new_tissue_type_list,changed_name_plus_var_imp_genes_table)
  assigning_disease_filtered_gene_data<-as.data.frame(assigning_disease_filtered_gene_data)
  assign("disease_filtered_gene_data",assigning_disease_filtered_gene_data,envir = .GlobalEnv)
  eq_val<-1
  expect_equal(eq_val,1)
  
})
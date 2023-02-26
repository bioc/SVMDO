#' @title SVMDO
#' @name globals
#' @return Including script files and global variables of GUI required to be 
#' initiated at the runApp file execution


source(file.path('R/table_server.R'))
source(file.path('R/table_ui.R'))

source(file.path('R/disc_gene_dw_server.R'))
source(file.path('R/disc_gene_dw_ui.R'))

source(file.path('R/directory_selection_server.R'))
source(file.path('R/directory_selection_ui.R'))

source(file.path('R/expression_dataset_input_server.R'))
source(file.path('R/expression_dataset_input_ui.R'))

source(file.path('R/test_data_selection_server.R'))
source(file.path('R/test_data_selection_ui.R'))

source(file.path('R/clinic_data_input_server.R'))
source(file.path('R/clinic_data_input_ui.R'))

source(file.path('R/disc_gene_dw_server.R'))
source(file.path('R/disc_gene_dw_ui.R'))

source(file.path('R/plot_list_server.R'))
source(file.path('R/plot_list_ui.R'))

source(file.path('R/plot_push_server.R'))
source(file.path('R/plot_push_ui.R'))

source(file.path('R/plot_show_server.R'))
source(file.path('R/plot_show_ui.R'))

source(file.path('R/surv_plot_dw_server.R'))
source(file.path('R/surv_plot_dw_ui.R'))

source(file.path('R/gene_list_table_visualization_ui.R'))
source(file.path('R/package_req_list.R'))

source(file.path('R/deg_server.R'))
source(file.path('R/deg_ui.R'))

source(file.path('R/top_val_server.R'))
source(file.path('R/top_val_ui.R'))

source(file.path('R/top_val_based_deg_filtration_server.R'))
source(file.path('R/top_val_based_deg_filtration_ui.R'))

source(file.path('R/do_based_gene_filtration_server.R'))
source(file.path('R/do_based_gene_filtration_ui.R'))

source(file.path('R/classification_server.R'))
source(file.path('R/classification_ui.R'))

source(file.path('R/survival_analysis_server.R'))
source(file.path('R/survival_analysis_ui.R'))

source(file.path('R/gui_obj_removal_server.R'))
source(file.path('R/gui_obj_removal_ui.R'))

globalVariables(c("final_discriminative_gene_set","top_genes","tcga_id_list", 
                  "sorted_new_bound_form_A","sorted_new_bound_form_B", "hr_list",
                  "fit1_val","max_plots","fit_list","modulename_list","p_list",
                  "..elected_val","..num_data", 
                  "..type_num_data","tissue_type_list","disease_filtered_gene_data",
                  "tcga_sample_comb","new_tissue_type_list","tissue_type",
                  "top_genes_test","total_exp_dataset"))
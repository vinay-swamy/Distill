library(tidyverse)
library(arrow)
library(glue)
load('clean_data/model_data_2018_08_01.Rdata')
load('data/master/raw_data_2018_08_08.Rdata')
load('clean_data/VPaC_12mtry_v11.Rdata')
load('clean_data/ogvfb_exome_cohort_2018_08_07.Rdata')
load('clean_data/colombia.Rdata')


general_train_data <-  model_data$ML_set__general_TT$train_set %>% as.data.frame
general_test_data <-  model_data$ML_set__general_TT$test_set %>% as.data.frame
eye_train_data <- model_data$ML_set__eye_TT$train_set %>% as.data.frame
eye_test_data <- model_data$ML_set__eye_TT$test_set %>% as.data.frame
other_train_data <- model_data$ML_set__other_TT$train_set
other_test_data <- model_data$ML_set__other_TT$test_set
write_feather(general_train_data, 'feather_data/raw/ml_set_general_train.ftr')
write_feather(general_test_data, 'feather_data/raw/ml_set_general_test.ftr')
write_feather(eye_train_data, 'feather_data/raw/ml_set_eye_train.ftr')
write_feather(eye_test_data, 'feather_data/raw/ml_set_eye_test.ftr')
write_feather(other_train_data, 'feather_data/raw/ml_set_other_train.ftr')
write_feather(other_test_data, 'feather_data/raw/ml_set_other_test.ftr')
write_feather(raw_data, 'feather_data/raw/raw_data.ftr')
write_feather(colombia_out, 'feather_data/raw/colombia.ftr')
write_feather( ogvfb_ML_set, 'feather_data/raw/ogvfb.ftr')
write_feather(model_data$Test_set__UK10K, 'feather_data/raw/uk10k.ftr')

panel <- readxl::read_excel('data/NISC100_Variant_Interpretation_June01_2018.xlsx')
ddl_path_cdot <- panel %>% filter(grepl('Path', `Interpretation Summary`, ignore.case = T)) %>% select(`#Chr`, End, Ref, Alt, avsnp147) %>% mutate(pos_id=paste0(`#Chr`, ':', End, '_', Ref, '_', Alt))
##############
### UK10K ####
##############
metadata <- readxl::read_excel(path='EGA_EGAD00001002656_NGS_reanalyze/data/1-s2.0-S0002929716305274-mmc3.xlsx') %>% mutate(Sample=Patient)
sample_gene_comp_het <- metadata %>% filter(Status=='Solved'  & Variant_HGVSc!='NA' & GT=='0/1') %>% group_by(Sample, Gene) %>% summarise(Count=n()) %>% filter(Count>1) 
metadata <- left_join(metadata, sample_gene_comp_het) %>% 
  mutate(Comp_Het_Path = case_when(Count >= 2 ~ 'CH', 
                                   TRUE ~ 'No')) %>% 
  select(-Count)
eyeIntegration <- c('eyeintegration_rnaseq_adipose_subcutaneous','eyeintegration_rnaseq_tpm_adipose_visceral_omentum','eyeintegration_rnaseq_tpm_adrenalgland','eyeintegration_rnaseq_tpm_artery_aorta','eyeintegration_rnaseq_tpm_artery_coronary','eyeintegration_rnaseq_tpm_artery_tibial','eyeintegration_rnaseq_tpm_brain_amygdala','eyeintegration_rnaseq_tpm_brain_anteriorcingulatecortex_ba24','eyeintegration_rnaseq_tpm_brain_caudate_basalganglia','eyeintegration_rnaseq_tpm_brain_cerebellarhemisphere','eyeintegration_rnaseq_tpm_brain_cerebellum','eyeintegration_rnaseq_tpm_brain_cortex','eyeintegration_rnaseq_tpm_brain_frontalcortex_ba9','eyeintegration_rnaseq_tpm_brain_hippocampus','eyeintegration_rnaseq_tpm_brain_hypothalamus','eyeintegration_rnaseq_tpm_brain_nucleusaccumbens_basalganglia','eyeintegration_rnaseq_tpm_brain_putamen_basalganglia','eyeintegration_rnaseq_tpm_brain_spinalcord_cervicalc_1','eyeintegration_rnaseq_tpm_brain_substantianigra','eyeintegration_rnaseq_tpm_breast_mammarytissue','eyeintegration_rnaseq_tpm_cells_ebv_transformedlymphocytes','eyeintegration_rnaseq_tpm_cells_transformedfibroblasts','eyeintegration_rnaseq_tpm_colon_sigmoid','eyeintegration_rnaseq_tpm_colon_transverse','eyeintegration_rnaseq_tpm_esc_stemcellline','eyeintegration_rnaseq_tpm_esophagus_gastroesophagealjunction','eyeintegration_rnaseq_tpm_esophagus_mucosa','eyeintegration_rnaseq_tpm_esophagus_muscularis','eyeintegration_rnaseq_tpm_heart_atrialappendage','eyeintegration_rnaseq_tpm_heart_leftventricle','eyeintegration_rnaseq_tpm_kidney_cortex','eyeintegration_rnaseq_tpm_liver','eyeintegration_rnaseq_tpm_lung','eyeintegration_rnaseq_tpm_minorsalivarygland','eyeintegration_rnaseq_tpm_muscle_skeletal','eyeintegration_rnaseq_tpm_nerve_tibial','eyeintegration_rnaseq_tpm_pancreas','eyeintegration_rnaseq_tpm_pituitary','eyeintegration_rnaseq_tpm_skin_notsunexposed_suprapubic','eyeintegration_rnaseq_tpm_skin_sunexposed_lowerleg','eyeintegration_rnaseq_tpm_smallintestine_terminalileum','eyeintegration_rnaseq_tpm_spleen','eyeintegration_rnaseq_tpm_stomach','eyeintegration_rnaseq_tpm_thyroid','eyeintegration_rnaseq_tpm_wholeblood')
numeric_predictors <- unique(c('ccr_pct_v1','cadd_raw','vest3_rankscore','cadd_phred','mis_z','pli','lof_z','phylop_100way','revel','hapmap2','hapmap1','n_mis','epilogos_quies','n_lof','precessive','pnull','adj_exp_lof','adj_exp_syn','dann','adj_exp_mis','syn_z','n_syn','epilogos_txwk','fitcons','m_cap_score','m_cap_rankscore','eigen_phred','eigen_raw','epilogos_tx','is_lof','eigen_pc_raw_rankscore','epilogos_reprpcwk','fathmm_mkl_coding_rankscore','metalr_score','fathmm_mkl_coding_score','metalr_rankscore','impact_severity','metasvm_rankscore','metasvm_score','epilogos_enh','genocanyon_score','fathmm_converted_rankscore','mpc','epilogos_enhg','af_exac_all','epilogos_reprpc','max_aaf_all','mutationassessor_score','gerp','polyphen_score','gerp_elements','mutationassessor_score_rankscore','stam_mean','an_exac_all','af_exac_nfe','provean_converted_rankscore','an_exac_nfe','lrt_score','lrt_omega','grantham','lrt_converted_rankscore','genocanyon_score_rankscore','an_exac_afr','an_exac_amr','an_exac_sas','epilogos_het','ac_exac_all','linsight','gno_an_popmax','exac_num_het','an_exac_eas','gno_an_all','ac_exac_nfe','mutationtaster_converted_rankscore','an_exac_oth','an_exac_fin','gno_an_nfe','gno_af_all','gno_an_afr','epilogos_tssaflnk','gno_af_popmax','epilogos_znf','segway_sum_score','aaf_esp_ea','epilogos_txflnk','provean_score','segway_mean_score','epilogos_tss','aaf_esp_all','af_exac_amr','gno_af_nfe','epilogos_enhbiv','af_exac_sas','sift_score','fathmm_score','ac_exac_amr','aaf_esp_aa','gno_ac_all','gno_af_afr','ac_exac_sas','af_exac_eas','gno_an_fin','af_exac_afr','gno_an_eas','gno_an_oth','gno_ac_nfe','gno_ac_popmax','ac_exac_eas','ac_exac_afr','epilogos_tssbiv','gno_ac_afr','vest3_score','sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01', row.names(VPaC_12mtry_v11$importance)))
nn_predictors <- c('ccr_pct_v1','cadd_raw','vest3_rankscore','cadd_phred','mis_z','pli','lof_z','phylop_100way','revel','hapmap2','hapmap1','n_mis','epilogos_quies','n_lof','precessive','pnull','adj_exp_lof','adj_exp_syn','dann','adj_exp_mis','syn_z','n_syn','epilogos_txwk','fitcons','m_cap_score','m_cap_rankscore','eigen_phred','eigen_raw','epilogos_tx','is_lof','eigen_pc_raw_rankscore','epilogos_reprpcwk','fathmm_mkl_coding_rankscore','metalr_score','fathmm_mkl_coding_score','metalr_rankscore','impact_severity','metasvm_rankscore','metasvm_score','epilogos_enh','genocanyon_score','fathmm_converted_rankscore','mpc','epilogos_enhg','af_exac_all','epilogos_reprpc','max_aaf_all','mutationassessor_score','gerp','polyphen_score','gerp_elements','mutationassessor_score_rankscore','stam_mean','an_exac_all','af_exac_nfe','provean_converted_rankscore','an_exac_nfe','lrt_score','lrt_omega','grantham','lrt_converted_rankscore','genocanyon_score_rankscore','an_exac_afr','an_exac_amr','an_exac_sas','epilogos_het','ac_exac_all','linsight','gno_an_popmax','exac_num_het','an_exac_eas','gno_an_all','ac_exac_nfe','mutationtaster_converted_rankscore','an_exac_oth','an_exac_fin','gno_an_nfe','gno_af_all','gno_an_afr','epilogos_tssaflnk','gno_af_popmax','epilogos_znf','segway_sum_score','aaf_esp_ea','epilogos_txflnk','provean_score','segway_mean_score','epilogos_tss','aaf_esp_all','af_exac_amr','gno_af_nfe','epilogos_enhbiv','af_exac_sas','sift_score','fathmm_score','ac_exac_amr','aaf_esp_aa','gno_ac_all','gno_af_afr','ac_exac_sas','af_exac_eas','gno_an_fin','af_exac_afr','gno_an_eas','gno_an_oth','gno_ac_nfe','gno_ac_popmax','ac_exac_eas','ac_exac_afr','epilogos_tssbiv','gno_ac_afr','vest3_score','sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01')
# random forest was initially trained on eyeIntegration and nuemric predictors
# nn pred is a subset of numeric preds
# first, identify duplicated variants 
id_cols <- c("variant_id", 'chrom', 'start', 'end', 'vcf_id', 'ref', 'alt', 'type', 'sub_type', 'Status')
model_data_ids <- lapply(c("ML_set__eye_TT", "ML_set__general_TT", "ML_set__other_TT"), 
                         function(x) bind_rows( select(model_data[[x]][['train_set']], any_of(id_cols)) %>% 
                                                  mutate(origin = glue('{x}_trainset') )  ,  
                                                select(model_data[[x]][['test_set']], any_of(id_cols)) %>% 
                                                  mutate(origin = glue('{x}_testset')) 
                         ) ) %>% bind_rows 

all_variants <- bind_rows(model_data_ids %>% mutate_all(as.character), 
                             select(model_data$Test_set__UK10K, any_of(id_cols)) %>% mutate_all(as.character) %>%
                               mutate(origin = 'uk10k'),
                             colombia_out %>% select(any_of(id_cols)) %>% mutate_all(as.character) %>% 
                               mutate(origin = 'colombia'), 
                             raw_data %>% select(any_of(id_cols)) %>% mutate_all(as.character) %>% 
                               mutate(origin = 'raw'),
                             ogvfb_ML_set %>% select(any_of(id_cols)) %>% mutate_all(as.character) %>%
                               mutate(origin = 'ogvfb')
) %>% distinct 

all_variants_distinct <-  all_variants %>% select(-origin) %>% distinct %>% mutate(var_id = paste('var_', 1:nrow(.)))
all_variant_ids <- left_join(all_variants, all_variants_distinct)

all_variant_duplicated <- all_variant_ids %>% 
  filter(duplicated(var_id)) %>% 
  select(var_id) %>% distinct %>%
  inner_join(all_variant_ids)
all_variant_duplicated_wide <- all_variant_duplicated %>% mutate(val =T) %>% 
  spread(key = origin, value = val, fill = F)
# no duplicated in main ML train and test, thats good, 
# and the ML eye set is a subset of the main one, so ignore those, so only real dup issues are in uk10k, ogvfb and raw
sum(all_variant_duplicated_wide$ML_set__general_TT_trainset & all_variant_duplicated_wide$ML_set__general_TT_testset) 

#th entirety of ml_set_other is contained within raw_data, so going to drop that 

main_set <- filter(all_variants, origin %in% c("ML_set__general_TT_trainset","ML_set__general_TT_testset")) %>% 
  select(-origin)

raw_data_no_dups <-  filter(all_variants, origin == 'raw') %>% 
  select(-Status) %>% 
  anti_join(main_set) %>% 
  select(-origin) %>% 
  distinct 

raw_data_clean <- raw_data_no_dups %>%   inner_join(raw_data, .)
uk10k_no_dups <- filter(all_variants, origin == 'uk10k') %>% 
  select(-Status) %>%
  anti_join(main_set) %>% 
  anti_join(raw_data_no_dups) %>% 
  select(-origin) %>% 
  mutate(variant_id = as.numeric(variant_id), 
         start = as.numeric(start), 
         end = as.numeric(end)) 
uk10k_clean <- uk10k_no_dups %>% 
  inner_join(model_data$Test_set__UK10K)
#only 15 uk10k samples are not dups; dropping 

ogvfb_nodups <- filter(all_variants, origin == 'ogvfb') %>% 
  select(-Status) %>%
  anti_join(main_set) %>% 
  anti_join(raw_data_no_dups) %>% 
  mutate(variant_id = as.numeric(variant_id), 
         start = as.numeric(start), 
         end = as.numeric(end)) %>% 
  anti_join(uk10k_no_dups) %>% 
  select(-origin) %>% 
  distinct
#all ogvfb candiates are dups
  
##correct cleaned data
allX_with_dups <- raw_data %>% 
  mutate(DataSet_o=DataSet) %>%
  mutate(start = as.numeric(as.character(start))) %>%
  mutate(Variant_genomic = paste0(chrom, ':', start + 1, ref, '>', alt)) %>%
  mutate(DataSet = case_when(DataSet_o == 'ddl_nisc_100_panel' ~ 'DDL NISC RD Cohort',
                             DataSet_o == 'clinvar' & status != 'PATHOGENIC_OTHER' ~ 'ClinVar HC',
                             DataSet_o == 'clinvar' & status == 'PATHOGENIC_OTHER' ~ 'ClinVar LC Path',
                             DataSet_o == 'grimm' & source == 'humvar' ~ 'Grimm HumVar',
                             DataSet_o == 'grimm' & source == 'exovar' ~ 'Grimm ExoVar',
                             DataSet_o == 'grimm' & source == 'exovar__humvar' ~ 'Grimm ExoVar/HumVar',
                             DataSet_o == 'grimm' & source == 'predictSNP' ~ 'Grimm PredictSNP',
                             DataSet_o == 'grimm' & source == 'swissvar' ~ 'Grimm SwissVar',
                             DataSet_o == 'grimm' & source == 'varibench' ~ 'Grimm VariBench',
                             DataSet_o == 'wellderly' ~ 'Wellderly',
                             grepl('homsy', DataSet_o) ~ 'Homsy',
                             grepl('unifun', DataSet_o) ~ 'UniFun',
                             grepl('samocha', DataSet_o) ~ 'Samocha',
                             DataSet_o == 'gnomad' & (pos_id %in% model_data$ML_set__general_TT$train_set$pos_id ||
                                                        pos_id %in% model_data$ML_set__general_TT$test_set$pos_id ||
                                                        pos_id %in% model_data$ML_set__other_TT$train_set$pos_id ||
                                                        pos_id %in% model_data$ML_set__other_TT$test_set$pos_id) ~ 'gnomAD Benign',
                             DataSet_o == 'UK10K' ~ 'UK10K',
                             TRUE ~ 'Other')) %>%
  mutate(Status = case_when((DataSet == 'DDL NISC RD Cohort' & pos_id %in% ddl_path_cdot$pos_id) |
                              (DataSet == 'DDL NISC RD Cohort' & end %in% ddl_path_cdot$End)  ~ 'Pathogenic',
                            DataSet == 'UK10K' & (Variant_genomic %in% (metadata %>% filter(Status == 'Solved') %>% pull(Variant_genomic))) ~ 'Pathogenic',
                            DataSet == 'UK10K' & (Variant_genomic %in% (metadata %>% filter(Status == 'Partially solved') %>% pull(Variant_genomic))) ~ 'Maybe Pathogenic',
                            grepl('pathogenic', DataSet_o, ignore.case = T) | grepl('pathogenic', status, ignore.case = T) ~ 'Pathogenic',
                            TRUE ~ 'NotPathogenic')) %>%
  mutate(Status = case_when(Status = grepl('Grimm', DataSet) ~ status,
                            TRUE ~ Status)) %>%
  #mutate_at(vars(one_of(numeric_predictors)), funs(as.numeric(.))) %>%
  #select(-status) %>%
  filter(DataSet != 'Other') %>%
  filter(Status != 'Maybe Pathogenic') %>% 
  select(all_of(id_cols),DataSet, Status)
allX_dups_collapsed <- allX_with_dups %>% 
  group_by_at(all_of(id_cols)) %>% 
  summarise(all_datasets = paste(DataSet, collapse = ';'),
            all_status = paste(Status, collapse = ';'))

allX_clean <- raw_data_clean %>% 
  mutate(DataSet_o=DataSet) %>%
  mutate(start = as.numeric(as.character(start))) %>%
  mutate(Variant_genomic = paste0(chrom, ':', start + 1, ref, '>', alt)) %>%
  mutate_at(vars(one_of(numeric_predictors)), funs(as.numeric(.))) %>% 
  inner_join(allX_dups_collapsed)

### VS: this fails  because of new tidyverse type coercion rules, need to add as.data.frame after all tidyverse code


allX_clean$fitcons_float <- allX_clean$fitcons

allX_clean <- as.data.frame(allX_clean)
allX_clean[is.na(allX_clean)] <- -1



##### next, subset to all common predictors 
id_cols <- c("variant_id", 'chrom', 'start', 'end', 'vcf_id', 'ref', 'alt', 'type', 'sub_type', 'Status')
####################
common_non_numeric_cols <- reduce( list(colnames(general_test_data), 
         colnames(colombia_out),
         colnames(allX_clean)), intersect ) %>%
  setdiff(numeric_predictors ) %>% 
  setdiff(eyeIntegration) %>% 
  setdiff(id_cols) %>% 
  .[!grepl('eyeintegration',.)]
count_missing <- function(df){
  df [is.na(df)] <- ''
  df %>% select(all_of(common_non_numeric_cols)) %>% apply(2, function(x) sum(x=='')) %>% {./nrow(df)}
}

missing_df <- tibble(cols= common_non_numeric_cols,
                     dtype = apply(allX_clean %>% select(all_of(common_non_numeric_cols)),2,class ),
                     nuniq = apply(general_train_data %>% select(all_of(common_non_numeric_cols)),2,n_distinct),
                     train = count_missing(general_train_data), 
                     test = count_missing(general_test_data), 
                     allx = count_missing(allX_clean), 
                     colombia = count_missing(colombia_out))
target_categorical_preds <- missing_df %>% filter(nuniq >1, nuniq <15, 
                                                  train<=.01, test <=.01, allx <=.01, colombia <=.01) %>% 
  pull(cols)

general_train_data_all <-  model_data$ML_set__general_TT$train_set %>% as.data.frame %>% 
  select(all_of(id_cols), all_of(numeric_predictors), all_of(eyeIntegration), all_of(target_categorical_preds)) %>% 
  mutate_at(vars(all_of(numeric_predictors)), funs(as.numeric(.))) %>% as.data.frame
general_train_data_all[is.na(general_train_data_all)] <- -1

general_test_data_all  <-  model_data$ML_set__general_TT$test_set %>% as.data.frame %>% 
  select(all_of(id_cols), all_of(numeric_predictors), all_of(eyeIntegration), all_of(target_categorical_preds)) %>% 
  mutate_at(vars(all_of(numeric_predictors)), funs(as.numeric(.)))
general_test_data_all[is.na(general_test_data_all)] <- -1

allX_all <- allX_clean %>% 
  select(all_of(id_cols), all_of(numeric_predictors), all_datasets, all_status,  all_of(eyeIntegration), all_of(target_categorical_preds))%>% 
  mutate_at(vars(all_of(numeric_predictors)), funs(as.numeric(.)))
allX_all[is.na(allX_all)] <- -1

colombia_all <- colombia_out %>% 
  select(all_of(id_cols), all_of(numeric_predictors), all_of(eyeIntegration), all_of(target_categorical_preds))%>% 
  mutate_at(vars(all_of(numeric_predictors)), funs(as.numeric(.)))
colombia_all[is.na(colombia_all)] <- -1

write_feather(general_train_data_all, 'feather_data/clean/ml_set_general_train.ftr')
write_feather(general_test_data_all, 'feather_data/clean/ml_set_general_test.ftr')
write_feather(allX_all, 'feather_data/clean/allX.ftr')
write_feather(colombia_all, 'feather_data/clean/colombia.ftr')


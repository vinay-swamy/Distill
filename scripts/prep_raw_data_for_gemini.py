'''
The goal of this script is to take all the Rmd's david made for prepping each sample, and merge everything into one step

'''
#%%
import pandas as pd 
import os 
import numpy as np
os.chdir('/data/swamyvs/Distill')
#%%
##########################
# clinvar_processing.Rmd #
##########################

clinvar_single = (pd.read_csv('data/raw_prep/clinvar_alleles.single.b37.tsv.gz', sep = '\t')
                  .assign(conflicted = lambda x:  x.conflicted.astype(str),
                          gold_stars = lambda x:  x.gold_stars.replace('-', '0').astype(float),
                          all_traits = lambda x:  x.all_traits.str.lower().fillna(''), 
                          xrefs = lambda x:  x.xrefs.str.lower())
                  )

pattern = r'stargardt|retina|leber|usher|cone-rod|rod-cone|macula|retinitis|eye|cornea'
path_eye = (clinvar_single
            .query('pathogenic > 0 | likely_pathogenic > 0')
            .query('conflicted  == "0" & gold_stars > 0')
            .pipe( lambda x: x[ (x.all_traits.str.contains(pattern, case = False)) |
                                (x.all_traits.str.contains(pattern, case = False))])
            .pipe(lambda x: x[~(x.all_traits.str.contains(r'cardio|cancer|carci|lynch', case = False))])
            .assign(ID= '.', 
                    QUAL = '.', 
                    FILTER = '.', 
                    INFO='STATUS=PATHOGENIC_EYE', 
                    FORMAT='GT', 
                    FAUX='0/1')
            )
#NOTE: 13 of samples that pass in the R code fail here, but looks like most of them also have a cancer diagnosis that should have failed the cancer regex



benign_eye = (clinvar_single
              .query('benign >0')
              .query('conflicted == "0"')
              .pipe(lambda x: x[~x.all_traits.str.contains(r'cardio|cancer|carci|lynch', 
                                                             case = False)])
              .pipe( lambda x: x[ (x.all_traits.str.contains(pattern, case = False)) |
                                  (x.all_traits.str.contains(pattern, case = False))])
              .assign(ID= '.', 
                    QUAL = '.', 
                    FILTER = '.', 
                    INFO='STATUS=PATHOGENIC_EYE', 
                    FORMAT='GT', 
                    FAUX='0/1')             
             )
#NOTE: 13 of samples that pass in the R code fail here, but looks like most of them also have a cancer diagnosis that should have failed the cancer regex(yes i copy and pasted this )

benign_all = (clinvar_single
              .pipe(lambda x: x[~x.all_traits.str.contains(r'cardio|cancer|carci|lynch', 
                                                             case = False)])
              .query('benign >0')
              .pipe(lambda x: x[~x.variation_id.isin(benign_eye['variation_id'])])
              .query('conflicted == "0"')
              .assign(ID= '.', 
                    QUAL = '.', 
                    FILTER = '.', 
                    INFO='STATUS=PATHOGENIC_EYE', 
                    FORMAT='GT', 
                    FAUX='0/1')  
            )
#NOTE: you already know 



out = pd.concat([path_eye, benign_eye, benign_all])
out.to_csv('data/metadata/clinvar_metadata.csv.gz', index =False)



header = '''##fileformat=VCFv4.2
##source=ClinVar
##INFO=<ID=STATUS,Number=1,Type=String,Description=\"Non-conflicting RD pathogenic or benign in ClinVar\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFAUX
'''
with open('data/vcfs/clinvar_RD.vcf', 'w+') as vcf:
        vcf.write(header)
(out
.filter(['chrom', 'pos', 'ID', 'ref', 'alt', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'FAUX'])
.sort_values(['chrom', 'pos'])
.to_csv('data/vcfs/clinvar_RD.vcf', mode = 'a', index=False, sep = '\t', header = False)
)
# %%       
###################
# build_grimm.Rmd #
###################

urls = ['http://structure.bmc.lu.se/VariBench/humvar_filtered_tool_scores.csv',
'http://structure.bmc.lu.se/VariBench/exovar_filtered_tool_scores.csv',
'http://structure.bmc.lu.se/VariBench/varibench_selected_tool_scores.csv',
'http://structure.bmc.lu.se/VariBench/predictSNP_selected_tool_scores.csv',
'http://structure.bmc.lu.se/VariBench/swissvar_selected_tool_scores.csv']
grimm_dfs = [pd.read_csv(u) for u in urls]

# %%
grimm_data = (pd.concat(grimm_dfs)
              .assign(Status = lambda x: np.where(x['True Label'] > 0, 'Pathogenic','NotPathogenic' ),
                      CHR = lambda x: x.CHR.str.replace('chr', '')
              .pipe(lambda x: x[x.CHR.str.contains(r'PATCH|MHC|LRC|GL')])
              )


# %%
# grimm_data_collapse <- grimm_data %>% 
#   group_by(CHR, `Nuc-Pos`, `REF-Nuc`, `ALT-Nuc`) %>% 
#   summarise(Status = paste(paste0('Status=',paste(unique(Status), collapse='__'))),
#             Source = paste(paste0('Source=',paste(unique(Source), collapse='__')))) %>% 
#   mutate(INFO = paste(Status, Source, sep=';')) %>% 
#   filter(!grepl('__', Status)) %>%   # remove conflicting
#   ungroup()
def join_col(s):
        return 'Status=' + '__'.join(s.unique())
# grimm_data_collapse = (grimm_data
#                         .groupby('CHR', 'Nuc-Pos', 'REF-Nuc', 'ALT-Nuc')
#                         .
 
#                      )


# %%

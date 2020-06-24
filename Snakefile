
rule all:
    input:
        'clean_data/assess_2018_08_27.Rdata' 
    shell:
        '''
        echo
        '''


rule build_grimm:
    output:
        'data/grimm/grimm.PED_faux.gemini.db', 
        'clean_data/grimm_ML.Rdata'
    shell:
        '''
        echo
        '''



rule ddl_process:
    output:
        'clean_data/ddl_nisc_panel_variants.Rdata',
        'data/ddl_nisc_100_panel/DDL_NISC_targeted_panel.PED_ddl_nisc.gemini.db',
    shell:
        '''
        echo
        '''


rule build_homsy_and_samocha:
    output:
        'clean_data/homsy_ML.Rdata',
        'clean_data/samocha_ML.Rdata',
        'data/homsy_benign/homsy.benign.PED_faux.gemini.db',
        'data/homsy_pathogenic/homsy.pathogenic.PED_faux.gemini.db',
        'data/samocha_benign/samocha.benign.PED_faux.gemini.db',
        'data/samocha_pathogenic/samocha.pathogenic.PED_faux.gemini.db'
    shell:
        '''
        echo
        '''

       
rule ogvfb_exomes_cohort:
    output:
        'clean_data/ogvfb_exome_cohort_2018_08_07.Rdata',
        'data/ogvfb/VCFs.GATK.PED_master.gemini.db'
    shell:
        '''
        echo
        '''


rule clinvar_processing:
    output:
        'data/clinvar/clinvar_RD.PED_faux.gemini.db',
         'data/clinvar/clinvar.gemini.tsv.gz' 
    shell:
        '''
        echo
        '''


rule unifun_process:
    output:
        'clean_data/unifun_ML.Rdata',
        'data/unifun_benign/UniFun_benign.fakeS.PED_faux.gemini.db',
        'data/unifun_pathogenic/UniFun_deleterious.fakeS.PED_faux.gemini.db'
    shell:
        '''
        echo
        '''

rule filter_gnomad:
    output:
        'data/gnomad_rare_benign_ish/gnomad.exomes.r2.0.2.sites.maxAF_0.01_20percent.PED_faux.gemini.db',
        'data/gnomad_rare_benign_ish/gnomad.exomes.r2.0.2.sites.maxAF_0.01_20percent.PED_faux.gemini.tsv.gz'
    shell:
        '''
        echo
        '''


rule query_and_process_gemini_database:
    output:
        'data/UK10K/UK10K_EGAD.hets.test.gz'
        'data/UK10K/UK10K_EGAD.homs.test.gz',
        'data/UK10K/UK10K_EGAD.hets.train.gz',
        'data/UK10K/UK10K_EGAD.homs.train.gz'
    shell:
        '''
        echo
        '''

rule build_UK10K:
    input:
        'data/UK10K/UK10K_EGAD.hets.test.gz',
        'data/UK10K/UK10K_EGAD.homs.test.gz',
        'data/UK10K/UK10K_EGAD.hets.train.gz',
        'data/UK10K/UK10K_EGAD.homs.train.gz'
    output:
        'clean_data/uk10k_gemini_rare_variants_2018_07_31.Rdata'
    shell:
        '''
        echo
        '''

rule Build_Data:
    input:
        uk10k_data = 'clean_data/uk10k_gemini_rare_variants_2018_07_31.Rdata',
        clinvar_file = 'data/clinvar/clinvar.gemini.tsv.gz' ,
        gnomad_file = 'data/gnomad_rare_benign_ish/gnomad.exomes.r2.0.2.sites.maxAF_0.01_20percent.PED_faux.gemini.tsv.gz'
    output:
        'clean_data/model_data_2018_08_01.Rdata'
    shell:
        '''
        echo
        '''

filler = ['VPaC__15mtry_v14_961734',  'VPaC__12mtry_v14_652018']

rule Build_Vpac_piece:
    input:
         'clean_data/model_data_2018_08_01.Rdata'
    output:
        expand('clean_data/VPaC_pieces/{rand}.Rdata', rand=filler)
    shell:
        '''
        echo
        '''


rule merge_models:
    input: 
        expand('clean_data/VPaC_pieces/{rand}.Rdata', rand=filler)
    output: 
        'clean_data/VPaC_12mtry_v11.Rdata'
    shell:
        '''
        echo
        '''


rule build_rawR:
    input:
        'data/grimm/grimm.PED_faux.gemini.db',
        'data/homsy_benign/homsy.benign.PED_faux.gemini.db',
        'data/homsy_pathogenic/homsy.pathogenic.PED_faux.gemini.db',
        #data/'ogvfb/VCFs.GATK.PED_master.gemini.db',
        'data/samocha_benign/samocha.benign.PED_faux.gemini.db',
        'data/samocha_pathogenic/samocha.pathogenic.PED_faux.gemini.db',
        'data/unifun_benign/UniFun_benign.fakeS.PED_faux.gemini.db',
        'data/unifun_pathogenic/UniFun_deleterious.fakeS.PED_faux.gemini.db',
        'data/wellderly/wellderly.coding.benign.fakeSample.PED_faux.gemini.db',#*** missing
        'data/ddl_nisc_100_panel/DDL_NISC_targeted_panel.PED_ddl_nisc.gemini.db',
        'data/gnomad_rare_benign_ish/gnomad.exomes.r2.0.2.sites.maxAF_0.01_20percent.PED_faux.gemini.db',
        'data/ogvfb/VCFs.GATK.PED_master.gemini.db',
        'data/clinvar/clinvar_RD.PED_faux.gemini.db'
    output:
        'data/master/raw_data_2018_08_08.Rdata'
    shell:
        '''
        echo
        '''

rule colombia:
    input:
        #/Volumes/data/projects/nei/hufnagel/colombia_cohort/columbia.gemini.het.tsv.gz
        'data/colombia_cohort/columbia.gemini.het.tsv.gz'
    output:
        'clean_data/colombia.Rdata'
    shell:
        '''
        echo
        '''

rule model_assess_builder:
    input:
        'clean_data/model_data_2018_08_01.Rdata', 
        'clean_data/VPaC_12mtry_v11.Rdata',
        'data/master/raw_data_2018_08_08.Rdata', 
        'clean_data/colombia.Rdata',
        'clean_data/ogvfb_exome_cohort_2018_08_07.Rdata'
    output: 
        'clean_data/assess_2018_08_27.Rdata' 
    shell:
        '''
        echo
        '''


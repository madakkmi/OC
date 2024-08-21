/****************************************************************************************************************************************
Basic adjusted models for overall OC and for 4 subtypes for testing menopause interaction

Notes:
1. 
2. 
3. 
*****************************************************************************************************************************************/

log using "OC\reports\basic_models_for_menopause_interaction.smcl", replace
use "OC\data\processed\ukb_OC_data_for_epi_analysis_preprocessed.dta", replace

/****************************************************************************************************************************************
Check and Create Outcome Variables for Subtypes
*****************************************************************************************************************************************/

capture confirm variable subtype_serous_oc
if !_rc {
    display "Outcome variables for subtypes exist, no need to create"
} else {
    display "Outcome variables for subtypes do not exist, now creating them."
    * Generate numeric outcome variables for subtypes
    gen subtype_serous_oc = 0 if ovarian_cancer == 0
    replace subtype_serous_oc = 1 if serous_oc == "incident"

    gen subtype_endometrioid_oc = 0 if ovarian_cancer == 0
    replace subtype_endometrioid_oc = 1 if Endometrioid_oc == "incident"

    gen subtype_clearcell_oc = 0 if ovarian_cancer == 0
    replace subtype_clearcell_oc = 1 if clearcell_oc == "incident"

    gen subtype_mucinous_oc = 0 if ovarian_cancer == 0
    replace subtype_mucinous_oc = 1 if mucinous_oc == "incident"
}

/****************************************************************************************************************************************
Define Variables
*****************************************************************************************************************************************/

/* Adjustment Variables */
local basic_adjust_vars i._x54_centre c._x53_yratcentre c.x21022_age c._x189_townsend i._x21000_ethn
local adjust_vars `basic_adjust_vars' // Later you can add additional adjustments

/* Binary Features */
local binary_features x20004__1355_bo x6141__2_sd x6142__2_retired x2724_mp x1835_mlive x2784_ocp x2674_bcscreen x6142__1_employed x2814_hrt x1797_flive x6179__100_nosuppl x6153__5_ocp x1960_fedup x6164__3_lightdiy x1677_bf x6149__6_denture x20110__101_none x6145__3_deathclose x2100_psych x1920_mood x1940_irrit x2000_worry x20004__1480_wisdom x20004__1507_ep x6162__1_car x6179__1_fishoil

/* Categorical/Ordinal Features */
local categorical_features x1369_beef x1628_alc x1647_cob x1697_hgt10 x2050_depfq x738_hhinc x981_plwalkdur

/* Continuous Features */
local continuous_features x20011_ageop x21022_age x2794_ageocpstart x30770_igf1 x130_pobeast x2734_nliveb x2744_bwtfirst x3064_pef x1807_fdage x30210_ephpct x2804_ageocplast x50_standhgt x30050_mch x20015_sithgt x30040_mcv x21002_wt x20023_mtimematch x129_pobnorth x1070_tv x3137_nmeasure x30870_tg x30720_cycc x48_wc x47_hgsright x30290_hlsrpct x30070_rdw x102_pr x20009_agenonca x23099_bcpct x30650_ast x4079_dbp x30270_mscv x30740_glc x30280_irf x4080_sbp x30710_crp x30120_lymct x46_hgsleft x30620_alt x1050_outtime x30520_uk x30690_cl x49_hip x30750_HbA1c x894_modactdur x699_lenaddr x709_nhh x136_nop x30150_ephct x30190_monopct x30200_nphpct x3062_fvc x30760_hdl x30860_tp

/* Outcome Variables */
local outcome_variables subtype_serous_oc subtype_endometrioid_oc subtype_clearcell_oc subtype_mucinous_oc

/* Interaction Variable */
local interact_var x2724_mp

/****************************************************************************************************************************************
Run Logistic Regression Models
*****************************************************************************************************************************************/

foreach outcome_var in `outcome_variables' {
    set more off
    tempname results
    postfile `results' str30 model str20 exposure str20 var_type n n_cancer beta se lci uci double p beta_interact se_interact lci_interact uci_interact double p_interact using "OC\reports\basic_binary_and_continuous_for_`outcome_var'_for_`interact_var'_interact.dta", replace

    /* Binary Features */
    foreach col in `binary_features' {
        if ("`col'" == "`interact_var'") continue

        logit `outcome_var' `col' `interact_var'#`col' `adjust_vars'
        count if e(sample)
        local n = r(N)
        count if e(sample) & `outcome_var' == 1
        local n_cancer = r(N)
        local beta = _b[`col']
        local se = _se[`col']
        local lci = _b[`col'] - (invnormal(0.975)*_se[`col'])
        local uci = _b[`col'] + (invnormal(0.975)*_se[`col'])
        test `col'
        local p = r(p)

        local b00 = _b[0.`interact_var'#0.`col']
        local b01 = _b[0.`interact_var'#1.`col']
        local b10 = _b[1.`interact_var'#0.`col']
        local b11 = _b[1.`interact_var'#1.`col']

        if (`b00' != 0) {
            local term 0.`interact_var'#0.`col'
        } 
        if (`b01' != 0) {
            local term 0.`interact_var'#1.`col'
        } 
        else if (`b10' != 0) {
            local term 1.`interact_var'#0.`col'
        } 
        else {
            local term 1.`interact_var'#1.`col'
        }

        local beta_interact = _b[`term']
        local se_interact = _se[`term']
        local lci_interact = _b[`term'] - (invnormal(0.975)*_se[`term'])
        local uci_interact = _b[`term'] + (invnormal(0.975)*_se[`term'])
        test `term'
        local p_interact = r(p)
        display `p_interact'

        post `results' ("basic_`interact_var'_interact") ("`col'") ("binary") (`n') (`n_cancer') (`beta') (`se') (`lci') (`uci') (`p') (`beta_interact') (`se_interact') (`lci_interact') (`uci_interact') (`p_interact')
    }

    /* Continuous Features */
    foreach col in `continuous_features' {
        logit `outcome_var' c.`col' `interact_var'#c.`col' `adjust_vars'
        count if e(sample)
        local n = r(N)
        count if e(sample) & `outcome_var' == 1
        local n_cancer = r(N)
        local beta = _b[`col']
        local se = _se[`col']
        local lci = _b[`col'] - (invnormal(0.975)*_se[`col'])
        local uci = _b[`col'] + (invnormal(0.975)*_se[`col'])
        test `col'
        local p = r(p)

        local beta_interact = _b[1.`interact_var'#c.`col']
        local se_interact = _se[1.`interact_var'#c.`col']
        local lci_interact = _b[1.`interact_var'#c.`col'] - (invnormal(0.975)*_se[1.`interact_var'#c.`col'])
        local uci_interact = _b[1.`interact_var'#c.`col'] + (invnormal(0.975)*_se[1.`interact_var'#c.`col'])
        test 1.`interact_var'#c.`col'
        local p_interact = r(p)

        post `results' ("basic_`interact_var'_interact") ("`col'") ("continuous") (`n') (`n_cancer') (`beta') (`se') (`lci') (`uci') (`p') (`beta_interact') (`se_interact') (`lci_interact') (`uci_interact') (`p_interact')
    }

    postclose `results'

    set more off
    tempname results
    postfile `results' str30 model str20 exposure str20 var_type n n_cancer double df double chi2 double p double df_interact double chi2_interact double p_interact using "OC\reports\basic_categorical_for_`outcome_var'_for_`interact_var'_interact.dta", replace

    /* Categorical/Ordinal Features */
    foreach col in `categorical_features' {
        logit `outcome_var' i.`col' `interact_var'#i.`col' `adjust_vars'
        estimates store m1

        logit `outcome_var' `adjust_vars' if e(sample)
        estimates store m0
        lrtest m1 m0, force
        local p = r(p)
        local df = r(df)
        local chi2 = r(chi2)
        local n = e(N)
        count if e(sample) & `outcome_var' == 1
        local n_cancer = r(N)

        logit `outcome_var' i.`col' `adjust_vars' if e(sample)
        estimates store m0
        lrtest m1 m0, force
        local p_interact = r(p)
        local df_interact = r(df)
        local chi2_interact = r(chi2)
        local n = e(N)
        count if e(sample) & `outcome_var' == 1
        local n_cancer = r(N)

        post `results' ("basic_`interact_var'_interact") ("`col'") ("categorical") (`n') (`n_cancer') (`df') (`chi2') (`p') (`df_interact') (`chi2_interact') (`p_interact')
    }

    postclose `results'

    set more off
    tempname results
    postfile `results' str30 model str20 exposure n n_cancer beta se lci uci double p using "OC\reports\basic_categorical_levels_for_`outcome_var'_for_`interact_var'_interact.dta", replace

    /* Categorical Levels */
    foreach col in `categorical_features' {
        logit `outcome_var' i.`col' `interact_var'#i.`col' `adjust_vars'
        levelsof `col' if e(sample), local(levels)

        foreach lev of local levels {
            local beta = _b[`lev'.`col']
            local se = _se[`lev'.`col']
            if `beta' == 0 & `se' == 0 {
                count if `col' == `lev' & e(sample)
                local n = r(N)
                count if `outcome_var' == 1 & `col' == `lev' & e(sample)
                local n_cancer = r(N)
                post `results' ("basic_`interact_var'_interact") ("`col'=`lev'") (`n') (`n_cancer') (0) (0) (0) (0) (.)
            } else {
                count if `col' == `lev' & e(sample)
                local n = r(N)
                count if `outcome_var' == 1 & `col' == `lev' & e(sample)
                local n_cancer = r(N)
                local beta = _b[`lev'.`col']
                local se = _se[`lev'.`col']
                local lci = _b[`lev'.`col'] - (invnormal(0.975)*_se[`lev'.`col'])
                local uci = _b[`lev'.`col'] + (invnormal(0.975)*_se[`lev'.`col'])
                test `lev'.`col'
                local p = r(p)
                post `results' ("basic_`interact_var'_interact") ("`col'=`lev'") (`n') (`n_cancer') (`beta') (`se') (`lci') (`uci') (`p')
            }
        }
    }

    postclose `results'
}

log close

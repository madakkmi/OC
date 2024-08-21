/****************************************************************************************************************************************
Descriptive statistics of the features to ensure coding is correct

Notes:
1. Check the total of the important features matches the total number of important features from the ML models
2. Features are segregated into binary, categorical/ordinal, and continuous. Tabulation for binary and categorical/ordinal; summary for continuous
3. Check for missing %, min and max values (esp. negative values), two values for binary features, and unlabeled values in categorical/ordinal features
*****************************************************************************************************************************************/

log using "reports\descriptive_statistics.smcl", replace
use "data\processed\ukb_OC_data_for_epi_analysis_preprocessed.dta", replace

/****************************************************************************************************************************************
Define important and other features
*****************************************************************************************************************************************/

/* Binary Features */
local binary_imp_features x20004__1355_bo x6141__2_sd x6142__2_retired x2724_mp x1835_mlive x2784_ocp x2674_bcscreen x6142__1_employed x2814_hrt x1797_flive x6179__100_nosuppl x6153__5_ocp x1960_fedup x6164__3_lightdiy x1677_bf x6149__6_denture x20110__101_none x6145__3_deathclose x2100_psych x1920_mood x1940_irrit x2000_worry x20004__1480_wisdom x20004__1507_ep x6162__1_car x6179__1_fishoil
local binary_other_features serous_oc Endometrioid_oc clearcell_oc mucinous_oc x20009_agenoncabin x20011_ageopbin

/* Categorical/Ordinal Features */
local categorical_imp_features x1369_beef x1628_alc x1647_cob x1697_hgt10 x20009_agenoncacat x2050_depfq x738_hhinc x981_plwalkdur
local categorical_other_features _x1717_skin _x6138_edu x1070_tvcat x129_pobnorthcat x130_pobeastcat x20011_ageopcat x20015_sithgtcat x20023_mtimematchcat _x20116_smokestat _x21000_ethn x2734_nlivebcat x2744_bwtfirstcat x2794_ageocpstartcat x2804_ageocplastcat x4079_dbpcat x4080_sbpcat x46_hgsleftcat x47_hgsrightcat _x53_yratcentre _x54_centre x699_lenaddrcat x709_nhhcat

/* Continuous Features */
local continuous_imp_features x20011_ageop x21022_age x2794_ageocpstart x30770_igf1 x130_pobeast x2734_nliveb x2744_bwtfirst x3064_pef x1807_fdage x30210_ephpct x2804_ageocplast x50_standhgt x30050_mch x20015_sithgt x30040_mcv x21002_wt x20023_mtimematch x129_pobnorth x1070_tv x30870_tg x3137_nmeasure x30720_cycc x48_wc x47_hgsright x30290_hlsrpct x30070_rdw x102_pr x20009_agenonca x23099_bcpct x30650_ast x4079_dbp x30270_mscv x30740_glc x30280_irf x4080_sbp x30710_crp x30120_lymct x46_hgsleft x30620_alt x1050_outtime x30520_uk x30690_cl x49_hip x30750_HbA1c x894_modactdur x699_lenaddr x709_nhh x136_nop x30150_ephct x30190_monopct x30200_nphpct x3062_fvc x30760_hdl x30860_tp
local continuous_other_features _x21001_bmi _x2714_period _x3581_agemp

/****************************************************************************************************************************************
Generate Descriptive Statistics
*****************************************************************************************************************************************/

/* Binary Features */
foreach col in `binary_imp_features' `binary_other_features' {
    tab `col' ovarian_cancer, miss col
    summ `col', det
    by ovarian_cancer, sort: summ `col', det
    mdesc `col', abb(25)
}

/* Categorical/Ordinal Features */
foreach col in `categorical_imp_features' `categorical_other_features' {
    tab `col' ovarian_cancer, miss col
    summ `col', det
    by ovarian_cancer, sort: summ `col', det
    mdesc `col', abb(25)
}

/* Continuous Features */
foreach col in `continuous_imp_features' `continuous_other_features' {
    summ `col', det
    by ovarian_cancer, sort: summ `col', det
    mdesc `col', abb(25)
}

log close

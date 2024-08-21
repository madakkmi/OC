/****************************************************************************************************************************************
Pre-process the data from ML models

Notes:
1. Add features in alphabetical order
2. For each feature, add shortname, define values, create relevant categorical and binary variables, and convert to SD values if required
*****************************************************************************************************************************************/

use "OC\data\processed\ukb_OC_data_for_epi_analysis.dta", replace

/****************************************************************************************************************************************
Common Value Labels
*****************************************************************************************************************************************/

* Binary variables
label define yes_no_lab 0 "No" 1 "Yes"

* Generic quarters variables
label define qtr_labels 1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4"

* Age categories
label define age_cat_lab1  0 "<50" 1 "50-59.9" 2 "60-69.9"  3 "70+"
label define age_cat_lab2  0 "<25" 1 "25-34.9" 2 "35-44.99" 3 "45+"
label define age_cat_lab3  0 "<50" 1 "50-59.9" 2 "60+"

/****************************************************************************************************************************************
Feature Engineering
*****************************************************************************************************************************************/

* Pulse rate
rename x102 x102_pr
label variable x102_pr "Pulse rate"

* Outdoors time
rename x1050 x1050_outtime
label variable x1050_outtime "Time spent outdoors in summer"
recode x1050_outtime -10=0.5 -3=. -1=.

* TV time
rename x1070 x1070_tv
label variable x1070_tv "Time spent watching television (TV)"
recode x1070_tv -10=0.5 -3=. -1=.

gen x1070_tvcat = x1070_tv
recode x1070_tvcat (0.5/2 = 1) (3/5 = 2) (6/max = 3) 
label define x1070_tvcat_lab  0 "0 hours/day" 1 "<1 to 2 hours/day" 2 "3 to 5 hours/day" 3 "Above 5 hours/day"
label values x1070_tvcat x1070_tvcat_lab
label variable x1070_tvcat "Time spent watching television (TV) (categorical)"

* Place of birth - north
rename x129 x129_pobnorth
label variable x129_pobnorth "Place of birth - north coordinate"
xtile x129_pobnorthcat = x129_pobnorth, nq(4)
label variable x129_pobnorthcat "Place of birth - north coordinate (categorical)"
label values x129_pobnorthcat qtr_labels

* Place of birth - east
rename x130 x130_pobeast
label variable x130_pobeast "Place of birth - east coordinate"
xtile x130_pobeastcat = x130_pobeast, nq(4)
label variable x130_pobeastcat "Place of birth - east coordinate (categorical)"
label values x130_pobeastcat qtr_labels

* Number of operations
rename x136 x136_nop
label variable x136_nop "Number of operations"

* Beef intake
rename x1369 x1369_beef
label variable x1369_beef "Beef intake"
label define mphe_100377 0 "Never" 1 "Less than once" 2 "Once a week" 3 "2-4 times a week or more" 4 "5-6 times or more a week" -1 "Do not know" -3 "Prefer not to answer" 
recode x1369_beef 4=3 5=3
label values x1369_beef mphe_100377

* Tea intake 
rename _x1488 _x1488_tea
label variable _x1488_tea "Tea intake"
recode _x1488_tea -10=0.5 -3=. -1=.

* Alcohol intake versus 10 years previously
rename x1628 x1628_alc
label variable x1628_alc "Alcohol intake versus 10 years previously"
label define x1628_alc_lab 1 "More nowadays" 2 "About the same" 3 "Less nowadays" -1 "Do not know" -3 "Prefer not to answer"
label values x1628_alc x1628_alc_lab

* Country of birth
rename x1647 x1647_cob
label variable x1647_cob "Country of birth (UK/elsewhere)"
label define mphe_100420 1 "England" 2 "Wales" 3 "Scotland" 4 "Northern Ireland/Republic of Ireland" 5 "Elsewhere" -1 "Do not know" -3 "Prefer not to answer"
recode x1647 5=4
recode x1647 6=5
label values x1647_cob mphe_100420

* Breastfed as a baby
rename x1677 x1677_bf
label variable x1677_bf "Breastfed as a baby"
label values x1677_bf yes_no_lab

* Comparative height at 10
rename x1697 x1697_hgt10
label variable x1697_hgt10 "Comparative height size at age 10"
label define x1697_lab 1 "Shorter" 2 "About average" 3 "Taller"
label values x1697_hgt10 x1697_lab

* Skin color
rename _x1717 _x1717_skin
label variable _x1717_skin "Skin color"
recode _x1717_skin 2=1 4=3
recode _x1717_skin -3=. -1=.
label define mp_100431 1 "Very fair/fair" 3 "Olive" 5 "Brown" 6 "Black" -1 "Do not know" -3 "Prefer not to answer"
label values _x1717_skin mp_100431

* Father alive?
rename x1797 x1797_flive
label variable x1797_flive "Father still alive?"
label values x1797_flive yes_no_lab

* Father's age at death
rename x1807 x1807_fdage
label variable x1807_fdage "Father's age at death"

* Mother alive?
rename x1835 x1835_mlive
label variable x1835_mlive "Mother still alive?"
label values x1835_mlive yes_no_lab

* Townsend deprivation index
xtile _x189_townsendcat = _x189, nq(4)
label variable _x189_townsendcat "Townsend deprivation index (category)"
label values _x189_townsendcat qtr_labels
rename _x189 _x189_townsend
label variable _x189_townsend "Townsend deprivation index"

* Mood swings
rename x1920 x1920_mood
label variable x1920_mood "Mood swings?"
label values x1920_mood yes_no_lab

* Irritability
rename x1940 x1940_irrit
label variable x1940_irrit "Irritability?"
label values x1940_irrit yes_no_lab

* Fed-up feelings
rename x1960 x1960_fedup
label variable x1960_fedup "Fed-up feelings?"
label values x1960_fedup yes_no_lab

* Worry too much
rename x2000 x2000_worry
label variable x2000_worry "Worry too long after embarrassment?"
label values x2000_worry yes_no_lab

* Bilateral oophorectomy
rename x20004__1355 x20004__1355_bo
label variable x20004__1355_bo "Operation - bilateral oophorectomy"
label values x20004__1355_bo yes_no_lab

* Wisdom teeth surgery
rename x20004__1480 x20004__1480_wisdom
label variable x20004__1480_wisdom "Operation - wisdom teeth surgery"

* Ectopic pregnancy surgery
rename x20004__1507 x20004__1507_ep
label variable x20004__1507_ep "Operation - ectopic pregnancy surgery"

* Age of first non-cancer illness
rename x20009 x20009_agenonca
label variable x20009_agenonca "Age when first non-cancer illness diagnosed"

gen x20009_agenoncabin = 1
replace x20009_agenoncabin = 0 if x20009_agenonca == .
label variable x20009_agenoncabin "Age when first non-cancer illness diagnosed available?"
label values x20009_agenoncabin yes_no_lab

gen x20009_agenoncacat = x20009_agenonca
recode x20009_agenoncacat (min/24.999999 = 0) (25/34.999999 = 1) (35/44.999999 = 2) (45/max = 3)
label values x20009_agenoncacat age_cat_lab2
label variable x20009_agenoncacat "Age when first non-cancer illness diagnosed (categorical)"

* Age when operation took place
rename x20011 x20011_ageop
label variable x20011_ageop "Age when operation took place"

gen x20011_ageopbin = 1
replace x20011_ageopbin = 0 if x20011_ageop == .
label variable x20011_ageopbin "Age when operation took place available?"
label values x20011_ageopbin yes_no_lab

xtile x20011_ageopcat = x20011_ageop, nq(4)
label variable x20011_ageopcat "Age when operation took place (categorical)"
label values x20011_ageopcat qtr_labels

* Sitting height
rename x20015 x20015_sithgt
label variable x20015_sithgt "Sitting height"
xtile x20015_sithgtcat = x20015_sithgt, nq(4)
label variable x20015_sithgtcat "Sitting height (categorical)"
label values x20015_sithgtcat qtr_labels

* Mean time to identify matches
rename x20023 x20023_mtimematch
label variable x20023_mtimematch "Mean time to correctly identify matches"
xtile x20023_mtimematchcat = x20023_mtimematch, nq(4)
label variable x20023_mtimematchcat "Mean time to correctly identify matches (categorical)"
label values x20023_mtimematchcat qtr_labels

* Illness of mother - none?
rename x20110__101 x20110__101_none
label variable x20110__101_none "Illness of mother - none (group 2 diseases)"
label values x20110__101_none yes_no_lab

* Smoking status
rename _x20116 _x20116_smokestat
label variable _x20116_smokestat "Smoking status"
recode _x20116_smokestat -1=. -3=.
label define _x20116_lab 0 "Never" 1 "Previous" 2 "Current"
label values _x20116_smokestat _x20116_lab

* Depressed mood frequency
rename x2050 x2050_depfq
label variable x2050_depfq "Frequency of depressed mood in last 2 weeks"
label define m_100484 1 "Not at all" 2 "Several days" 3 "More than half the days" 4 "Nearly every day" -1 "Do not know" -3 "Prefer not to answer"
label values x2050_depfq m_100484

* Seen a psychiatrist for nerves, anxiety, etc.
rename x2100 x2100_psych
label variable x2100_psych "Seen a psychiatrist for nerves, anxiety, tension or depression?"
label values x2100_psych yes_no_lab

* Ethnicity
generate _x21000_ethn = _x21000
label variable _x21000_ethn "Ethnicity"
replace _x21000_ethn = 1 if inlist(_x21000, 1, 1001, 1002, 1003)
replace _x21000_ethn = 2 if inlist(_x21000, 3, 3001, 3002, 3003, 3004)
replace _x21000_ethn = 2 if inlist(_x21000, 5)
replace _x21000_ethn = 4 if inlist(_x21000, 4, 4001, 4002, 4003)
replace _x21000_ethn = 5 if inlist(_x21000, 2, 6, 2001, 2002, 2003, 2004)
replace _x21000_ethn = 5 if _x21000 == .
replace _x21000_ethn = 5 if inlist(_x21000, -1, -3)
label define ethnicity_lab 1 "White European" 2 "Asian" 4 "Black African" 5 "Other/mixed/unknown"
label values _x21000_ethn ethnicity_lab

* BMI
rename _x21001 _x21001_bmi
label variable _x21001_bmi "Body mass index (BMI)"

* Weight
rename x21002 x21002_wt
label variable x21002_wt "Weight"

* Age
rename x21022 x21022_age
label variable x21022_age "Age"

gen x21022_agebin = 0
replace x21022_agebin = 1 if x21022_age >= 65
label variable x21022_agebin "Age above 64?"
label values x21022_agebin yes_no_lab

gen x21022_agecat = x21022_age
recode x21022_agecat (min/49.9 = 0) (50/59.9 = 1) (60/69.9 = 2) (70/max = 3)
label define age_cat 0 "<50" 1 "50-59.9" 2 "60-69.9" 3 "70+"
label values x21022_agecat age_cat
label variable x21022_agecat "Age (category)"

* Body fat percentage
rename x23099 x23099_bcpct
label variable x23099_bcpct "Body fat percentage"

* Had BC screening?
rename x2674 x2674_bcscreen
label variable x2674_bcscreen "Ever had breast cancer screening / mammogram?"
label values x2674_bcscreen yes_no_lab

* Age when period started
rename _x2714 _x2714_period
label variable _x2714_period "Age when periods started (menarche)"
recode _x2714_period -3=. -1=.

* Had menopause?
rename x2724 x2724_mp
label variable x2724_mp "Had menopause?"
label values x2724_mp yes_no_lab

* Number of live births
rename x2734 x2734_nliveb
label variable x2734_nliveb "Number of live births"

gen x2734_nlivebcat = x2734_nliveb
recode x2734_nlivebcat (2/3 = 2) (4/max = 3)
label define x2734_nlivebcat_lab 0 "None" 1 "One" 2 "Two/three" 3 "Above three"
label values x2734_nlivebcat x2734_nlivebcat_lab
label variable x2734_nlivebcat "Number of live births (categorical)"

* Birth weight of first child
rename x2744 x2744_bwtfirst
label variable x2744_bwtfirst "Birth weight of first child"
xtile x2744_bwtfirstcat = x2744_bwtfirst, nq(4)
label variable x2744_bwtfirstcat "Birth weight of first child (categorical)"
label values x2744_bwtfirstcat qtr_labels

* Ever taken oral contraceptives?
rename x2784 x2784_ocp
label variable x2784_ocp "Ever taken oral contraceptive pill?"
label values x2784_ocp yes_no_lab

* Age when started oral contraceptives?
rename x2794 x2794_ageocpstart
label variable x2794_ageocpstart "Age started oral contraceptive pill"

gen x2794_ageocpstartcat = x2794_ageocpstart
recode x2794_ageocpstartcat (min/24.999999 = 0) (25/34.999999 = 1) (35/44.999999 = 2) (45/max = 3)
label values x2794_ageocpstartcat age_cat_lab2
label variable x2794_ageocpstartcat "Age when started oral contraceptive pill (categorical)"

* Age when last used oral contraceptives?
rename x2804 x2804_ageocplast
label variable x2804_ageocplast "Age when last used oral contraceptive pill"
replace x2804_ageocplast = x21022_age if x2804_ageocplast == -11

gen x2804_ageocplastcat = x2804_ageocplast
recode x2804_ageocplastcat (min/24.999999 = 0) (25/34.999999 = 1) (35/44.999999 = 2) (45/max = 3)
label values x2804_ageocplastcat age_cat_lab2
label variable x2804_ageocplastcat "Age when last used oral contraceptive pill (categorical)"

* Had HRT?
rename x2814 x2814_hrt
label variable x2814_hrt "Ever used hormone-replacement therapy (HRT)?"
label values x2814_hrt yes_no_lab

* Biomarkers

gen x30040_orig = x30040
gen x30050_orig = x30050
gen x30070_orig = x30070
gen x30120_orig = x30120
gen x30150_orig = x30150
gen x30190_orig = x30190
gen x30200_orig = x30200
gen x30210_orig = x30210
gen x30270_orig = x30270
gen x30280_orig = x30280
gen x30290_orig = x30290
gen x30520_orig = x30520
gen x30620_orig = x30620
gen x30650_orig = x30650
gen x30690_orig = x30690
gen x30710_orig = x30710
gen x30720_orig = x30720
gen x30740_orig = x30740
gen x30750_orig = x30750
gen x30760_orig = x30760
gen x30770_orig = x30770
gen x30860_orig = x30860
gen x30870_orig = x30870

 
rename x30040 x30040_mcv
label variable x30040_mcv "Mean corpuscular volume"

rename x30050 x30050_mch
label variable x30050_mch "Mean corpuscular haemoglobin"

rename x30070 x30070_rdw
label variable x30070_rdw "Red blood cell (erythrocyte) distribution width"

rename x30120 x30120_lymct
label variable x30120_lymct    "Lymphocyte count"

rename x30150  x30150_ephct
label variable x30150_ephct    "Eosinophill count"

rename x30190 x30190_monopct
label variable x30190_monopct   "Monocyte percentage"

rename x30200 x30200_nphpct
label variable x30200_nphpct    "Neutrophill percentage"

rename x30210 x30210_ephpct
label variable x30210_ephpct    "Eosinophill percentage"

rename x30270  x30270_mscv
label variable x30270_mscv     "Mean sphered cell volume"

rename x30280 x30280_irf
label variable x30280_irf "Immature reticulocyte fraction"

rename x30290 x30290_hlsrpct
label variable x30290_hlsrpct "High light scatter reticulocyte percentage"


rename x30520 x30520_uk
label variable x30520_uk "Potassium in urine"

rename x30620 x30620_alt
label variable x30620_alt "Alanine aminotransferase"

rename x30650 x30650_ast
label variable x30650_ast "Aspartate aminotransferase"

rename x30690  x30690_cl                
label variable x30690_cl  "Cholesterol"

rename x30710 x30710_crp
label variable x30710_crp "C-reactive protein"


rename x30720 x30720_cycc
label variable x30720_cycc "Cystatin C"

rename x30740 x30740_glc
label variable x30740_glc "Glucose"

rename x30750  x30750_HbA1c             
label variable x30750_HbA1c  "Glycated haemoglobin (HbA1c)"

rename x30760 x30760_hdl
label variable x30760_hdl "HDL cholesterol"

rename x30770 x30770_igf1
label variable x30770_igf1 "IGF-1"

rename x30860 x30860_tp
label variable x30860_tp "Total protein"

rename x30870 x30870_tg
label variable x30870_tg "Triglycerides"



local biomarkers_features x30040_mcv x30050_mch x30070_rdw x30120_lymct x30150_ephct x30190_monopct x30200_nphpct x30210_ephpct x30270_mscv x30280_irf x30290_hlsrpct x30520_uk x30620_alt x30650_ast x30690_cl x30710_crp x30720_cycc x30740_glc x30750_HbA1c x30760_hdl x30770_igf1 x30860_tp x30870_tg 

foreach col in `biomarkers_features' {
    summ `col'
    replace `col' = (`col' - r(mean)) / r(sd)
}

* FVC
rename x3062 x3062_fvc
label variable x3062_fvc "Forced vital capacity (FVC)"

* PEF
rename x3064 x3064_pef
label variable x3064_pef "Peak expiratory flow (PEF)"

* Number of measurements
rename x3137 x3137_nmeasure
label variable x3137_nmeasure "Number of measurements made"

* Age at menopause
rename _x3581 _x3581_agemp
label variable _x3581_agemp "Age at menopause (last menstrual period)"
recode _x3581_agemp -3=. -1=.

* Diastolic BP
rename x4079 x4079_dbp
label variable x4079_dbp "Diastolic blood pressure"
xtile x4079_dbpcat = x4079_dbp, nq(4)
label variable x4079_dbpcat "Diastolic blood pressure (categorical)"
label values x4079_dbpcat qtr_labels

* Systolic BP
rename x4080 x4080_sbp
label variable x4080_sbp "Systolic blood pressure"
xtile x4080_sbpcat = x4080_sbp, nq(4)
label variable x4080_sbpcat "Systolic blood pressure (categorical)"
label values x4080_sbpcat qtr_labels

* Hand grip strength - left
rename x46 x46_hgsleft
label variable x46_hgsleft "Hand grip strength (left)"
xtile x46_hgsleftcat = x46_hgsleft, nq(4)
label variable x46_hgsleftcat "Hand grip strength (left) (categorical)"
label values x46_hgsleftcat qtr_labels

* Hand grip strength - right
rename x47 x47_hgsright
label variable x47_hgsright "Hand grip strength (right)"
xtile x47_hgsrightcat = x47_hgsright, nq(4)
label variable x47_hgsrightcat "Hand grip strength (right) (categorical)"
label values x47_hgsrightcat qtr_labels

* Waist circumference
rename x48 x48_wc
label variable x48_wc "Waist circumference"

* Hip circumference
rename x49 x49_hip
label variable x49_hip "Hip circumference"

* Standing height
rename x50 x50_standhgt
label variable x50_standhgt "Standing height"

* Assessment centre visit
gen temp1 = substr(_x53, 1, 4)
encode temp1, generate(_x53_yratcentre)
label variable _x53_yratcentre "Year of attending assessment centre"
rename _x53 _x53_dtatcentre
label variable _x53_dtatcentre "Date of attending assessment centre"

* Assessment centre
generate _x54_centre = _x54
label variable _x54_centre "Assessment centre"
rename _x54 _x54_centre_notused
label variable _x54_centre_notused "UK Biobank assessment centre"
label define centre_short_lab 0 "Scotland" 1 "England - north and middle" 2 "England - south" 3 "Wales"
recode _x54_centre (11004 11005=0) (10003 11001 11006 11008 11009 11010 11013 11014 11016 11017 11021 11023=1) (11002 11007 11011 11012 11018 11020=2) (11003 11022=3)
label values _x54_centre centre_short_lab

* Education
rename _x6138 _x6138_edu
label variable _x6138_edu "Education"
recode _x6138_edu -7=0 -3=. 1=20 6=20 2=10 3=10 4=10 5=10
recode _x6138_edu 10=1 20=2
label define _x6138_edu_lab 0 "None" 1 "NVQ/CSE/A-levels" 2 "Degree/professional"
label values _x6138_edu _x6138_edu_lab

* HH - son and/or daughter
rename x6141__2 x6141__2_sd
label variable x6141__2_sd "People in household - son and/or daughter"
label values x6141__2_sd yes_no_lab

* Employed?
rename x6142__1 x6142__1_employed
label variable x6142__1_employed "In paid employment or self-employed"
label values x6142__1_employed yes_no_lab

* Retired
rename x6142__2 x6142__2_retired
label variable x6142__2_retired "Retired"
label values x6142__2_retired yes_no_lab

* Death of close relative?
rename x6145__3 x6145__3_deathclose
label variable x6145__3_deathclose "Illness, injury, bereavement, stress in last 2 years - death of a close relative"
label values x6145__3_deathclose yes_no_lab

* Denture?
rename x6149__6 x6149__6_denture
label variable x6149__6_denture "Dental problems - denture"
label values x6149__6_denture yes_no_lab

* Medication - oral contraceptive pills
rename x6153__5 x6153__5_ocp
label variable x6153__5_ocp "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones - oral contraceptive pill or minipill"
label values x6153__5_ocp yes_no_lab

* Transport - car/motor vehicle
rename x6162__1 x6162__1_car
label variable x6162__1_car "Types of transport used - car/motor vehicle"
label values x6162__1_car yes_no_lab

* Doing light DIY
rename x6164__3 x6164__3_lightdiy
label variable x6164__3_lightdiy "Types of physical activity in last 4 weeks - light DIY (e.g., pruning, watering the lawn)"
label values x6164__3_lightdiy yes_no_lab

* No mineral or other supplements?
rename x6179__100 x6179__100_nosuppl
label variable x6179__100_nosuppl "Mineral and other dietary supplements - none"

* Minerals/other supplements - fish oil
rename x6179__1 x6179__1_fishoil
label variable x6179__1_fishoil "Mineral and other dietary supplements - fish oil"

* Length at current address
rename x699 x699_lenaddr
label variable x699_lenaddr "Length of time at current address"
recode x699_lenaddr -1 -3 = . -10 = 0.5
xtile x699_lenaddrcat = x699_lenaddr, nq(4)
label variable x699_lenaddrcat "Length of time at current address (categorical)"
label values x699_lenaddrcat qtr_labels

* Number of people in the household
rename x709 x709_nhh
label variable x709_nhh "Number in household"
recode x709_nhh -1 -3 = .
gen x709_nhhcat = x709_nhh
recode x709_nhhcat (4/max = 4)
label define x709_nhhcat_lab 1 "1" 2 "2" 3 "3" 4 "4+"
label values x709_nhhcat x709_nhhcat_lab
label variable x709_nhhcat "Number in household (categorical)"

* Household income
rename x738 x738_hhinc
label variable x738_hhinc "Average total household income before tax"
label define m_100294 1 "Less than 18,000" 2 "18,000 to 30,999" 3 "31,000 to 51,999" 4 "52,000 to 100,000" 5 "Greater than 100,000" -1 "Do not know" -3 "Prefer not to answer"
label values x738_hhinc m_100294

* Duration of walk
rename _x874 _x874_walkdur
label variable _x874_walkdur "Duration of walks"
recode _x874_walkdur -3=. -1=.

* Duration of moderate activity
rename x894 x894_modactdur
label variable x894_modactdur "Duration of moderate activity"
recode x894_modactdur -3=. -1=.

* Duration walking for pleasure
rename x981 x981_plwalkdur
label variable x981_plwalkdur "Duration walking for pleasure"
recode x981_plwalkdur 3=2 4=2 5=3 6=3 7=4
label define x981_plwalkdurlab 1 "Less than 15 minutes" 2 "Between 15 and 1.5 hours" 3 "Between 1.5 hours and 3 hours" 4 "Over 3 hours"
label values x981_plwalkdur x981_plwalkdurlab

* Save the processed data
save "OC\data\processed\ukb_OC_data_for_epi_analysis_preprocessed.dta", replace

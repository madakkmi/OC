/****************************************************************************************************************************************
Test factors identified from literature for their use in risk factor adjusted models

Notes:
1. 
2. 
3. 
*****************************************************************************************************************************************/

/* Merge with additional dataset */
merge 1:1 userID using "OC\data\raw\phy_act_alcohol_sleep_from_ang_07oct2022.dta"
drop if _merge != 3

/* Recode variables */
recode n_1160_0_0 -3=. -1=.
recode n_1558_0_0 -3=. -1=.

/* Define adjustment variables and variables to test */
local basic_adjust_vars i._x54_centre c._x53_yratcentre c.x21022_age c._x189_townsend i._x21000_ethn
local vars_to_test i._x6138_edu family_bc_pc c._x21001_bmi c._x874_walkdur c.pa_met c._x1488_tea i._x20116_smokestat i.n_1558_0_0 c.x2734_nliveb x2724_mp c._x3581_agemp c._x2714_period x2814_hrt c.n_1160_0_0 x2784_ocp x20004__1355_bo i.x738_hhinc

/* Start logging */
log using "OC\reports\testing_factors_for_risk_factor_adjusted_models.smcl", replace

/* Run logistic regression for each variable to test */
foreach col in `vars_to_test' {
    logit ovarian_cancer `col' `basic_adjust_vars', or
    estimates store m1
    logit ovarian_cancer `basic_adjust_vars' if e(sample)
    estimates store m0
    lrtest m1 m0
    local p = r(p)
    display `p'
}

log close

/*
Results:

Education                          _x6138_edu              - NOT SIGNIFICANT
Family history of OC, BC, PC       family_bc_pc            - SIGNIFICANT
BMI                                _x21001_bmi             - NOT SIGNIFICANT
Daily walking time                 _x874_walkdur           - NOT SIGNIFICANT
Physical activity                  pa_met                  - NOT SIGNIFICANT
Tea intake                         _x1488_tea              - NOT SIGNIFICANT
Smoking status                     _x20116_smokestat       - NOT SIGNIFICANT
Alcohol consumption                n_1558_0_0              - NOT SIGNIFICANT
Number of live births              x2734_nliveb            - SIGNIFICANT
Menopausal status                  x2724_mp                - NOT SIGNIFICANT
Age at menopause                   _x3581_agemp            - NOT SIGNIFICANT
Age at menarche                    _x2714_period           - NOT SIGNIFICANT
Hormone therapy                    x2814_hrt               - NOT SIGNIFICANT
Daily sleeping time                n_1160_0_0              - NOT SIGNIFICANT
Use of oral contraceptives         x2784_ocp               - SIGNIFICANT
Tubal ligation                     x20004__1355_bo         - SIGNIFICANT
Family income                      x738_hhinc              - SIGNIFICANT
*/

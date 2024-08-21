/****************************************************************************************************************************************
Create SHAP Pie Chart (Category-wise Summary)

Notes:
1. Make sure every feature has a category.
2. 
3. 
*****************************************************************************************************************************************/

* Import the CSV file
import delimited "OC\reports\imp_features_model_shap_feature_imp_avg_for_150_runs_with_display_category.csv"

* Generate and set display order for categories
gen display_order = 0

replace display_order = 1 if displaycategory == "A - Baseline & personal characteristics"
replace display_order = 2 if displaycategory == "B - Sociodemographics"
replace display_order = 3 if displaycategory == "C - Lifestyle and environment"
replace display_order = 4 if displaycategory == "D - Physical measurements"
replace display_order = 5 if displaycategory == "E - Cognitive function"
replace display_order = 6 if displaycategory == "F - Psychosocial factors"
replace display_order = 7 if displaycategory == "G - Self-reported diseases"
replace display_order = 8 if displaycategory == "H - Medications & Operations"
replace display_order = 9 if displaycategory == "I - Health and medical history"
replace display_order = 11 if displaycategory == "J - Biomarkers"
replace display_order = 12 if displaycategory == "K - Sex-specific factors"


* Rename the existing category field to numeric category and not used in the project
rename category category_num  
gen category = ""

replace category = "Baseline & personal characteristics" if displaycategory == "A - Baseline & personal characteristics"
replace category = "Sociodemographics" if displaycategory == "B - Sociodemographics"
replace category = "Lifestyle and environment" if displaycategory == "C - Lifestyle and environment"
replace category = "Physical measurements" if displaycategory == "D - Physical measurements"
replace category = "Cognitive function" if displaycategory == "E - Cognitive function"
replace category = "Psychosocial factors" if displaycategory == "F - Psychosocial factors"
replace category = "Self-reported diseases" if displaycategory == "G - Self-reported diseases"
replace category = "Medications & Operations" if displaycategory == "H - Medications & Operations"
replace category = "Health and medical history" if displaycategory == "I - Health and medical history"
replace category = "Biomarkers" if displaycategory == "J - Biomarkers"


* Rename shap_importance to shap for simplicity
rename shap_importance shap

* Create the pie chart
graph pie shap, over(category) sort(display_order) angle(178) plabel(_all percent, size(relative1p5) format(%2.0f)) line(lalign(outside)) legend(colfirst size(small)) xsize(10) ysize(10)

* Save the graph
graph save "OC\reports\shap_pie_chart_by_category.gph", replace

* Export the graph to PDF and SVG formats
graph export "OC\reports\shap_pie_chart_by_category.pdf", as(pdf) name("Graph") replace
graph export "OC\reports\shap_pie_chart_by_category.svg", as(svg) name("Graph") replace

* Save the dataset
save "OC\data\interim\shap_pie_chart_data.dta", replace

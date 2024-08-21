/****************************************************************************************************************************************
Table 1

Notes:
1. Generate a descriptive statistics table for specified variables.
2. Use Python to capture and format the results.
3. Export the table to a CSV file.
*****************************************************************************************************************************************/

version 17.0
use "OC\data\processed\ukb_OC_data_for_epi_analysis_preprocessed.dta", replace

* Variables to be included in Table 1
local table1_vars x21022_agecat _x21000_ethn _x189_townsendcat _x6138_edu family_bc_pc x2734_nlivebcat x2784_ocp x20004__1355_bo x738_hhinc

* Using Python to capture the results
python:
from sfi import Matrix, Macro, Scalar, Data, ValueLabel
import pandas as pd
end

* Create a pandas dataframe to hold the results
python:
cols1 =  ['overall_control', 'overall_case']
cols2 =  ['variable', 'level', 'overall_control', 'overall_case']
df = pd.DataFrame(columns=cols2)
end

* Loop through the variables and generate the table
foreach var in `table1_vars' {
    tab `var' ovarian_cancer, missing matcell(x)
    
    python:
    cols = Matrix.getColNames('x')
    index = Matrix.getRowNames('x')
    data = Matrix.get('x')
    df_overall = pd.DataFrame(data=data, columns=cols, index=index)
    df_temp = df_overall.copy()
    df_temp.columns = ['overall_control', 'overall_case']
    df_temp['level'] = df_temp.index
    df_temp['variable'] = "`var'"
    df = df.append(df_temp)
    print(df)
    end
}

* Finalize the dataframe and save to CSV
python:
df.reset_index(inplace=True, drop=True)
df.to_csv("OC/reports/table_1.csv", index=False)
end

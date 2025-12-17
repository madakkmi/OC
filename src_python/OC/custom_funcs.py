# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 2023

@author: madakkmi

"""

#%% packages

import OC.config as config
import numpy as np
import pandas as pd
import copy
import random
from collections import Counter
from catboost.utils import get_roc_curve
from sklearn.metrics import confusion_matrix, classification_report, auc
from statsmodels.stats.multitest import fdrcorrection

from pathlib import Path

#%% constants
TOP_CATS_EXCLUDE = [2, 100089, 100088, 100314, 100091, 100093]
N_ITERATIONS = 5

#%% functions

def log_df_changes(changes_df, df, df_name, comment):
    shape = str(df.shape)
    row = df.shape[0]
    col = df.shape[1]
    
    temp = pd.DataFrame(columns=['dataframe', 'comments', 'shape', 'rows', 'columns'], 
                        data=[[df_name, comment, shape, row, col]])
    return pd.concat([changes_df, temp], ignore_index=True)

def get_baseline_features(phesant_df, uk_bio_dd_df, cats_hier_df, notna_counts_df, 
                          notna_cutoff_percent=config.notna_cutoff_percent):
    '''
    For this analysis we need only baseline features. This procedure
    restricts features to baseline features.
       
    Parameters
    ----------
    phesant_df : pandas.DataFrame
        PHESANT processed data frame
    
    uk_bio_dd_df : pandas.DataFrame
        UK biobank data dictionary
    
    cats_hier_df : pandas.DataFrame
        UK biobank category hierarchy tree    
    
    notna_counts_df : pandas.DataFrame
        Data frame containing phesant_df columns with non-missing information
    
    notna_cutoff_percent : int
        Cutoff value (minimum) for including potential baseline features, default = 95%
    
    Returns
    -------
    list : columns_needed_phesant
        Baseline and health-related outcome columns needed 
    
    Example
    -------
    columns_needed_phesant = restrict_to_baseline_health_outcomes(phesant_df)
    '''
    populated_cols_bio = set(notna_counts_df[notna_counts_df.notna_percent >= notna_cutoff_percent].field_id)
    populated_cols_bio.discard('serID')
    populated_cols = []

    s = pd.Series(data=phesant_df.columns)

    for col_to_be_in in populated_cols_bio:
        populated_cols.extend(list(s[s.str[:] == f'x{str(col_to_be_in)}']))   
        populated_cols.extend(list(s[s.str.startswith(f'x{str(col_to_be_in)}__')]))

    print(f'Number of features populated for at least {notna_cutoff_percent}% non-missing: {len(populated_cols)}')

    all_sub_cats = TOP_CATS_EXCLUDE.copy()
    current_sub = TOP_CATS_EXCLUDE.copy()

    for _ in range(N_ITERATIONS):
        new_parents = [cats_hier_df[cats_hier_df.parent_id == e].child_id.tolist() for e in current_sub]
        new_parents = [e for e in new_parents if e]
        
        if not new_parents:
            break
        
        current_sub = []
        for element in new_parents:
            current_sub.extend(element)
        
        all_sub_cats.extend(current_sub)

    columns_exclude = []
    for cat_out in all_sub_cats:
        columns_exclude.extend(list(uk_bio_dd_df[uk_bio_dd_df['Category'] == cat_out]['FieldID']))

    columns_exclude_phesant = []
    s = pd.Series(data=phesant_df.columns)

    for col_to_be_out in columns_exclude:
        columns_exclude_phesant.extend(list(s[s.str[:] == f'x{str(col_to_be_out)}']))     
        columns_exclude_phesant.extend(list(s[s.str.startswith(f'x{str(col_to_be_out)}__')]))

    columns_needed_phesant = set(populated_cols) - set(columns_exclude_phesant) - {'userID'}
    columns_needed_phesant = list(columns_needed_phesant)

    return columns_needed_phesant

def get_high_corr_features_matrix(corr_matrix, threshold, notna_counts_df, uk_bio_dd_df):
    '''
    Get a reduced size correlation matrix with potential features to be dropped
    and merged with notna_count and notna_percent and UK biobank dictionary to decide
    which features to be dropped.
   
    Parameters
    ----------
    corr_matrix : pandas.DataFrame
        Spearman correlation matrix with absolute values and nulls replaced with zeros

    threshold : float
        Threshold value to be considered for dropping

    notna_counts_df : pandas.DataFrame
        Data frame with index as PHESANT processed features containing columns notna_count and notna_percent

    uk_bio_dd_df : pandas.DataFrame
        UK biobank data dictionary

    Returns
    -------
    pandas.DataFrame : corr_matrix_analyse
        Data frame containing potential features to be dropped with additional columns for 
        making decisions on dropping

    Example
    -------
    corr_matrix_analyse = get_high_corr_features_matrix(corr_matrix, 0.90, notna_counts_df, uk_bio_dd_df)
    '''    
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
    to_be_dropped_upper = [col for col in upper.columns if any(upper[col] > threshold)]
     
    lower = corr_matrix.where(np.tril(np.ones(corr_matrix.shape), k=-1).astype(bool))
    to_be_dropped_lower = [col for col in lower.columns if any(lower[col] > threshold)]

    to_be_dropped_potential = set(to_be_dropped_upper) | set(to_be_dropped_lower)

    corr_matrix_analyse = corr_matrix.loc[to_be_dropped_potential, to_be_dropped_potential]

    corr_matrix_analyse['FieldID'] = corr_matrix_analyse.index
    corr_matrix_analyse['Orig_index'] = corr_matrix_analyse.index
    corr_matrix_analyse['FieldID'] = corr_matrix_analyse.FieldID.apply(
        lambda s: s[1:] if '__' not in s else s[1:s.find('__')]).astype(int)

    corr_matrix_analyse = corr_matrix_analyse.merge(uk_bio_dd_df[['FieldID', 'Field', 'Participants']], 
                                                    left_on='FieldID', right_on='FieldID', how='left')
    corr_matrix_analyse.index = corr_matrix_analyse['Orig_index']
    corr_matrix_analyse = corr_matrix_analyse.merge(notna_counts_df, how='left', left_index=True, right_index=True)
    corr_matrix_analyse.drop(columns=['Orig_index'], inplace=True)

    return corr_matrix_analyse

def get_features_to_be_dropped_based_on_correlation(corr_matrix_analyse, threshold, manual_keep=set(), manual_drop=set()):
    '''
    Identify features to be dropped based on correlation analysis.
    
    Parameters
    ----------
    corr_matrix_analyse : pandas.DataFrame
        Correlation matrix for analysis

    threshold : float
        Correlation threshold for dropping features

    manual_keep : set
        Manually specified features to keep

    manual_drop : set
        Manually specified features to drop

    Returns
    -------
    tuple : (final_to_drop, final_to_keep)
        Sets of features to drop and keep

    Example
    -------
    final_to_drop, final_to_keep = get_features_to_be_dropped_based_on_correlation(corr_matrix_analyse, 0.90)
    '''
    cma = corr_matrix_analyse.copy(deep=True)

    for ind in cma.index:
        cma.at[ind, ind] = 0.0
    
    global_count_dict = {}
    global_keep_dict = {}
    global_drop_dict = {}
    global_manual_dict = {}

    def make_decision(ind_count, col_count):
        if ind_count > col_count:
            return 'DROP'
        elif ind_count < col_count:
            return 'KEEP'
        else:
            return 'MANUAL'
    
    for ind in cma.index:
        for_each_count_dict = {}
        for_each_keep_dict = {}
        for_each_drop_dict = {}
        for_each_manual_dict = {}
        
        for col in cma.index:
            if cma.at[ind, col] >= threshold:
                ind_count = cma.at[ind, 'notna_count']
                col_count = cma.at[col, 'notna_count']
                for_each_count_dict[col] = col_count

                decision = make_decision(ind_count, col_count)
                if decision == 'KEEP':
                    for_each_keep_dict[col] = 'KEEP'
                elif decision == 'DROP':
                    for_each_drop_dict[col] = 'DROP'
                else:
                    for_each_manual_dict[col] = 'MANUAL'

        global_count_dict[ind] = for_each_count_dict
        global_keep_dict[ind] = for_each_keep_dict
        global_drop_dict[ind] = for_each_drop_dict
        global_manual_dict[ind] = for_each_manual_dict

    global_drop_copy = copy.deepcopy(global_drop_dict)
    global_keep_copy = copy.deepcopy(global_keep_dict)
    global_manual_copy = copy.deepcopy(global_manual_dict)

    for dropkey in global_drop_dict:
        for dropcol in global_drop_dict[dropkey]:
            for keepkey in global_keep_dict:
                if dropcol in global_keep_dict[keepkey]:
                    global_keep_copy[keepkey].pop(dropcol, None)
                    print(f'Popped column {dropcol} from global_keep_copy[{keepkey}]')
                    
            for manualkey in global_manual_dict:
                if dropcol in global_manual_dict[manualkey]:
                    global_manual_copy[manualkey].pop(dropcol, None)
                    print(f'Popped column {dropcol} from global_manual_copy[{manualkey}]')

    features_to_drop = set().union(*global_drop_copy.values())
    features_to_keep = set().union(*global_keep_copy.values())

    global_manual_copy = {k: v for k, v in global_manual_copy.items() if v}

    for key in global_manual_copy:
        if key not in manual_keep and key not in manual_drop and key not in features_to_drop:
            manual_keep.add(key)
        for mkey in global_manual_copy[key]:
            if mkey not in manual_keep and mkey not in manual_drop and mkey not in features_to_keep:
                manual_drop.add(mkey)

    final_to_drop = features_to_drop | manual_drop
    final_to_keep = features_to_keep | manual_keep

    assert len(final_to_drop) + len(final_to_keep) == cma.shape[0]
    return final_to_drop, final_to_keep

def split_dataset(phesant_df, n_train, n_test, n_val, seed):
    '''
    Split dataset into training, testing, and validation sets.
   
    Parameters
    ----------
    phesant_df : pandas.DataFrame
        PHESANT processed data frame

    n_train : int
        Number of training samples required

    n_test : int
        Number of testing samples required

    n_val : int
        Number of validation samples required
        
    seed : int
        Seed for random number generator

    Returns
    -------
    tuple : (X, y, X_train, y_train, X_test, y_test, X_val, y_val)
        Resulting datasets

    Example
    -------
    X, y, X_train, y_train, X_test, y_test, X_val, y_val = split_dataset(phesant_df, 10000, 2000, 2000, seed)
    '''
    random.seed(seed)
    X = phesant_df.iloc[:, :-1]
    y = phesant_df.iloc[:, -1]
    
    n_samples_reqd = n_train + n_test + n_val
    random_indices = random.sample(range(len(X)), n_samples_reqd)
    
    train_indices = random_indices[:n_train]
    test_indices = random_indices[n_train: n_train + n_test]
    val_indices = random_indices[n_train + n_test:]
    
    X_train = X.iloc[train_indices, :]
    y_train = y.iloc[train_indices]
    
    X_test = X.iloc[test_indices, :]
    y_test = y.iloc[test_indices]
    
    X_val = X.iloc[val_indices, :]
    y_val = y.iloc[val_indices]

    y_train = y_train.astype(int)
    y_test = y_test.astype(int)
    y_val = y_val.astype(int)
    
    return X, y, X_train, y_train, X_test, y_test, X_val, y_val

def do_sampling(option, X, y, params, seed):
    '''
    Perform sampling on dataset based on the specified option.
    
    Parameters
    ----------
    option : str
        Sampling option ('weighted', 'oversample', 'undersample', or 'none')
    
    X : pandas.DataFrame
        Feature matrix
    
    y : pandas.Series
        Target vector
    
    params : dict
        Model parameters
    
    seed : int
        Seed for random number generator

    Returns
    -------
    tuple : (X, y, params)
        Sampled feature matrix, target vector, and updated parameters
    '''
    if option == 'weighted':
        weights = dict(Counter(y))
        positive_class_wt = round(weights[0] / weights[1], 4)
        params['scale_pos_weight'] = positive_class_wt
        return X, y, params
    
    elif option == 'oversample':
        over = RandomOverSampler(sampling_strategy='minority', random_state=seed)
        X, y = over.fit_resample(X, y)
        return X, y, params
    
    elif option == 'undersample':
        under = RandomUnderSampler(sampling_strategy='majority', random_state=seed)
        X, y = under.fit_resample(X, y)
        return X, y, params
    
    else:
        return X, y, params



def print_cb_model_params(model):
    """Print all parameters of the CatBoost model."""
    print("-----------------------------------------------")
    print(model.get_all_params())
    print("-----------------------------------------------")

def record_model_performance(model, model_name, test_pool, y_test, output_path, return_auc=False):
    """
    Record the performance of a CatBoost model on the test pool. It also records 
    confusion matrix, classification report, y_test, y_pred, and y_pred_proba. 
    ROC AUC is calculated, and FPR, TPR, and thresholds are recorded as well.

    Parameters
    ----------
    model : catboost.core.CatBoostClassifier
        A trained CatBoost classifier

    model_name : str
        Model name to be recorded

    test_pool : catboost.core.Pool
        Test pool to be used for assessing performance

    y_test : array-like
        Ground truth labels for the test set

    output_path : str
        Path to save the output files

    return_auc : bool, optional
        Whether to return the AUC value, by default False

    Returns
    -------
    float, optional
        ROC AUC value if return_auc is True

    Example
    -------
    record_model_performance(cb_model, 'catboost_model1', test_pool, y_test, output_path)
    """
    print_cb_model_params(model)
    
    # Evaluate the model
    metric_values = model.eval_metrics(test_pool, config.eval_metrics, plot=False)
    pd.DataFrame(metric_values).to_csv(f'{output_path}/{model_name}_test_pool_metric_values.csv')

    # Confusion matrix
    y_pred = model.predict(test_pool)
    y_pred_proba = model.predict_proba(test_pool)[:, 1]
    y_pred_raw = model.predict(test_pool, prediction_type='RawFormulaVal')

    conf_matrix = pd.DataFrame(confusion_matrix(y_test, y_pred), 
                               index=["Actual_0", "Actual_1"], 
                               columns=["Predicted_0", "Predicted_1"])
    conf_matrix.to_csv(f'{output_path}/{model_name}_conf_matrix.csv')
    print(conf_matrix)

    # Classification report
    report = classification_report(y_test, y_pred, digits=4)
    print(report)
    with open(f'{output_path}/{model_name}_classification_report.txt', "w") as text_file:
        text_file.write(report)

    # ROC AUC
    curve = get_roc_curve(model, test_pool)
    fpr, tpr, thresholds = curve

    roc_auc = auc(fpr, tpr)
    print(f"ROC AUC for {model_name} is {roc_auc}")
    with open(f'{output_path}/{model_name}_auc_value.txt', "w") as text_file:
        text_file.write(f"ROC AUC for {model_name} is {roc_auc}")

    if return_auc:
        return roc_auc

def get_cb_builtin_feature_importance(model, X_df, ukb_dd_df):
    """
    Returns feature importance identified using CatBoost's built-in default feature importance method.

    Parameters
    ----------
    model : catboost.core.CatBoostClassifier
        A trained CatBoost classifier

    X_df : pandas.DataFrame
        Dataset to get feature names

    ukb_dd_df : pandas.DataFrame
        UK Biobank data dictionary

    Returns
    -------
    pandas.DataFrame : feature_imp_df
        A dataframe containing PHESANT processed feature, original UK Biobank feature, feature importance, and other
        fields from the UK Biobank data dictionary

    Example
    -------
    feature_imp_df = get_cb_builtin_feature_importance(cb_model, X_train, ukb_dd_df)
    """
    feature_imp_df = pd.DataFrame(model.feature_importances_, index=X_df.columns, columns=["Importance"])
    feature_imp_df = feature_imp_df.sort_values(by="Importance", ascending=False)
    feature_imp_df.Importance = feature_imp_df.Importance / feature_imp_df.Importance.sum() * 100

    feature_imp_df["PHESANT_feature"] = feature_imp_df.index
    feature_imp_df["FieldID"] = feature_imp_df.index

    feature_imp_df.loc["userID", "FieldID"] = str(-9999999)
    feature_imp_df["FieldID"] = feature_imp_df.FieldID.apply(lambda s: s[1:] if '__' not in s else s[1:s.find('__')]).astype(int)

    feature_imp_df = feature_imp_df.merge(ukb_dd_df, left_on="FieldID", right_on="FieldID", how="left")
    fields_to_return = ["PHESANT_feature", "Importance", "FieldID", "Field", "Participants", "Category"]

    return feature_imp_df.loc[:, fields_to_return]

def get_cb_shap_feature_importance(model, pool, df_columns, ukb_dd_df):
    """
    Returns feature importance identified using CatBoost's SHAP values.

    Parameters
    ----------
    model : catboost.core.CatBoostClassifier
        A trained CatBoost classifier

    pool : catboost.core.Pool
        Dataset to calculate SHAP values on

    df_columns : iterable
        Column names

    ukb_dd_df : pandas.DataFrame
        UK Biobank data dictionary

    Returns
    -------
    pandas.DataFrame : feature_imp_df
        A dataframe containing PHESANT processed feature, original UK Biobank feature, SHAP feature importance, and other
        fields from the UK Biobank data dictionary

    Example
    -------
    feature_imp_df = get_cb_shap_feature_importance(cb_model, cb_train_pool, X_train.columns, ukb_dd_df)
    """
    shap_values = model.get_feature_importance(pool, type="ShapValues")
    shap_sum = np.abs(shap_values).mean(axis=0)[:-1]

    importance_df = pd.DataFrame([df_columns.tolist(), shap_sum.tolist()]).T
    importance_df.columns = ["column_name", "shap_importance"]
    importance_df = importance_df.sort_values("shap_importance", ascending=False)
    importance_df.shap_importance = importance_df.shap_importance / importance_df.shap_importance.sum() * 100

    importance_df.index = importance_df.column_name
    importance_df["PHESANT_feature"] = importance_df.column_name
    importance_df["FieldID"] = importance_df.column_name

    if "userID" in importance_df.index:
        importance_df.loc["userID", "FieldID"] = str(-9999999)

    importance_df["FieldID"] = importance_df.FieldID.apply(lambda s: s[1:] if '__' not in s else s[1:s.find('__')]).astype(int)

    importance_df = importance_df.merge(ukb_dd_df, left_on="FieldID", right_on="FieldID", how="left")
    fields_to_return = ["PHESANT_feature", "column_name", "shap_importance", "FieldID", "Field", "Participants", "Category", "Path"]

    return importance_df.loc[:, fields_to_return]

def get_average_shap_value(n_runs, begin_seed, files_begin_with):
    """
    Calculate the average SHAP values across multiple runs.

    Parameters
    ----------
    n_runs : int
        Number of runs to average

    begin_seed : int
        Starting seed value

    files_begin_with : str
        Prefix of the SHAP values files

    Returns
    -------
    pandas.DataFrame : shap_avg
        Dataframe with averaged SHAP values

    Example
    -------
    shap_avg = get_average_shap_value(10, 0, 'model_shap_values')
    """
    dfs = [pd.read_csv(f'{files_begin_with}_{begin_seed + i}_shap_feature_imp.csv', index_col=0) for i in range(n_runs)]
    for df in dfs:
        df.index = df.PHESANT_feature

    shap_avg = dfs[0].copy(deep=True)
    shap_avg.drop(columns='shap_importance', inplace=True)

    shap_imp = dfs[0].shap_importance
    for i in range(1, n_runs):
        shap_imp += dfs[i].shap_importance

    shap_imp /= n_runs
    shap_avg = shap_avg.join(shap_imp, how='left')
    shap_avg = shap_avg.sort_values(by='shap_importance', ascending=False)

    return shap_avg

def sim_score(Si, Sj):
    """
    Calculate the similarity score between two sets.

    Parameters
    ----------
    Si : set
        First set

    Sj : set
        Second set

    Returns
    -------
    float : similarity_score
        Similarity score between the sets

    Example
    -------
    score = sim_score(set1, set2)
    """
    return len(Si & Sj) / len(Si | Sj)

def get_average_auc_value(n_runs, begin_seed, files_begin_with):
    """
    Calculate the average AUC value across multiple runs.

    Parameters
    ----------
    n_runs : int
        Number of runs to average

    begin_seed : int
        Starting seed value

    files_begin_with : str
        Prefix of the AUC value files

    Returns
    -------
    float : auc_avg
        Averaged AUC value

    Example
    -------
    auc_avg = get_average_auc_value(10, 0, 'model_auc_value')
    """
    auc_avg = 0.0

    for i in range(n_runs):
        with open(f'{files_begin_with}_{begin_seed + i}_auc_value.txt', 'r') as text_file:
            line = text_file.read()
            start = line.rfind('0.')
            auc_val = float(line[start:])
            auc_avg += auc_val

    auc_avg /= n_runs

    return auc_avg


def add_custom_display_categories(custom_categories_df, input_csvfilename, output_dir):
    """
    Add custom display categories to the input CSV file and save the result.

    Parameters
    ----------
    custom_categories_df : pandas.DataFrame
        DataFrame containing custom categories to be added

    input_csvfilename : str
        Path to the input CSV file

    output_dir : pathlib.Path
        Directory to save the output files

    Returns
    -------
    None

    Example
    -------
    add_custom_display_categories(custom_categories_df, 'input.csv', output_dir)
    """
    input_df = pd.read_csv(input_csvfilename, index_col=0)
    output_df = input_df.merge(custom_categories_df, left_on='FieldID', right_on='FieldID', how='left')
    
    feature_ids_without_categories = set(output_df[output_df.DisplayCategory.isna()].column_name)
    print(f'Features without display categories: {feature_ids_without_categories}')  
    
    output_df.to_csv(output_dir / f'{input_csvfilename.stem}_with_display_category.csv')
    output_df[output_df.DisplayCategory.isna()].to_csv(output_dir / f'{input_csvfilename.stem}_without_display_category.csv')

def create_stata_dta_files_epi_analysis(phe_imp_copy_df, outcome_df, raw_df, additional_raw_fields, replace_with_raw_fields, output_dir):
    """
    Create STATA DTA files for logistic regression analysis.

    Parameters
    ----------
    phe_imp_copy_df : pandas.DataFrame
        DataFrame with processed phenotype data

    outcome_df : pandas.DataFrame
        DataFrame with outcome data

    raw_df : pandas.DataFrame
        DataFrame with raw data

    additional_raw_fields : list
        List of additional raw fields to include

    replace_with_raw_fields : list
        List of fields to replace with raw fields

    output_dir : pathlib.Path
        Directory to save the output files

    Returns
    -------
    None

    Example
    -------
    create_stata_dta_files_epi_analysis(phe_imp_copy_df, outcome_df, raw_df, additional_raw_fields, replace_with_raw_fields, output_dir)
    """
    phe_imp_copy_df.drop(columns=replace_with_raw_fields, inplace=True)
    replace_with_raw_fields = [int(c[1:]) for c in replace_with_raw_fields]

    additional_raw_fields_new_names = dict(zip([f'x{c}' for c in additional_raw_fields], [f'_x{c}' for c in additional_raw_fields]))
    additional_raw_fields_copy = copy.deepcopy(additional_raw_fields)
    
    additional_raw_fields.extend(replace_with_raw_fields)
    additional_raw_fields = list(set([f'x{c}_0_0' for c in additional_raw_fields]))
    additional_raw_fields.append('userId')

    raw_features_needed_df = raw_df.filter(items=additional_raw_fields, axis='columns')
    raw_features_needed_df.index = raw_features_needed_df.userId.astype(int)
    raw_features_needed_df.rename({'userId': 'userID'}, axis='columns', inplace=True)
    
    new_and_old = dict(zip(raw_features_needed_df.columns, [s if not s.startswith('x') else f"{s[:s.find('_0_0')]}" for s in raw_features_needed_df.columns]))
    raw_features_needed_df.rename(new_and_old, axis='columns', inplace=True)
    raw_features_needed_df.rename(columns={'x399_0_': 'x399'}, inplace=True)
    
    mother_ill_df = raw_df.filter(like='x20110_0_', axis='columns')
    father_ill_df = raw_df.filter(like='x20107_0_', axis='columns')
    sibling_ill_df = raw_df.filter(like='x20111_0_', axis='columns')
    
    mother_bc_df = (mother_ill_df == 5).sum(axis='columns') > 0
    mother_pc_df = (mother_ill_df == 13).sum(axis='columns') > 0
    
    father_bc_df = (father_ill_df == 5).sum(axis='columns') > 0
    father_pc_df = (father_ill_df == 13).sum(axis='columns') > 0
    
    sibling_bc_df = (sibling_ill_df == 5).sum(axis='columns') > 0
    sibling_pc_df = (sibling_ill_df == 13).sum(axis='columns') > 0
    
    family_bc_df = mother_bc_df + father_bc_df + sibling_bc_df
    family_pc_df = mother_pc_df + father_pc_df + sibling_pc_df
    family_bc_pc_df = family_bc_df + family_pc_df
    
    family_bc_pc_df.name = 'family_bc_pc'
    phe_imp_copy_df['family_bc_pc'] = family_bc_pc_df
    
    raw_features_needed_df.drop(columns=['userID'], inplace=True)
    phe_imp_copy_df = phe_imp_copy_df.merge(raw_features_needed_df, how='left', left_index=True, right_index=True)
    phe_imp_copy_df = phe_imp_copy_df.merge(outcome_df, how='left', left_index=True, right_index=True)
    
    cols = phe_imp_copy_df.select_dtypes(np.float16).columns
    change_dtypes = dict(zip(cols, [np.float32] * len(cols)))
    phe_imp_copy_df = phe_imp_copy_df.astype(change_dtypes)
    
    phe_imp_copy_df = phe_imp_copy_df.rename(columns=additional_raw_fields_new_names)
    phe_imp_copy_df.to_stata(output_dir / 'ukb_OC_data_for_epi_analysis.dta', write_index=False)
    phe_imp_copy_df.index = phe_imp_copy_df.userID.astype(int)

def get_fdr_corrected_significant_features(variable_value_label_link, epi_continuous_bin_dta_file, 
                                           epi_categorical_dta_file, epi_category_levels_dta_file,
                                           pre_processed_epi_dta_file, all_features_shap_file, imp_features_shap_file):
    """
    Identify significant features after FDR correction.

    Parameters
    ----------
    variable_value_label_link : dict
        Dictionary linking variable values to labels

    epi_continuous_bin_dta_file : str
        Path to the STATA DTA file for continuous and binary variables

    epi_categorical_dta_file : str
        Path to the STATA DTA file for categorical variables

    epi_category_levels_dta_file : str
        Path to the STATA DTA file for category levels

    pre_processed_epi_dta_file : str
        Path to the pre-processed STATA DTA file

    all_features_shap_file : str
        Path to the CSV file containing SHAP values for all features

    imp_features_shap_file : str
        Path to the CSV file containing SHAP values for important features

    Returns
    -------
    None

    Example
    -------
    get_fdr_corrected_significant_features(variable_value_label_link, 'epi_continuous_bin.dta', 
                                           'epi_categorical.dta', 'epi_category_levels.dta',
                                           'pre_processed_epi.dta', 'all_features_shap.csv', 'imp_features_shap.csv')
    """
    def has_digits(index):
        return any(char.isdigit() for char in index)    
    
    epi_data = pd.read_stata(pre_processed_epi_dta_file, iterator=True)
    variable_labels = epi_data.variable_labels()
    value_labels = epi_data.value_labels()
    
    epi_results_df = pd.DataFrame.from_dict(variable_labels, orient='index', columns=['feature_desc'])
    epi_results_df['feature_id'] = epi_results_df.index
    
    def get_part_before_last_underscore(string):
        return string.rsplit('_', 1)[0]
    
    for ind in epi_results_df.index:
        if ind.endswith('_notused') or ind.endswith('cat') or ind.endswith('bin') or ind.endswith('orig') or ind.endswith('_oc') \
        or ind.startswith('ca_') or ind.startswith('n_') or ind.startswith('_est') \
        or ind.startswith('_x') or ind.startswith('oc_') or '_' not in ind or not has_digits(ind):
            print(f'Dropping {ind}')
            epi_results_df.drop(index=ind, inplace=True)
    
    epi_results_df['feature_id_shorter'] = epi_results_df.feature_id.apply(get_part_before_last_underscore)
    
    all_features_shap_df = pd.read_csv(all_features_shap_file, index_col=1)
    imp_features_shap_df = pd.read_csv(imp_features_shap_file, index_col=1)
    
    epi_results_df = epi_results_df.merge(imp_features_shap_df[['DisplayCategory']], left_on='feature_id_shorter',
                                          right_index=True, how='left')
                                          
    all_features_shap_df.rename(mapper={'shap_importance': 'all_features_shap_avg'}, axis='columns', inplace=True)
    epi_results_df = epi_results_df.merge(all_features_shap_df[['all_features_shap_avg']],
                                          left_on='feature_id_shorter', right_index=True, how='left')
    
    epi_results_df.sort_values(by='all_features_shap_avg', ascending=False, inplace=True)
    epi_results_df.insert(epi_results_df.shape[1], 'all_features_shap_rank', range(1, epi_results_df.shape[0] + 1))
    
    imp_features_shap_df.rename(mapper={'shap_importance': 'imp_features_shap_avg'}, axis='columns', inplace=True)
    epi_results_df = epi_results_df.merge(imp_features_shap_df[['imp_features_shap_avg']],
                                          left_on='feature_id_shorter', right_index=True, how='left')
    
    epi_results_df.sort_values(by='imp_features_shap_avg', ascending=False, inplace=True)
    epi_results_df.insert(epi_results_df.shape[1], 'imp_features_shap_rank', range(1, epi_results_df.shape[0] + 1))
    
    epi_results_just_shap_df = epi_results_df.copy(deep=True)
    
    continuous_binary_results_df = pd.read_stata(epi_continuous_bin_dta_file)
    continuous_binary_results_df.index = continuous_binary_results_df[['exposure', 'model']].apply(lambda x: f'{x[0]}_{x[1]}', axis='columns')
    
    categorical_results_df = pd.read_stata(epi_categorical_dta_file)
    categorical_results_df.index = categorical_results_df[['exposure', 'model']].apply(lambda x: f'{x[0]}_{x[1]}', axis='columns')
    
    columns_main_effects = ['exposure', 'var_type', 'model', 'p']
    
    p_values_df = pd.DataFrame(columns=['exposure', 'var_type', 'model', 'p_value', 'term'])
    
    temp = continuous_binary_results_df[columns_main_effects]
    temp.rename(columns={'p': 'p_value'}, inplace=True)
    temp = temp.assign(term='main')
    p_values_df = pd.concat([p_values_df, temp])
    
    temp = categorical_results_df[columns_main_effects]
    temp.rename(columns={'p': 'p_value'}, inplace=True)
    temp = temp.assign(term='main')
    p_values_df = pd.concat([p_values_df, temp])
    
    p_values_df = p_values_df[p_values_df.p_value.notna()]
    
    rejected, q_value = fdrcorrection(p_values_df.p_value, alpha=config.alpha, method='indep', is_sorted=False)
    p_values_df['rejected'] = rejected
    p_values_df['q_value'] = q_value
    
    def get_counts(f, df):
        n = df.loc[df[df['exposure'] == f].index, 'n'].values[0]
        n_cancer = df.loc[df[df['exposure'] == f].index, 'n_cancer'].values[0]
        return n, n_cancer
    
    def get_first_coeff(f, df):
        b = df.loc[df[df['exposure'] == f].index, 'beta'].values[0]
        se = df.loc[df[df['exposure'] == f].index, 'se'].values[0]
        lci = df.loc[df[df['exposure'] == f].index, 'lci'].values[0]
        uci = df.loc[df[df['exposure'] == f].index, 'uci'].values[0]
        p = df.loc[df[df['exposure'] == f].index, 'p'].values[0]
        return b, se, lci, uci, p
    
    def get_exps(b, lci, uci):
        return np.exp(b), np.exp(lci), np.exp(uci)
    
    def cust_round(num, deci=2):
        if np.isnan(num):
            return num
        
        if num == int(num):
            return num
        else:
            while True:
                new_num = round(num, deci)
                if new_num != int(new_num):
                    return new_num
                else:
                    deci += 1
    
    epi_results_df[['n', 'n_cancer', 'beta_exp', 'se', 'lci_exp', 'uci_exp', 'p']] = pd.DataFrame([[np.nan] * 7])
    
    continuous_binary_features = p_values_df[p_values_df.var_type != 'categorical'].exposure
    significant_continuous_binary_features = p_values_df[(p_values_df.q_value < config.alpha) & (p_values_df.var_type != 'categorical')].exposure
    
    for feature in continuous_binary_features:
        n, n_cancer = get_counts(feature, continuous_binary_results_df)
        beta, se, lci, uci, p = get_first_coeff(feature, continuous_binary_results_df)
        beta_exp, lci_exp, uci_exp = get_exps(beta, lci, uci)
        
        epi_results_df.loc[feature, ['n', 'n_cancer', 'beta_exp', 'se', 'lci_exp', 'uci_exp', 'p']] = [n, n_cancer, beta_exp, se, lci_exp, uci_exp, p]
    
    epi_results_df['significant'] = 0
    epi_results_df.loc[significant_continuous_binary_features, 'significant'] = 1
    
    two_dec_cols = ['all_features_shap_avg', 'imp_features_shap_avg', 'beta_exp', 'lci_exp', 'uci_exp']
    needed_cols = ['feature_desc', 'DisplayCategory', 'all_features_shap_avg', 'all_features_shap_rank', 
                   'imp_features_shap_avg', 'imp_features_shap_rank', 'n', 'n_cancer', 'beta_exp', 'lci_exp', 
                   'uci_exp', 'p_value', 'q_value', 'significant']
    modified_col_names = ['Feature Description', 'Category', 'SHAP value (all)', 'SHAP Rank (all)', 
                          'SHAP value (imp)', 'SHAP Rank (imp)', 'N', 'Cases', 'OR', 'LCI', 'UCI', 
                          'P_value', 'Q_value', 'Significant']
    
    epi_results_shorter_df = epi_results_df[epi_results_df.n.notna()]
    
    for col in two_dec_cols:
        epi_results_shorter_df[col] = epi_results_shorter_df[col].apply(lambda s: cust_round(s, deci=2))
        
    epi_results_shorter_df = epi_results_shorter_df.merge(p_values_df[['exposure', 'p_value', 'q_value']], how='left', left_index=True, right_on='exposure')
    epi_results_shorter_df.index = epi_results_shorter_df.feature_id
    epi_results_shorter_df = epi_results_shorter_df[needed_cols]
    epi_results_shorter_df.columns = modified_col_names
    
    epi_results_shorter_df.to_csv(config.reports_dir / f'{epi_continuous_bin_dta_file.stem}_fdr_corrected.csv')
    epi_results_df.to_csv(config.reports_dir / f'{epi_continuous_bin_dta_file.stem}_fdr_corrected_additional_columns.csv')
    
    categorical_levels_results_df = pd.read_stata(epi_category_levels_dta_file)
    categorical_levels_results_df.rename(mapper={'p': 'p_for_level'}, axis='columns', inplace=True)
    categorical_levels_results_df.index = categorical_levels_results_df.exposure.apply(lambda s: s[:s.rfind('=')])
    categorical_levels_results_df['feature_id'] = categorical_levels_results_df.index
    categorical_levels_results_df['categorical_level_desc'] = ''    
    
    categorical_levels_results_df[['beta_exp', 'lci_exp', 'uci_exp']] = np.exp(categorical_levels_results_df[['beta', 'lci', 'uci']])
    categorical_levels_results_df['categorical_level'] = categorical_levels_results_df.exposure.apply(lambda s: s[s.rfind('=')+1:])
    
    cat_q_values = p_values_df[['p_value', 'q_value', 'exposure']].copy(deep=True)
    cat_q_values.index = cat_q_values.exposure
    cat_q_values.drop(columns=['exposure'], inplace=True)
    categorical_levels_results_df = categorical_levels_results_df.merge(cat_q_values, left_index=True, right_index=True, how='left')
    
    categorical_levels_results_df[['feature_desc', 'feature_id_shorter', 'DisplayCategory',
        'all_features_shap_avg', 'all_features_shap_rank', 'imp_features_shap_avg', 'imp_features_shap_rank']] = \
        epi_results_df[['feature_desc', 'feature_id_shorter', 'DisplayCategory', 'all_features_shap_avg', 'all_features_shap_rank', 
                        'imp_features_shap_avg', 'imp_features_shap_rank']]
    
    col_pos = dict(zip(categorical_levels_results_df.columns, range(len(categorical_levels_results_df.columns))))
    for i in range(len(categorical_levels_results_df)):
        feature_id = categorical_levels_results_df.iat[i, col_pos['feature_id']]
        categorical_level = int(categorical_levels_results_df.iat[i, col_pos['categorical_level']])
        value_def = variable_value_label_link[feature_id]
        value_label = value_labels[value_def][categorical_level]
        categorical_levels_results_df.iat[i, col_pos['categorical_level_desc']] = value_label
    
    epi_results_df = categorical_levels_results_df.copy(deep=True)
    epi_results_df.drop(index=epi_results_df.index, inplace=True)
    
    categorical_features = p_values_df[p_values_df.var_type == 'categorical'].exposure
    significant_categorical_features = p_values_df[(p_values_df.q_value < config.alpha) & (p_values_df.var_type == 'categorical')].exposure
    
    for feature in categorical_features:
        epi_results_df = pd.concat([epi_results_df, categorical_levels_results_df.loc[feature, :]])
        
    epi_results_df['significant'] = 0
    epi_results_df.loc[significant_categorical_features, 'significant'] = 1
    
    two_dec_cols = ['all_features_shap_avg', 'imp_features_shap_avg', 'beta_exp', 'lci_exp', 'uci_exp']
    needed_cols = ['feature_desc', 'DisplayCategory', 'all_features_shap_avg', 'all_features_shap_rank', 
                   'imp_features_shap_avg', 'imp_features_shap_rank', 'categorical_level', 'categorical_level_desc', 
                   'n', 'n_cancer', 'beta_exp', 'lci_exp', 'uci_exp', 'p_for_level', 'p_value', 'q_value', 'significant']
    modified_col_names = ['Feature Description', 'Category', 'SHAP value (all)', 'SHAP Rank (all)', 
                          'SHAP value (imp)', 'SHAP Rank (imp)', 'Level', 'Level Desc', 'N', 'Cases', 'OR', 
                          'OR (LCI)', 'OR (UCI)', 'P value (level)', 'P value', 'Q value', 'Significant']
    
    epi_results_df_shorter = epi_results_df[needed_cols]
    
    for col in two_dec_cols:
        epi_results_df_shorter[col] = epi_results_df_shorter[col].apply(lambda s: cust_round(s, deci=2))
    
    epi_results_df_shorter.columns = modified_col_names
    epi_results_df_shorter.to_csv(config.reports_dir / f'{epi_categorical_dta_file.stem}_fdr_corrected.csv')
    epi_results_df.to_csv(config.reports_dir / f'{epi_categorical_dta_file.stem}_fdr_corrected_additional_columns.csv')

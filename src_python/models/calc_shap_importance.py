# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 2023

Author: madakkmi
"""

import sys
import os
import time
import copy
from pathlib import Path
import pandas as pd
import catboost as cb
from collections import Counter
import OC.config as config
import OC.custom_funcs as funcs
import pickle
import lzma


def ensure_correlation_matrix_exists():
    if not os.path.exists(config.interim_data_dir / 'corr_matrix.csv'):
        print('Cannot find correlation analysis results, please do not continue')
        sys.exit(1)
    else:
        print('corr_matrix.csv exists, so continue')


def load_data_files():
    ukb_dd_df = pd.read_csv(config.raw_data_dir / 'ukb_data_dictionary_showcase.csv')
    ukb_dd_df.index = ukb_dd_df.FieldID
    ukb_dd_df.index.rename('FieldID_asindex', inplace=True)

    phe_df = pd.read_feather(config.processed_data_dir / 'ukb_OC_data_for_ML_models.feather')
    phe_df.index = phe_df.userID.astype(int)
    print(f'Shape of phe_df: {phe_df.shape}')

    corr_matrix = pd.read_csv(config.interim_data_dir / 'corr_matrix.csv', index_col=0)
    notna_counts_df = pd.read_feather(config.interim_data_dir / 'notna_counts.feather')
    notna_counts_df.index = notna_counts_df.feature_id

    phe_orig_df = phe_df.copy(deep=True)

    return ukb_dd_df, phe_df, corr_matrix, notna_counts_df, phe_orig_df


def drop_highly_correlated_features(phe_df, corr_matrix, notna_counts_df, ukb_dd_df):
    corr_matrix_processed = funcs.get_high_corr_features_matrix(
        corr_matrix, config.high_corr_threshold, notna_counts_df, ukb_dd_df)
    corr_matrix_processed.to_csv(config.interim_data_dir / 'corr_matrix_processed.csv')

    features_to_drop, features_to_keep = funcs.get_features_to_be_dropped_based_on_correlation(
        corr_matrix_processed, config.high_corr_threshold)

    features_to_drop.update({'x34', 'x21003'})  # Manually add these age proxy features

    pd.Series(list(features_to_drop)).to_csv(config.interim_data_dir / 'corr_matrix_features_dropped.csv')
    pd.Series(list(features_to_keep)).to_csv(config.interim_data_dir / 'corr_matrix_features_retained.csv')

    phe_df.drop(columns=features_to_drop, inplace=True)
    print(f'Shape of phe_df after dropping highly correlated features: {phe_df.shape}')
    return phe_df


def calculate_sample_splits(phe_df):
    n_test_samples = int(phe_df.shape[0] * config.test_percent)
    n_val_samples = int(phe_df.shape[0] * config.val_percent)
    n_train_samples = phe_df.shape[0] - (n_test_samples + n_val_samples)
    assert phe_df.shape[0] == n_train_samples + n_test_samples + n_val_samples
    return n_train_samples, n_test_samples, n_val_samples


def run_model_and_calculate_shap_importance(phe_df, ukb_dd_df, n_train_samples, n_test_samples, n_val_samples):
    options = ['none', 'oversample', 'undersample', 'weighted']
    opt = 3
    params = config.params

    for i in range(config.n_runs):
        X, y, X_train, y_train, X_test, y_test, X_val, y_val = funcs.split_dataset(
            phe_df, n_train_samples, n_test_samples, n_val_samples, config.seed + i)

        print(f'Currently chosen option is {options[opt]}')
        X_train, y_train, params = funcs.do_sampling(options[opt], X_train, y_train, params, config.seed + i)

        cb_train_pool = cb.Pool(X_train, y_train)
        cb_test_pool = cb.Pool(X_test, y_test)
        cb_val_pool = cb.Pool(X_val, y_val)

        cb_model = cb.CatBoostClassifier(**params)
        cb_model.fit(X=cb_train_pool, eval_set=cb_val_pool, verbose=config.n_verbose)
        funcs.record_model_performance(cb_model, f'all_features_model_seed_{config.seed + i}', cb_test_pool, y_test, config.model_output_dir)

        all_features_model_default_imp = funcs.get_cb_builtin_feature_importance(cb_model, X_train, ukb_dd_df)
        all_features_model_default_imp.to_csv(config.model_output_dir / f'all_features_model_seed_{config.seed + i}_default_feature_imp.csv')

        all_features_model_shap_imp = funcs.get_cb_shap_feature_importance(cb_model, cb_val_pool, X_val.columns, ukb_dd_df, shap_calc_type='Regular')
        all_features_model_shap_imp.to_csv(config.model_output_dir / f'all_features_model_seed_{config.seed + i}_shap_feature_imp.csv')


def average_shap_importance_and_auc():
    shap_files_begin_with = f'{config.model_output_dir}/all_features_model_seed'
    shap_avg = funcs.get_average_shap_value(n_runs=config.n_runs, begin_seed=config.seed, files_begin_with=shap_files_begin_with)
    shap_avg.to_csv(config.reports_dir / f'all_features_model_shap_feature_imp_avg_for_{config.n_runs}_runs.csv')

    auc_files_begin_with = f'{config.model_output_dir}/all_features_model_seed'
    auc_avg = funcs.get_average_auc_value(n_runs=config.n_runs, begin_seed=config.seed, files_begin_with=auc_files_begin_with)
    print(f'Average AUC value for {config.n_runs} runs = {auc_avg}')
    with open(config.reports_dir / f'all_features_model_auc_avg_for_{config.n_runs}_runs.txt', 'w') as f:
        f.write(f'Average AUC value for {config.n_runs} runs = {auc_avg}')

    return shap_avg


def calculate_feature_similarity(shap_avg, phe_df):
    top_n = int((phe_df.shape[1] - 1) * config.feature_imp_top_percent)
    main_set = set(shap_avg.iloc[:top_n, :].index)
    similarity_df = pd.DataFrame(columns=['main_set', 'other_set', 'sim_score'])

    for n in range(1, config.n_runs, 1):
        shap_avg_incr = funcs.get_average_shap_value(n_runs=n, begin_seed=config.seed, files_begin_with=shap_files_begin_with)
        other_set = set(shap_avg_incr.iloc[:top_n, :].index)
        s_score = funcs.sim_score(main_set, other_set)
        similarity_df = similarity_df.append(pd.DataFrame(columns=['main_set', 'other_set', 'sim_score'], data=[[f'{config.n_runs}th', f'{n}th', s_score]]))

    similarity_df.to_csv(config.reports_dir / f'all_features_model_similarity_scores_for_{config.n_runs}_runs.csv')


def run_important_features_model(phe_df, n_train_samples, n_test_samples, n_val_samples):
    shap_avg = pd.read_csv(config.reports_dir / f'all_features_model_shap_feature_imp_avg_for_{config.n_runs}_runs.csv', index_col=0)
    shap_avg = shap_avg.sort_values(by='shap_importance', ascending=False)
    top_n = int((phe_df.shape[1] - 1) * config.feature_imp_top_percent)
    imp_features = list(shap_avg.index[:top_n])
    imp_features_and_outcome = imp_features + ['ovarian_cancer']

    phe_imp_df = phe_df.filter(items=imp_features_and_outcome, axis='columns')
    phe_imp_df['userID'] = phe_imp_df.index
    phe_imp_df.reset_index(drop=True, inplace=True)
    phe_imp_df.to_feather(config.interim_data_dir / 'ukb_OC_imp_features.feather')
    phe_imp_df.index = phe_imp_df.userID
    phe_imp_df.drop(columns=['userID'], inplace=True)

    params_imp = config.params
    options = ['none', 'oversample', 'undersample', 'weighted']
    opt = 3

    for i in range(config.n_runs):
        X_imp, y_imp, X_train_imp, y_train_imp, X_test_imp, y_test_imp, X_val_imp, y_val_imp = funcs.split_dataset(
            phe_imp_df, n_train_samples, n_test_samples, n_val_samples, config.seed + i)

        print(f'Currently chosen option is {options[opt]}')
        X_train_imp, y_train_imp, params_imp = funcs.do_sampling(options[opt], X_train_imp, y_train_imp, params_imp, config.seed)

        cb_train_imp_pool = cb.Pool(X_train_imp, y_train_imp)
        cb_test_imp_pool = cb.Pool(X_test_imp, y_test_imp)
        cb_val_imp_pool = cb.Pool(X_val_imp, y_val_imp)

        cb_model_imp = cb.CatBoostClassifier(**params_imp)
        cb_model_imp.fit(X=cb_train_imp_pool, eval_set=cb_val_imp_pool, verbose=config.n_verbose)
        funcs.record_model_performance(cb_model_imp, f'imp_features_model_seed_{config.seed + i}', cb_test_imp_pool, y_test_imp, config.model_output_dir)

        imp_features_model_default_imp = funcs.get_cb_builtin_feature_importance(cb_model_imp, X_train_imp, ukb_dd_df)
        imp_features_model_default_imp.to_csv(config.model_output_dir / f'imp_features_model_seed_{config.seed + i}_default_feature_imp.csv')

        imp_features_model_shap_imp = funcs.get_cb_shap_feature_importance(cb_model_imp, cb_val_imp_pool, X_val_imp.columns, ukb_dd_df, shap_calc_type='Regular')
        imp_features_model_shap_imp.to_csv(config.model_output_dir / f'imp_features_model_seed_{config.seed + i}_shap_feature_imp.csv')


def average_important_features_shap_and_auc():
    shap_files_begin_with = f'{config.model_output_dir}/imp_features_model_seed'
    shap_avg = funcs.get_average_shap_value(n_runs=config.n_runs, begin_seed=config.seed, files_begin_with=shap_files_begin_with)
    shap_avg.to_csv(config.reports_dir / f'imp_features_model_shap_feature_imp_avg_for_{config.n_runs}_runs.csv')

    auc_files_begin_with = f'{config.model_output_dir}/imp_features_model_seed'
    auc_avg = funcs.get_average_auc_value(n_runs=config.n_runs, begin_seed=config.seed, files_begin_with=auc_files_begin_with)
    print(f'Average AUC value for {config.n_runs} runs = {auc_avg}')
    with open(config.reports_dir / f'imp_features_model_auc_avg_for_{config.n_runs}_runs.txt', 'w') as f:
        f.write(f'Average AUC value for {config.n_runs} runs = {auc_avg}')


def add_feature_categories_to_files():
    custom_categories_df = pd.read_csv(config.raw_data_dir / 'ukb_field_custom_categories.csv')
    custom_categories_df.index = custom_categories_df.FieldID
    custom_categories_df.index.name = 'FieldID_idx'

    input_csvfilename = config.reports_dir / 'imp_features_model_shap_feature_imp_avg_for_150_runs.csv'
    funcs.add_custom_display_categories(custom_categories_df, input_csvfilename, config.reports_dir)

    input_csvfilename = config.reports_dir / 'all_features_model_shap_feature_imp_avg_for_150_runs.csv'
    funcs.add_custom_display_categories(custom_categories_df, input_csvfilename, config.reports_dir)


def create_stata_dta_files():
    with lzma.open(config.raw_data_dir / 'ukb_dataset_not_PHESANT_preprocessed.xz') as f:
        raw_df = pd.read_pickle(f)
    raw_df.index = raw_df.userId

    phe_imp_copy_df = pd.read_feather(config.interim_data_dir / 'ukb_OC_imp_features.feather')
    phe_imp_copy_df.index = phe_imp_copy_df.userID

    outcome_df = pd.read_csv(config.raw_data_dir / 'ukb_oc_outcome.txt', sep='\t')
    outcome_df.index = outcome_df.n_eid

    additional_raw_fields = [53, 54, 1558, 189, 6138, 21000, 20116, 21001, 1717, 1160, 1488, 3581, 2714, 874, 22040]
    replace_with_raw_fields = ['x1050', 'x1070', 'x136', 'x709', 'x699', 'x894']
    funcs.create_stata_dta_files_epi_analysis(phe_imp_copy_df, outcome_df, raw_df, additional_raw_fields, replace_with_raw_fields, config.processed_data_dir)


def adjust_fdr_corrected_p_values(variable_value_label_link, outcome_vars, shap_file_1, shap_file_2):
    for outcome_var in outcome_vars:
        epi_continuous_bin_dta_file = Path(config.reports_dir / f'basic_binary_and_continuous_for_{outcome_var}.dta')
        epi_categorical_dta_file = Path(config.reports_dir / f'basic_categorical_for_{outcome_var}.dta')
        epi_category_levels_dta_file = Path(config.reports_dir / f'basic_categorical_levels_for_{outcome_var}.dta')
        funcs.get_fdr_corrected_significant_features(variable_value_label_link, epi_continuous_bin_dta_file,
                                                     epi_categorical_dta_file, epi_category_levels_dta_file,
                                                     pre_processed_epi_dta_file, shap_file_1, shap_file_2)
        time.sleep(2)


def main():
    ensure_correlation_matrix_exists()
    ukb_dd_df, phe_df, corr_matrix, notna_counts_df, phe_orig_df = load_data_files()
    phe_df = drop_highly_correlated_features(phe_df, corr_matrix, notna_counts_df, ukb_dd_df)
    n_train_samples, n_test_samples, n_val_samples = calculate_sample_splits(phe_df)
    run_model_and_calculate_shap_importance(phe_df, ukb_dd_df, n_train_samples, n_test_samples, n_val_samples)
    shap_avg = average_shap_importance_and_auc()
    calculate_feature_similarity(shap_avg, phe_df)
    run_important_features_model(phe_df, n_train_samples, n_test_samples, n_val_samples)
    average_important_features_shap_and_auc()
    add_feature_categories_to_files()
    create_stata_dta_files()

    variable_value_label_link = {
        'x1070_tvcat': 'x1070_tvcat_lab',
        'x1369_beef': 'mphe_100377',
        'x1628_alc': 'x1628_alc_lab',
        'x1647_cob': 'mphe_100420',
        'x1697_hgt10': 'x1697_lab',
        'x2050_depfq': 'm_100484',
        'x3137_nmeasure': 'x3137_nmeasurelab',
        'x738_hhinc': 'm_100294',
        'x981_plwalkdur': 'x981_plwalkdurlab'
    }

    outcome_vars = ['ovarian_cancer', 'subtype_serous_oc', 'subtype_endometrioid_oc', 'subtype_clearcell_oc', 'subtype_mucinous_oc']
    pre_processed_epi_dta_file = Path(config.processed_data_dir / 'ukb_OC_data_for_epi_analysis_preprocessed.dta')
    all_features_shap_file = Path(config.reports_dir / 'all_features_model_shap_feature_imp_avg_for_150_runs_with_display_category.csv')
    imp_features_shap_file = Path(config.reports_dir / 'imp_features_model_shap_feature_imp_avg_for_150_runs_with_display_category.csv')

    adjust_fdr_corrected_p_values(variable_value_label_link, outcome_vars, all_features_shap_file, imp_features_shap_file)


if __name__ == "__main__":
    main()

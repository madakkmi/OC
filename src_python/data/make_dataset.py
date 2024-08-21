# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 2023

@author: madakkmi
"""

import os
import sys
import copy
import time
import random
import numpy as np
import pandas as pd
from collections import Counter
import catboost as cb
import csv
import OC.config as config
import OC.custom_funcs as funcs
import pickle
import lzma

def load_phesant_data():
    """Load PHESANT pre-processed data."""
    phe_df = pd.read_feather(config.interim_data_dir / 'ukb_PHESANT_preprocessed.feather')
    phe_df.index = phe_df.userID.astype(int)
    phe_df.drop(columns=['userID'], inplace=True)
    print(f'Shape of PHESANT pre-processed UKB data: {phe_df.shape}')
    return phe_df

def load_raw_data():
    """Load raw UKB data."""
    with lzma.open(config.raw_data_dir / 'ukb_dataset_not_PHESANT_preprocessed.xz') as f:
        raw_df = pd.read_pickle(f)
    return raw_df

def remove_withdrawn_participants(phe_df):
    """Remove withdrawn participants from the dataset."""
    withdrawn_df = pd.read_csv(config.raw_data_dir / 'ukb_withdrawn_app20175_dt20220222.csv')
    withdrawn_remove = set(phe_df.index) & set(withdrawn_df.eid)
    phe_df.drop(index=withdrawn_remove, inplace=True)
    return phe_df

def save_column_names(phe_df):
    """Save column names for reference."""
    phe_cols_df = pd.DataFrame(data=phe_df.columns, columns=["ukb_col_name"])
    phe_cols_df.to_csv(config.interim_data_dir / 'ukb_col_names.csv', index=False)

def load_data_dictionary():
    """Load data dictionary obtained from UKB."""
    ukb_dd_df = pd.read_csv(config.raw_data_dir / 'ukb_data_dictionary_showcase.csv')
    ukb_dd_df.index = ukb_dd_df.FieldID
    ukb_dd_df.index.rename('FieldID_asindex', inplace=True)
    return ukb_dd_df

def load_feature_hierarchy():
    """Load hierarchical relationships of feature categories."""
    cat_hier_df = pd.read_csv(config.raw_data_dir / 'ukb_catbrowse.csv')
    return cat_hier_df

def calculate_missing_percentage(phe_df):
    """Calculate the percentage of missing values for each feature."""
    notna_counts_df = pd.DataFrame(data=phe_df.notna().sum(axis="index"), index=phe_df.columns, columns=["notna_count"])
    notna_counts_df["notna_percent"] = (notna_counts_df.notna_count / phe_df.shape[0] * 100.0)
    notna_counts_df.sort_values(by="notna_percent", ascending=False, inplace=True)
    
    # Get field ID for dummy columns (the part before '__')
    notna_counts_df['field_id'] = notna_counts_df.index
    notna_counts_df['field_id'] = notna_counts_df.field_id.apply(lambda s: s[1:] if s.find('__') < 0 else s[1:s.find('__')])
    
    # Save missing percentage for reference
    notna_counts_df['feature_id'] = notna_counts_df.index
    notna_counts_df.reset_index(drop=True, inplace=True)
    notna_counts_df.to_feather(config.interim_data_dir / 'notna_counts.feather')
    notna_counts_df.index = notna_counts_df.feature_id
    notna_counts_df.to_csv(config.interim_data_dir / 'notna_counts.csv')
    
    return notna_counts_df

def restrict_data_to_cases_controls(phe_df, outcome_df):
    """Restrict data file to cases and controls."""
    phe_df = phe_df.loc[outcome_df.index, :]
    assert phe_df.shape[0] == outcome_df.shape[0]
    return phe_df

def preprocess_data(phe_df, raw_df, outcome_df, ukb_dd_df, cat_hier_df):
    """Preprocess and merge data."""
    # Calculate missing percentages and restrict features
    notna_counts_df = calculate_missing_percentage(phe_df)
    
    # Get baseline features satisfying missing percentage condition
    baseline_features = funcs.get_baseline_features(phe_df, ukb_dd_df, cat_hier_df, notna_counts_df, notna_cutoff_percent=config.notna_cutoff_percent)
    print(f'Number of baseline features with at least {config.notna_cutoff_percent}% non-missing: {len(baseline_features)}')
    
    # Add essential features
    baseline_features.extend(['userID', 'ovarian_cancer'])
    
    # Append the outcome variable
    outcome_df['ovarian_cancer'] = 0
    outcome_df.loc[outcome_df.ovarianc == 'incident', 'ovarian_cancer'] = 1
    phe_df['ovarian_cancer'] = outcome_df['ovarian_cancer']
    
    # Drop non-selected features
    phe_df = phe_df.filter(items=baseline_features, axis="columns")
    phe_df.reset_index(drop=True, inplace=True)
    phe_df.to_feather(config.interim_data_dir / 'ukb_PHESANT_preprocess_baseline.feather')
    phe_df.index = phe_df.userID.astype(int)
    
    # Replace PHESANT pre-processed continuous variables with raw ones
    initial_shape = phe_df.shape
    
    uniq_counts = phe_df.nunique()
    cont_features = uniq_counts[uniq_counts > 20]
    biomark_features = cont_features.filter(regex='x30[0-9]+[0-9]+[0-9]+').index
    cont_features = cont_features.drop(labels=biomark_features)
    
    cont_features_list = list(cont_features.index)
    cont_features_list.remove('userID')
    cont_features_list.remove('x399')  # Number of incorrect matches in round
    cont_features_list.remove('x400')  # Time to complete round
    
    ukb_raw_features_needed = [f'{c}_0_0' for c in cont_features_list]
    ukb_raw_features_needed.append('userId')
    
    ukb_raw_subset_df = raw_df.filter(items=ukb_raw_features_needed, axis='columns')
    ukb_raw_subset_df.index = ukb_raw_subset_df.userId.astype(int)
    ukb_raw_subset_df.rename({'userId': 'userID'}, axis='columns', inplace=True)
    ukb_raw_subset_df.drop(columns=['userID'], inplace=True)
    
    # Replace negative values with np.nan except for TDI
    townsend = raw_df.x189_0_0
    ukb_raw_subset_df.drop(columns=['x189_0_0'], inplace=True)
    ukb_raw_subset_df.replace([-1, -3, -10], np.nan, inplace=True)
    ukb_raw_subset_df['x189_0_0'] = townsend
    
    ukb_raw_subset_df.rename(dict(zip(ukb_raw_features_needed, cont_features_list)), axis='columns', inplace=True)
    
    # Remove the PHESANT pre-processed columns
    phe_df.drop(columns=cont_features_list, inplace=True)
    
    # Store ovarian_cancer outcome variable to add as the last column
    oc_outcome = phe_df.ovarian_cancer
    phe_df.drop(columns=['ovarian_cancer'], inplace=True)
    
    # Merge with the raw data and append the outcome variable
    phe_df = phe_df.merge(ukb_raw_subset_df, how='left', left_index=True, right_index=True)
    phe_df = pd.concat([phe_df, oc_outcome], axis=1)
    
    assert phe_df.shape == initial_shape
    
    # Remove age proxy variables before the correlation analysis
    exclusions = ['x34', 'x21003']
    phe_df = phe_df.drop(columns=exclusions)
    
    return phe_df

def save_processed_data(phe_df):
    """Save processed data."""
    phe_df['userID'] = phe_df.index
    phe_df.reset_index(drop=True, inplace=True)
    phe_df.to_feather(config.processed_data_dir / 'ukb_OC_data_for_ML_models.feather')
    phe_df.index = phe_df.userID
    phe_df.drop(columns=['userID'], inplace=True)
    print(f'Final shape of phe_data for ML models: {phe_df.shape}')

def main():
    phe_df = load_phesant_data()
    raw_df = load_raw_data()
    phe_df = remove_withdrawn_participants(phe_df)
    save_column_names(phe_df)
    ukb_dd_df = load_data_dictionary()
    cat_hier_df = load_feature_hierarchy()
    
    outcome_df = pd.read_csv(config.raw_data_dir / 'outcome_data.csv')  # Assumed to be the outcome file
    phe_df = restrict_data_to_cases_controls(phe_df, outcome_df)
    
    phe_df = preprocess_data(phe_df, raw_df, outcome_df, ukb_dd_df, cat_hier_df)
    save_processed_data(phe_df)

if __name__ == "__main__":
    main()

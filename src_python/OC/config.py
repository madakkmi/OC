# -*- coding: utf-8 -*-
"""
Configuration settings for the project.

Created on Wed Jan 18 2023

Author: madakkmi
"""

from pathlib import Path

# directories for data and reports
interim_data_dir = Path('/OC/data/interim')
processed_data_dir = Path('/OC/data/processed')
raw_data_dir = Path('/OC/data/raw')
reports_dir = Path('/OC/reports')
model_output_dir = Path('/OC/reports/model_output')
figures_dir = Path('/OC/reports/figures')

# catBoost training configuration
early_stopping = 100
n_iterations = 2000 
n_verbose = 50
task_type = 'CPU'
boost_type = 'Plain'
eval_metrics = ['AUC', 'BalancedAccuracy', 'Logloss', 'Accuracy', 'Precision', 'Recall']

params = {
    'loss_function': 'Logloss',  # Objective function
    'custom_loss': ['AUC'],  
    'eval_metric': 'AUC',  # For early stopping
    'use_best_model': True,
    'od_type': 'Iter',
    'od_wait': early_stopping,
    'task_type': task_type,
    'iterations': n_iterations,
    'boosting_type': boost_type,
}

# Dataset split ratios
val_percent = 0.2
test_percent = 0.2

# Feature selection parameters
feature_imp_top_percent = 0.03
notna_cutoff_percent = 70
high_corr_threshold = 0.9

# Random seed for reproducibility
seed = 123
n_runs = 150  # Number of times to calculate feature importance

# Plot settings
transparent = False
dpi = 300
small_size = 8
medium_size = 11
bigger_size = 12

# Alpha value for statistical tests
alpha = 0.05

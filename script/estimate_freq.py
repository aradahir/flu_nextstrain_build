#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from augur.dates import get_numerical_dates, numeric_date_type
from augur.frequencies import format_frequencies
# KdeFrequencies is the correct class for this method
from augur.frequency_estimators import get_pivots, KdeFrequencies 
import json
from augur.io import read_metadata
from augur.utils import write_json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Estimate sequence frequencies from metadata using the KDE method, stratified by clade.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # REQUIRED INPUTS
    parser.add_argument("--metadata", required=True, help="TSV file of metadata with at least 'strain', 'date', and specified clade columns")
    parser.add_argument("--tree-json", required=True, help="Annotated tree JSON file (Used only for tip names/clade aggregation)")
    
    # KDE-SPECIFIC PARAMETERS (REPLACING DIFFUSION PARAMETERS)
    parser.add_argument("--narrow-bandwidth", required=True, type=float, help="Narrow bandwidth (sigma) for KDE calculation")
    parser.add_argument("--proportion-wide", type=float, default=0.0, help="Proportion of wide bandwidth to use for smoothing")
    parser.add_argument("--wide-bandwidth", type=float, default=0.25)
    # GENERAL PARAMETERS
    parser.add_argument("--pivot-interval-units", choices=["weeks", "months"], default="months")
    parser.add_argument("--pivot-interval", type=int, default=2, help="Interval between pivots in weeks")
    parser.add_argument("--min-date", type=numeric_date_type, help="Minimum date to estimate frequencies for")
    parser.add_argument("--max-date", type=numeric_date_type, help="Maximum date to estimate frequencies for")
    parser.add_argument("--clade-column", nargs='+', default=['clade'], help="Metadata columns to group strains by (e.g., 'clade subclade')")
    parser.add_argument("--output", required=True, help="JSON file in tip-frequencies format")
    args = parser.parse_args()

    # --- 1. Load Data ---
    # Load Tree (required to map tips to clades later)
    with open(args.tree_json, 'r') as f:
        tree_data = json.load(f)
    
    # Load Metadata (must load all clade columns for aggregation)
    columns_to_load = ["strain", "date"] + args.clade_column # Load date and all grouping columns
    metadata = read_metadata(
        args.metadata,
        columns=columns_to_load,
        dtype="string",
    )
    dates = get_numerical_dates(metadata, fmt='%Y-%m-%d')
    
    # --- 2. Collect Valid Strains and Observations Globally ---
    valid_strains_df = metadata.copy()
    valid_strains_df['num_date'] = pd.Series({strain: np.mean(dates[strain]) for strain in dates})
    
    # Filter on date AND all specified clade columns
    valid_strains_df.dropna(subset=['num_date'] + args.clade_column, inplace=True)
    
    all_strains = valid_strains_df.index.tolist()
    all_observations = valid_strains_df['num_date'].values
    
    # --- 3. Calculate Global Pivots and Tip Frequency Matrix (KDE) ---
    pivots = get_pivots(
        all_observations,
        args.pivot_interval,
        args.min_date,
        args.max_date,
        args.pivot_interval_units,
    )
    
    # Initialize the KDE Estimator (no 'method' parameter needed here)
    frequencies_estimator = KdeFrequencies(
        sigma_narrow=args.narrow_bandwidth,
        sigma_wide=args.wide_bandwidth,
        proportion_wide=args.proportion_wide,
        pivot_frequency=args.pivot_interval/ (12 if args.pivot_interval_units == "months" else 52),
        start_date=args.min_date,
        end_date=args.max_date,
    )
    
    # CRITICAL: Calculate ONE global frequency matrix for ALL valid strains
    # This matrix contains the KDE time-series for every sequence (tip)
    global_frequency_matrix = frequencies_estimator.estimate_frequencies(
        all_observations,
        pivots,
    )

    # --- INSERT THIS BLOCK AFTER CALCULATING global_frequency_matrix ---

    # 1. Calculate the sum of all tip frequencies at each pivot point (column sums)
    # The result is a 1D array where each element is the total frequency at that time point.
    pivot_sums = np.sum(global_frequency_matrix, axis=0)

    # 2. Renormalize the entire matrix to ensure each column sums exactly to 1.0.
    # This addresses floating point inaccuracies inherent in KDE estimation.
    # np.divide safely handles division by zero (if a pivot sum is 0, output is 0).
    normalized_matrix = np.divide(
        global_frequency_matrix,
        pivot_sums,
        where=pivot_sums!=0, # Only divide where the sum is not zero
        out=np.zeros_like(global_frequency_matrix)
    )

    # 3. Replace the original matrix with the normalized one
    global_frequency_matrix = normalized_matrix

    # --- CONTINUE with the rest of your script (Steps 4A, 4B, 5) ---
    # The subsequent steps (saving tip frequencies and aggregating clade frequencies)
    # will now use the perfectly normalized matrix.
    frequency_dict = {"pivots": list(pivots)}

    # --- 4. Part A: Save Tip Frequencies ---
    for index, strain in enumerate(all_strains):
        tip_freq_array = global_frequency_matrix[index]
        if tip_freq_array.sum() > 0:
            frequency_dict[strain] = {
                "frequencies": format_frequencies(tip_freq_array)
            }

    # --- 4.5. Part B: Aggregate and Save Group Frequencies (using the matrix) ---
    for group_column in args.clade_column:
        
        unique_group_values = valid_strains_df[group_column].unique()
        
        for group_value in unique_group_values:
            # Find the index position in the global_frequency_matrix for all strains in this group
            indices_in_matrix = [
                all_strains.index(s) 
                for s in valid_strains_df.index[valid_strains_df[group_column] == group_value].tolist()
            ]
            
            if indices_in_matrix:
                # Sum the tip frequencies for all strains belonging to this group
                aggregate_freq = np.sum(global_frequency_matrix[indices_in_matrix], axis=0)

                # -----------------------------------------------------------------
                # 🎯 NEW: Explicitly check and normalize the aggregated array
                # This corrects for potential residual float errors during aggregation
                # and prevents errors if the sum somehow exceeds 1.0 (though unlikely 
                # after the initial matrix normalization).
                aggregate_freq = np.clip(aggregate_freq, 0.0, 1.0)
                # -----------------------------------------------------------------
                
                # Create the required "clade|value" key
                node_name = f"{group_column}|{group_value}"
                
                if aggregate_freq.sum() > 0:
                    frequency_dict[node_name] = {
                        "frequencies": format_frequencies(aggregate_freq)
                    }

    # --- 5. Part C: Add Explicit Global Frequency ---
    global_freq_array = np.array([1.0] * len(pivots))
    frequency_dict["global"] = {
        "frequencies": format_frequencies(global_freq_array)
    }

    write_json(frequency_dict, args.output)

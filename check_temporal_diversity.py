#!/usr/bin/env python3
import pandas as pd

# Check L segment temporal diversity
print("=== L SEGMENT TEMPORAL ANALYSIS ===")
df_L = pd.read_csv('results/segments/L/filtered/rvf_L_metadata.tsv', sep='\t')
print(f"Total L segment records: {len(df_L)}")
print(f"Unique dates: {df_L['date'].nunique()}")
print(f"Date range: {df_L['date'].min()} to {df_L['date'].max()}")
print("Date distribution:")
print(df_L['date'].value_counts().sort_index())

print("\n=== M SEGMENT TEMPORAL ANALYSIS ===")
df_M = pd.read_csv('results/segments/M/filtered/rvf_M_metadata.tsv', sep='\t')
print(f"Total M segment records: {len(df_M)}")
print(f"Unique dates: {df_M['date'].nunique()}")
print(f"Date range: {df_M['date'].min()} to {df_M['date'].max()}")
print("Date distribution:")
print(df_M['date'].value_counts().sort_index())

print("\n=== S SEGMENT TEMPORAL ANALYSIS ===")
df_S = pd.read_csv('results/segments/S/filtered/rvf_S_metadata.tsv', sep='\t')
print(f"Total S segment records: {len(df_S)}")
print(f"Unique dates: {df_S['date'].nunique()}")
print(f"Date range: {df_S['date'].min()} to {df_S['date'].max()}")
print("Date distribution:")
print(df_S['date'].value_counts().sort_index())

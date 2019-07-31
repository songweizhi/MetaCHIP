import pandas as pd
from copy import deepcopy


def drop_small_value_rows(df_in, min_value):

    df_in_copy = deepcopy(df_in)
    df_in_copy['max_in_row'] = df_in_copy.max(axis=1)
    df_out_tmp = df_in_copy[df_in_copy['max_in_row'] >= min_value]
    df_out = df_out_tmp.drop(labels='max_in_row', axis=1)

    return df_out


def drop_small_value_cols_and_rows(input_file, min_value, output_file):
    input_df = pd.read_csv(open(input_file), delimiter='\t', header=0, index_col=0)

    input_df_filtered_row = drop_small_value_rows(input_df, min_value)
    input_df_filtered_row_t = input_df_filtered_row.transpose()
    input_df_filtered_row_col_t = drop_small_value_rows(input_df_filtered_row_t, min_value)
    input_df_filtered_row_col = input_df_filtered_row_col_t.transpose()

    input_df_filtered_row_col.to_csv(output_file, sep='\t')


input_file = 'Test.csv'
min_value = 0
output_file = 'Test2.csv'


drop_small_value_cols_and_rows(input_file, min_value, output_file)

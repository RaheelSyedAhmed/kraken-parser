#!/usr/bin/env python3
import argparse
import pathlib
import kraken_parsers
import pandas as pd


# set the parser variable equal to an instance of the ArgumentParser class of the argparse module
# allows parser to add the arguments we supplied in the parentheses
# sets the args variable to the parsed version of arguments that the user supplies
parser = argparse.ArgumentParser() # foundation setup
parser.add_argument('-i', dest='input_files', help='Supply input files of kraken bin output.', nargs='+', required=True) # establshes way to interpret shell expression and account for flag -i
args = parser.parse_args() # parses arguments supplied by the user and puts them in args

# iterate all the args.input_files which are the kraken bin output files
patient_dict: dict[str, list[kraken_parsers.PatientData.BinData]] = {}
for file in args.input_files: # note: input_files are a list of strings
    full_path = pathlib.PosixPath(file).resolve()
    patient_name = full_path.parent.as_posix().split("/")[-1]
    bin_name = full_path.stem

    current_bin = kraken_parsers.PatientData.BinData(bin_name, "/path/to/bins", patient_name)
    current_bin.set_kraken_output_by_file(full_path)
    patient_dict[patient_name] = patient_dict.get(patient_name, []) + [current_bin]

final_contaminants_df = pd.read_excel("/path/to/contaminants_list.xlsx", sheet_name="final_contaminants_list")
final_contaminants = final_contaminants_df['IDs'].to_list()

for patient_name, bins in patient_dict.items():
    for bin_data in bins:
        bin_data.set_kraken_output(bin_data.filter_kraken_output_by_taxa(final_contaminants))
        bin_data.replace_node_names("output/path/for/filtered_bins")
        bin_data.generate_tax_group_summary("/output/path/for/filtered_group_summaries")

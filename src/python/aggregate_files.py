import sys
import pandas as pd
import glob

prefix = sys.argv[1]
suffix = sys.argv[2]
output_file = sys.argv[3]

# Get a list of all the CSV files in the current directory
csv_files = glob.glob(f"{prefix}*.{suffix}")

# Initialize an empty DataFrame
result = pd.DataFrame()

# Iterate over the CSV files
for file in csv_files:
    factor_name = file.split(prefix)[1].split(suffix)[0].replace(".", "")
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file, header=None, names=["Gene", factor_name])

    # Merge the DataFrames, using the "Gene" column as the key
    if result.empty:
        result = df
    else:
        result = result.merge(df, on="Gene", how="outer")

# Sort the DataFrame by the "Gene" column
result = result.sort_values(by="Gene")

# Reset the index and drop the old index
result = result.reset_index(drop=True)

# Save the aggregated DataFrame to a new CSV file
result.to_csv(output_file, index=False)

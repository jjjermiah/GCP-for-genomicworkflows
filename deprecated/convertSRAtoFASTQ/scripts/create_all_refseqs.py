import pandas as pd

# Initialize an empty list to store all refseqs
all_refseqs = []

print(snakemake.input)

for i in snakemake.input:
    try:
        refseqs = pd.read_csv(i, header=None).values.tolist()
        for j in refseqs:
            all_refseqs.append(j[0])
    except Exception as e:
        print(f"Error reading file {i}: {e}")
        continue
    

unique_refseqs = list(set(all_refseqs))


# Write the unique refseqs to the output file
output_file = "RNA/refseqlists/all_refseqs.csv"
output_file = snakemake.output[0]
with open(output_file, "w") as f:
    for refseq in unique_refseqs:
        f.write(f"{refseq}\n")

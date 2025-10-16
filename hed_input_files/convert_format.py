#!/usr/bin/env python3

import os

# Input and output file paths
input_file = "/data/t1k_hed_input_run19.txt"
output_file = "/data/haplomat_input.txt"

# Create a proper MAC header line
with open(output_file, 'w') as outfile:
    # Write the header line: ID and loci names
    outfile.write("ID\tA\tA\tB\tB\tC\tC\n")
    
    # Process the input data
    with open(input_file, 'r') as infile:
        # Skip the header line
        next(infile)
        
        for line in infile:
            fields = line.strip().split('\t')
            
            if len(fields) < 7:
                continue  # Skip malformed lines
                
            sample_id = fields[0]
            a1 = fields[1].replace("A*", "")
            a2 = fields[2].replace("A*", "")
            b1 = fields[3].replace("B*", "")
            b2 = fields[4].replace("B*", "")
            c1 = fields[5].replace("C*", "")
            c2 = fields[6].replace("C*", "")
            
            # Format MAC format with tabs (no loci names in data lines)
            haplomat_line = f"{sample_id}\t{a1}\t{a2}\t{b1}\t{b2}\t{c1}\t{c2}"
            outfile.write(haplomat_line + "\n")

print(f"Conversion complete! Output file: {output_file}")

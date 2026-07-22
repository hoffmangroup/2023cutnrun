import pysam
import numpy as np
import sys
import os
import math
from collections import defaultdict

def bootstrap_bam(original_bam, shuffled_bam, output_dir):
    original_bam = pysam.AlignmentFile(original_bam, "rb")
    shuffled_bam = pysam.AlignmentFile(shuffled_bam, "rb")
    
    # Get all reads from original BAM
    print("------- read original BAM -------")
    original_read_pairs = {}

    for read in original_bam:
        read_name_prefix = read.query_name.rsplit('.', 1)[0]  # Get the name prefix
        if read_name_prefix not in original_read_pairs:
            original_read_pairs[read_name_prefix] = []
        original_read_pairs[read_name_prefix].append(read)

    original_read_pairs_list = list(original_read_pairs.values())

    # Get all reads from shuffle BAM
    print("------- read shuffle BAM -------")
    boot_read_pairs = {}

    for read in shuffled_bam.fetch():
        read_name_prefix = read.query_name.rsplit('.', 1)[0]  # Get the name prefix
        if read_name_prefix not in boot_read_pairs:
            boot_read_pairs[read_name_prefix] = []
        boot_read_pairs[read_name_prefix].append(read)

    boot_read_pairs_list = list(boot_read_pairs.values())

    total_pairs = original_bam.mapped // 2  # Divide by 2 because each pair consists of two reads

    # Sample read pairs randomly without replacement
    for fraction in [0.05, 0.1, 0.25, 0.5]:
        print(f"------- sample {int(fraction * 100)}% -------")
        boot_pairs = math.ceil(total_pairs * fraction)
        original_pairs = total_pairs - boot_pairs
        print(f"Total pairs: {total_pairs}; original pairs: {original_pairs}; boot pairs: {boot_pairs}.")

        rng = np.random.default_rng()
        selected_original_pairs = rng.choice(original_read_pairs_list, size=(original_pairs, ), replace=False)
        selected_boot_pairs = rng.choice(boot_read_pairs_list, size=(boot_pairs, ), replace=False)
        selected_pairs = np.concatenate((selected_original_pairs, selected_boot_pairs))

        print("------- write -------")
        # Open output BAM file for writing
        output_bam = os.path.join(output_dir + f"/CTCF_50000_main_{int(fraction * 100)}percent_noise.bam")
        print(f"write to: {output_bam}")
        output = pysam.AlignmentFile(output_bam, "wb", template=original_bam)

        # Write selected read pairs to output BAM file
        for read_pair in selected_pairs:
            for read in read_pair:
                output.write(read)
        
        output.close()

    # Close BAM files
    original_bam.close()
    shuffled_bam.close()
    

if __name__ == "__main__":
    original_bam_file = sys.argv[1]
    shuffled_bam_file = sys.argv[2]
    output_dir = sys.argv[3]
    bootstrap_bam(original_bam_file, shuffled_bam_file, output_dir)

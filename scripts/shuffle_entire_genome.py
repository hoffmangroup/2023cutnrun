import pysam
import sys
import random

def get_chr_lengths(bamfile):
    """Get chromosome lengths from the BAM file header."""
    return {ref: length for ref, length in zip(bamfile.references, bamfile.lengths)}

def shuffle_bam(input_bam, output_bam):
    # Open input BAM file
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    shuffled_bamfile = pysam.AlignmentFile(output_bam, "wb", header=bamfile.header)

    # Get chromosome lengths
    chr_lengths = get_chr_lengths(bamfile)

    print("------- read BAM -------")
    read_pairs = {}
    for read in bamfile:
        qname, mate = read.query_name.rsplit('.', 1)
        if mate[0] == "1":
            if qname in read_pairs:
                read_pairs[qname][0] = read
            else:
                read_pairs[qname] = [read, None]
        elif mate[0] == "2":
            if qname in read_pairs:
                read_pairs[qname][1] = read
            else:
                read_pairs[qname] = [None, read]
        

    print("------- shuffle BAM -------")
    for read1, read2 in read_pairs.values():
        # Get the chromosome length
        chr_name = read1.reference_name
        chr_length = chr_lengths[chr_name]

        # Calculate the fragment length
        fragment_length = max(read1.reference_end, read2.reference_end) - min(read1.reference_start, read2.reference_start)
        # Generate a random start position for the fragment
        new_start = random.randint(0, chr_length - fragment_length)
        new_end = new_start + fragment_length

        # Adjust the read positions
        if read1.reference_start < read2.reference_start:
            read1_start = new_start
            read2_start = new_start + (read2.reference_start - read1.reference_start)
        else:
            read2_start = new_start
            read1_start = new_start + (read1.reference_start - read2.reference_start)
        
        # Update read positions
        read1.reference_start = read1_start
        read2.reference_start = read2_start
        read1.next_reference_start = read2.reference_start
        read2.next_reference_start = read1.reference_start
        
        shuffled_bamfile.write(read1)
        shuffled_bamfile.write(read2)

    bamfile.close()
    shuffled_bamfile.close()
    

if __name__ == "__main__":
    input_bam = sys.argv[1]
    output_bam = sys.argv[2]

    shuffle_bam(input_bam, output_bam)

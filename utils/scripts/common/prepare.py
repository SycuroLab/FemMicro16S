import os
import gzip
import multiprocessing
import sys
import re

def count_reads(file_path):
    """Check if a FASTQ file contains any reads, handling both compressed (.gz) and uncompressed files."""
    if file_path.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open
    with open_func(file_path, 'rt') as f:
        try:
            next(f)  # Try reading the first line to determine if there's any content
            return True
        except StopIteration:
            return False

def process_file(file_path):
    """Process a single file to determine if it has reads, supporting both .fastq and .fastq.gz file extensions."""
    # Adjust this condition to check for both possible file extensions
    if file_path.endswith('.fastq.gz') or file_path.endswith('.fastq'):
        if count_reads(file_path):
            return file_path
        else:
            print(f"No reads in {file_path}")
            return None
    else:
        return None  # Skip files that do not match the expected extensions

def main(path):
    """The main processing function to handle file checks, creation of samples directory, and file processing."""
    if not os.path.isdir(path):
        print(f"Input directory '{path}' does not exist. Exiting.")
        return

    # Check and possibly create a 'samples' directory
    samples_dir = 'samples'
    if not os.path.exists(samples_dir):
        os.makedirs(samples_dir)
    else:
        print(f"'{samples_dir}' directory already exists.")

    # Set up parallel processing of files
    pool = multiprocessing.Pool()
    files = [os.path.join(path, file) for file in os.listdir(path) if file.endswith(('.fastq.gz', '.fastq'))]
    read_files = [file for file in pool.map(process_file, files) if file]

    # Extracting unique sample IDs
    unique_sample_ids = set()
    for file in read_files:
        base_name = os.path.basename(file)
        sample_id = re.split(r"_R[1-2]_|_R[1-2]\.|_[1-2]\.", base_name)[0]
        unique_sample_ids.add(sample_id)

    # Separating files into R1 and R2
    R1_files = [file for file in read_files if '_R1_' in file or '_1.' in file]
    R2_files = [file for file in read_files if '_R2_' in file or '_2.' in file]

    # Creating the final tab-delimited file
    with open(os.path.join(samples_dir, 'samples.tsv'), 'w') as output_file:
        output_file.write("Sample_ID\tR1\tR2\n")
        for id in unique_sample_ids:
            R1_file = next((file for file in R1_files if id in file), None)
            R2_file = next((file for file in R2_files if id in file), None)
            output_file.write(f"{id}\t{R1_file}\t{R2_file}\n")

if __name__ == "__main__":
    main(sys.argv[1])

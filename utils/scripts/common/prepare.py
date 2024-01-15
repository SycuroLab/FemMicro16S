
import os
import gzip
import multiprocessing
import sys

def count_reads(file_path):
    """Check if a FASTQ file contains any reads."""
    with gzip.open(file_path, 'rt') as f:
        try:
            next(f)
            return True
        except StopIteration:
            return False

def process_file(file_path):
    """Process a single file to determine if it has reads."""
    if file_path.endswith('.fastq.gz'):
        if count_reads(file_path):
            return file_path
        else:
            print(f"No reads in {file_path}")
            return None

def main(path):
    if not os.path.isdir(path):
        print(f"Input directory '{path}' does not exist. Exiting.")
        return

def main(path):
    if os.path.isdir(path) and not os.path.exists('samples'):
        os.makedirs('samples')
    else:
        print(f"samples directory already exists.")
        
    # Parallel processing of files
    pool = multiprocessing.Pool()
    files = [os.path.join(path, file) for file in os.listdir(path)]
    read_files = [file for file in pool.map(process_file, files) if file]

    # Extracting unique sample IDs
    unique_sample_ids = set()
    for file in read_files:
        base_name = os.path.basename(file)
        sample_id = base_name.split('_R')[0]
        unique_sample_ids.add(sample_id)

    # Separating files into R1 and R2
    R1_files = [file for file in read_files if '_R1_' in file or '_R1' in file]
    R2_files = [file for file in read_files if '_R2_' in file or '_R2' in file]

    
    # Creating the final tab-delimited file
    with open('samples/samples.tsv', 'w') as output_file:
        output_file.write("Sample_ID\tR1\tR2\n")
        for id in unique_sample_ids:
            R1_file = next((file for file in R1_files if id in file), None)
            R2_file = next((file for file in R2_files if id in file), None)
            output_file.write(f"{id}\t{R1_file}\t{R2_file}\n")

if __name__ == "__main__":
    main(sys.argv[1])

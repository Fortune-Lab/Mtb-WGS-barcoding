
def process_fastq_files(file_path1, file_path2):
    try:
        with open(file_path1, 'r') as file1, open(file_path2, 'r') as file2:
            # Process the files line by line or in chunks
            for line1, line2 in zip(file1, file2):
                # Perform operations using data from both files
                print(f"File 1: {line1.strip()}, File 2: {line2.strip()}")
    except FileNotFoundError:
        print("Error: One or both files not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage
file_path_1 = 'file1.fastq'
file_path_2 = 'file2.fastq'
process_fastq_files(file_path_1, file_path_2)


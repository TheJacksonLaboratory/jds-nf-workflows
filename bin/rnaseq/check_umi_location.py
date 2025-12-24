import gzip
import sys
import re

def is_dna(sequence):
    # Checks if a string contains only DNA characters (and N)
    return bool(re.match(r'^[ACGTN]+$', sequence))

def check_fastq_header(file_path):
    try:
        # Handle both gzipped and plain text files
        if file_path.endswith('.gz'):
            opener = gzip.open
        else:
            opener = open

        print(f"Checking file: {file_path}...")
        
        with opener(file_path, 'rt') as f:
            # Check the first few reads to be safe
            check_count = 0
            has_umi = True
            
            for line in f:
                if check_count >= 10: # Check first 10 reads
                    break
                
                # Only process header lines (start with @)
                if line.startswith('@'):
                    check_count += 1
                    
                    # Split line by space to get the ID (remove comments)
                    read_id = line.strip().split()[0]
                    
                    # Split by colon
                    parts = read_id.split(':')
                    
                    # Standard Illumina ID has 7 parts. 
                    # If UMI is injected, it usually has 8 parts.
                    if len(parts) == 8:
                        potential_umi = parts[-1]
                        
                        # Verify the last part is actually DNA
                        if not is_dna(potential_umi):
                            print(f"  [Warn] Read {check_count} has 8 fields, but last field '{potential_umi}' is not DNA.")
                            has_umi = False
                    else:
                        has_umi = False
            
            if has_umi and check_count > 0:
                print(f"VERIFIED: UMI tags detected in headers (8 fields, last field is DNA).")
                print(f"Example UMI from last checked read: {parts[-1]}")
                return True
            else:
                print(f"NO UMI DETECTED: Headers appear to be standard Illumina format (or other non-UMI format).")
                return False

    except Exception as e:
        print(f"Error reading file: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python check_umi.py <path_to_fastq>")
        sys.exit(1)
        
    check_fastq_header(sys.argv[1])
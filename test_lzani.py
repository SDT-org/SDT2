import subprocess
import sys
import os

# Test lzani directly
lzani_path = r"C:\Users\mclund2\Desktop\SDT2\backend\bin\lz-ani.exe"
fasta_path = r"C:\Users\mclund2\Desktop\CodePen\ANI\anello_ORF1.fasta"

# First, try to get version/help
print("Testing lz-ani executable...")
print(f"Executable exists: {os.path.exists(lzani_path)}")
print(f"FASTA exists: {os.path.exists(fasta_path)}")

# Try running with --help
print("\n--- Testing --help ---")
try:
    result = subprocess.run([lzani_path, "--help"], capture_output=True, text=True, timeout=5)
    print(f"Return code: {result.returncode}")
    print(f"Stdout: {result.stdout[:500]}")
    print(f"Stderr: {result.stderr[:500]}")
except subprocess.TimeoutExpired:
    print("Timed out!")
except Exception as e:
    print(f"Error: {e}")

# Try running with minimal command
print("\n--- Testing minimal all2all command ---")
try:
    cmd = [
        lzani_path,
        "all2all",
        "--in-fasta", fasta_path,
        "--out", "test_output.tsv",
        "--out-ids", "test_ids.tsv"
    ]
    print(f"Command: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    print(f"Return code: {result.returncode}")
    print(f"Stdout: {result.stdout[:500]}")
    print(f"Stderr: {result.stderr[:500]}")
    
    # Check if output files were created
    if os.path.exists("test_output.tsv"):
        print("\nOutput file created!")
        with open("test_output.tsv", "r") as f:
            print(f"First 200 chars: {f.read(200)}")
    else:
        print("\nNo output file created")
        
except subprocess.TimeoutExpired:
    print("Timed out after 30 seconds!")
except Exception as e:
    print(f"Error: {e}")

# Try with shell=True
print("\n--- Testing with shell=True ---")
try:
    cmd_str = f'"{lzani_path}" all2all --in-fasta "{fasta_path}" --out test_output2.tsv --out-ids test_ids2.tsv'
    print(f"Command: {cmd_str}")
    
    result = subprocess.run(cmd_str, capture_output=True, text=True, shell=True, timeout=30)
    print(f"Return code: {result.returncode}")
    print(f"Stdout: {result.stdout[:500]}")
    print(f"Stderr: {result.stderr[:500]}")
except subprocess.TimeoutExpired:
    print("Timed out after 30 seconds!")
except Exception as e:
    print(f"Error: {e}")

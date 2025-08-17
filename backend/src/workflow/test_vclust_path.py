from pathlib import Path
import os

# This file is now correctly in the same directory as vclust_pipeline.py
current_dir = os.path.dirname(os.path.abspath(__file__))
print(f"Current directory: {current_dir}")

# Method 1: Using os.path
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(current_dir)))
print(f"Root directory (Method 1): {root_dir}")
vclust_dir = os.path.join(root_dir, "vclust")
vclust_script = os.path.join(vclust_dir, "vclust.py")
print(f"Vclust script path: {vclust_script}")
print(f"Exists: {os.path.exists(vclust_script)}")

# Method 2: Using Path
file_path = Path(__file__).resolve()
print(f"Current file: {file_path}")
root_dir_alt = file_path.parent.parent.parent.parent
print(f"Root directory (Method 2): {root_dir_alt}")
vclust_script_alt = root_dir_alt / "vclust" / "vclust.py"
print(f"Alternative vclust script path: {vclust_script_alt}")
print(f"Exists (alt): {os.path.exists(vclust_script_alt)}")

# Method 3: Working backward from current working directory
cwd = os.getcwd()
print(f"Current working directory: {cwd}")
vclust_dir_cwd = os.path.join(cwd, "vclust")
vclust_script_cwd = os.path.join(vclust_dir_cwd, "vclust.py")
print(f"Vclust script path (from CWD): {vclust_script_cwd}")
print(f"Exists (CWD): {os.path.exists(vclust_script_cwd)}")

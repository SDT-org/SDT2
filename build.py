import os
import sys
import shutil
import subprocess
import platform

venv_path = os.path.join(".", "venv-pywebview")
build_path = os.path.join(".", "build")
os.makedirs(build_path, exist_ok=True)
assets_path = "./assets"

bio_data_path = ["Bio", "Align", "substitution_matrices", "data"]
bio_data_src_path = os.path.join(
    venv_path,
    "Lib" if os.name == "nt" else os.path.join("lib", "python3.11"),
    "site-packages",
    *bio_data_path,
)
bio_data_dest_path = os.path.join(*bio_data_path)
bio_data_option = f"--include-data-dir={bio_data_src_path}={bio_data_dest_path}"

python_executable = (
    os.path.join(venv_path, "Scripts", "python")
    if os.name == "nt"
    else os.path.join(venv_path, "bin", "python")
)

if platform.system() == "Windows":
    icon_option = f"--windows-icon-from-ico={os.path.join(assets_path, 'app.ico')}"
elif platform.system() == "Darwin":
    icon_option = f"--macos-app-icon={os.path.join(assets_path, 'app.icns')}"
else:
    icon_option = ""


def create_command(path, args):
    file_without_extension = os.path.splitext(os.path.basename(path))[0]
    report_path = os.path.join(
        "build", f"{file_without_extension}-compilation-report.xml"
    )
    command = [
        python_executable,
        "-m",
        "nuitka",
        bio_data_option,
        f"--report={report_path}",
    ] + args
    if os.name == "posix" and sys.platform == "darwin":
        command.append("--macos-create-app-bundle")
    command.append(icon_option)
    command.append(path)
    return command


def run_command(command):
    result = subprocess.Popen(command, stdout=subprocess.PIPE)
    result.communicate()
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    return result.returncode == 0


def main():
    command = create_command(
        os.path.join("backend", "src", "app.py"),
        [
            "--standalone",
            "--onefile",
            "--include-data-dir=gui=gui",
            "--nofollow-import-to=matplotlib",
            "--nofollow-import-to=doctest",
            "--output-filename=SDT2",
            "--windows-disable-console",
            f"--output-dir={build_path}",
        ],
    )

    run_command(command)


if __name__ == "__main__":
    main()

import os
import sys
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
    print(command)
    match platform.system():
        case "Windows":
            command.extend(
                [
                    f"--windows-icon-from-ico={os.path.join(assets_path, 'app.ico')}",
                    "--windows-disable-console",
                ]
            )
        case "Darwin":
            command.extend(
                [
                    "--macos-create-app-bundle",
                    f"--macos-app-icon={os.path.join(assets_path, 'app.icns')}",
                ]
            )
        case _:
            pass

    if not platform.system() == "Darwin":
        command.append("--onefile")

    command.append(path)

    print(command)

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
            "--include-data-dir=gui=gui",
            "--include-data-dir=docs=docs",
            "--nofollow-import-to=matplotlib",
            "--nofollow-import-to=doctest",
            "--output-filename=SDT2",
            f"--output-dir={build_path}",
        ],
    )

    if run_command(command) and platform.system() == "Darwin":
        os.rename(
            os.path.join(".", "build", "app.app"),
            os.path.join(".", "build", "SDT2.app"),
        )


if __name__ == "__main__":
    main()

import os
import sys
import subprocess
import platform
import argparse

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

python_executable = (
    os.path.join(venv_path, "Scripts", "python")
    if os.name == "nt"
    else os.path.join(venv_path, "bin", "python")
)


def get_parasail_library_name():
    if sys.platform == "win32":
        return "parasail.dll"
    elif sys.platform == "darwin":
        return "libparasail.dylib"
    elif sys.platform.startswith("linux"):
        return "libparasail.so"
    else:
        raise RuntimeError("Unsupported platform")


def get_parasail_library_path(library_name):
    venv_path = os.environ.get("VIRTUAL_ENV")
    if not venv_path:
        raise RuntimeError("VIRTUAL_ENV not set. Activate your virtual environment.")

    for root, _, files in os.walk(venv_path):
        if library_name in files:
            return os.path.join(root, library_name)

    raise RuntimeError(f"{library_name} not found in the virtual environment.")


def make_platform_build_command(settings):
    path = os.path.join("backend", "src", "app.py")
    file_without_extension = os.path.splitext(os.path.basename(path))[0]
    report_path = os.path.join(
        "build", f"{file_without_extension}-compilation-report.xml"
    )
    parasail_library_name = get_parasail_library_name()
    parasail_library_path = get_parasail_library_path(parasail_library_name)

    command = [
        python_executable,
        "-m",
        "nuitka",
        f"--include-data-dir={bio_data_src_path}={bio_data_dest_path}",
        f"--include-data-file={parasail_library_path}={os.path.join('parasail', parasail_library_name)}",
        f"--report={report_path}",
        "--standalone",
        "--include-data-dir=gui=gui",
        "--include-data-dir=docs=docs",
        "--nofollow-import-to=matplotlib",
        "--nofollow-import-to=doctest",
        "--output-filename=SDT2",
        f"--output-dir={build_path}",
        "--assume-yes-for-downloads",
    ]

    match platform.system():
        case "Windows":
            command.extend(
                [
                    f"--windows-icon-from-ico={os.path.join(assets_path, 'app.ico')}",
                    "--windows-console-mode=disable",
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

    if not platform.system() == "Darwin" and not settings.disable_onefile:
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
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--disable-onefile",
        action="store_true",
        help="Disable --onefile on Windows and Linux",
    )
    settings = parser.parse_args()

    command = make_platform_build_command(settings)

    if run_command(command) and platform.system() == "Darwin":
        os.rename(
            os.path.join(".", "build", "app.app"),
            os.path.join(".", "build", "SDT2.app"),
        )


if __name__ == "__main__":
    main()

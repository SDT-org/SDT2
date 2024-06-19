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


def copy_script_binaries():
    if sys.platform == "darwin":
        sdt_path = os.path.join(build_path, "SDT2.app", "Contents", "MacOS", "SDT2")
        cluster_path = os.path.join(
            build_path, "cluster.app", "Contents", "MacOS", "cluster"
        )
    elif sys.platform == "linux":
        sdt_path = os.path.join(build_path, "SDT2.bin")
        cluster_path = os.path.join(build_path, "cluster.bin")
    elif os.name == "nt":
        sdt_path = os.path.join(build_path, "SDT2.exe")
        cluster_path = os.path.join(build_path, "cluster.exe")
        print(sdt_path)

    app_bin_path = os.path.join(
        build_path,
        (
            os.path.join("app.app", "Contents", "MacOS")
            if sys.platform == "darwin"
            else "app.dist"
        ),
        "bin",
    )

    if not os.path.exists(app_bin_path):
        os.makedirs(app_bin_path)
    shutil.copy2(sdt_path, app_bin_path)
    shutil.copy2(cluster_path, app_bin_path)
    print(f"Moved SDT and Cluster binaries into {app_bin_path}")


def main():
    commands = {
        "sdt": create_command(
            os.path.join("backend", "scripts", "SDT2.py"),
            [
                "--onefile",
                "--standalone",
                "--nofollow-import-to=matplotlib",
                "--nofollow-import-to=doctest",
                "--nofollow-import-to=pywebview",
                f"--output-dir={build_path}",
            ],
        ),
        "cluster": create_command(
            os.path.join("backend", "scripts", "cluster.py"),
            [
                "--onefile",
                "--standalone",
                "--nofollow-import-to=matplotlib",
                f"--output-dir={build_path}",
            ],
        ),
        "app": create_command(
            os.path.join("backend", "src", "app.py"),
            [
                "--standalone",
                "--include-data-dir=gui=gui",
                "--output-filename=SDT2",
                "--windows-disable-console",
                f"--output-dir={build_path}",
            ],
        ),
    }

    build_thing = None
    if len(sys.argv) > 1:
        build_thing = sys.argv[1]

    if build_thing in commands:
        run_command(commands[build_thing])
        if build_thing == "app":
            copy_script_binaries()
    elif build_thing != None and len(build_thing):
        print(
            "Specify what to build: 'sdt' or 'cluster' scripts, or 'app' to build everything."
        )
    else:
        results = [
            run_command(commands["sdt"]),
            run_command(commands["cluster"]),
            run_command(commands["app"]),
        ]

        if all(result == True for result in results):
            copy_script_binaries()


if __name__ == "__main__":
    main()

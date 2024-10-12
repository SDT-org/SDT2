## Developing

### Prerequisites

- Install Bun https://bun.sh/docs/installation
- (Optional but recommended) [Install pyenv](https://github.com/pyenv/pyenv?tab=readme-ov-file#installation)
- (Linux) Make sure the following are installed:
  - libffi-dev
  - autoconf
  - automake
  - libtools
- Install Python 3.11
  - With pyenv you can do this with `pyenv install 3.11` then `pyenv local 3.11`
- Install virtualenv: `pip install virtualenv`

### Setup

- Run `bun run init`

### Run

- Run `bun run start`
  - To enable dev mode on MacOS/Linux, run with the DEBUG flag: `DEBUG=true bun run start`
  - To enable dev mode in Windows powershell, run `$env:DEBUG="true"` once before starting

## License

This project is licensed under the terms of the GNU General Public License v3.0 (GPL-3.0).

## Acknowledgements

This project uses several open-source libraries:

1. **pywebview** - [BSD License](https://github.com/r0x0r/pywebview/blob/master/LICENSE)
2. **biopython** - [Biopython License Agreement](https://github.com/biopython/biopython/blob/master/LICENSE.rst)
3. **numpy** - [BSD License](https://github.com/numpy/numpy/blob/main/LICENSE.txt)
4. **scipy** - [BSD License](https://github.com/scipy/scipy/blob/main/LICENSE.txt)
5. **psutil** - [BSD License](https://github.com/giampaolo/psutil/blob/master/LICENSE)
6. **pandas** - [BSD License](https://github.com/pandas-dev/pandas/blob/main/LICENSE)
7. **nuitka** - [Apache License 2.0](https://github.com/Nuitka/Nuitka/blob/main/LICENSE.txt)
8. **react** - [MIT License](https://github.com/facebook/react/blob/main/LICENSE)

Please see the respective repositories for more information on their licensing.

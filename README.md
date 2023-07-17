# vesta-espresso

[VESTA](https://jp-minerals.org/vesta/en/) is an extremely powerful and useful tool heavily used in the condensed matter physics, quantum chemistry and materials science communities. However, it does not natively support viewing Quantum ESPRESSO (QE) input files (which contain the structural information). Unfortunately, since VESTA is closed-source, third-party developers cannot add this functionality to VESTA.

This repository contains an extremely simple shell script that allows you to open Quantum ESPRESSO input files (which we consider to be files with `.pwi` or `.in` extensions) in VESTA. For example, by aliasing the script in your `.bashrc` (or equivalent) to, e.g., `vesta`, you can open Quantum ESPRESSO input files in VESTA by simply typing `vesta <input file>` in your terminal. There is also a simple Mac App that allows you to drag Quantum ESPRESSO input files onto it to open them in VESTA, or open files from within finder as you would with natively supported file types (e.g., `POSCAR`, `.cif`, etc.)

The script can accept both Quantum ESPRESSO input files and other VESTA-supported formats (which will be passed to VESTA as is without any processing), or a combination thereof. QE files are converted to a hidden temporary CIF file using [pymatgen-io-espresso](https://github.com/oashour/pymatgen-io-espresso) and then opened in VESTA.

Since the script uses `pymatgen-io-espresso` under the hood, it can only convert the QE input files that `pmg-espresso` can read. At the time of this writing, `pmg-espresso` supports all input files, including `ibrav != 0` (although it is not thoroughly tested), except for `crystal_sg` (i.e., specifying the space group and inequivalent atoms only). See the `pmg-espresso` documentation for more details.

## Requirements
1. [VESTA](https://jp-minerals.org/vesta/en/). I tested with version 3.5.6, but any decently modern version should work.
2. [pymatgen-io-espresso](https://github.com/oashour/pymatgen-io-espresso) in a virtual environment (see below)
3. A *nix operating system. I tested on macOS 13.2, but I see no reason the script wouldn't work on Linux. The App (for Finder integration) is naturally macOS only, but you should be able to get similar functionality quite easily with most Linux window managers just using the script. You might be able to get this working on Windows using [WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10), but I haven't tested this.

## Simple Installation (Mac)
These instructions assume you have already downloaded [VESTA](https://jp-minerals.org/vesta/en/) and installed it in the default location (`/Applications/VESTA.app`). 
1. Install [pymatgen-io-espresso](https://github.com/oashour/pymatgen-io-espresso) in a venv located at `$HOME/.vesta_espresso_venv`
2. Save the `vesta_espresso.sh` script from this repository somewhere convenient. You can download it from the [releases](https://github.com/oashour/vesta-espresso/releases) page.
3. In your `$HOME/.bashrc` file (or equivalent for your shell), add the following line
```bash
alias vesta='path/to/vesta_espresso.sh"
```
4. Run `source $HOME/.bashrc` (or equivalent) to load the changes. You can now use the `vesta` command to open Quantum ESPRESSO input files (and all other natively supported files) in VESTA. If you didn't add the alias, you have to call the script with its full path.
```bash
vesta pw.in test.cif POSCAR.vasp
```
5. Download the `.dmg` file from the [releases](https://github.com/oashour/vesta-espresso/releases) page. This contains a `.app` file that you can drag into your `/Applications/VESTA` folder. If you add this to your Dock, you'll be able to drag Quantum ESPRESSO input files with the extensions `.in` or `.pwi` onto the icon to open them in VESTA. You can also change the default app of those files in Finder to `VESTA (QE).app` to open them with VESTA.
6. If you would like to be able to right click `.in/.pwi` files and select "Open with -> VESTA (QE).app", then you will need to run the following command in your terminal:
```bash
/System/Library/Frameworks/CoreServices.framework/Versions/A/Frameworks/LaunchServices.framework/Versions/A/Support/lsregister -f /Applications/VESTA/VESTA\ \(QE\).app
```

Note that you can't use the `open` dialogue within VESTA to open Quantum ESPRESSO files. The script and App simply act as wrappers around VESTA. Either use Finder (whether by setting `VESTA (QE).app` as the default file or right click -> open with), or drag the file onto the App icon in the Dock.

## Customized Installation (Mac/Linux)
If you're not using a Mac or you just want to use the script without the App (i.e., no Finder integration), you can put your virtual environment and VESTA in any location you want. You just need to set the following environment variables in your `.bashrc` (or equivalent) file:

```bash
export PMG_VENV="path/to/your/pmg/venv/" # Not including bin/activate
export VESTA="path/to/vesta" # Vesta itself not the script
alias vesta='path/to/vesta_espresso.sh"
```

You can also repackage the Mac app using [Platypus](https://sveinbjorn.org/platypus) with the included `.platypus` profile to use a custom location for VESTA and the virtual environment.

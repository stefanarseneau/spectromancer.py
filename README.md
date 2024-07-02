# spectromancer.py
spectromancer.py is a wizard for managing MAGE spectrum observations of white dwarfs. It expects a directory structure in the form:
```
observations/
├── target1/
|   ├── exposure1.fits
│   ├── exposure2.fits
│   └── exposure1.fits
├── target2/
│   ├── exposure1.fits
│   ├── exposure2.fits
│   └── exposure3.fits
└── obsfile.csv
```
Then, spectromancer.py can be used to automatically coadd the spectra and measure parameters. Each target should have its own folder, which contains all of its individual exposures. The file `obsfile.csv` can contain any information, but must have a column called `name` that corresponds to the name of each target. An example observation and obsfile are provided in the repository. Then, spectromancer.py can be called from the command line or imported as a python module.

```
usage: spectromancer.py [-h] [--measure-rvs] [-o [OUTFILE]] [-v] path

spectromancer.py is a wizard for managing MAGE spectrum observations of white dwarfs

positional arguments:
  path                  path to the observation folder

options:
  -h, --help            show this help message and exit
  --measure-rvs         append RVs to the observation table?
  -o [OUTFILE], --outfile [OUTFILE]
                        file to which to write the observation table
  -v, --verbose         save plots?
```
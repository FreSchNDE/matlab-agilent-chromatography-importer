# MATLAB Agilent Chromatography Importer

Three self-contained MATLAB functions that read the binary data files an Agilent ChemStation / OpenLab LC(-MS) or GC(-MS) instrument writes into a `.D` data folder, and return their contents (signal data + metadata) as a plain MATLAB struct, plus a small helper to turn that metadata into readable text.

| Function | Reads | Returns |
|---|---|---|
| `importAgilentUV` | `*.UV` (diode-array 3D UV/Vis) | absorbance matrix (retention time × wavelength) + metadata |
| `importAgilentCH` | `*.CH` (single detector channel: DAD/MWD/ELSD/FID) | one signal trace over retention time + metadata |
| `importAgilentMS` | `*.MS` (mass spectrometer) | TIC + ion-abundance matrix (retention time × m/z, optionally sparse) + metadata |
| `describeFileContent` | the struct any of the above returns | a multiline `"structName.fieldName: value"` text summary |

Each function is a single, dependency-free `.m` file, copy the one you need.

## Requirements

- MATLAB **R2021a or newer** (uses `arguments` blocks and `name=value` call syntax). Written and tested in R2025b.
- (Development) The optional reference unit tests need R with the [chromConverter](https://github.com/ethanbass/chromConverter) package.

## Usage

```matlab
fileContent = importAgilentUV("path/to/DAD1.UV");        % scaled by default (see signal.unit)
fileContent = importAgilentUV("path/to/DAD1.UV", ApplyScaling=false);  % raw counts

fileContent = importAgilentCH("path/to/DAD1A.ch");
fileContent = importAgilentMS("path/to/MSD1.MS");
fileContent = importAgilentMS("path/to/MSD1.MS", Precision=0);  % unit-mass m/z axis
fileContent = importAgilentMS("path/to/MSD1.MS", Sparse=true);  % sparse ion-abundance matrix
```

A mass spectrum only records the ions detected in each scan, so for full-scan data the ion-abundance matrix is usually well over 90% zeros.
`Sparse=true` builds it directly as a sparse matrix, which takes a fraction of the memory (e.g. 9 MB instead of 107 MB); use `full()` for functions that don't accept sparse input.
When `Precision` merges several ions of one scan into the same m/z, their abundances are summed.

Every function returns a scalar struct with grouped metadata (`file`, `sample`, `method`, `instrument`, ...) and a `signal` group holding the data.
All structs returned by the same function share one fixed field layout independent of file version, so you can build an array of them.
Fields a file's version does not provide are set to a missing / `NaN`.

`describeFileContent` turns that metadata into a readable multiline string, e.g. to populate a free-text description field. Fields that aren't a single value (like `signal`, which holds vectors/matrices) are skipped automatically; use `ExcludeFields` to additionally drop scalar fields you don't want:

```matlab
text = describeFileContent(fileContent, ExcludeFields="method.datetime");
```

See the full field list and per-version notes in each function's help:

```matlab
help importAgilentUV
help importAgilentCH
help importAgilentMS
help describeFileContent
```

## Supported file versions

Only the versions produced by the instrument these were developed against are verified; the rest are implemented from format references but **untested** since no sample files were available.

| File | Versions | Tested |
|---|---|---|
| `.UV` | `131` (LC) | ✅ |
|       | `131` (OL), `31` | untested |
| `.CH` | `130` | ✅ |
|       | `30`, `8`, `81`, `179`, `181` | untested |
| `.MS` | `2` LC-MS | ✅ |
|       | `2` GC-MS | untested |

An unrecognized version raises an error.

## Development & Testing

The test class `tests/TestImportAgilent.m` runs over every `.D` folder found
under a dataset directory. It has two groups of checks:

- **Self-contained** checks that need no external reference and run on any `.D` folder:

  Each file imports without error, all structs from one importer share a fixed schema (and concatenate into an array), field types are stable, the `.CH` header min/max match the decoded signal, scaling is consistent, the sparse `.MS` import matches the full one, and coarser m/z rounding keeps each scan's total abundance.
- **Reference** checks that compare the decoded data bit-for-bit against chromConverter.
  These are automatically skipped (marked *incomplete*) if the reference files are not present.

### Run on your own `.D` folders

Point the test at your own data (no sample datasets are shipped with this repo):

```matlab
setenv("AGILENT_TEST_DATA", "C:\path\to\folder\of\D\folders");
runtests("tests")
```

`AGILENT_TEST_DATA` is any folder containing `.D` directories (searched recursively).
If it is unset, the test looks for a `datasets` folder next to the repo.
The self-contained checks will run against whatever files are found.

### (Optional) enable the bit-exact reference checks

Generate reference decodings with chromConverter, then run the tests:

```sh
# from the tests/ folder; needs R + chromConverter installed
Rscript generate_groundtruth.R  [datasetRoot]  [outputDir]
```

- `datasetRoot` defaults to `../datasets`; pass your own folder to match `AGILENT_TEST_DATA`.
- `outputDir` defaults to `<system temp>/agilent_importer_truth`, which the test reads automatically.
  Override with the `AGILENT_TEST_TRUTH` environment variable.

Files whose version chromConverter cannot read are simply skipped.

## Acknowledgments

The binary-format knowledge and several decoder routines were derived from these open-source projects:

- **chromConverter** (Ethan Bass): R package, GPL-3. <https://github.com/ethanbass/chromConverter>
- **rainbow** (Evan Shi): Python library + format documentation, GPL-3. <https://github.com/evanyeyeye/rainbow>, <https://rainbow-api.readthedocs.io>
- **Chromatography Toolbox** (James Dillon): MATLAB, MIT. <https://github.com/chemplexity/chromatography>
- **OpenChrom / ChemClipse** (Lablicate): used to cross-check the `.UV`/`.CH` format. <https://github.com/OpenChrom/openchrom>

## License

Parts of the `.UV`, `.CH` and `.MS` decoders were adapted from **chromConverter** and **rainbow**, both licensed **GPL-3.0**.
Because this code derives from GPL-3.0 sources, it is distributed under **GPL-3.0** as well.

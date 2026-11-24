# ASKAP-NAASHI spectral line validation report

`askapnaashi_validation.py` builds the NAASHI validation report: an HTML page (`index.html`) plus a `Figures/` folder. It is based on the ASKAP-FLASH script (`askapflash_validation_17sep2021.py`) and the report has almost the same structure. If you know the FLASH script, this README lists what you need to know to run this one, and what is different.


## 1. Requirements

- Python 3
- numpy, scipy, astropy, matplotlib, Pillow
- astroquery: only needed if you re-enable the RACS-low catalogue query (see section 6).
- Linux or macOS: the script calls `grep`, `sed`, `tail` and `rm`, as the FLASH script did.


## 2. How to run

```bash
# from inside the SBID directory (same as FLASH)
python askapnaashi_validation.py -s <SBID> -c <cal_SBID>

# from anywhere, pointing at the SBID directory (new)
python askapnaashi_validation.py -s <SBID> -c <cal_SBID> -d /path/to/<SBID>
```

| Option | Meaning |
|---|---|
| `-s`, `--sbid` | Science SBID (required) |
| `-c`, `--cal_sbid` | Calibrator SBID (required) |
| `-i`, `--imagebase` | Base string for image names (default `i.SB%s.cube`), as in FLASH |
| `-d`, `--data_dir` | **New.** SBID data directory. The script moves into it with `chdir`, so every relative path, and the outputs, are inside it. Default: `DATA_DIR` in the script (`'.'`) |

The outputs go into the SBID directory: `index.html`, and `Figures/`, with the noise-statistics outputs in `Figures/noise/`.


## 3. Where to edit settings

All the settings are at the top of the script, in two blocks.

### `>>> ADMIN: EDIT PATHS`

Every server-specific path and file name is here. Each line is marked `# >>> ADMIN: EDIT PATH`. All paths are relative to the SBID directory.

| Setting | Default | Used for |
|---|---|---|
| `DATA_DIR` | `'.'` | SBID directory. Can contain `{sbid}`, e.g. `'/data/naashi/{sbid}'` |
| `METADATA_DIR` | `metadata` | mslist, footprint, fieldlist files |
| `SLURM_DIR` | `slurmOutput` | `*.sh` (ASKAPsoft version) |
| `DIAGNOSTICS_DIR` | `diagnostics` | cubestats dirs, flagging summaries, beam plots |
| `CUBESTATS_DIR_GLOB` | `cubestats-*` | Every matching dir is searched, e.g. `cubestats-NAASHI_143`. The name does not have to match the field name |
| `FLAG_SUMMARY_GLOB` | `Flagging_Summaries/*averaged.ms.flagSummary` | Flagging statistics |
| `BEAMLOG_DIR` | `SpectralCube_BeamLogs` | Beam logs |
| `SOURCE_SPECTRA_DIR` / `NOISE_SPECTRA_DIR` | `SourceSpectra` / `NoiseSpectra` | Extracted spectra |
| `SOURCE_SPECTRA_FILE` / `NOISE_SPECTRA_FILE` | `spec_SB{sbid}_component_{cid}.fits` / `noise_SB{sbid}_component_{cid}.fits` | Spectrum file names |
| `SELAVY_CATALOGUE_GLOBS` | `casda/catalogues/selavy-image.i*.SB{sbid}.cont.*taylor.0.restored*.components.xml` (and 3 fallbacks) | Selavy component catalogue. The first match is used |
| `CONT_IMAGE_GLOBS` | `image.i.*SB{sbid}.cont.taylor.0.restored.conv.fits`, … | Continuum image |
| `FOOTPRINT_GLOBS`, `FIELDLIST_GLOB` | `metadata/footprintOutput-sb{sbid}-*.txt`, `metadata/fieldlist*.txt` | Beam footprint and field centre |
| `FIG_DIR`, `HTML_NAME` | `Figures`, `index.html` | Outputs |


## 4. Expected input layout (inside the SBID directory)

```
<SBID>/
├── metadata/
│   ├── mslist-*.txt                    (incl. mslist-scienceData*.txt)
│   ├── footprintOutput-sb<SBID>-*.txt
│   └── fieldlist*.txt                  (field centre; falls back to the mslist position)
├── slurmOutput/*.sh
├── diagnostics/
│   ├── cubestats-*/                    (cubeStats*linmos.contsub.txt, per-beam cubeStats, *.cube*linmos*.png)
│   ├── beamNoise*.png, beamMinMax*.png (may also be inside cubestats-*/)
│   └── Flagging_Summaries/*averaged.ms.flagSummary
├── SpectralCube_BeamLogs/beamlog.image.restored.*.txt
├── SourceSpectra/spec_SB<SBID>_component_<id>.fits
├── NoiseSpectra/noise_SB<SBID>_component_<id>.fits
├── casda/catalogues/selavy-image.i*.SB<SBID>.cont.*taylor.0.restored*.components.xml
├── image.i.*SB<SBID>.cont.taylor.0.restored.conv.fits
└── image.restored.<imagebase>.contsub.fits   (optional; otherwise the beam size comes from the beam logs)
```

The new sections are optional. If a file is missing, for example the spectra or the continuum image, the report is still built. It shows "not available" or "(not found)" in that place, or skips the section.

If a required file is missing (mslist, slurm `.sh`, linmos cubeStats, flagging summaries), the script stops with a message that names the missing pattern. It no longer fails with an `IndexError`.


## 5. What changed in the report compared with FLASH

- **Title:** "ASKAP NAASHI Spectral Line Data Validation Report", in blue. The bandpass-smoothing footnote was removed.
- **Mosaic Statistics:**
  - The top row shows only the Image Cube and the Continuum Subtracted Cube. The Residual Cube was removed.
  - The "Bandpass Smoothing Issue" cell was replaced by **Catalogue Statistics**, which shows:
    - number of Selavy components;
    - number of components brighter than 30 mJy (`flux_int`);
    - number of extracted spectra (files in `SourceSpectra/`).
- **Beams Statistics:** the Residual Beam Cube column was removed.
- **New section, Spectra Extracted from 15 Brightest Sources:**
  - Components 1a–15a, in a 5×3 grid.
  - Plotted from the FITS spectra in the style of E. Mahony's `plot_spec_local.py`: 8 stacked panels, a ±5σ noise band, and beam-forming intervals.
  - Name and flux in the plot title come from the Selavy catalogue.
- **New section, Spectral and Continuum Noise Statistics:** a port of H. Yoon's scripts 01–06, with five plots:
  1. Spectral Noise Inside vs Outside 3.2 deg
  2. RMS Noise vs Distance From Field Centre
  3. Spectral Noise Sky Distribution
  4. Continuum Noise Sky Distribution (Selavy `rms_image`)
  5. Distribution of Bright Components (continuum image + Selavy components + beam footprint)

  A "Click here" button under plot 1 downloads `Figures/noise/SB<SBID>_combined_all.txt`. Its columns are: component, RA, Dec, flux_peak, mid-third noise, continuum noise, distance.

  The median text files are also written to `Figures/noise/`, but they are not linked from the report.
- **RACS-low reference catalogue:** this replaces NVSS, but it is currently **disabled**. Every related line in the main program is commented out with a `# [RACS] ` prefix:
  - the `--catalogue` / `--no_catalogue` options;
  - the query block;
  - the report section.

  To re-enable it:
  ```bash
  sed -i 's/^# \[RACS\] //' askapnaashi_validation.py   # macOS: sed -i '' ...
  ```
  Then run with `--catalogue`, or set `QUERY_CATALOGUE = True`. The VizieR query is slow.


## 6. Calculation changes

**Changes to the FLASH code:**
- **Missing data channels:** a channel counts as missing when its Max (9th column of the linmos contsub cubeStats) is 0. This rule is unchanged from FLASH.
- **The n=1 bug:** when there were no missing or bad channels, FLASH displayed "n= 1" (the "none" placeholder was counted). It now shows "n= 0".
- **Other bugs fixed:**
  - a crash when a beam had fewer than 2 fully flagged baselines;
  - the flagged-antenna list showing "none" when exactly one antenna was flagged;
  - the conversion of negative Dec values.
- **Diagnostic PNGs** are matched by name (image/contsub), not by count or order. Every `cubestats-*` directory is searched.
- **Code clean-up:**
  - removed two leftover `print('yes')` debug lines;
  - removed two functions that were never called (`cal_binnedAvg`, `html_button`);
  - removed two unused imports (`import PIL`, `OptimizeWarning`).

  None of these changes the results.

**Changes to the original scripts 01–05:**
- **Channel ranges:** the first, mid and last ranges are each a third of the channels (`nchan // 3`), instead of a fixed 5184. For the full 15552-channel cube this gives the same 5184.
- **Distance from the field centre:** this is now a true angular separation, and RA offsets include cos(Dec). Inside/outside membership and the medians can therefore differ slightly from the original scripts. To reproduce the old flat-sky values, set `DISTANCE_FLAT_SKY = True`.
- **Script 04:** the plot title said "mid", but the code read the last third. It now uses the mid third.
- **Bright-component count:** uses `flux_int > 30 mJy`, matching how the pipeline selects sources for spectrum extraction.


### 7. Survey configuration

These are the survey values. NAASHI uses the same band as FLASH (711.5–999.5 MHz), so the FLASH values were kept:
- Tsys 60 K, beam-forming intervals every 5 MHz from 713.5 MHz;
- QC thresholds (beam logs, bad channels, noise rank);
- report title and contact.

Other settings for the new sections are in the same area:

| Setting | Default | Meaning |
|---|---|---|
| `SOURCE_COMPONENTS` | 1a–15a | Components shown in the spectra section |
| `CAT_BRIGHT_MJY` | 30 | Threshold in mJy for the bright-component count (uses `flux_int`) |
| `INNER_RADIUS_DEG` | 3.2 | Inside/outside split, in degrees from the field centre |
| `DISTANCE_FLAT_SKY` | `False` | `True` uses the original flat `sqrt(dRA² + dDec²)` distance |
| `CONT_NOISE_COLUMNS` | `['rms_image']` | Selavy column used for the continuum noise |
| `CONT_IMAGE_STEP` | 4 | Plot every Nth pixel of the continuum image, for speed |


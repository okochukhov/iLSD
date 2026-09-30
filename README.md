# iLSD: improved least-squares deconvolution

[iLSD](https://ui.adsabs.harvard.edu/abs/2010A%26A...524A...5K/abstract) is a Fortran 90 code for least-squares deconvolution (LSD) of stellar intensity and polarization spectra.
LSD describes an observation as a superposition of identical profiles, scaled by the weights of a line mask and shifted to the position of every line, and recovers the mean profile common to all lines by solving the corresponding linear least-squares problem.
The mean profile reaches a signal-to-noise ratio far above that of the individual lines, which is what makes LSD the standard tool for detecting weak stellar magnetic fields and for deriving the Stokes profiles used in Doppler and Zeeman-Doppler imaging.

Beyond the [classical single-profile deconvolution](https://ui.adsabs.harvard.edu/abs/1997MNRAS.291..658D/abstract), iLSD implements

- **multiprofile LSD** -- several mean profiles reconstructed simultaneously from different columns of line weights, for instance to separate lines of different chemical elements, strength, excitation or magnetic sensitivity;
- **logarithmic deconvolution** -- an intensity spectrum can be deconvolved as log(I/Ic) instead of the line depth 1-I/Ic, [a formulation](https://ui.adsabs.harvard.edu/abs/2024MNRAS.529.2071D/abstract) which treats superposition of blended and saturated lines more accurately;
- **regularization** -- first-order Tikhonov smoothing of the profile in velocity, which raises the S/N of the mean profile at a controlled cost in velocity resolution;
- **line weight adjustment** -- an iterative correction of the individual mask weights so that the LSD model fits the observation better;
- **full error information** -- error bars scaled by the quality of the fit and the complete covariance matrix of the reconstructed profiles.

## Requirements

A Fortran compiler. The Makefile uses `gfortran`; any other standard-conforming compiler will do. The BLAS/LAPACK routines needed to solve the LSD problem are bundled in `src/lineq4.f`.

## Compiling

Run `make` in `src/`; it builds the executable `ilsd`.

## Running

    ilsd <configuration_file> [<Niter>]

The names of all output files are formed from the name of the configuration file with its extension replaced, so `ilsd test1.cfg` writes `test1.lsd`, `test1.mod` and so on.
A path is allowed and is preserved: `ilsd runs/star1.cfg` writes into the `runs` directory.
The optional second argument `<Niter>` is described under [Line weight adjustment](#line-weight-adjustment) below.

The `examples/` directory holds two ready-to-run cases, together with the output they produce.

## Configuration file

The configuration file (see `examples/test1.cfg`) holds one item per line:

1. name of the file with the observed spectrum, in quotes
2. name of the file with the line mask, in quotes
3. minimum velocity, maximum velocity and velocity step of the LSD profile grid, in km/s
4. interpretation of the input spectrum:
   - `0` -- intensity, deconvolved as the line depth 1-I/Ic
   - `1` -- polarization, V/Ic, Q/Ic or U/Ic
   - `2` -- intensity, deconvolved as log(I/Ic)
5. regularization parameter *(optional)*
6. output flag *(optional)*

**Velocity grid.** The number of velocity bins is (Vmax-Vmin)/Vstep rounded to the nearest whole number, plus one. The maximum velocity is then adjusted to fall on the grid, so it can come out slightly larger than requested; the code reports the grid it has actually set up.

**Input spectrum flag.** Mode 2 requires all intensities to be positive and all mask weights to be smaller than one, and converts the mask weights to the logarithmic scale as well. In every mode the profile in the output file is written on the same scale as the input spectrum, so a mode 0 or 2 profile has its continuum at 1 and a mode 1 profile has its continuum at 0.

**Regularization.** Omitting line 5, or setting it to 0, gives the unregularized solution. A positive value applies first-order Tikhonov regularization, which improves the S/N of the mean profile at the cost of smoothing in velocity. The parameter is dimensionless and has a comparable effect for different spectra and masks, but not for different velocity steps: its effect scales as 1/Vstep^2. There is no universally best value, but values of order 0.1 to 1 are a useful starting point for a grid step close to the spectral pixel.

A negative value instead introduces a light smoothing of the reconstructed  profile with a three-point kernel whose neighbour weight is read off the sub-diagonal of the autocorrelation matrix rather than chosen by the user.

**Output flag.** Line 6 controls how much is written; each level adds files to those of the level below:

| flag | files written |
|:----:|---------------|
| `0` (or absent) | `<prefix>.lsd` only |
| `1` | adds `<prefix>.wln` and `<prefix>.mod` |
| `2` | adds `<prefix>.cov` |

The `<prefix>.cor` file is written whenever the line weights are iterated, independently of this flag.
Because the two optional lines are read in order, an output flag can only be given if a regularization parameter is given before it; use 0 on line 5 if no regularization is wanted.

## Input files

- **Observed spectrum.** A 2 or 3 column table (see `examples/test1.obs`). The columns are wavelength, continuum normalized observation, and an optional error bar. When the third column is absent the code assumes the same error bar, 0.01, for all spectral points. Only those points that fall within the velocity range of at least one mask line are used; the remainder are discarded.

- **Line mask.** A table with at least 2 columns (see `examples/test1.lin`). The first is the line central wavelength, the second the line weight. Further weights may be given in the third, fourth and following columns for multiprofile LSD, in which case one mean profile is reconstructed per weight column. A weight of exactly zero excludes the line from calculation of the corresponding profile.

## Output files

- `<prefix>.lsd` -- the LSD profile(s). The first record gives the number of
  velocity bins, the number of profiles, an estimate of the S/N of the first
  profile, and the reduced chi-square of the fit to the observations. It is
  followed by one record per velocity bin, giving velocity in km/s, the profile,
  and its error bar. For multiprofile LSD each further profile is preceded by one
  record with its own S/N estimate and the same reduced chi-square.

- `<prefix>.mod` -- the LSD approximation of the retained parts of the input
  spectrum. The columns give wavelength, LSD model, observed spectrum and
  uncertainty of the observed spectrum, all on the scale of the input observations.

- `<prefix>.wln` -- the total observational weight, the sum of 1/sigma^2 over all
  spectral points a line contributes to, of every line in the mask. The columns
  give line wavelength and weight. This output is useful for computing weighted
  averages of different line parameters.

- `<prefix>.cov` -- the full covariance matrix of the reconstructed profiles. The
  matrix is NOT scaled by the reduced chi-square, whereas the error bars in the
  `.lsd` file are.

- `<prefix>.cor` -- the accumulated line weight corrections, when the weights have
  been iterated. The columns give line wavelength and the total additive
  correction applied to its weight, so that the final weight is the one in the
  mask plus this number. For a mask with several weight columns the correction
  belongs to the first column in which the line has a non-zero weight, and in
  mode 2 it is additive on the logarithmic scale in which that mode works.

## Line weight adjustment

With a second argument the code repeats the deconvolution that many times, each pass correcting every mask weight to fit the observations better and recomputing the profile from the corrected mask.
All output refers to the last pass, the accumulated corrections go to `<prefix>.cor`, and most of the gain comes in the first two or three passes.
A warning is printed if a pass makes the fit clearly worse, which means the adjustment is not converging and its result should not be used.

The observations constrain the summed weight of a group of blended lines far better than the individual values, so the corrections are held near the input mask by an internal prior.
Even so, corrections of heavily blended lines are not robust individually and should not be used as a measure of individual line strengths.

## Important notes

- It is entirely the responsibility of the user to make sure that LSD profile(s)
  are fitted to sensible observed spectra. The program does not perform outlier
  rejection, removal of zero or NaN pixels, or removal of spectral regions
  affected by telluric lines.

- It is entirely the responsibility of the user to provide appropriate line
  weights in the input file. The program does not modify or renormalize line
  weights and does not renormalize, shift or continuum-correct the calculated
  LSD profiles.

- The error bars of the LSD profile(s) are scaled by the square root of the
  reduced chi-square of the LSD model to account for possible systematic errors
  of the LSD approximation and/or a wrong estimate of the input spectrum
  uncertainties. The chi-square is saved together with the LSD profiles, allowing
  one to recover the unscaled error bars.

- The amplitude of the LSD profile(s) is fixed by the normalization of the mask
  weights applied by the user, not by the observations or by the code. The
  normalization used should be stated when LSD profiles or quantities derived
  from them are published (see the recommendations in Sect. 2.5 of Kochukhov et al. 2010).

## Author and citation

iLSD is developed by [Oleg Kochukhov](https://www.astro.uu.se/~oleg) (Department of Physics and Astronomy, Uppsala University, Sweden).

If you use iLSD in your research, please cite [Kochukhov et al. (2010)](https://ui.adsabs.harvard.edu/abs/2010A%26A...524A...5K/abstract).

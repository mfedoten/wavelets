# Continuous wavelet transform

This is a MATLAB script I'm using to obtain continuous wavelet transform (CWT).  It uses built-in MATLAB functions to calculate the transform (cwt.m and cwtft.m), the main interest here is how to chose scales/frequency and how to compute cone of influence (COI).

This function allows two ways of computing CWT:
- straightforward, based on convolution;
- more computationally efficient, based on FFT.

The choice of scales can be done also using two approaches:
- construct linearly spaced *frequency* vector and then convert it to scales;
- geometrical spacing of *scales* as described in [1].

Cone of influence is the region of the wavelet transform, which is influenced by edge effect. COI is computed according to the method used for calculation of CWT (convolution- or FFT-based):
- convolution-based: COI is defined as regions which are large than wavelet support at each scale.
- FFT-based: COI is defined as e-folding time for the autocorrelation of wavelet power at each scale [2].

This function estimates maximum frequency (minimal scale) as half of Nyquist frequency or takes the one provided by user. Minimum frequency (maximal scale) is chosen such as COI at this scale/frequency would affect half of time points.

This functions returns scalogram, percentage energy for each coefficient of CWT. It also plots CWT (if such option is specified), all the values on the plot are **linear**.
Plot function displays COI as hatched regions, to do so an additional function is required. [Hatchfill function](http://www.mathworks.com/matlabcentral/fileexchange/30733-hatchfill) was used for that. I modifies this function slightly in order to control color of hatch lines and added to the repo for convenience. Otherwise, instead of using hatched regions, COI can be indicated by using MATLAB `patch` function with alpha set to a value less than 1.


## Syntax
```matlab
[tfr,t,f,scales,coi] = wt(sig,property_structure)
[tfr,t,f,scales,coi] = wt(sig,'propertyname',propertyvalue)
```

|Input  |                   |
|-------|-------------------|
|sig    | signal to analyse |
|opt    | options for wavelet transform, can be defined either as structure or as name-value pairs |

**Note:** for now, only analytical Morlet wavelet is supported for FFT-based CWT.

|Output |                   |
|-------|-------------------|
|tfr    | matrix with scalogram values (energy per coefficient of CWT), where rows correspond to frequency, columns to time |
|f      | frequency vector |
|t      | time vector |
|coi    | cone of influence of specified wavelet function.  Coefficients within the cone should be treated with care. Returns matrix of size `Nsc\*2`, where `Nsc` is number of scales, the columns correspond to left and right borders of COI |
|scales | vector of scales |

|Options  | Possible values      |                   |
|---------|----------------------|-------------------|
|wname    | string               | contains name of the wavelet to be used. For FFT-based only analytical Morlet (`morl`) is supported. For convolution-based CWT wavelet name can be any supported by MATLAB's `cwt`-function (see `waveinfo` for more info). |
|fs       | float                | sampling frequency (Hz) |
|type     | **'fft'**, 'conv'    | which type of CWT to use: convolution of FFT-based |
|sampling | **'freq'**, 'scales' | which sampling to use to construct vector of scales: linear sampling of frequencies or geometrical sampling of scales |
|fmax     | float                | CWT will be computed up to this frequency (should be less or equal to Nyquist frequency) |
| fmin    | float                | minimal frequency. If nothing is provided, then it is chosen such as only 50% of coefficients are affected by the COI. |
|fstep    | float                | if sampling was chosen as 'freq', specifies frequency resolution |
|nscales  | int                  | if sampling was chosen as 'scales', specifies desired number of scales |
|chi      | int (default is 85)  | if sampling was chosen as 'scales', specifies the percentage of overlap between wavelet basis, see [1] for details. |
|F        | vector of floats     | can be used to create vector of scales, by conversion instead of creating vector inside the function (substitute 'fstep' option) |
|norm     | boolean (**False**)  | if `True`, then each row of resulting CWT will be normalized by the corresponding scale. |
|plot     | boolean              | if `True` plots CWT with contour plots, values and scales/frequencies are linearly scaled |


## Example
Here is a simple example of the function's use. More examples are in the examples folder.

```matlab
% create chirp signal
N  = 2048;
fs = 1024;
t  = (0:N-1)/fs;
x  = chirp(t,0,t(end),fs/2,'quadratic')';
x  = x + 0.5*randn(size(x));
figure; plot(t,x);
% get cwt from min frequency to fs/2 with default step and display results
[WT,f_wt,t_wt,coi,scales] = wt(x,'type','fft','fs',fs,'plot',1);
% or you could specify options in the form of structure
opt = struct('type','fft','fs',fs,'plot',1);
[WT,f,t,coi,scales] = wt(x,opt);
```
![Example signal](/examples/example_sig.png)
![Example CWT](/examples/example.png)


## References
1. Jordan, D., Miksad, R. W. & Powers, E. J. Implementation of the continuous wavelet transform for digital time series analysis. Review of Scientific Instruments 68, 1484 (1997).
2. Torrence, C. & Compo, G. P. A Practical Guide to Wavelet Analysis.  Bulletin of the American Meteorological Society 79, 61-78 (1998).

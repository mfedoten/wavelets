# Continuous wavelet transform

This is MATLAB script I'm using for continuous wavelet transform (CWT).
It uses built-in MATLAB functions to calculate the transform (cwt.m and cwtft.m), the main interest here is how to chose scales/frequency and how to compute cone of influence (COI).

This function allows two ways of computing CWT:
- straightforward, based on convolution;
- more computationally efficient, based on FFT;

The choice of scales can be done also using two approaches:
- construct linearly spaced *frequency* vector and then convert it to scales;
- geometrical spacing of *scales* as described in [1].

Cone of influence is the region of the wavelet transform, which is influenced by edge effect. COI is computed according to the method used for calculation of CWT (convolution- or FFT-based):
- convolution-based: COI is defined as regions which are large than wavelet support at each scale.
- FFT-based: COI is defined as e-folding time for the autocorrelation of wavelet power at each scale [2].

## Syntax

<<<<<<< HEAD
### Input:

---
|sig    | signal to analyse |
|fs     | sampling frequency |
|wname  | name of the wavelet function |
|opt    | options for wavelet transform, can be defined either as structure or as name-value pairs |
---

OUTPUT: tfr - matrix with coefficients of wavelet transform, where rows 
              correspond to frequency, columns to time.
        f   - frequency vector;
        t   - time vector;
        coi - cone of influence of specified wavelet function.
              Coefficients within the cone should be treated with care.
              Returns matrix of size Nsc\*2. Where Nsc is number of
              scales,columns correspond to left and right borders of COI.
OPTIONS:
type:
sampling:
fmax:
fstep:
F:
nscales:
plot:

REFERENCES:
1. Jordan, D., Miksad, R. W. & Powers, E. J. Implementation of the 
   continuous wavelet transform for digital time series analysis. Review 
   of Scientific Instruments 68, 1484 (1997).
2. Torrence, C. & Compo, G. P. A Practical Guide to Wavelet Analysis. 
   Bulletin of the American Meteorological Society 79, 61?78 (1998).
=======
% INPUT: sig    - signal to analyse;
%        fs     - sampling frequency;
%        wname  - name of the wavelet function;
%        opt    - options for wavelet transform, can be defined either as
%                 structure or as name-value pairs.
% 
% OUTPUT: tfr - matrix with coefficients of wavelet transform, where rows 
%               correspond to frequency, columns to time.
%         f   - frequency vector;
%         t   - time vector;
%         coi - cone of influence of specified wavelet function.
%               Coefficients within the cone should be treated with care.
%               Returns matrix of size Nsc\*2. Where Nsc is number of
%               scales,columns correspond to left and right borders of COI.
% OPTIONS:
% type:
% sampling:
% fmax:
% fstep:
% F:
% nscales:
% plot:
%
% REFERENCES:
% 1. Jordan, D., Miksad, R. W. & Powers, E. J. Implementation of the 
%    continuous wavelet transform for digital time series analysis. Review 
%    of Scientific Instruments 68, 1484 (1997).
% 2. Torrence, C. & Compo, G. P. A Practical Guide to Wavelet Analysis. 
%    Bulletin of the American Meteorological Society 79, 61?78 (1998).
>>>>>>> ae8fb3d9e60c21d4c2a48551250a03959e3404eb

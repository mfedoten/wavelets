# Continuous wavelet transform

Here is MATLAB script I'm using for continuous wavelet transform (CWT).
It uses built-in MATLAB functions to calculate the transform (cwt.m and cwtft.m), the main interest here is how to chose scales/frequency and how to compute cone of influence (COI).

This function allows two ways of computing CWT:
- straightforward, based on convolution;
- more computationally efficient, based on FFT;

The choice of scales can be done also using two approaches:
- construct linearly spaced *frequency* vector and then convert it to scales;
- geometrical spacing of *scales* as described in [1].

Cone of influence is the region of the wavelet transform, which is influenced by edge effect. COI is computed according to the method used for calculation of CWT (convolution- or FFT-based):
- convolution-based: COI is defined as:
$ | t - u | \leq sB $

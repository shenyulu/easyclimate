easyclimate.filter.waveletFunctions
===================================

.. py:module:: easyclimate.filter.waveletFunctions


Attributes
----------

.. autoapisummary::

   easyclimate.filter.waveletFunctions.__author__


Functions
---------

.. autoapisummary::

   easyclimate.filter.waveletFunctions.wavelet
   easyclimate.filter.waveletFunctions.wave_bases
   easyclimate.filter.waveletFunctions.wave_signif
   easyclimate.filter.waveletFunctions.chisquare_inv
   easyclimate.filter.waveletFunctions.chisquare_solve


Module Contents
---------------

.. py:data:: __author__
   :value: 'Evgeniya Predybaylo, Michael von Papen'


.. py:function:: wavelet(Y, dt, pad=0, dj=-1, s0=-1, J1=-1, mother=-1, param=-1, freq=None)

   WAVELET  1D Wavelet transform with optional significance testing
   wave, period, scale, coi = wavelet(Y, dt, pad, dj, s0, J1, mother, param)

   Computes the wavelet transform of the vector Y (length N),
   with sampling rate DT.

   By default, the Morlet wavelet (k0=6) is used.
   The wavelet basis is normalized to have total energy=1 at all scales.

   INPUTS:

   Y = the time series of length N.
   DT = amount of time between each Y value, i.e. the sampling time.

   OUTPUTS:

   WAVE is the WAVELET transform of Y. This is a complex array
   of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
   ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
   The WAVELET power spectrum is ABS(WAVE)**2.
   Its units are sigma**2 (the time series variance).

   OPTIONAL INPUTS:

   *** Note *** if none of the optional variables is set up, then the program
   uses default values of -1.

   PAD = if set to 1 (default is 0), pad time series with zeroes to get
           N up to the next higher power of 2. This prevents wraparound
           from the end of the time series to the beginning, and also
           speeds up the FFT's used to do the wavelet transform.
           This will not eliminate all edge effects (see COI below).

   DJ = the spacing between discrete scales. Default is 0.25.
           A smaller # will give better scale resolution, but be slower to plot.

   S0 = the smallest scale of the wavelet.  Default is 2*DT.

   J1 = the # of scales minus one. Scales range from S0 up to S0*2**(J1*DJ),
       to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.

   MOTHER = the mother wavelet function.
               The choices are 'MORLET', 'PAUL', or 'DOG'

   PARAM = the mother wavelet parameter.
           For 'MORLET' this is k0 (wavenumber), default is 6.
           For 'PAUL' this is m (order), default is 4.
           For 'DOG' this is m (m-th derivative), default is 2.


   OPTIONAL OUTPUTS:

   PERIOD = the vector of "Fourier" periods (in time units) that corresponds
           to the SCALEs.

   SCALE = the vector of scale indices, given by S0*2**(j*DJ), j=0...J1
           where J1+1 is the total # of scales.

   COI = if specified, then return the Cone-of-Influence, which is a vector
       of N points that contains the maximum period of useful information
       at that particular time.
       Periods greater than this are subject to edge effects.


.. py:function:: wave_bases(mother, k, scale, param)

   WAVE_BASES  1D Wavelet functions Morlet, Paul, or DOG

   DAUGHTER,FOURIER_FACTOR,COI,DOFMIN = wave_bases(MOTHER,K,SCALE,PARAM)

   Computes the wavelet function as a function of Fourier frequency,
   used for the wavelet transform in Fourier space.
   (This program is called automatically by WAVELET)

   INPUTS:

   MOTHER = a string, equal to 'MORLET' or 'PAUL' or 'DOG'
   K = a vector, the Fourier frequencies at which to calculate the wavelet
   SCALE = a number, the wavelet scale
   PARAM = the nondimensional parameter for the wavelet function

   OUTPUTS:

   DAUGHTER = a vector, the wavelet function
   FOURIER_FACTOR = the ratio of Fourier period to scale
   COI = a number, the cone-of-influence size at the scale
   DOFMIN = a number, degrees of freedom for each point in the wavelet power
               (either 2 for Morlet and Paul, or 1 for the DOG)


.. py:function:: wave_signif(Y, dt, scale, sigtest=0, lag1=0.0, siglvl=0.95, dof=None, mother='MORLET', param=None, gws=None)

   WAVE_SIGNIF  Significance testing for the 1D Wavelet transform WAVELET

   SIGNIF = wave_signif(Y,DT,SCALE,SIGTEST,LAG1,SIGLVL,DOF,MOTHER,PARAM)

   INPUTS:

   Y = the time series, or, the VARIANCE of the time series.
       (If this is a single number, it is assumed to be the variance...)
   DT = amount of time between each Y value, i.e. the sampling time.
   SCALE = the vector of scale indices, from previous call to WAVELET.


   OUTPUTS:

   SIGNIF = significance levels as a function of SCALE
   FFT_THEOR = output theoretical red-noise spectrum as fn of PERIOD


   OPTIONAL INPUTS:
   SIGTEST = 0, 1, or 2.    If omitted, then assume 0.

           If 0 (the default), then just do a regular chi-square test,
               i.e. Eqn (18) from Torrence & Compo.
           If 1, then do a "time-average" test, i.e. Eqn (23).
               In this case, DOF should be set to NA, the number
               of local wavelet spectra that were averaged together.
               For the Global Wavelet Spectrum, this would be NA=N,
               where N is the number of points in your time series.
           If 2, then do a "scale-average" test, i.e. Eqns (25)-(28).
               In this case, DOF should be set to a
               two-element vector [S1,S2], which gives the scale
               range that was averaged together.
               e.g. if one scale-averaged scales between 2 and 8,
               then DOF=[2,8].

   LAG1 = LAG 1 Autocorrelation, used for SIGNIF levels. Default is 0.0

   SIGLVL = significance level to use. Default is 0.95

   DOF = degrees-of-freedom for signif test.
           IF SIGTEST=0, then (automatically) DOF = 2 (or 1 for MOTHER='DOG')
           IF SIGTEST=1, then DOF = NA, the number of times averaged together.
           IF SIGTEST=2, then DOF = [S1,S2], the range of scales averaged.

       Note: IF SIGTEST=1, then DOF can be a vector (same length as SCALEs),
           in which case NA is assumed to vary with SCALE.
           This allows one to average different numbers of times
           together at different scales, or to take into account
           things like the Cone of Influence.
           See discussion following Eqn (23) in Torrence & Compo.

   GWS = global wavelet spectrum, a vector of the same length as scale.
           If input then this is used as the theoretical background spectrum,
           rather than white or red noise.


.. py:function:: chisquare_inv(P, V)

   CHISQUARE_INV  Inverse of chi-square cumulative distribution function (cdf).

     X = chisquare_inv(P,V) returns the inverse of chi-square cdf with V
     degrees of freedom at fraction P.
     This means that P*100 percent of the distribution lies between 0 and X.

     To check, the answer should satisfy:   P==gammainc(X/2,V/2)

   Uses FMIN and CHISQUARE_SOLVE


.. py:function:: chisquare_solve(XGUESS, P, V)

   CHISQUARE_SOLVE  Internal function used by CHISQUARE_INV

     PDIFF=chisquare_solve(XGUESS,P,V)  Given XGUESS, a percentile P,
     and degrees-of-freedom V, return the difference between
     calculated percentile and P.

   Uses GAMMAINC

   Written January 1998 by C. Torrence

   extra factor of V is necessary because X is Normalized



This repository contains C# code which implements the DCTR spectral calculation method:

"A discrete cosine transform-based rapid algorithm for high-resolution, full-spectrum calculations over inhomogeneous gas paths"
Sudarshan P. Bharadwaj
Journal of Quantitative Spectroscopy and Radiative Transfer
Volume 316, April 2024, 108895

(or alternatively)

https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4549617

The DCTR method is a fast line-by-line spectral calculation method, which convolves Lorentz and Doppler line profiles in the frequency domain, using the Discrete Cosine Transform (DCT). The method can handle inhomogeneous gas paths (specified as a set of homogeneous slabs) in a single shot, unlike other methods which require an individual pass per homogeneous slab. The method was compared against the fast spectral solver RADIS. The benchmark for both DCTR and RADIS was the HAPI (HITRAN Application Programming Interface) line-by-line solver.

The advantages of DCTR over RADIS are:

1. DCTR is 10 times faster than RADIS for the same problem class (speed comparisons reflect the "full computation time" including overheads such as file reads and writes, and are made without the use of GPUs)
2. DCTR is 10 to 20 times more accurate than RADIS, as compared to the benchmark HAPI
3. DCTR consumes the same, or somewhat less memory than RADIS
4. DCTR can handle ~3X longer spectral ranges per run than RADIS
5. DCTR is able to perform computations at spectral resolutions of up to 0.0001 cm-1, as compared to 0.001 cm-1 for RADIS, for the same amount of computer RAM (32 GB)
6. DCTR can handle inhomogeneous gas paths, consisting of any number of homogeneous slabs, in a single pass, as compared to RADIS, which needs multiple passes (one pass per slab)
7. DCTR uses the Discrete Cosine Transform (as compared to RADIS, which uses the Fourier Transform) - as a result:
   1. There is a small additional speed gain for DCTR (on top of speed gains from other aspects of the DCTR algorithm)
   2. DCT calculations only involve real numbers, as opposed to complex numbers for FFTs
   3. Only one half of the Lorentz profile needs to be computed for each spectral line, which doubles the effective line wing  
8. RADIS utilizes a third-party Python library (SYMPY) to estimate line intensity split weights between bins, whereas DCTR calculates these weights analytically, thereby saving significant CPU time

More details may be found in the paper cited above.

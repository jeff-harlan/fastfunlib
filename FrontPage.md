A C library for fast multiprecision evaluation of transcendental functions, using fixed-point arithmetic directly on top of GMP/MPIR. The goal is to provide optimized base implementations of transcendental functions, with minimal overhead at low precision as well as asymptotic speed. It should be usable directly for computations that work well in fixed point (i.e. involving both inputs and outputs of order of magnitude close to unity), or as a base library for implementation of floating-point arithmetic with bells and whistles.

# Current status #

The current target is to test various algorithms, mainly for elementary functions. There is no code with a usable library interface yet. The svn trunk contains test implementations for exp (the same code may also compute cos or sin), log, and gamma.

## Exp ##

The test program computes exp(0.37) and compares both accuracy and speed to MPFR-2.4.1. This implementation is uniformly faster than MPFR up to at least several hundred bits of precision on one AMD64 system. Accuracy is within 1-2 bits of the correctly rounded value (as computed by MPFR). Times (5th and 6th columns) are in nanoseconds.

```
fredrik@airy:~/src/fastfunlib/fastfunlib/exptest$ ./exptest
 prec   acc   J   r     mpfr     this   faster
   53    52   1   2     5070     1240   4.089
   66    65   1   6     5200     1390   3.741
   82    81   1   8     5590     1460   3.829
  102   101   1   6     5810     1670   3.479
  127   126   2   8     7170     2280   3.145
  158   157   2   7     8390     2500   3.356
  197   196   1  11     8610     3030   2.842
  246   246   1   7     9870     3900   2.531
  307   307   2  12    12260     4440   2.761
  383   381   3   9    12860     5220   2.464
  478   477   2  11    16260     6580   2.471
  597   596   2  14    19220     8039   2.391
  746   746   3  13    23400    10500   2.229
  932   931   3   8    27900    14100   1.979
 1165  1164   3  14    35400    18700   1.893
 1456  1455   3  14    48500    27000   1.796
 1820  1819   3  12    67000    38000   1.763
 2275  2274   4  11    93000    56500   1.646
 2843  2842   4  17   130500    84000   1.554
 3553  3553   5  14   190500   124000   1.536
 4441  4440   6  14   278000   187000   1.487
 5551  5550   6  17   400500   279500   1.433
 6938  6937   6  19   599000   429500   1.395
 8672  8670   7  18   891000   645000   1.381
10840 10839   7  25  1310000   973000   1.346
13550 13549   8  23  2024500  1472000   1.375
16937 16935   9  25  3121000  2299000   1.358
21171 21170   9  27  4647000  3446000   1.349
26463 26461   9  36  5335000  5209500   1.024
```

The J and r parameters are algorithmic tuning parameters. The test program performs an exhaustive search and reports the best combination. Eventually a program for pre-tuning will have to be written.

This is still a naive implementation and optimistically the speed at low precision could be improved by perhaps an additional factor two.

## Log ##

Logarithm using Taylor series and argument reduction is twice as fast as MPFR up to a few hundred bits:

```
fredrik@airy:~/src/fastfunlib/fastfunlib/logtest$ ./logtest
 prec   acc   J   r     mpfr     this   faster
   53    53   1   1     4760     1750   2.720
   66    66   2   1     4800     1960   2.449
   82    82   1   3     5370     2390   2.247
  102   102   2   3     6630     2560   2.590
  127   127   3   2     8150     3310   2.462
  158   158   2   2     8160     3960   2.061
  197   197   2   2     9530     4760   2.002
  246   246   4   2    13430     6150   2.184
  307   307   4   2    15180     7300   2.079
  383   383   4   3    18120     9180   1.974
  478   478   4   3    23040    11340   2.032
  597   597   5   3    27760    14820   1.873
  746   746   6   4    30900    18900   1.635
  932   932   5   4    39800    25100   1.586
 1165  1165   5   5    49400    34600   1.428
 1456  1456   5   5    59500    46000   1.293
 1820  1820   5   7    79500    65000   1.223
 2275  2275   6   7   108000    92000   1.174
 2843  2843   7   9   133000   131500   1.011
 3553  3552   6  11   189500   188000   1.008
 4441  4441   7  11   253500   278500   0.910
 5551  5551   7  12   350500   409000   0.857
```

Logarithm using a 256KB (2^9 values x 4096 bits) precomputed lookup table buys another 2x speedup at low precision:

```
fredrik@airy:~/src/fastfunlib/fastfunlib/logtest2$ ./logtest2
 prec   acc   J   r     mpfr     this   faster
   53    53   1   0     4710      900   5.233
   66    66   1   0     4770      930   5.129
   82    82   1   0     5700     1040   5.481
  102   102   1   0     6660     1200   5.550
  127   127   1   0     7800     1430   5.455
  158   158   2   0     8660     1750   4.949
  197   197   1   0     9280     2170   4.276
  246   246   2   0    13480     2670   5.049
  307   307   2   0    15020     3120   4.814
  383   383   3   0    18100     4080   4.436
  478   478   3   0    23640     5120   4.617
  597   597   3   0    28620     6880   4.160
  746   746   4   0    31000     9400   3.298
  932   932   4   0    39800    12200   3.262
 1165  1165   4   0    50200    17600   2.852
 1456  1456   4   0    59000    25500   2.314
 1820  1820   4   0    80000    38500   2.078
 2275  2275   6   0   109000    56500   1.929
 2843  2843   6   0   134500    87500   1.537
 3553  3553   6   0   191000   133000   1.436
 4441  4441   7  11   254500   278500   0.914
 5551  5550   7  12   352000   410500   0.857
```

## Gamma function ##

The gammatest program implements the gamma function using Taylor series, designed for small arguments, e.g. x < 100 (for very large x, a second implementation based on Stirling's series should be used). This approach pays off greatly, giving 20-100x speedups over MPFR.

Times are for computing gamma(5.7) (Columns: prec, accurate, MPFR time, this time, speedup.)

```
fredrik@airy:~/src/fastfunlib/fastfunlib/gammatest$ ./gammatest
   53    53      76350       3150   24.238
   66    66      83210       3530   23.572
   82    82      90030       3770   23.881
  102   102     108540       4220   25.720
  127   127     121530       5850   20.774
  158   158     140270       6920   20.270
  197   197     174500       9090   19.197
  246   246     238250      12100   19.690
  307   307     309100      15800   19.563
  383   383     409400      19720   20.761
  478   478     604360      25260   23.926
  597   597     847160      33280   25.456
  746   746    1292500      48000   26.927
  932   932    1872700      66200   28.289
 1165  1165    3168700     111200   28.496
 1456  1456    5332000     171000   31.181
 1820  1820    9475000     288000   32.899
 2275  2275   16993000     474000   35.850
 2843  2843   30980500     809000   38.295
 3553  3553   56078000    1362000   41.173
 4441  4441  109397000    2362500   46.306
 5551  5551  219725000    3975500   55.270
 6938  6938  416183500    6976500   59.655
 8672  8672  835161500   11983500   69.693
10840 10840 1652142000   20552500   80.386
13550 13550 3450440500   34704000   99.425
```


Generating the coefficients for 1000-digit precision from scratch should take about a second (a small fraction of a second for lower precisions). Currently, the test code must read pre-generated coefficients from a file. The svn repository includes a data file generated with mpmath that works up to about 15000 bits.

There should be an interface for loading and saving coefficient data to files between sessions. The fastfunlib project could also provide a downloadable repository of expensive data, so that not every user who needs many extremely high-precision function evaluations has to generate the same data from scratch (which could take hours).

## Other functions ##

Initially exp, cosh/sinh, cos/sin, log, atan, zeta constants/Bernoulli numbers, gamma and digamma functions should be added; more functions may be added later on.
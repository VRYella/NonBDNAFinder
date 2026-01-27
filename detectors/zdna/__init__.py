"""Z-DNA detector module for detecting left-handed double helix structures.

This module implements Z-DNA detection using:
1. 10-mer scoring table from Ho et al. (1986)
2. eGZ-motif pattern recognition (Herbert 1997)

References:
    - Ho PS, Ellison MJ, Quigley GJ, Rich A. (1986). A computer aided thermodynamic 
      approach for predicting the formation of Z-DNA in naturally occurring sequences.
      EMBO J. 5(10):2737-44.
    - Herbert A, Alfken J, Kim YG, Mian IS, Nishikura K, Rich A. (1997). 
      A Z-DNA binding domain present in the human editing enzyme, double-stranded 
      RNA adenosine deaminase. Proc Natl Acad Sci U S A. 94(16):8421-6.
"""

from .detector import ZDNADetector

__all__ = ['ZDNADetector']

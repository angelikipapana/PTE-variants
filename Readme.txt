This directory includes an ensemble of matlab and R codes. 

The codes are developed for the computation of the following causality measures: 

1. transfer entropy (TE):  TE.m

2. partial transfer entropy (PTE): PTE.m

3. eight (8) connectivity-based variants of PTE (1A - 4B): connectivity_basedPTEvariants.m

4. three (3) RF-based variants of PTE (5A - 5C): RF_basedPTEvarinats.m
   (where RF stands for Random Forests)

The p-values of the null hypothesis H0: no causality, are extracted from an one-sided test on 'nsur' time shifted 
surrogates (randomizing the driving time series).

All measures are computed using the k-nearest neighbors' (KNN) estimator 
(Kraskov, A., Stögbauer, H., & Grassberger, P. (2004). Estimating mutual information. Physical review E, 69(6), 066138).


The auxiliary codes are necessary to run the above matlab functions.

The descripition of the PTE variants can be found in;

A. Papana, A. Papana-Dagiasis, E. Siggiridou. Shortcomings of transfer entropy and partial transfer 
entropy: Extending them to escape the curse of dimensionality. 
https://arxiv.org/abs/2004.11760
 
Please use the above reference if these codes are used for reported results.
 
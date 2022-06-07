# Structure-based virtual screening for PDL1 dimerizers: evaluating generic scoring functions

An innovative mechanism to inhibit the PD1/PDL1 interaction is PDL1 dimerization induced by small-molecule PDL1 binders. Structure-based virtual screening is a promising approach to discovering such small-molecule PD1/PDL1 inhibitors. Here we investigate which type of generic scoring functions is most suitable to tackle this problem. We consider CNN-Score, an ensemble of convolutional neural networks, as the representative of machine-learning scoring functions. We also evaluate Smina, a commonly used classical scoring function, and IFP, a top structural fingerprint similarity scoring function. These three types of scoring functions were evaluated on two test sets sharing the same set of small-molecule PD1/PDL1 inhibitors, but using different types of inactives: either true inactives (molecules with no in vitro PD1/PDL1 inhibition activity) or assumed inactives (property-matched decoy molecules generated from each active). On both test sets, CNN-Score performed much better than Smina, which in turn strongly outperformed IFP. The fact that the latter was the case, despite precluding any possibility of exploiting decoy bias, demonstrates the predictive value of CNN-Score for PDL1. These results suggest that re-scoring Smina-docked molecules with CNN-Score is a promising structure-based virtual screening method to discover new small-molecule inhibitors of this important therapeutic target.

# Paper Preprint 
The preprint version of the file can be found in [ChemRxiv](https://chemrxiv.org/engage/chemrxiv/article-details/623f1d1d8ab37367e372b017).

For reproducibility, we have created the binder for this github repo.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/sawsimeon/PDL1_Generic/HEAD)

# Introduction

This repository contains Julia files for some chapters in Bruce E. Hansen's *Econometrics*. They have (with the author's permission) been ported from the matlab files for chapters 3, 4, 8 and 10. Those original files and the data are found at [Hansen's website](https://users.ssc.wisc.edu/~bhansen/econometrics/).

The "translations" are pretty much line for line, with a few small tweaks when crucial (like avoiding some of the allocations in bootstrap simulations). That could be a start, at least since it probably mirrors the author’s intended approach. On the other hand, this could be rewritten in a much more Julian way. Maybe that could be round 2. Please let me know (or make a PR) if you want to contribute.

# On the Files
1. printmat.jl contains some functions for pretty-printing matrices.
2. ExtraFunctions.jl contains a few functions for reshuffling matrices and for emulating matlab's `quantile()` function.
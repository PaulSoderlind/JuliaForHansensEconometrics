# Introduction

This repository contains Julia files for some chapters in Bruce E. Hansen's *Econometrics*. They have (with the author's permission) been ported from the matlab files for some chapters. Those original files and the data are found at [Hansen's website](https://users.ssc.wisc.edu/~bhansen/econometrics/). The code here assumes that the data files are in the subfolder "Data" and that a subfolder "Results" is writable.

The "translations" are pretty much line for line, with a few small tweaks when crucial (like avoiding some of the allocations in bootstrap simulations). That could be a start, at least since it probably mirrors the author’s intended approach. On the other hand, this could be rewritten in a much more Julian way. Maybe that could be round 2. Please let me know (or make a PR) if you want to contribute.

The Julia files use either the .txt or the .dta data files. To load the latter, the [`ReadStatTables.jl`](https://github.com/junyuan-chen/ReadStatTables.jl) package is used.

# On the function files (in the subfolder "jlFiles"):
- printmat.jl contains some functions for pretty-printing matrices.
- ExtraFunctions.jl and UtilityFunctions.jl contain some functions for reshuffling data etc.

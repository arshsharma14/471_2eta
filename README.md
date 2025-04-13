# Growth curve analysis for E. coli with BrkA expression

This GitHub compiles the data and R scripts used to perform growth curve analysis and Western blot densitometry for the purpose of reproducibility. See R script for detailed comments.

## Abbreviations
- WT: E. coli BW25113
- Oliver: E. coli JW0052 from Oliver Micb 471 lab
- Tropini: E. coli JW0052 from Tropini lab
- WTBrkA, OliverBrkA, TropiniBrkA: respective strains transformed with pPALMC1, expressing BrkA

## Plate reader runs
Biological replicates and optimization performed as follows:
- 27Feb2025: WT, Oliver, WTBrkA, OliverBrkA
  - 3 technical replicates
  - 0.05% SDS and gradient EDTA
  - discarded due to low dynamic range from high SDS concentration
- 14Mar2025: WT, Tropini at gradient SDS for optimization, 3 technical replicates
- 21Mar2025: WT, Tropini, WTBrkA, TropiniBrkA
  - 3 technical replicates
  - 0.01% SDS and gradient EDTA
- 27Mar2025: WT, Oliver, Tropini, and all 3 with BrkA
  - 2 technical replicates each
  - 0.01% SDS and gradient EDTA
- 04Apr2025: same as above
- 10Apr2025: same as above

## Data wrangling
Biological replicates were compiled in Excel as follows:
- WT:
  1) 21Mar2025
  2) 27Mar2025
  3) 04Apr2025
  4) 10Apr2025
- Oliver:
  1) 27Mar2025
  2) 04Apr2025
  3) 10Apr2025
- Tropini:
  1) 21Mar2025
  2) 27Mar2025
  3) 04Apr2025
  4) 10Apr2025
- WTBrkA:
  1) 21Mar2025
  2) 27Mar2025
  3) 04Apr2025
  4) 10Apr2025
- OliverBrkA:
  1) 27Mar2025
  2) 04Apr2025
  3) 10Apr2025
- TropiniBrkA: (discard 04Apr2025 as outlier)
  1) 21Mar2025
  2) 27Mar2025
  3) 10Apr2025
Data was then loaded into R as .csv files and had subsequent analyses done.

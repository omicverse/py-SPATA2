# py-SPATA2 Benchmark

Dataset: 6402 observations x 100 genes.

| Operation | Seconds | Result shape |
|---|---:|---|
| `getCoordsDf` | 0.000418 | `(6402, 3)` |
| `getFeatureDf` | 0.001804 | `(6402, 6)` |
| `identifyTissueOutline` | 0.001991 | `(5, 2)` |
| `identifySpatialOutliers` | 0.010231 | `(6402,)` |
| `removeSpatialOutliers` | 0.005265 | `(6400, 100)` |

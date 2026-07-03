This is an AI generated R package. The package is currently Under Development and not validated.


# Example
```r
library(devtools)
install_github("PrabinDahal/parasiteclearance")
library(parasiteclearance)

# Example data
dat <- data.frame(
  id = c(
    1, 1, 1,
    5, 5, 5, 5, 5, 5, 5,
    8, 8, 8, 8, 8, 8, 8
  ),
  time = c(
    0, 6, 12,
    0, 6, 12, 18, 24, 36, 48,
    0, 6, 12, 18, 24, 36, 48
  ),
  parasitaemia = c(
    57500, 640, 48,
    2416, 1680, 960, 304, 128, 64, 16,
    1760, 1104, 560, 320, 288, 128, 32
  )
)

# Estimate Slope HL
estimate_clearance_batch(dat, detection_limit = 16)

# Output 
| id | status    | reason | model_type | tlag | clearance_rate_constant | slope_half_life |
|----|----------|--------|------------|------|--------------------------|------------------|
| 1  | estimated | ok     | linear     | 0    | 0.591                    | 1.17             |
| 5  | estimated | ok     | cubic      | 0    | 0.126                    | 5.48             |
| 8  | estimated | ok     | cubic      | 0    | 0.0789                   | 8.78             |


#==========================================
# How does it compare with WWARN PCE tool?
#==========================================
https://www.iddo.org/wwarn/parasite-clearance-estimator-pce

# Matches for id 1 and 8, not quite for id 5

| id | Clearance_rate_constant | Slope_half_life |
|----|--------------------------|------------------|
| 1  | 0.59                     | 1.17             |
| 5  | 0.11                     | 6.43             |
| 8  | 0.08                     | 8.78             |

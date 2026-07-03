# Example

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

estimate_clearance_batch(dat, detection_limit = 16)

# Output
> estimate_clearance_batch(dat, detection_limit = 16)
# A tibble: 3 × 7
     id status    reason model_type  tlag clearance_rate_constant slope_half_life
  <dbl> <chr>     <chr>  <chr>      <dbl>                   <dbl>           <dbl>
1     1 estimated ok     linear         0                  0.591             1.17
2     5 estimated ok     cubic          0                  0.126             5.48
3     8 estimated ok     cubic          0                  0.0789            8.78

## How does it compare with WWARN PCE tool?
# Matches for id 1 and 5.

https://www.iddo.org/wwarn/parasite-clearance-estimator-pce

| id | Clearance_rate_constant | Slope_half_life |
|----|--------------------------|------------------|
| 1  | 0.59                     | 1.17             |
| 5  | 0.11                     | 6.43             |
| 8  | 0.08                     | 8.78             |

# Aortic aneurysm progression data

This dataset contains longitudinal measurements of grades of aortic
aneurysms, measured by ultrasound examination of the diameter of the
aorta.

## Format

A data frame containing 4337 rows, with each row corresponding to an
ultrasound scan from one of 838 men over 65 years of age.

|         |           |                                      |
|---------|-----------|--------------------------------------|
| `ptnum` | (numeric) | Patient identification number        |
| `age`   | (numeric) | Recipient age at examination (years) |
| `diam`  | (numeric) | Aortic diameter                      |
| `state` | (numeric) | State of aneurysm.                   |

The states represent successive degrees of aneurysm severity, as
indicated by the aortic diameter.

|         |                   |          |
|---------|-------------------|----------|
| State 1 | Aneurysm-free     | \< 30 cm |
| State 2 | Mild aneurysm     | 30-44 cm |
| State 3 | Moderate aneurysm | 45-54 cm |
| State 4 | Severe aneurysm   | \> 55 cm |

683 of these men were aneurysm-free at age 65 and were re-screened every
two years. The remaining men were aneurysmal at entry and had successive
screens with frequency depending on the state of the aneurysm. Severe
aneurysms are repaired by surgery.

## Source

The Chichester, U.K. randomised controlled trial of screening for
abdominal aortic aneurysms by ultrasonography.

## References

Jackson, C.H., Sharples, L.D., Thompson, S.G. and Duffy, S.W. and Couto,
E. Multi-state Markov models for disease progression with classification
error. *The Statistician*, 52(2): 193–209 (2003)

Couto, E. and Duffy, S. W. and Ashton, H. A. and Walker, N. M. and
Myles, J. P. and Scott, R. A. P. and Thompson, S. G. (2002)
*Probabilities of progression of aortic aneurysms: estimates and
implications for screening policy* Journal of Medical Screening
9(1):40–42

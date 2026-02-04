# FEV1 measurements from lung transplant recipients

A series of measurements of the forced expiratory volume in one second
(FEV1) from lung transplant recipients, from six months onwards after
their transplant.

## Format

A data frame containing 5896 rows. There are 204 patients, the rows are
grouped by patient number and ordered by days after transplant. Each row
represents an examination and containing an additional covariate.

|         |           |                                                                                   |
|---------|-----------|-----------------------------------------------------------------------------------|
| `ptnum` | (numeric) | Patient identification number.                                                    |
| `days`  | (numeric) | Examination time (days after transplant).                                         |
| `fev`   | (numeric) | Percentage of baseline FEV1. A code of 999 indicates the patient's date of death. |
| `acute` | (numeric) | 0/1 indicator for whether the patient suffered an acute infection or rejection    |
|         |           | within 14 days of the visit.                                                      |

## Source

Papworth Hospital, U.K.

## Details

A baseline "normal" FEV1 for each individual is calculated using
measurements from the first six months after transplant. After six
months, as presented in this dataset, FEV1 is expressed as a percentage
of the baseline value.

FEV1 is monitored to diagnose bronchiolitis obliterans syndrome (BOS), a
long-term lung function decline, thought to be a form of chronic
rejection. Acute rejections and infections also affect the lung function
in the short term.

## References

Jackson, C.H. and Sharples, L.D. Hidden Markov models for the onset and
progression of bronchiolitis obliterans syndrome in lung transplant
recipients *Statistics in Medicine*, 21(1): 113–128 (2002).

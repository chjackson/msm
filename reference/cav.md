# Heart transplant monitoring data

A series of approximately yearly angiographic examinations of heart
transplant recipients. The state at each time is a grade of cardiac
allograft vasculopathy (CAV), a deterioration of the arterial walls.

## Format

A data frame containing 2846 rows. There are 622 patients, the rows are
grouped by patient number and ordered by years after transplant, with
each row representing an examination and containing additional
covariates.

|            |           |                                                                         |
|------------|-----------|-------------------------------------------------------------------------|
| `PTNUM`    | (numeric) | Patient identification number                                           |
| `age`      | (numeric) | Recipient age at examination (years)                                    |
| `years`    | (numeric) | Examination time (years after transplant)                               |
| `dage`     | (numeric) | Age of heart donor (years)                                              |
| `sex`      | (numeric) | sex (0=male, 1=female)                                                  |
| `pdiag`    | (factor)  | Primary diagnosis (reason for transplant)                               |
|            |           | IHD=ischaemic heart disease, IDC=idiopathic dilated cardiomyopathy.     |
| `cumrej`   | (numeric) | Cumulative number of acute rejection episodes                           |
| `state`    | (numeric) | State at the examination.                                               |
|            |           | State 1 represents no CAV, state 2 is mild/moderate CAV                 |
|            |           | and state 3 is severe CAV. State 4 indicates death.                     |
| `firstobs` | (numeric) | 0 = record represents an angiogram or date of death.                    |
|            |           | 1 = record represents transplant (patient's first observation)          |
| `statemax` | (numeric) | Maximum observed state so far for this patient (added in version 1.5.1) |

## Source

Papworth Hospital, U.K.

## References

Sharples, L.D. and Jackson, C.H. and Parameshwar, J. and Wallwork, J.
and Large, S.R. (2003). Diagnostic accuracy of coronary angiopathy and
risk factors for post-heart-transplant cardiac allograft vasculopathy.
Transplantation 76(4):679-82

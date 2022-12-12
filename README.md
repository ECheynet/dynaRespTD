# dynaRespTD
The coupled dynamic response of a suspension bridge to wind turbulence is computed in the time domain

[![View Buffeting response of a suspension bridge (time domain) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/73436-buffeting-response-of-a-suspension-bridge-time-domain)
[![DOI](https://zenodo.org/badge/249145922.svg)](https://zenodo.org/badge/latestdoi/249145922)


[![Donation](https://camo.githubusercontent.com/a37ab2f2f19af23730565736fb8621eea275aad02f649c8f96959f78388edf45/68747470733a2f2f77617265686f7573652d63616d6f2e636d68312e707366686f737465642e6f72672f316339333962613132323739393662383762623033636630323963313438323165616239616439312f3638373437343730373333613266326636393664363732653733363836393635366336343733326536393666326636323631363436373635326634343666366536313734363532643432373537393235333233303664363532353332333036313235333233303633366636363636363536353264373936353663366336663737363737323635363536653265373337363637)](https://www.buymeacoffee.com/echeynet)


## Summary

The lateral, vertical and torsional response of a suspension is computed in the time domain using computed turbulent velocity time histories. The simplified Bridge model of the Lysefjord Bridge [1] is considered for the modelling of the structure. Turbulence is modelled using the Kaimal model [2]. The quasi-steady theory, as well as the strip assumption, are used. Modal coupling between the lateral, vertical and torsional motions are accounted for in the time-domain model.

## Content

The present submission contains:
- The function dynaRespTD.m that computes the bridge response in the time domain (Non-linear  load + modal coupling)
- The function dynaResp_noCouplingTD.m that computes the bridge response in the time domain (linearised  load + no modal coupling)
- An example file  that compares the time-domain approach with the frequency domain approach 
- Various functions used for the example file, including simulation of correlated wind histories [3], computation of the bridge modal parameters [4], computation of the bridge response in the frequency domain [5].

Any suggestion, comment or question is welcomed.


## References

[1] Cheynet, E., Jakobsen, J. B., & Snæbjörnsson, J. (2016). Buffeting response of a suspension bridge in complex terrain. Engineering Structures, 128, 474-487.

[2] Kaimal, J. C., Wyngaard, J. C. J., Izumi, Y., & Coté, O. R. (1972). Spectral characteristics of surface‐layer turbulence. Quarterly Journal of the Royal Meteorological Society, 98(417), 563-589.

[3]https://se.mathworks.com/matlabcentral/fileexchange/68632-wind-field-simulation-the-fast-version 

[4] https://se.mathworks.com/matlabcentral/fileexchange/51815-calculation-of-the-modal-parameters-of-a-suspension-bridge 

[5] https://se.mathworks.com/matlabcentral/fileexchange/51970-buffeting-response-of-a-suspension-bridge-frequency-domain 

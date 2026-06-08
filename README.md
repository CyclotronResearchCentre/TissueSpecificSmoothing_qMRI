# TissueSpecificSmoothing_qMRI

## Overview

This repository contains the MATLAB and Python pipelines used to evaluate and compare several smoothing strategies for quantitative MRI (qMRI) data, including:

* **Tissue-Weighted Smoothing (TWS)**
* **gTSPOON (generalized Tissue SPecific smOOthing compeNsated)**
* **SUSANs (Smallest Univalue Segment Assimilating Nucleus smoothing)**

The objective is to investigate the impact of different smoothing approaches on quantitative MRI parameter maps while preserving tissue specificity and anatomical boundaries.

The framework was developed for whole-brain voxel-wise analyses of multiparameter mapping (MPM) data, including R1, R2*, MTsat and PD maps.

---

## Smoothing Methods

### Tissue-Weighted Smoothing (TWS)

TWS (Draganski et al. 2014) performs tissue-specific Gaussian smoothing using tissue probability maps (TPMs) as weighting functions. For each voxel, a normalized weighted average is computed:

$$
TWS(x) =
\frac{(g \ast (\omega\, s(\phi)))(x)}{(g \ast \omega)(x)}
\quad\text{with}\quad
\begin{cases}
TPM(x) > 0.05 \\
(g \ast \omega)(x) > 0.05
\end{cases}
$$

where:

* s(\phi) is the quantitative warped MRI map into standard space,
* \omega is the modulated tissue weights,
* g is a Gaussian kernel,
* * denotes convolution.

This approach minimizes signal contamination across tissue boundaries by restricting smoothing to voxels with similar tissue probabilities.

---

### gTSPOON

gTSPOON (TSPOON originally introduced by Lee et al. 2009) extends tissue-specific smoothing by incorporating probabilistic tissue information and Gaussian weighting while explicitly addressing partial-volume effects. The method aims to improve tissue specificity compared with conventional Gaussian filtering while maintaining robust noise reduction.

$$
\underbrace{
\frac{g * (M_{WM}\, s(\phi))(x)}{g * M_{WM}(x)}
}_{=\mathrm{TSPOON}(x)}
\;\xrightarrow[\;\;]{}\;
\underbrace{
\frac{g * (M_{TC}\, s(\phi))(x)}{g * M_{TC}(x)}
}_{=\mathrm{gTSPOON}(x)}
\quad\text{with}\quad
(g * M_{TC})(x) > 0.05
$$

where:
* s(\phi) is the quantitative maps warped into standard space,
* M_{TC} is the "majority and greater than 20%" tissue mask,
* g is a Gaussian kernel,
* * denotes convolution.

---

### SUSAN

The Smallest Univalue Segment Assimilating Nucleus (SUSAN) (Smith and Brady, 1997) is a non-linear edge-preserving smoothing algorithm implemented in FSL.

For a voxel (x), SUSANs can be expressed:

$$
SUSANs(x) =
\; \text{susan} \ast (M_{TC}\, s(\phi))(x)
$$

where:
* s(\phi) is the quantitative maps warped into standard space,
* M_{TC} is the "majority and greater than 20%" tissue mask,
* susan* denotes the SUSAN smoothing operator.


For qMRI applications, the brightness threshold ((bt)) was computed as:

[
bt = 0.75 \times \mathrm{median}(I_{valid})
]

where (I_{valid}) denotes positive finite voxels within the masked image. Spatial smoothing was specified by a target FWHM and converted to the SUSAN parameter (dt) using:

[
dt = \frac{\mathrm{FWHM}}{\sqrt{8\ln(2)}}.
]

---

## Comparison Framework

The different smoothing approaches were compared using voxel-wise analyses including voxel-wise log-likekihood, Bland-Altman plots performed on quantitative MRI maps. All smoothing methods were applied using matched smoothing scales to ensure fair comparisons.

---

## Statistical Analysis

Voxel-wise statistical analyses were performed using:

### SPM12

Statistical Parametric Mapping (SPM12) was used to build and estimate General Linear Models (GLMs). Parametric statistical inference was performed using Random Field Theory (RFT) under stationary and non stationary assumption procedures implemented in SPM.

### SnPM13

Statistical NonParametric Mapping (SnPM13) was used to perform permutation-based inference. SnPM provides robust family-wise error (FWE) control without relying on Gaussian random field assumptions and is particularly well suited for qMRI studies with moderate sample sizes or non-Gaussian residual distributions.

Both parametric (SPM) and non-parametric (SnPM) approaches were used to evaluate the robustness of smoothing-dependent findings.

---

## Software Requirements

* MATLAB (R2019b or later recommended)
* SPM12
* SnPM13
* FSL (for SUSAN smoothing)
* Python 3.x



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

[
TWS(x)=\frac{(g * (\omega s))(x)}
{(g * \omega)(x)}
]

where:

* (s) is the quantitative MRI map,
* (\omega) is the tissue probability map,
* (g) is a Gaussian kernel,
* (*) denotes convolution.

This approach minimizes signal contamination across tissue boundaries by restricting smoothing to voxels with similar tissue probabilities.

---

### gTSPOON

gTSPOON (TSPOON originally introduced by Lee et al. 2009) extends tissue-specific smoothing by incorporating probabilistic tissue information and Gaussian weighting while explicitly addressing partial-volume effects. The method aims to improve tissue specificity compared with conventional Gaussian filtering while maintaining robust noise reduction.

---

### SUSAN

SUSAN (Smith and Brady, 1997) is a non-linear edge-preserving smoothing algorithm implemented in FSL.

For a voxel (p) and neighboring voxel (x), SUSAN combines:

* a spatial weighting term:

[
c_s(p,x)=
\exp\left(
-\frac{|x-p|^2}{2dt^2}
\right)
]

* an intensity similarity term:

[
c_b(p,x)=
\exp\left(
-\left(
\frac{I(x)-I(p)}{bt}
\right)^2
\right)
]

The final weight is:

[
w(p,x)=c_s(p,x)c_b(p,x)
]

and the smoothed intensity is obtained as a normalized weighted average.

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

The different smoothing approaches were compared using voxel-wise analyses performed on quantitative MRI maps.

Comparisons focused on:

* preservation of tissue specificity,
* sensitivity to partial-volume effects,
* noise reduction performance,
* impact on statistical inference,
* spatial extent and localization of significant effects.

All smoothing methods were applied using matched smoothing scales to ensure fair comparisons.

---

## Statistical Analysis

Voxel-wise statistical analyses were performed using:

### SPM12

Statistical Parametric Mapping (SPM12) was used to build and estimate General Linear Models (GLMs), including:

* multiple regression analyses,
* group comparisons,
* covariate-based models.

Parametric statistical inference was performed using Random Field Theory (RFT) correction procedures implemented in SPM.

### SnPM13

Statistical NonParametric Mapping (SnPM13) was used to perform permutation-based inference.

SnPM provides robust family-wise error (FWE) control without relying on Gaussian random field assumptions and is particularly well suited for qMRI studies with moderate sample sizes or non-Gaussian residual distributions.

Both parametric (SPM) and non-parametric (SnPM) approaches were used to evaluate the robustness of smoothing-dependent findings.

---

## Software Requirements

* MATLAB (R2021a or later recommended)
* SPM12
* SnPM13
* FSL (for SUSAN smoothing)
* Python 3.x

  * NumPy
  * Pandas
  * NiBabel

---

## References

Smith SM, Brady JM. *SUSAN—A New Approach to Low Level Image Processing*. International Journal of Computer Vision. 1997;23(1):45–78. doi:10.1023/A:1007963824710

Chen Y, Hopp FR, Malik M, Wang PT, Woodman K, Youk S, Weber R. *Reproducing FSL's fMRI Data Analysis via Nipype: Relevance, Challenges, and Solutions*. Frontiers in Neuroimaging. 2022;1:953215. doi:10.3389/fnimg.2022.953215

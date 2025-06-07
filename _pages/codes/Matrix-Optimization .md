---
layout: single
title: "Codes under the Matrix Optimization Project"
permalink: /Codes/Matrix-Optimization/
author_profile: true
---

- [QSDPNAL (version 1.0): a MATLAB software for solving convex quadratic  semidefinite programming (QSDP)](https://github.com/MatOpt/QSDPNAL) ([click here for an introduction on how to use the package](https://blog.nus.edu.sg/mattohkc/softwares/qsdpnal/))  [<span style="color:red">CAUTION</span>: this software is for research purpose. It is neither intended nor designed to be a general purpose software at the moment.] For the details of the software, please check the following papers:
 

  [[Xudong Li](https://www.lixudong.info/), Defeng Sun, and [Kim Chuan Toh](https://blog.nus.edu.sg/mattohkc/),  “[QSDPNAL: A two-phase augmented Lagrangian method for convex quadratic semidefinite programming](https://www.polyu.edu.hk/ama/profile/dfsun/Li_et_al-2018-Mathematical_Programming_Computation.pdf)”, **Mathematical Programming Computation**, 10 (2018) 703--743.]

  [[Xudong Li](https://www.lixudong.info/), Defeng Sun, and [Kim Chuan Toh](https://blog.nus.edu.sg/mattohkc/), “[A block symmetric Gauss-Seidel decomposition theorem for convex composite quadratic programming and its applications](https://www.polyu.edu.hk/ama/profile/dfsun/Li2019_Article_ABlockSymmetricGaussSeidelDeco.pdf)”, **Mathematical Programming** 175 (2019) 395--418. [arXiv:1703.06629](https://arxiv.org/abs/1703.06629)]

- <a href="{{ '/files/SDPNAL+v1.0.zip' | relative_url }}" download>SDPNAL+</a>: a MATLAB software for solving large scale semidefinite programming with bound constraints ([click here for an introduction on how to use the package](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)) [awarded the triennial [Beale–Orchard-Hays Prize](https://www.mathopt.org/?nav=boh) for Excellence in Computational Mathematical Programming by the [Mathematical Optimization Society](https://www.mathopt.org/) at Bordeaux, France, July 2-6, 2018. See [Picture 1]({{ '/files/beale-orchard_hays-award2018.jpg' | relative_url }}), [Picture 2]({{ '/files/Ceremony_BOH.jpeg' | relative_url }}), and [Picture 3]({{ '/files/BOH_MedalSunDF.jpeg' | relative_url }}).]  [CAUTION: this software is NOT designed for solving small to medium sized SDP problems, for which interior point methods based software such as [SDPT3](https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/) is a better option.] For the details of the software, please check the following papers:

  [Defeng Sun, [Kim Chuan Toh](https://blog.nus.edu.sg/mattohkc/), [Y.C. Yuan](https://www.polyu.edu.hk/ama/people/academic-staff/dr-yuan-yancheng/?sc_lang=en), [Xinyuan Zhao](https://scholar.google.com/citations?user=nFG8lEYAAAAJ&hl=en), [SDPNAL+: A Matlab software for semidefinite programming with bound constraints (version 1.0)]({{ '/files/SDPNALplus-OMS-revision-2.pdf' | relative_url }}), to appear in **Optimization Methods and Software** (2019).]

  [Liuqin Yang, Defeng Sun, and [Kim Chuan Toh](https://blog.nus.edu.sg/mattohkc/), [SDPNAL+: a majorized semismooth Newton-CG augmented Lagrangian method for semidefinite programming with nonnegative constraints]({{ '/files/SDPNAL+.pdf' | relative_url }}), **Mathematical Programming Computation**, 7 (2015), pp. 331-366.]

  [Defeng Sun, [Kim Chuan Toh](https://blog.nus.edu.sg/mattohkc/), and Liuqin Yang, “[A convergent 3-block semi-proximal alternating direction method of multipliers for conic programming with 4-type constraints]({{ '/files/A%20CONVERGENT%203-BLOCK%20SEMIPROXIMALADMM2015.pdf' | relative_url }})”, **SIAM Journal on Optimization** Vol. 25, No. 2 (2015) 882–915. [Detailed computational results for over 400 problems tested in the paper]({{ '/files/PADMM3c-full-tables.pdf' | relative_url }}). You may also find [a supplementary note here]({{ '/files/Comparing-different-ADMMs.pdf' | relative_url }}) on more detailed comparisons between the performance of our proposed algorithm and various variants of ADMMs.]

  [[Xinyuan Zhao](https://scholar.google.com/citations?user=nFG8lEYAAAAJ&hl=en), D.F. Sun, and [Kim Chuan Toh](https://blog.nus.edu.sg/mattohkc/), [A Newton-CG augmented Lagrangian method for semidefinite programming]({{ '/files/NewtonCGAugLag.pdf' | relative_url }}), **SIAM Journal on Optimization**, 20 (2010), pp. 1737--1765.]

- **"Solving log-determinant optimization problems by a Newton-CG proximal point algorithm"**. See the brief user's guide [logdet-0-guide.pdf]({{ '/files/logdet-0-guide.pdf' | relative_url }}).

- <a href="{{ '/files/CorMatHdm_general.m' | relative_url }}" download>CorMatHdm_general.m</a> Computing the **H-weighted Nearest Correlation Matrix with fixed elements and lower and upper bounds** [H should not have too many zero elements for better numerical performance; otherwise, see CaliMatHdm] Testing example: <a href="{{ '/files/testCorMatHdm_general.m' | relative_url }}" download>testCorMatHdm_general.m</a> (uploaded on September 14, 2009).

- <a href="{{ '/files/CaliMatHdm.zip' | relative_url }}" download>CaliMatHdm.zip</a> Calibrating the **H-weighted Nearest Covariance Matrix** [H is allowed to have a large number of zero elements] (uploaded in April 2010).
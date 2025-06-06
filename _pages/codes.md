---
layout: single
title: "Codes in Matlab and others"
permalink: /Codes/
author_profile: true
---



- [Codes for nearest (covariance) correlation matrix problems](/Codes/correlation-matrix/)  

- [Codes under the Matrix Optimization (<span style="color:red">MatOpt</span>) Project](/Codes/Matrix-Optimization/)  

- [Codes under the Statistical Optimization (<span style="color:red">StaOpt</span>) Project](/Codes/Statistical-Optimization/)  

- [Codes for rank constrained problems](/Codes/rank-constrained/)  

- [Codes for other problems](/Codes/other/) 


<!-- ## Codes for nearest (covariance) correlation matrix problems
- Codes for the Nearest Correlation Matrix problem (the problem was initially introduced by Prof. [Nick Higham](https://www.maths.manchester.ac.uk/~higham/)):  <a href="../files/CorrelationMatrix.m" download>CorrelationMatrix.m</a> is a Matlab code written for computing the nearest correlation matrix problem (first uploaded in August 2006; last updated on August 30, 2019). This code should be good enough for most Matlab users.  If your Matlab version is very low and you really need a faster code, you can download <a href="../files/mexeig.mexw64" download>mexeig.mexw64</a> (for win64 operating system) and if use win32 or Linux system, you need to download the installmex file <a href="../files/installmex.m" download>installmex.m</a> and the c-file <a href="../files/mexeig.c" download>mexeig.c</a> by running the installmex.m first. For a randomly generated  <span style="color:blue">3,000 by 3,000</span> pseudo correlation matrix (the code is insensitive to input data), the code needs <span style="color:red">24</span> seconds to reach a solution with the relative duality gap less than 1.0e-3 after 3 iterations and 43 seconds  with the relative duality gap less than 1.0e-10 after 6 iterations in my Dell Desktop with Intel (R) Core i7 processor and for an invalid <span style="color:red">10,000 by 10,000</span> pseudo correlation matrix, the code needs 15 minutes to reach a solution with the relative duality gap less than 1.0e-4 after 4 iterations and 24 minutes with the relative duality gap less than 1.0e-12 after 7 iterations. For practitioners, you may set the stopping criterion (relative duality gap) to stay between 1.0e-1 and 1.0e-3 to run the code (typically, 1 to 3 iterations). If you need a C/C++ code, download <a href="../files/main.c" download>main.c</a> and <a href="../files/main.h" download>main.h</a>, which were written by [Pawel Zaczkowski](https://www.linkedin.com/in/pawel-zaczkowski-13a6a233/?originalSubdomain=uk) under a summer research project. If you are a client to [The Numerical Algorithms Group](https://nag.com/) (NAG), you may also enjoy [their commercialized implementations](https://nag.com/IndustryArticles/Nearest_Correlation_Matrix.pdf). The code in R <a href="../files/CorrelationMatrix.R" download>CorrelationMatrix.R</a> was written by [Ying Cui](https://sites.google.com/site/optyingcui/) (last updated on August 31, 2019; for efficiency, please use [Microsoft R open](https://mran.microsoft.com/open)) and the code in Python <a href="../files/CorrelationMatrix.py" download>CorrelationMatrix.py</a> was written by [Yancheng Yuan](https://www.polyu.edu.hk/ama/people/academic-staff/dr-yuan-yancheng/?sc_lang=en) (last updated on May 11, 2017), <span style="color:blue">respectively.</span>
-  <a href="../files/CorNewton3.m" download>CorNewton3.m</a> Computing the **Nearest Correlation <span style="color:red">Matrix</span> with fixed diagonal and off diagonal elements** (uploaded on September 14, 2009). The code in **R** <a href="../files/CorNewton3.R" download>CorNewton3.R</a> was provided by Professor Luca Passalacqua ([luca.passalacqua@uniroma1.it](mailto:luca.passalacqua@uniroma1.it)) (uploaded on **October 7, 2016**; for efficiency, please use [Microsoft R open](https://mran.microsoft.com/open)).
- <a href="../files/CorNewton3_Wnorm.m" download>CorNewton3_Wnorm.m</a> Computing the **W-norm Nearest Correlation Matrix with fixed diagonal and off diagonal elements** Testing example: <a href="../files/testCorMatWnorm.m" download>testCorMatWnorm.m</a>(uploaded on September 14, 2009).
- <a href="../files/CorMatHdm.m" download>CorMatHdm.m</a> Calibrating the **H-weighted Nearest Correlation Matrix** Testing example: <a href="../files/testCorMatHdm.m" download>testCorMatHdm.m</a> (uploaded in June 2008; last updated on September 10, 2009).
- <a href="../files/CorMatHdm_general.m" download>CorMatHdm_general.m</a> Computing the **H-weighted Nearest Correlation Matrix with fixed elements and lower and upper bounds** [H should not have too many zero elements for better numerical performance; otherwise, see CaliMatHdm] Testing example:  <a href="../files/testCorMatHdm_general.m" download>testCorMatHdm_general.m</a> (uploaded on September 14, 2009).
-  <a href="../files/LagDualNewton.m" download>LagDualNewton.m</a> (this is superseded by CorNewton3.m) Testing example: <a href="../files/testLagDualNewton.m" download>testLagDualNewton.m</a> (LagDualNewton method for the **Band Correlation Stress Testing**, "CorNewton1.m" will be called). 
- <a href="../files/CorNewtonSchur.m" download>CorNewtonSchur.m</a> Testing example: <a href="../files/testCorNewtonSchur.m" download>testCorNewtonSchur.m</a> (Schur decomposition based method for the **Local Correlation Stress Testing**, "CorNewton1.m" will be called).
- <a href="../files/AugLagNewton.m" download>AugLagNewton.m</a> (this is superseded by CorMatHdm_general.m) Testing example: <a href="../files/testAugLagNewton.m" download>testAugLagNewton.m</a> (AugLagNewton method for the **Band Correlation Stress Testing**, "CorNewton1.m" will be called). (uploaded in March 2007).
- <a href="../files/CaliMat1Mex.zip" download>CaliMat1Mex.zip</a> (Codes and testing example for) Calibrating **Covariance Matrix Problems with Inequality and/or Equality Constraints** (uploaded in April 2010).
- <a href="../files/CaliMatHdm.zip" download>CaliMatHdm.zip</a> Calibrating the **H-weighted Nearest Covariance Matrix** [H is allowed to have a large number of zero elements] (uploaded in April 2010).
- <a href="../files/Rank_CaliMat.zip" download>Rank_CaliMat.zip</a> Calibrating the **Nearest Correlation Matrix with Rank Constraints** (uploaded in April 2010).
- <a href="../files/Rank_CaliMatHdm.zip" download>Rank_CaliMatHdm.zip</a> Calibrating the **H-weighted Nearest Correlation Matrix with Rank Constraints** (uploaded in April 2010; last updated in October 2010 by including the refined Major codes).

---

## Codes under the Matrix Optimization (<span style="color:red">MatOpt</span>) Project
- [QSDPNAL (version 1.0): a MATLAB software for solving convex quadratic  semidefinite programming (QSDP)](https://github.com/MatOpt/QSDPNAL) ([click here for an introduction on how to use the package](https://blog.nus.edu.sg/mattohkc/softwares/qsdpnal/))  [<span style="color:red">CAUTION</span>: this software is for research purpose. It is neither intended nor designed to be a general purpose software at the moment.] For the details of the software, please check the following papers:
 

  [[Xudong Li](https://www.lixudong.info/), Defeng Sun, and [Kim Chuan Toh](https://blog.nus.edu.sg/mattohkc/),  “[QSDPNAL: A two-phase augmented Lagrangian method for convex quadratic semidefinite programming](https://www.polyu.edu.hk/ama/profile/dfsun/Li_et_al-2018-Mathematical_Programming_Computation.pdf)”, **Mathematical Programming Computation**, 10 (2018) 703--743.]

  [[Xudong Li](https://www.lixudong.info/), Defeng Sun, and [Kim Chuan Toh](https://blog.nus.edu.sg/mattohkc/), “[A block symmetric Gauss-Seidel decomposition theorem for convex composite quadratic programming and its applications](https://www.polyu.edu.hk/ama/profile/dfsun/Li2019_Article_ABlockSymmetricGaussSeidelDeco.pdf)”, **Mathematical Programming** 175 (2019) 395--418. [arXiv:1703.06629](https://arxiv.org/abs/1703.06629)]

- <a href="../files/SDPNAL+v1.0.zip" download>SDPNAL+</a>: a MATLAB software for solving large scale semidefinite programming with bound constraints ([click here for an introduction on how to use the package](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)) [awarded the triennial [Beale–Orchard-Hays Prize](https://www.mathopt.org/?nav=boh) for Excellence in Computational Mathematical Programming by the [Mathematical Optimization Society](https://www.mathopt.org/) at Bordeaux, France, July 2-6, 2018. See [Picture 1](../files/beale-orchard_hays-award2018.jpg), [Picture 2](../files/Ceremony_BOH.jpeg), and [Picture 3](../files/BOH_MedalSunDF.jpeg).]  [CAUTION: this software is NOT designed for solving small to medium sized SDP problems, for which interior point methods based software such as [SDPT3](https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/) is a better option.] For the details of the software, please check the following papers:

  [Defeng Sun, [Kim Chuan Toh](https://blog.nus.edu.sg/mattohkc/), [Y.C. Yuan](https://www.polyu.edu.hk/ama/people/academic-staff/dr-yuan-yancheng/?sc_lang=en), [Xinyuan Zhao](https://scholar.google.com/citations?user=nFG8lEYAAAAJ&hl=en), [SDPNAL+: A Matlab software for semidefinite programming with bound constraints (version 1.0)](../files/SDPNALplus-OMS-revision-2.pdf), to appear in **Optimization Methods and Software** (2019).]

  [Liuqin Yang, Defeng Sun, and [Kim Chuan Toh](https://blog.nus.edu.sg/mattohkc/), [SDPNAL+: a majorized semismooth Newton-CG augmented Lagrangian method for semidefinite programming with nonnegative constraints](../files/SDPNAL+.pdf), **Mathematical Programming Computation**, 7 (2015), pp. 331-366.]

  [Defeng Sun, [Kim Chuan Toh](https://blog.nus.edu.sg/mattohkc/), and Liuqin Yang, “[A convergent 3-block semi-proximal alternating direction method of multipliers for conic programming with 4-type constraints](../files/A%20CONVERGENT%203-BLOCK%20SEMIPROXIMALADMM2015.pdf)”, **SIAM Journal on Optimization** Vol. 25, No. 2 (2015) 882–915. [Detailed computational results for over 400 problems tested in the paper](../files/PADMM3c-full-tables.pdf). You may also find [a supplementary note here](../files/Comparing-different-ADMMs.pdf) on more detailed comparisons between the performance of our proposed algorithm and various variants of ADMMs.]

  [[Xinyuan Zhao](https://scholar.google.com/citations?user=nFG8lEYAAAAJ&hl=en), D.F. Sun, and [Kim Chuan Toh](https://blog.nus.edu.sg/mattohkc/), [A Newton-CG augmented Lagrangian method for semidefinite programming](../files/NewtonCGAugLag.pdf), **SIAM Journal on Optimization**, 20 (2010), pp. 1737--1765.]

- **"Solving log-determinant optimization problems by a Newton-CG proximal point algorithm"**. See the brief user's guide [logdet-0-guide.pdf](../files/logdet-0-guide.pdf).

- <a href="../files/CorMatHdm_general.m" download>CorMatHdm_general.m</a> Computing the **H-weighted Nearest Correlation Matrix with fixed elements and lower and upper bounds** [H should not have too many zero elements for better numerical performance; otherwise, see CaliMatHdm] Testing example: <a href="../files/testCorMatHdm_general.m" download>testCorMatHdm_general.m</a> (uploaded on September 14, 2009).

- <a href="../files/CaliMatHdm.zip" download>CaliMatHdm.zip</a> Calibrating the **H-weighted Nearest Covariance Matrix** [H is allowed to have a large number of zero elements] (uploaded in April 2010).

---

## Codes under the Statistical Optimization (<span style="color:red">StaOpt</span>) Project

- [**SuiteLasso**: a MATLAB suite for regression problems with generalized Lasso regularizers (GitHub)](https://github.com/MatOpt/SuiteLasso) [last updated in April 2021 with all source codes available]. [See the introduction on how to use it](https://github.com/MatOpt/SuiteLasso/blob/main/README.txt).

- [Square_Root_PMM](https://github.com/StatisticsLearningOPT/square_root_PMM): [A MATLAB software for square-root regression problems (GitHub)](https://github.com/StatisticsLearningOPT/square_root_PMM/blob/main/README.txt) [Last updated in January 2021]. Copyright (c) 2021 by Peipei Tang, Chengjing Wang, Defeng Sun, and Kim-Chuan Toh. This is a software package for solving the square-root regression problem:         min{ \|X \beta - b \|_2+\lambda p(\beta) - q(\beta)}. <span style="color:blue">For the details of the software, please check the following paper:</span>

  [Peipei Tang, Chengjing Wang, Defeng Sun, and [Kim Chuan Toh](https://blog.nus.edu.sg/mattohkc/),  “[A sparse semismooth Newton based proximal majorization-minimization algorithm for nonconvex square-root-loss regression problems](../files/19-247_Published.pdf)”, [Journal of Machine Learning Research](https://jmlr.org/papers/v21/19-247.html) 21(226):1--38, 2020.]

- <a href="../files/ConvexClustering.zip" download>**ConvexClustering**</a>: [a MATLAB package for convex clustering](https://blog.nus.edu.sg/mattohkc/softwares/convexclustering/) [last updated in June 2021]. [See the introduction on how to use it](https://blog.nus.edu.sg/mattohkc/softwares/convexclustering/).

---

## Codes for rank constrained problems
- <a href="../files/Rank_CaliMat.zip" download>Rank_CaliMat.zip</a> Calibrating the **Nearest Correlation Matrix with Rank Constraints** (uploaded in April 2010).

- <a href="../files/Rank_CaliMatHdm.zip" download>Rank_CaliMatHdm.zip</a> Calibrating the **H-weighted Nearest Correlation Matrix with Rank Constraints** (uploaded in April 2010; last updated in October 2010 by including the refined Major codes).

---

## Codes for other problems

- <a href="../files/IQEP_Newton.m" download>IQEP_Newton.m</a> Computing the **Inverse Quadratic Eigenvalue Problems** Testing example: <a href="../files/testIQEP_Newton.m" download>testIQEP_Newton.m</a> (uploaded in March 2008; last updated on July 15, 2016 by Ying Cui ([cuiying@u.nus.edu](mailto:cuiying@u.nus.edu))). -->
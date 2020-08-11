# HERMIT: A Benchmark Suite for the Internet of Medical Things [![GitHub release](https://img.shields.io/github/release-pre/ankurlimaye/HERMIT-BenchmarkSuite.svg)](https://github.com/ankurlimaye/HERMIT-BenchmarkSuite/releases) 


The HERMIT benchmark suite presents a collection of open-source and portable applications from the 
Internet of Medical Things (IoMT) domain. HERMIT's primary goal is to facilitate research into novel 
microarchitecture optimizations and right-provisioned system architectures for the IoT domain by 
providing a set of use-case applications. HERMIT applications span a variety of medical fields, 
including medical image processing algorithms, inverse Radon transform, implantable heart monitoring
 algorithms, activity monitoring, etc. HERMIT also include two supplementary applications for 
security and data compression that are essential for IoT. HERMIT applications are:
 
 | Medical Application | Description |
 |:-----------:|:-----------|
 | _activity_ | Physical Activity Estimation |
 | _apdet_ | Sleep Apnea Detection |
 | _hrv_ | Heart-rate Variability |
 | _imghist_ | Histogram Equalization |
 | _iradon_ | Inverse Radon Transform |
 | _kmeans_ | k-means Clustering |
 | _sqrs_ | QRS Detection in ECG |
 | _wabp_ | Arterial Blood Pressure Monitor |
 
 | IoT Application | Description |
 |:---------------:|:------------|
 | _aes_ | Advanced Encryption Standard |
 | _lzw_ | Lempel-Ziw-Welch Compression |
 


## How to cite

Ankur Limaye and Tosiron Adegbija, "[HERMIT: A Benchmark Suite for the Internet of Medical Things](https://ieeexplore.ieee.org/document/8392676)," in *IEEE Internet of Things Journal*, vol. 5, no. 5, pp. 4212 - 4222, October 2018.

[![HERMIT_DOI_badge](https://img.shields.io/badge/DOI-https%3A%2F%2Fdoi.org%2F10.1109%2FJIOT.2018.2849859-blue.svg)](https://doi.org/10.1109/JIOT.2018.2849859)
[![HERMIT_BibTexdownload](https://img.shields.io/badge/BibTex-download-blue.svg)](https://github.com/ankurlimaye/HERMIT-BenchmarkSuite/blob/master/CITATION.bib)

## Build Instructions

### Install dependencies

(Insert text here)

### Build HERMIT

*To build all kernels*:

    make all
    
*To build a subset of kernels*: (e.g., only build activity and imghist)

    make activity imghist
    
*Cross compilation of kernels*:

    make CC=/path/to/compiler

## Folder structure

    \bin (HERMIT binaries are copied to this folder after building)
    \include (Includes files necessary for building HERMIT applications)
    \inputs (Contains test-inputs)
    \prebuilt-bin-riscv (Contains HERMIT binaries cross-compiled for RISC-V)
    \prebuilt-bin-x86 (Contains HERMIT binaries cross-compiled for x86)
    \scripts (Contains scripts to run HERMIT binaries)
    \src (HERMIT application source files)

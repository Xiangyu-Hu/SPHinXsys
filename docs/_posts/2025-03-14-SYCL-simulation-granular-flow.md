---
layout: post
title: "Experience of SYCL-Based Development on Granular Material Large Deformation"
date: 2025-03-14
categories: sycl examples
---

# Experience of SYCL-Based Development on Granular Material Large Deformation

**Authors:** Shuang Li, Xiangyu Hu  

## Introduction  

### Aim  

This project focuses on **heterogeneous SPH simulation of granular materials**, particularly leveraging **GPU acceleration**.  

### Development Route  

We rewrote the existing **plasticity-related functions** in **SYCL** to achieve efficient GPU computation.  

### What You Can Gain  

- Access to an **open-source heterogeneous SPH simulation** framework for granular materials.  
- Insights into **SYCL-based development**, including errors I met and solutions.  

---

## New Test Cases  

### 1. **test_2d_column_collapse_sycl**  and **test_2d_column_collapse_ck**

### 2. **test_3d_repose_angle_sycl**   and **test_3d_repose_angle_ck**

> **Note:** These cases are SYCL-based implementations of their original versions.  

---

## GPU Acceleration Efficiency  

For a simulation with **1.6 million particles** (dp = 0.0015 m), the computation time for **0.1s** is:  

Device | Time (with ParticleSort) | Time (without ParticleSort)  
--------|--------------------------|-----------------------------  
**CPU** (Intel Xeon Platinum 8377C) | ~40 min | N/A  
**GTX 2080Ti** | ~16 min | N/A  
**GTX 1030** | ~210 min | ~280 min  

> **Conclusion:** GPU acceleration significantly outperforms CPU computation, making it crucial for **large-scale field simulations**.  

---

## Development Experience  

### (1) Variable Initialization Order vs. Definition Order  

**Issue:** The initialization order of variables differed from their definition order, leading to potential bugs.  

**Solution:** Add the following flag in `CMAKEPRESETS.json` to enforce correct ordering:  

CMAKE_CXX_FLAGS: "-Werror=reorder"  

---

### (2) Build Error on Linux in GitHub CI  

**Issue:** The test passes on Windows but fails on Linux.  

**Solution:** The **GCC compiler on Linux** enforces stricter rules and treats warnings as errors. To ensure compatibility, first build the project locally using **g++** and resolve any warnings before pushing to GitHub.  

---

### (3) Math Function Linking Errors in SYCL  

**Issue:** Functions like `sin()` caused linking errors in SYCL computations.  

**Solution:** Use SYCL's built-in math functions explicitly.  

#### **Modification in `base_data_type.h`**

# if SPHINXSYS_USE_FLOAT  
using Real = float;  
using UnsignedInt = u_int32_t;  
namespace math = sycl;  
# else  
using Real = double;  
using UnsignedInt = size_t;  
namespace math = std;  
# endif  

#### **Usage in Code**

Real sin_x = math::sin(x);  

---

### (4) Regression Testing on GitHub  

**Issue:** GitHub CI/CD only supports CPU-based testing, but the project requires SYCL-based regression tests.  

**Solution:**  

1. Generate the regression test database using **CPU**.  
2. Modify `CMAKEPRESETS.json` by removing the following line:  

   "SPHINXSYS_SYCL_TARGETS": "nvptx64-nvidia-cuda"  

3. With `SPHINXSYS_USE_SYCL="ON"`, the SYCL functions **will still be executed**, but on the CPU instead of GPU.  

---

## Conclusion  

This work demonstrates the feasibility and advantages of **SYCL-based GPU acceleration** in **granular material simulations**.

# 🚀 PBA-V3
Convergence Analysis of Parallax Bundle Adjustment (PBA) vs. Sparse Bundle Adjustment (SBA) Across 100+ BA Datasets

---

## 📌 Overview

PBA-V3 is a research-oriented project for:

🔍 Systematic convergence analysis of Bundle Adjustment (BA) algorithms, focusing on  
Parallax Bundle Adjustment (PBA) vs. Sparse Bundle Adjustment (SBA)

Unlike standard BA implementations, this project emphasizes:
- 📊 First-order and second-order optimization behavior
- 📈 Iteration-wise convergence diagnostics
- 🧪 Large-scale evaluation across 100+ datasets

---

## 🎯 Project Purpose

This project aims to answer key questions in BA optimization:

- Is PBA more stable than SBA?
- Does PBA converge faster or more robustly?
- How do their:
  - Gradient behaviors differ?
  - Hessian structures compare?
  - Numerical stability characteristics evolve?

To achieve this, we build a unified convergence analysis framework.
  

## 🖥 Tested Platforms
This project has been tested on:

- **Operating System:** Windows 11  
- **Compiler / IDE:** Visual Studio 2022 (Win32)  
- **CMake version:** ≥ 3.10  

> Note: Other platforms may work, but they have not been verified.


## 📦 Dependencies
- **Cholmod**: 1.5.0   
- **Eigen**: 3.3.5

---

## ⚙️ Installation

### 1️⃣ Clone the repository
```bash
git clone https://github.com/Polar-vision/PBA-V3.git
cd PBA-V3
```

### 2️⃣ Build the BA library (ba)
```bash
cd ba
mkdir build
cmake -B build -S . -G "Visual Studio 17 2022" -A Win32
cmake --build build --config Release
```

### 3️⃣ Build the demo (example)
```bash
cd ../example
mkdir build
cmake -B build -S . -G "Visual Studio 17 2022" -A Win32
cmake --build build --config Release
```

### 4️⃣ Run the Demo
```bash
build/Release/example.exe
```

💡 Note:
Make sure you have successfully built the ba library before running the demo.
The executable file (example.exe) will be generated in the build/Release directory.

# 🚀 Datasets

## 📥 Download
You can download the dataset package from the following link:  
[**Download Datasets**](https://drive.google.com/file/d/12bO4WTqCzckXtI5Bt97uV-gMgx5J_gyU/view?usp=sharing)

After downloading, **extract the contents into the `PBA-V3/` directory**.

The final folder structure should look like this:  
PBA-V3/  
├── datasets/  
│ ├── KD1-problem-83-48102/  
│ ├── ...  
├── 3rdparty/  
├── ba/    
└── example/    

> 💡 **Tip:**  
> Make sure the dataset folder is placed correctly before running the program, otherwise the relative paths may not resolve properly.  
> If you need to use a custom dataset path, you can modify it in  
> [`3rdparty/datapath.h`](./3rdparty/datapath.h).


---
## 🗂 Data Format

Each dataset is organized in the following structure:

- **Intrinsics (`cal.txt`)**  
  Contains the camera intrinsic parameters: `fx`, `fy`, `cx`, `cy` for each camera.

- **Extrinsics (`Cam.txt`)**  
  Contains the camera poses, including Euler angles (`ez`, `ey`, `ex`), the perspective center `(Xc, Yc, Zc)`, and the camera ID.  
  - `ez` : rotation around the z-axis  
  - `ey` : rotation around the y-axis  
  - `ex` : rotation around the x-axis

- **3D Points (`XYZ.txt`)**  
  Lists the 3D coordinates of the object points: `X`, `Y`, `Z`.

- **Feature Tracks (`Feature.txt`)**  
  Each line represents a feature track, including the number of views, the corresponding image indices, and the (u, v) coordinates in each image.

# 📊 Convergence Analysis Framework

We analyze BA optimization from **three perspectives**:

---

## 1️⃣ First-Order Metrics (Gradient Behavior)

### 🔹 Gradient Lipschitz Continuity
\[
L_k = \frac{\left\|\nabla f(x_k) - \nabla f(x_{k-1})\right\|}{\left\|x_k - x_{k-1}\right\|}
\]
- Measures gradient smoothness
- Large values indicate instability or oscillation

---

### 🔹 Descent Direction Quality
\[
\cos(\theta_k) = \frac{\nabla f(x_k)^\top \Delta x_k}{\left\|\nabla f(x_k)\right\| \left\|\Delta x_k\right\|}
\]
- $\approx -1 \rightarrow$ ideal descent
- $\approx 0 \rightarrow$ ineffective
- $> 0 \rightarrow$ wrong direction

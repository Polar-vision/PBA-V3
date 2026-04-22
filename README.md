# 🚀 PBA-V3
Convergence Analysis of Parallax Bundle Adjustment (PBA) vs. Sparse Bundle Adjustment (SBA) Across 100+ BA Datasets

---

## 📌 Overview

PBA-V3 is a research-oriented project for:

🔍 Systematic convergence analysis of Bundle Adjustment (BA) algorithms, focusing on  
PBA vs. SBA

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
$$
L_k = \frac{\left\|g(x_{next}) - g(x_{current})\right\|}{\left\|x_{next} - x_{current}\right\|}
$$

<div align="center">
  <img src="images/gradient_lipschitz.png" alt="Gradient Lipschitz Continuity Diagram" width="550">
  <p align="center"><em>Gradient Lipschitz continuity check in GN/LM optimization</em></p>
</div>

This metric quantifies the **local smoothness of the gradient**, measuring how rapidly the gradient changes with respect to the parameters. It directly reflects the curvature of the objective function, which is critical for analyzing optimization stability.

---

#### ❓ Is a large rate of change in the gradient between consecutive iterations necessarily bad?

Not necessarily. A large gradient update ratio is not inherently problematic, and can even be a sign of healthy progress:

- ✅ **Optimization early stage**: Large gradient changes often reflect rapid progress toward the solution, as the objective function’s curvature naturally produces steep gradients far from the minimum.
- ✅ **Valid large steps**: When the optimizer (e.g., Levenberg–Marquardt) takes large, valid steps, the gradient magnitude and direction can change significantly while still converging toward a minimum.

However, large gradient changes are concerning if they coincide with:

- ⚠️ Erratic oscillations in the residual or cost function value;
- ⚠️ A loss of descent direction quality (e.g., a positive dot product between the negative gradient and the update direction);
- ⚠️ No sustained reduction in the objective function, even over multiple iterations.

> **Summary**: A large gradient change is not a direct indicator of failure, but rather a signal that must be interpreted alongside descent direction quality, residual reduction, and the optimization stage.

---

### 🔹 Gradient Direction Quality
$$
\cos(\theta_k) = \frac{-\nabla f(x_k)^\top \Delta x_k}{\left\|\nabla f(x_k)\right\| \left\|\Delta x_k\right\|}
$$

This metric measures the alignment between the negative gradient (steepest descent direction) and the actual parameter update direction $\Delta x_k$. It quantifies how effectively the optimization step leverages the gradient to reduce the objective function.

- $\approx 1$ → **Ideal descent**: The update direction is nearly identical to the negative gradient, maximizing the first-order descent gain.
- $\approx 0$ → **Ineffective**: The update direction is nearly orthogonal to the gradient, resulting in negligible progress.
- $< 0$ → **Wrong direction**: The update direction has a positive dot product with the gradient, meaning it is aligned with the ascent direction and will increase the objective function value.

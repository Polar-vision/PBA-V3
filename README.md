# 🚀 Building from Source
**PBA-V3:** Convergence Analysis of PBA vs. SBA Across 100+ BA Datasets 

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

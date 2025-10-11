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

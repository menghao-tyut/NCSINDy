# NC-SINDy: Noise-Corrected Data-Driven Approach for Nonlinear Dynamical System Modeling in Frequency Domain

This repository contains the open-source code and supplementary materials for the paper **"A noise-corrected data-driven approach for nonlinear dynamical system modeling in frequency domain"**, including all numerical simulation scripts, predefined noise data, and reproducible results for Sections 4.1 and 4.2.

***

## 📄 Paper Information

* **Title**: A noise-corrected data-driven approach for nonlinear dynamical system modeling in frequency domain

* **Authors**: Jian Tang1, Jiahui Xu2, Hao Meng2,3,4,*, Zhihong Zhang1, Huidong Xu2,3,4, Xiaoyan Xiong5, Zhihua Wang2,3,4

* **Affiliations**:

  1. College of Mechanical Engineering, Taiyuan University of Technology, Taiyuan 030024, China

  2. College of Aeronautics and Astronautics, Taiyuan University of Technology, Taiyuan 030024, China

  3. Shanxi Key Laboratory of Material Strength & Structural Impact, Taiyuan University of Technology, Taiyuan 030024, China

  4. Shanxi Research Center of Basic Discipline of Mechanics, Taiyuan University of Technology, Taiyuan 030024, China

  5. College of Robotics Science and Engineering, Taiyuan University of Technology, Taiyuan 030024, China

* **Corresponding Author**: Hao Meng

  **Email**: `menghao@tyut.edu.cn`

  For questions about the code or paper, please contact Hao Meng at menghao@tyut.edu.cn

***

## 🛠️ Environment & Dependencies

This code is written for **Wolfram Mathematica 11** (and compatible with later versions).

* No additional Mathematica packages are required.

* `.xls` data files can be read directly in Mathematica.

***

## 📂 File Structure & Description

### &#x20;Detailed File Functions

1. **`NCSINDyPackage.m`** The core Mathematica package implementing the noise-corrected sparse identification of nonlinear dynamics (NC-SINDy) algorithm proposed in the paper, including frequency-domain processing, noise correction, and nonlinear sparse regression functions.&#x20;

2. **`Noise41.xls` & `Noise42.xls`** Pre-generated fixed noise datasets for Sections 4.1 and 4.2. Including these files ensures full reproducibility of the results in the paper by eliminating random variations.&#x20;

3. **`Section4.1.nb` / `Section4.1.pdf`** - `Section4.1.nb`: Mathematica notebook for reproducing Section 4.1 simulations. `Section4.1.pdf` : Static PDF export of the notebook, showing all code and outputs for quick reference without running Mathematica.&#x20;

4. **`Section4.2.nb` / `Section4.2.pdf`** - `Section4.2.nb`: Mathematica notebook for reproducing Section 4.2 simulations. `Section4.2.pdf`: Static PDF export of the notebook, showing all code and outputs for Section 4.2 without running Mathematica.

## 📂 How to Reproduce the Results&#x20;

1. Clone or download this repository to your local machine.

2.Move `NCSINDyPackage.m` to `...\Mathematica11\AddOns\Packages` .Move`Noise41.xls` and `Noise42.xls` in the proper directory for successful noise data loading.

3. Open the Mathematica notebook for the section you want to reproduce (e.g., `Section4.1.nb`).&#x20;

4. Run the notebook in Wolfram Mathematica 11. It will automatically load `NCSINDyPackage.m` and the corresponding noise data file.&#x20;

5. All figures, tables, and numerical results from the paper will be generated in the notebook.&#x20;

6. For quick reference, view the corresponding `.pdf` file to check expected outputs.

## 📜 Citation

If you use this code or the NC-SINDy method in your research, please cite our paper:

@article{tang2026ncsindy, title={A noise-corrected data-driven approach for nonlinear dynamical system modeling in frequency domain}, author={Tang, Jian and Xu, Jiahui and Meng, Hao and Zhang, Zhihong and Xu, Huidong and Xiong, Xiaoyan and Wang, Zhihua}, journal={[To be filled in after publication]}, year={2026}, volume={}, pages={}, doi={[To be filled in after publication]} }

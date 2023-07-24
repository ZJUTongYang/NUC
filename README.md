# NUC

This repository contains the supplementary materal (code and video) for the journal paper entitled ``Template-free Non-revisiting Uniform Coverage Path Planning on Curved Surfaces", IEEE/ASME Transactions on Mechatronics (T-Mech). The coverage skeleton generation in C++ has been wrapped to python for easy usage. 

## Video Link

[Supplementary Video](https://drive.google.com/file/d/1sYnp-nKgyRzVhqUaI8ly20HRpIq9SC3B/view?usp=sharing)

## Usage

1. Compile the C++ code into python
```
c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) nuc.cpp -o nuc$(python3-config --extension-suffix)
```

2. Run the demo python file
```
python3 nuc.py
```

3. Run the MATLAB script for visualisation
```
visualisation.m
```

## Cite This Paper
```
@article{Yang2023Template,
  title={Template-Free Nonrevisiting Uniform Coverage Path Planning on Curved Surfaces},
  author={Yang, Tong and Miro, Jaime Valls and Nguyen, Minh and Wang, Yue and Xiong, Rong},
  journal={IEEE/ASME Transactions on Mechatronics},
  year={2023},
  publisher={IEEE}
}
```


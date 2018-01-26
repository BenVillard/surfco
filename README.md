# SurFCo: Surface Mesh From Contours

<p align="center">
  <img src="/Ims/Framework.png" width="50%">
</p>



# About 

SurFCo is a Matlab toolbox with some C++ Scripts to compute a surface mesh from biological delineations, contours or segmentations. It has been developped at the Institute of Biomedical Engineering (IBME), at the University of Oxford, under the supervision of Professor Vicente Grau, and Dr. Ernesto Zacur. 

This work has been published in the following publications, if you use the code we would highly appreciate you citing them:

[1] B. Villard, V. Grau, and E. Zacur, Surface mesh reconstruction from cardiac MRI contours, J. Imaging, vol. 4(1), no. 16, 2018.

[2]  B. Villard, V. Carapella, R. Ariga,  V. Grau, and E. Zacur, Cardiac Mesh Reconstruction from Sparse, Heterogeneous Contours. In: Valdés Hernández M., González-Castro V. (Eds.) Medical Image Understanding and Analysis. MIUA 2017. Communications in Computer and Information Science, Vol. 723. Springer, Cham

# Download and Build

SurFCo uses VTK files. Currently this package only works on windows machines. If you have VTK enabled on MAC, please replace the relevant VTK files (MAC Compiled) inside the "Library/VTK" folder. 

The main file to run is "SurFCo.m". Please refer to [1] for help with parameter selection. Note: SurFCo.m can be run as is. 

An example dataset is provided to ensure SurFCo runs properly. Please run Example_1.m as a test run. 

<p align="center">
  <img src="/Ims/Mesh_2.png" width="40%">
</p>


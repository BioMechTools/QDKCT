# DKCT_Quatitative_Measurement
If you find this code useful in your research, please consider citing:
'A robust and semi-automatic quantitative measurement of patellofemoral instability based on four dimensional computed tomography'

# Table of Content

# Preparation
* Download the codes and data
* Move the folder "Data" to root directory, which is parent directory of "DKCT_Quatitative_Measurement"
* Download the CPD library from "Point Set Registration: Coherent Point Drift" https://arxiv.org/pdf/0905.2635.pdf, or use the version we used in this experiment "CPD_Registration" 

# First glance of the works
You can have a first glance of the final result by running the "Animation/DKCT_Animation.m";

# Workflow
* To redo the experiments, firstly you need to prepare the femur, tibia, patella shape models in stl format; 
* The first step is to determine the anatomical coordinate frame through running "AnatomicalCoordinateFrame\coordinate_main.m"
* The second step is to determine the transform information from the 3D CT to the 4D CT through "CPD_Registration\cpd_morph_static2dynamic_cut_version.m"
* The third step is to determine the transform information among 4D CT scans through "CPD_Registration\cpd_morph_dynamics_cut_version_5_1_11.m"
* The fourth step is to identify the landmarks from the shape models "Quantitative_measurement\cal_landmark_main.m"
* The fifth step is to identify the TT-TG/PC-TG distance through "Quantitative_measurement\cal_TT_TG_PC_TG_dynamic_distance_main.m"

# Others
* To make use of Yao's method to adjust the TT-TG distance by considering the disagreement between the axial CT scan and the anatomical axial scan, you can reference the code "Quantitative_measurement\cal_TT_TG_distance_main.m"

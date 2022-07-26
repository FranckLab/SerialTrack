# SerialTrack MATLAB code package

### SerialTrack: ScalE and Rotation Invariant Augmented Lagrangian Particle Tracking    

This repository contains the Matlab m-files to run our SerialTrack particle tracking algorithm. This code package includes both 2D and 3D particle tracking scripts for full-field 2D and 3D displacement fields, respectively. More details can be found in our SerialTrack paper (https://doi.org/10.48550/arXiv.2203.12573).
  
### Important pages
* [Download latest version v1.0!](https://github.com/FranckLab/SerialTrack/releases)
* [Example data] (https://minds.wisconsin.edu/handle/1793/82901) 
* [Manual](https://github.com/FranckLab/SerialTrack/blob/main/manual_v1.0.pdf) 
* [Questions/Issues](https://github.com/FranckLab/SerialTrack/issues)
* [Bug Fixes/history](https://github.com/FranckLab/SerialTrack/wiki/Bug-Fixes!)
* [Cite](https://github.com/FranckLab/SerialTrack#cite)
* [Franck Lab](https://www.franck.engr.wisc.edu/)
 
## Install and execute SerialTrack

### Code installation

To run SerialTrack, please download MATLAB version later than 2021a and install these code packages/toolboxes:
* Curve fitting toolbox
* Image processing toolbox
* Parallel computing toolbox
* Statistics and Machine Learning Toolbox
* Wavelets analysis toolbox


### Input Image Requirements

* 2D image sequences (at least two frames)
* 3D volumetric images (at least two volumetric stacks)

Images should have segment-able and localizable particles: limited or no particle clusters, particles distinguishable from noise, roughly spherically/cylindrically symmetric, size (roughly) 3 px to 20 px diameter. See examples of particles in the test cases. Custom pre-processing, segmentation, and localization codes can be straightforwardly incorporated (and are encouraged for challenging imaging or specimen conditions).

### Running included example case

1. Make sure that the main files and the supplemental m files (from file exchange) are added to the path in Matlab.
2. Download and save the [example data sets] in the "./imgFolder" folder. 
3. Run the example*.m files in the folder "./Example_main_files"
    1. Edit the paths for SerialTrack and images (as needed)
    2. Check for lines labeled "TODO" and follow instructions (as needed)
   

## Cite
If used please cite:
[1] Jin Yang, Yue Yin, Alexander K. Landauer, Selda Buyuktozturk, Jing Zhang, Luke Summey, Alexander McGhee, Matt K. Fu, John O. Dabiri, Christian Franck. SerialTrack: ScalE and Rotation Invariant Augmented Lagrangian Particle Tracking. (2022) https://doi.org/10.48550/arXiv.2203.12573
 
```
@article{yang2022serialtrack,
  title = {{SerialTrack: ScalE and Rotation Invariant Augmented Lagrangian Particle Tracking}},
  author = {Yang, Jin and Yin, Yue and Landauer, Alexander K. and Buyuktozturk, Selda and Zhang, Jing and Summey, Luke and McGhee, Alexander and Fu, Matt K. and Dabiri, John O. and Franck, Christian},
  doi = {10.48550/ARXIV.2203.12573},
  url = {https://arxiv.org/abs/2203.12573},
  year = {2022}
  publisher = {arXiv}}
```


## Contact and support
For questions, please first refer to [FAQ](https://github.com/FranckLab/SerialTrack#faq) and [Questions/Issues](https://github.com/FranckLab/SerialTrack/issues). Add a new question if similar issue hasn't been reported. The author's contact information can be found at [Franck Lab](https://www.franck.engr.wisc.edu/).

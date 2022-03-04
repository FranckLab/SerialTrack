# SerialTrack MATLAB code package

### SerialTrack: ScalE and Rotation Invariant Augmented Lagrangian Particle Tracking    

 
This repository contains the Matlab m-files to run our SerialTrack particle tracking algorithm. This code package includes both 2D and 3D particle tracking scripts for full-field 2D and 3D displacement fields, respectively. More details can be found in our SerialTrack paper (***Link to be updated***).
  
### Important pages
* [Download latest version v1.0!](https://github.com/FranckLab/SerialTrack/releases)
* [Example data] https://uwmadison.box.com/s/vohfmyoho82zymjiiwkyaqhpfwlv4hmr (***Datasets link to be updated***)
* [FAQ](https://github.com/FranckLab/SerialTrack#faq)
* [Questions/Issues](https://github.com/FranckLab/SerialTrack/issues)
* [Bug Fixes/history](https://github.com/FranckLab/SerialTrack/wiki/Bug-Fixes!)
* [Cite](https://github.com/FranckLab/SerialTrack#cite)
* [Franck Lab](https://www.franck.engr.wisc.edu/)
 
## Install and execute SerialTrack

### Code installation

To run SerialTrack, please download MATLAB (development is on 2021a, older versions may work but are not supported) and install these code packages/toolboxes:
* System Identification Toolbox
* Image Processing Toolbox
* Statistics and Machine Learning Toolbox
* Partial Differential Equation Toolbox
* Wavelet Toolbox
* Curve Fitting Toolbox
* Parallel Computing Toolbox
* MATLAB Parallel Server
* Polyspace Bug Finder

### Input Image Requirements

* 2D image sequences (at least two frames)
* 3D volumetric images (at least two volumetric stacks)

Images should have segment-able and localizable particles: limited or no particle clusters, particles distinguishable from noise, roughly spherically/cylindrically symmetric, size (roughly) 3 px to 20 px diameter. See examples of particles in the test cases. Custom pre-processing, segmentation, and localization codes can be straightforwardly incorporated (and are encouraged for challenging imaging or specimen conditions).

### Running included example case

1. Make sure that the main files and the supplemental m files (from file exchange) are added to the path in Matlab.
2. Download and save the [example data sets](***Datasets link to be updated***) in the "./imgFolder" folder. 
3. Run the example*.m files in the folder "./Example_main_files"
   

## Cite
If used please cite:
[](*** Paper link to be updated ***)

```bibtex
@article{ 
}
}
```

## Contact and support
For questions, please first refer to [FAQ](https://github.com/FranckLab/SerialTrack#faq) and [Questions/Issues](https://github.com/FranckLab/SerialTrack/issues). Add a new question if similar issue hasn't been reported. The author's contact information can be found at [Franck Lab](https://www.franck.engr.wisc.edu/).

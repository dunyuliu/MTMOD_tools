name: mtmod
channels:
  - conda-forge
  - anaconda  
  - defaults
dependencies:
dependencies:
 - python>=3
 - matplotlib
 - numpy
 - scipy
 - xarray
 - gmt
 - pygmt
 - pip:
    - gdal==3.5.0
    - okada-wrapper==18.12.7.3
    - pyqt5-sip==12.9.0
    - tectonic-utils==0.0.8
    - cutde

# Lazy Selector on ITK Showcase
Lazy Selector / visualisation of connected components. And itk showcase.

Dependencies
------------

* ITK 4.6.1
* VTK 6.2.0

Building
--------

Built with cmake.

You need a compiler that fully supports `C++11`. It has been tested with `g++-4.7.4`.

ITK has to be built with Review ON and vtkGlue.

    mkdir build
    cd build
    cmake ..
    make

Running
-------
    
Run with

    ./lazySelecThor pathToVolume [binary threshold]
    
License
-------
GPL v3 see LICENSE
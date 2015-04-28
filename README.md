# Graph Cuts

Dependencies
------------

* ITK 4.6.1
* VTK 6.2.0

The code is header only but the code main has an extra dependency:

* [tclap](http://tclap.sourceforge.net/)



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

    ./graphCutter -s ~/code/data/merged.tif -i ~/code/data/image.tif -g ~/code/data/seedSink.tif -v
    
You can find sample data at [cvlab](https://documents.epfl.ch/groups/c/cv/cvlab-unit/public/espina/sample_data/graphcuts/)
    
    
License
-------
GPL v3 see LICENSE

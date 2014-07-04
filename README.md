clFFT-D
=======

D bindings to clFFT (https://github.com/clMathLibraries/clFFT), an FFT library for opencl.

Usage
-----

```import clFFT;``` and the library is ready to go. However, unlike the clFFT C header, you must import your own choice of opencl headers to access the openCL API which is needed for using clFFT.  These bindings internally use DerelictCL (https://github.com/MeinMein/DerelictCL and http://code.dlang.org/packages/derelict_extras-opencl) but do not expose it publicly.

To build, simply add the relevant dependancy to your project's ```dub.json``` and let dub take care of it and DerelictCL.
If you are not using dub then you're on your own, but something like this should work to create a static library for you to link your project to:
```dmd clFFT.d -lib -I$DERELICTCL_SOURCE_DIR -L-L$PATH_TO_DERELICTCL_LIBS -L-lDerelictCL```

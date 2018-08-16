# Biotissue
multiscale biological tissue implementation using amsi

## Installation
1) Install [SCOREC/core](https://github.com/SCOREC/core)
   Make sure you compile with Simmetrix support. Follow directions in the [Scorec Core wiki](https://github.com/SCOREC/core/wiki/General-Build-instructions).

2) Install [tobinw/las](https://github.com/tobinw/las)

3) Install [tobinw/amsi](https://github.com/tobinw/amsi)

4) Install [yaml-cpp-0.3.0](https://github.com/jbeder/yaml-cpp/releases/tag/release-0.3.0)
   This package is used for some tests of the microscale, and setting up microscale only jobs. We use an old release of this package to enable building on the BGQ which has old versions of Cmake and does not have full C++11 support.
    ```bash
     git clone https://github.com/jbeder/yaml-cpp.git
     cd yaml-cpp
     git checkout tags/release-0.3.0
     mkdir build
     cd build
     cmake .. -DCMAKE_INSTALL_PREFIX=$DEVROOT/install/yaml-cpp/0.3.0
     make -j 16 install
     mkdir -p yaml-cpp-module $DEVROOT/module/yaml-cpp/
     cp yaml-cpp-module $DEVROOT/module/yaml-cpp/0.3.0
    ```

5) Build and Install Biotissue
   ````bash
   cd scripts
   ./config.sh
   cd ../build_debug
   make -j 16 install
   ```

6) [Optional] Run tests
  ```bash
  cd build_debug
  ctest
  ```

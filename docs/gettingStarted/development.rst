Setting up a development environment
===============================================

#. Download the scorec-mumfim repository.

   .. code-block:: bash

      # create parent directory for scorec-mumfim
      mkdir mumfim-develop && cd mumfim-develop

      # download scorec-mumfim
      https://github.com/SCOREC/mumfim.git
#. Create a spack environment in development folder.

   .. code-block:: bash

      cd mumfim-develop

      # add spack to the command line
      source $SPACK_ROOT/spack/share/spack/setup-env.sh

      # create environment
      spack env create -d 

      # activate environment
      spack env activate /path/to/mumfim-develop

Visual Studio Code (vscode)
```````````````````````````````````````````````
#. Setup vscode for C++. Follow directions given `here <https://code.visualstudio.com/docs/languages/cpp>`_.
#. Get CMake extension.
#. Create a .vscode folder in directory (.git level). Create c_cpp_properties.json and settings.json files.

   .. code-block:: yaml
      :caption: c_cpp_properties.json

      {
         "configurations": [
            {
               "name": "CMake",
               "compileCommands": "${config:cmake.buildDirectory}/compile_commands.json",
               "configurationProvider": "ms-vscode.cmake-tools"
            }
         ],
         "version": 4
      }

   .. code-block:: yaml
      :caption: settings.json

      {
         "cmake.configureArgs": [
            "-DBUILD_EXTERNAL=OFF",
            "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"
         ],
         "cmake.configureEnvironment": {
         "CMAKE_PREFIX_PATH":"/home/username/mumfim-develop/.spack-env/view"
         },
         "cmake.cmakePath": "/home/username/spack/opt/spack/linux-ubuntu20.04-skylake/gcc-9.4.0/cmake-3.23.1-76ep73jjuzdtvkgpcsgu2omgsddwqtq6/bin/cmake"
      }
#. Open file in vscode.

   .. code-block:: bash

      # directory with .vscode folder
      cd path/to/mumfim

      # replace . with any file name
      # include file extension
      code .


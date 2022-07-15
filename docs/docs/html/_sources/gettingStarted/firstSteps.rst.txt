First steps 
===============================================
This guide will cover basic installation of scorec-mumfim. 

Windows
```````````````````````````````````````````````
#. Install windows subsystem for linux (wsl). Follow directions given `here <https://docs.microsoft.com/en-us/windows/wsl/install>`_. 
#. Begin initial setup. Open the command line. Enter single commands as separated by comments. 

   .. code-block:: bash

      # create parent directory for all spack things
      mkdir spack && cd spack
      export SPACK_ROOT=`pwd`

      # download spack
      git clone https://github.com/spack/spack.git

      # make sure to clear any modules which you don't want to interfere with spack
      module purge

      # add spack to the command line
      source $SPACK_ROOT/spack/share/spack/setup-env.sh

      # make directory to store all external spack repos
      mkdir repos && cd repos

      # download the mumfim spack repo
      git clone git@github.com:jacobmerson/mumfim-spack.git

      # add repositories to spack
      spack repo add mumfim-spack/mumfim  
#. Install scorec-mumfim.

   .. code-block:: bash

      # install mumfim
      spack install mumfim
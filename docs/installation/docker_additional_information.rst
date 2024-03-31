.. _docker_additional_information:

CellOracle docker image quick start and additional notes
========================================================


In this documentation, We provide example commands for our celloracle docker image below. But we are not aiming to provide detailed information about general docker usage.
We highly recommend learning about docker if you are not familiar with it, and make sure you have adequate knowledge about docker prior to start celloracle analysis with docker.

Quick start
^^^^^^^^^^^
1. Download celloracle docker image from docker Hub.

::

    docker pull kenjikamimoto126/celloracle_ubuntu:0.18.0


2. Make docker container and start running it.

- As we recommend using jupyter notebook (or jupyter lab), we need to set up port. Here, we connect docker container port 8888 with local host port 8888.

- To access files smoothly, we can use `bind mount <https://docs.docker.com/storage/bind-mounts/>`_. Here, we will mount the `data_folder`, a directory in your local machine.


::

    mkdir data_folder # Create data folder in your local environment.
    docker run -dit \
      --name celloracle_container \
      -p 8888:8888 \
      -v $(pwd)/data_folder:/root/data_folder \
      kenjikamimoto126/celloracle_ubuntu:0.18.0



3. Enter the docker container.

::

    docker container exec -it celloracle_container bash



4. In the docker container environment, start jupyter notebook as follows.

::

    cd
    jupyter notebook --port=8888 --ip=0.0.0.0 --allow-root --no-browser


After starting jupyter, please open your browser and enter http://localhost:8888 to access jupyter notebook running in the docker container.
You need to enter a token to access jupyter. The token can be found in your terminal running jupyter notebook.

.. warning::
   - We found that the CellOracle calculations may be EXTREMELY SLOW in a Windows Subsystem for Linux (WSL). Please don't use docker on Windows OS.
   - You may have issue if enough computational resource is not assigned to your docker container. Please make sure you have enough memory.
   - If you found bug or error in the docker installation, please let us know through `GitHub issue page <https://github.com/morris-lab/CellOracle/issues>`_.



Docker image build information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We built our docker image using Dockerfile automatic build function.
The Dockerfile is available `here <https://github.com/morris-lab/CellOracle/blob/master/other_files/Dockerfile>`_.
You can modify it to create custom docker image by yourself.
If you make custom environment, please do so on your responsibility.

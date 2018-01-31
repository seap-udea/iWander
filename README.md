# iWander
### Dynamics of Interstellar Wanderers

[![arXiv](http://img.shields.io/badge/arXiv-1711.09397-orange.svg?style=flat)](http://arxiv.org/abs/1711.09397)

Presentation
------------

This package is an implementation of the general method devised by
Zuluaga et al. (2017) for assesing the origin of interstellar small
bodies (asteroids and comets).

The package include data and tools that can be used in general for
studying the dynamics of an interstellar vagabond object (small-body,
interstellar spaceship and even stars).

If you use or adapt the package for any scientific purpose please
refer to this paper:

> Zuluaga,Sanchez-Hernandez, Sucerquia & Ignacio Ferrin, A general
> method for assesing the procedence of an interstellar small body:
> the case of 1I/´Oumuamuaa (1I/2017 U1),
> [arXiv:1711.09397](https://arxiv.org/abs/1711.09397).

The paper is under review for publication in the Astronomical Journal.

> **NOTE** As the revision process is undergoing, un updated version
  of the paper will be available in [this
  link](doc/Zuluaga_et_al_2017-AssesingOriginProbability.pdf).

Progenitor Candidates
----------------------

A list of the progenitor candidates of the interstellar object
**1I/2017 U1 ('Oumuamua)**, identified using this method along with
their corresponding origin probabilities are available (and will be
updated) [in this file](CANDIDATES-Oumuamua-Past.md).

> **NOTE**: This list could change when better astrometric information
  or improvements on the methodology be available. Stay tuned!

Getting the package
-------------------

You may obtain the package in different ways:

- Cloning anonymously the github repository:

  ```  
  git clone http://github.com/seap-udea/iWander.git
  ```  

- Downloading the tarball of the latest release:

  ```  
  wget http://github.com/seap-udea/iWander/archive/master.zip
  ```  

- Cloning the package as a developer (credentials required):

  ```  
  git clone git@github.com:seap-udea/iWander.git
  ```  

The size of the package is large (several hundreds of MBs).  This is
mainly due to the data required to run some of the modules.  

Installation
----------------

For installing the package run:

    bash install.sh

The installation script performs several operations required to start
using the package:

1. Copy the template files in the ``util/conf`` directory to the
   root package directory.

2. Install the required linux and python packages.

3. Unpack large files.  Large files are splitted in 20MB chunks of
   data inside the `.store` directory.  

4. Attempt to compile a small test program.  If compilation fails
   comment/uncomment the proper line in the file ``compiler.in``:

   ```  
   ###################################################
   #CHOOSE YOUR ARCHITECTURE
   ###################################################
   #ARCH=32
   ARCH=64
   ```  

Quickstart
----------

1. Prepare the configuration file for the object and place it in the
   ``objects`` directoy. See the already available examples there
   (Oumuamua, Voyager 1 and Voyager 2).

1. Edit the general configuration files ``iwander.conf`` and
   ``wanderer.conf`` for setting up the name of the wanderer and its
   properties.

2. Compile the key programs:

   ```  
   make all
   ```  

3. Generate the surrogate objects and propagate them until the time of ingress.

   ```  
   bash bin/run.sh wanderer.exe
   ```  

4. Compute the minimum distance to all the stars in the input catalog
   and select the candidates.

   ```  
   bash bin/run.sh encounters.exe
   ```  

5. Find progenitor candidates and compute their origin probability:

   ```  
   bash bin/run.sh probability.exe
   ```  

The output of this process is the file ``progenitors-<wanderer>.csv``
having a list of the progenitor candidates with their respective
origin probability.

Optionally you can:

6. Generate the table of ingress properties:

   ```  
   python3 bin/ingress.py
   ```  

7. Sort out the candidates according to position probability or
   minimum distance:

   ```  
   python3 bin/progenitors.py
   ```  

Structure of the package
------------------------

The package is made of three type of components: programs, scripts and
databases.

Programs are used to compute the core functions (propagate wanderers,
find encounters, compute interstellar origin probabilities, etc.)

Scripts are used for pre and post processing of the information
required or produced by the package.

Databases contain the information required to run some of the
functionalities of the package.

Components
----------

- **wanderer**: This program integrate the orbit of a moving object
  inside the Solar System.

  * Function: 

    This program perform three different tasks:

    1) Calculate the time t_asy when the single conic approximation is
       good enough to predict the future position of the interstellar
       object.

    2) Calculate the time t_ing at which the object was at a distance
       equivalent to the truncation tidal radius of the Solar System.

    3) Predict the position and velocity of the surrogate objects at
       t_ing.

  * Input: None

  * Output: 

    * wanderer.csv: properties of all the surrogate objects.

    * ingress.dat: a summary of the ingress orbit properties including
      the epoch of asymptotic elements and their covariance matrix,
      the time of ingress, the radiant and velocity at ingress.

- **encounters**: This program integrate the orbit of a moving object
  inside the Solar System.
  
  * Function: 

    This program perform two different tasks:

    1) Compute the LMA minimum distance and time to all stars in the
       AstroRV catalogue...

    2) Select the progenitor candidates.

  * Input: 
    - wanderer.csv

    Output: 

    - encounters-<Wanderer>.csv: all the columns of the input catalog (AstroRV)
      plus additional information computed from the LMA approximation.

    - candidates-<Wanderer>.csv: list of objects fulfilling certain
      selection criteria that classify them as close encounters candidates.

- **probability**: This program integrate the orbit of a moving object
  inside the Solar System.

  * Function: calculate the IOP for a list of stellar candidates.

  * Input:
    - wanderer-<object>.csv
    - candidates-<object>.csv

  * Output: 
    - progenitors-<object>.csv

For the developer
-----------------

iWander uses GSL and Spice as backbone utility libraries.  The latest
precompiled version of both libraries, along witth the header files
are provided with the package in the `util` directory.

Input parameters are passed to the programs using a configuration file
`<programa.conf>`.  The configuration file has the structure of a C
program.  The declarations and actions in the program are included
directly into the `main` of the corresponding program.

Naming conventions:

* Configuration variables: Capitalized. Example: Wanderer,
  Npart.

* Macros and global variables: Fully capital. Example: FILENAME,
  REARTH.

* Routines: Umbrella style. Example: vectorAllocate, integrateEOM.

* Local variables: Free naming rules.

Acknowledgements
----------------

This package has been developed thanks to the incredible work made by
previous scientist and developers. Most of the work of those who make
this package possible has been cited in our papers.  Others are
mentioned in the software itself.

License
--------------
Copyright (C) 2017 Jorge I. Zuluaga, Oscar Sanchez-Hernandez, Mario Sucerquia & Ignacio Ferrin

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and the databases associated (the "Package"),
to deal in the Package without restriction, including without
limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Package, and to
permit persons to whom the Package is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Package.  A
reference to the Package shall be included in all scientific
publications that make use of the Package.

THE PACKAGE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE DATABASE OR THE USE OR OTHER DEALINGS IN THE DATABASE.
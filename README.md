# iWander
Dynamics of Interstellar Wanderers

Presentation
------------

This package is an implementation of the general method devised by
Zuluaga et al. (2017) for assesing the origin of interstellar small
bodies (asteroids and comets).

The package include a series of databases and tools that can be used
in general for studying the dynamics of an interstellar vagabond
object (small-body, interstellar spaceship and even stars).

If you use or adapt the package for any scientific purpose please
refer to [this paper](http://arxiv.org):

      Zuluaga,Sánchez-Hernández, Sucerquia & Ignacio Ferrin, A general
      method for assesing the procedence of an interstellar small
      body: the case of 1I/´Oumuamuaa (1I/2017 U1), arXiv:1120456.

Getting the package
-------------------

You may obtain the package in different ways:

- Cloning anonymously the github repository:
```  
git clone http://github.com/seap-udea/iWander.git
```  

- Download the tarball:
```  
wget http://github.com/seap-udea/iWander/archive/master.zip
```  

- Cloning the package as developer (permissions required):
```  
git clone git@github.com:seap-udea/iWander.git
```  

The size of the package is large (several hundreds of MBs).  This is
mainly due to the data required to run some of the modules (GAIA and
Radial velocities databases, SPICE Kernels, etc.).  

Before you start
----------------

Before you start using the package set up the `makefile`.  Create a
local copy of `compiler.in.temp`:

```  
cp compiler.in.temp compiler.in
```  

Edit the file to properly choose your system architecture.  You just
need to comment/uncomment the proper line in the file:

```  
###################################################
#CHOOSE YOUR ARCHITECTURE
###################################################
#ARCH=32
ARCH=64
```  

Test it:

```  
make
```  

Unpacking large files
---------------------

Large files are splitted in 20MB chunks of data inside the `.store`
directory.  Before you use some critical modules you must unpack those
large files:

```  
make unpack
```  

Quickstart
----------

The travel of an 'interstellar wanderer' starts in the Solar System.
The first step to follow the path of the body is to propagate their
initial conditions into the far realms of our planetary system.

This is achieved using the program `wanderer.exe`.  In order to run
the program execute:

```  
make wanderer.exe
```  

If you want to modify the options of the program edit the file
`wanderer.conf`.

This program start with the initial elements of the body (which are
fully specified in the configuration file) and propagate the orbit
into the distant future or past (this is specified with the variable
`duration`).

Since the initial position of the body is always uncertain, the
program allows you to integrate the path of many test particles with
orbital elements compatible to those of your body.

Once the integration is fully completed the results are stored in the
file `wanderer.csv`.

For a fully explanation about the input options and output information see
components section.

Once you have propagated your wanderer out of the Solar System is time
to check if it found or will find any star in their path.  This
achieved by running the program:

```  
make encounters.exe
```  

This will create a list of stars in the huge database included with
the package, to which the wanderer will flyby to a maximum distance
defined by the user in a desired time window.

Structure of the package
------------------------

The package is made of three components: programs, scripts and
databases.  

Programs are the components used to compute the core functions
(propagate wanderers, find encounters, compute procedence or capture
probabilities, etc.)

Scripts are used for pre and post processing of the information
required or produced by the package.  Here we mainly use python and
Ipython scripts.

Databases contain the information required to run some of the
functionalities of the package.

Components
----------

- **wanderer**: This program integrate the orbit of a moving object
  inside the Solar System.
  
  * Usage: 

    ``
    make wanderer.exe
    ``

  * Input information: wanderer.conf

  * Output: 

    + File: wanderer.csv.  Contains the solution of the
      backward/forward integration of the wanderer and/or the test
      particles created from their errors.

     ```
	  0:i
	  1-6:Initial elements, q,e,i,W,w,Mo,to,mu
	  7:tdb (terminal)
	  8-13:Position Ecliptic J2000
	  14-19:Position J2000
	  20-25:Position Galactic J2000
	  26:RA(h) (terminal)
	  27:DEC(deg)
	  28:l(deg)
	  29:b(deg)
	  30:d(AU)
	  31-36:Asymptotic elements, q,e,i,W,w,Mo,to,mu
    ```


For the developer
-----------------

iWander uses GSL and Spice as backbone utility libraries.  The latest
precompiled version of both libraries, along witth the header files
are provided with the package in the `util` directory.

Acknowledgements
----------------

This package has been developed thanks to the incredible work made by
previous scientist and developers. Most of the work of those who make
this package possible has been cited in our papers.  Others are
mentioned in the software itself.
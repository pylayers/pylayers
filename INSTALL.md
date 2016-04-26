# PyLayers Installation

* Installation notes for [Linux / OS X](https://github.com/pylayers/pylayers/blob/master/INSTALL.md#linux--os-x) and [Windows](https://github.com/pylayers/pylayers/blob/master/INSTALL.md#windows-install)
* [Source/Project Directories Set-up](https://github.com/pylayers/pylayers/blob/master/INSTALL.md#sourceproject-directories-set-up)
* [Usual Install Issues](https://github.com/pylayers/pylayers/blob/master/INSTALL.md#usual-install-issues)
* [Dependencies](https://github.com/pylayers/pylayers/blob/master/INSTALL.md#dependencies)


## Linux / OS X

The prefered way to install PyLayers is to first install the free python distribution [Anaconda](https://store.continuum.io/cshop/anaconda/) on your platform. 
By doing so, you're installing most of the required dependencies.

1. Download [Anaconda](https://store.continuum.io/cshop/anaconda/)
2. Install Anaconda with the given command ( $bash Anaconda-<version>-<plateform>.sh ). 
Important : Say yes to add Anaconda to your path
3. Clone PyLayers : **git clone https://github.com/pylayers/pylayers.git**
4. Run "./installer_unix" file from the PyLayers directory (YOUR FIRST NEED to add exe rights using the command: chmod +x ./installer_unix )
5. Done

## Windows Install

1. Download and Install [Anaconda](https://store.continuum.io/cshop/anaconda/) 
2. Download the Shapely package corresponding to your platform from 
[Unofficial Windows Binaries for Python Extension Packages #Shapely](http://www.lfd.uci.edu/~gohlke/pythonlibs/#shapely) *
3. Download the Basemap package corresponding to your platform from [Unofficial Windows Binaries for Python Extension Packages #Basemap](http://www.lfd.uci.edu/~gohlke/pythonlibs/#basemap)
4. Install Shapely and Basemap using the downloaded files
5. 3. Clone PyLayers : **git clone https://github.com/pylayers/pylayers.git**
6. Run installer.bat
7. Done


* Despite shapely is part of the Anaconda distribution, the Windows distribution doesn't contains the libgeos dependencies.


# Source/project Directories Set-up

The source and Project directory path are automatically filled into a ".pylayers" file during the install. The .pylayers file is stored at the root of user's home (~/).

* The source directory, gathers all the source files.
* The project directory , gathers configuration files.

The default ".pylayers" contains :
    
    ```
    source
    /the/path/where/installer/has/been/run/pylayers
    project
    ~/pylayers_project
    ```

The source path indicates where the PyLayers' source files are located. You generally don't have to modify this path

The project path indicates where the configurations files are located.
If you decide to move your project directory elsewhere, you must change the project path in consequence.


**Notes for previous users of PyLayers:**

The environment variables $PYLAYERS and $BASENAME are no longer required.
If you try to re-install PyLayers on an existing version:
-the source path of the .pylayers file will be set at the value contained into the $PYLAYERS environment variable
-the project path of the .pylayers file will be set at the value contained into the $BASENAME environment variable




# Usual Install Issues

Most of the issue encountered during the install is due to 
1. Multiple python installations and/or  
2. Multiple installation of some packages ( especially shapely/basemap under windows)

## 1. Multiple python installations

Often, multiple python installation co-exists on a same computer, and this can cause complex PyLayers installation issues. To avoid that kind of issue, please,, make sure that you're using the python installation provided by Anaconda ,and not another one.
To this end, you can check this on unix platform with :
$ which python

and under windows by opening 'cmd' and using :
c:\ for %i in (python.exe) do @echo.   %~$PATH:i

If the answer doesn't point to the Anaconda installation directory, you're using the wrong python version.
You have to update the default python installation. I cannot describe here the detail procedure because it 
is really depend to your configuration.

## 2. Multiple installation of some packages

Another source of issues comes from some package conflicts.
This is especially the case by using a WINDOWS OS, when the shapely installation has to be perfomed by using an external installer and not the one provided by Anaconda (because it does not contains required dependencies). 

If after the installation, the import of a PyLayers module cause an Error involving shapely and/or basemap, you're probably facing the package conflict issue. 

The easiest way to solve it is to remove the '.egg' of the incriminated package in :

$ <Path to Anaconda>/lib/python2.7/site-packages/   (UNIX)

or 
c:\<Path to Anaconda>\Lib\site-packages (Windows)




# Dependencies 
## heading<a name="headin"></a>

The following modules are required to use PyLayers 

numpy>=1.6.1  
scipy>=0.10.1  
networkx>=1.7  
matplotlib>=1.1.0  
shapely>=1.2.14  
descartes>=1.0  
SimPy>=2.2  
PIL>=1.1.5  
bitstring>=3.0.2  
pyintervall>=1.0  
osmapi>=0.3  
imposm  
basemap>=1.0  
pandas >=0.12.0  
mayavi >=4.3.1  


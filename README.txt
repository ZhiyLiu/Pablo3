Pablo2 

Using CMAKE

Product of MIDAG 2012

CMAKE use with help of CREATIS 2012 (c) 

lib will contains as many directories as you need
---

lib/library name is just ... a sample (example) library

appli/bin  contains the executables and tests

doc/to see documentation of pablo

data/ to see some data for the test for example


Compilation of Pablo in Windows Systems

Download and install ITK and VTK

Use CMAKE to configure the project, the libraries LAPACK, FLTK will be find automatically, they are included in the project. in the folder include_ext.

When the project is configured using CMAKE, you can compile it.



Compilation of Pablo in Linux Systems


External libraries needed:

LAPACK
FLTK

You can either download and install this libraries or install the distribution packages.

ITK
VTK

It is recommended to Download the sources, compile them and install them.


To compile pablo:

mkdir PabloBin
cd PabloBin
ccmake <path to pablo>

press 'c' configure until the generate 'g' appears.
	During the configuration step, cmake might ask to point to the path to the libraries FLTK, ITK and VTK
	If you have compiled and install FLTK, is recommended to manually set the path for each variable shown in cmake. 
	It is highly recommended to use the distribution version as this will pose no problem when configuring Pablo.

make

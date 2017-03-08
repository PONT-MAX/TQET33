Model-CS Toolbox v1.1 . 
Release date: Feb 20, 2012

—————————————
Changes in the latest version
—————————————

Fixes small bugs all around, particularly in the 2D-wavelet tree model. 


————————
Original README.
————————

------------------------------------------------------------------
Software for Model-based compressive sensing
------------------------------------------------------------------

Original release date : August 7, 2009
Reference       : "Model-based compressive sensing", submitted to IEEE Trans. Info. Theory
Authors           : Richard G. Baraniuk, Volkan Cevher, Marco F. Duarte and Chinmay Hegde
Download        : http://arxiv.org/abs/0808.3572

Questions/suggestions/comments: chinmay@rice.edu


--------
Notes:
---------
* This package contains Matlab implementations of various CS recovery algorithms that have been developed using the framework described in the above paper. 

* Each folder in the package essentially consists of a CS recovery algorithm based on a particular signal model, and a script that tests that recovery algorithm. The names of the scripts typically end with '_example.m'

* Most of the code is plain Matlab code
Each m-file starts with a header that contains a brief description of its I/O interface, details of extra software needed for its proper functioning, known problems with code usage and any additional references to newer, published papers. 
Thus, it is IMPORTANT that you carefully read the header of any function which you plan to use!

* Additional packages needed for some functions (and which are always useful) :
    - The Rice Wavelet Toolbox (http://dsp.rice.edu/software/rice-wavelet-toolbox) 
    - l1magic (http://www.acm.caltech.edu/l1magic/)
    - CVX (http://www.stanford.edu/~boyd/cvx/)
  
-------------------------------
Pest control
-------------------------------
 Most of the scripts, as they stand, have been tested on Windoze, 32-bit UNIX and Mac OS X platforms. 
 If you run into a bug and if Matlab complains, do one of the following:
   - make sure you've installed the necessary packages as listed above.
   - ensure said packages are included in your path.
   - look at the header of the buggy file and check for I/O errors. There are virtually zero checkpoints for wrong argument types, so 
     make sure you've got your matrix sizes etc. right when you pass a set of arguments to a function. 
  -  look at the Matlab source around the line at which the code crashed and look for comments/warnings/pitfalls that need to addressed.
  
   If all fails, send me an email and I will try and help you out!
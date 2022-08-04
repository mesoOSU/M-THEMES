# Implimentation of NFFT's NuFFT into M-THEMES

### *Arvind Rajagopalan, Summer 2022*



## Synopsis
---
the main branch of M-THEMES (Microstructure THErmo-MEchanical Simulator) uses fftw2 for its internal computations. We would like to replace this with NFFT, so we can utilize the non-uniform fast fourier transform (Nu-FFT) instead. This internship will be considered a success if we can either:
* 1. Compile M-THEMES using NFFT instead of fftw2 and replace all the fft calls with nu-fft calls, allowing for non-uniform calculations

* 2. Compile M-THEMES using NFFT instead of fftw2, but just use NFFT's uniform fft, which is mostly identical to fftw2 but allows for compiling with modern MPI libraries and easier implimentation of Nu-fft in the future

OR

* 3. a clear and concise documentation of what you tried, so the next person who works on this has a road map of what to do and not do ([here is a good example written by a previous Hish School intern](https://github.com/mesoOSU/M-THEMES/blob/hdf5/writeup.md) )


## Quick Background
---
### Basic Material Science for Metals
Almost all metals are, on a microscopic scale, a composite of individual single crystals. Think of a formica countertop, except the individual grains are between 1-100 microns (for reference, visible light is between 0.1 and 1 microns, and a human hair is 60).

If you have a single crystal, the stress-strain relationship (ie, how a material bends given a force) is very Anisotropic (ie, depends on a direction). This is like how bamboo is easy to slice lenghwise, but not accross. Most materials are polycrystalline, but if a bunch of your crystals are pointing in the same direction, that anisotropy starts showing through. easy example of this is folding metal for swords, which makes them able to bend instead of break (folding aligns the grians, don't worry about how though)

### M-THEMES
Now, modeling this for any process is hard, but modeling it for industrial processes with multi-axis twists and turns is REALLY hard, especially when modeling tiny things with only a few thousand grains (this is often called the Meso lengh scale). even if we know what happens, its a lot of math. It turns out (again, don't worry about why) the math is exponentially easier if solved in fourier space instead of euclidian space. this is a problem that comes up a lot in other fields, so Comp Sci people in the 90s got really good at writing efficient fft solvers, with fftw being one of the most famous.

M-THEMES is a Mesoscale microanalysis solver that uses a 3d fast fourier transform (3d-fft) to solve stress and strain equation. in simpler terms, it takes the stress strain math, converts it to a fourier series, passes that to an fft solver, then interprets the results.

### fft solvers

The problem right now is M-THEMES uses fftw2 for the fourier math, which has been superceeded by fftw3 since 1999, unsupported since mid-2010, and doesn't allow for non-uniform calculations. The lack of support forces us to use other outdated libraries and limits our ability to incorporate new capabilities written by other researchers, whereas the lack of a non-uniform solver puts a limit on the accuracy of our simulations.

Both of these problems can be solved by using NFFT instead of FFTW, which is an alternative fft library which wraps fftw3 and adds new capabilities. Even if we ONLY get M-THEMES working with NFFT, what eliminates several current problems other users are having. If we can switch out the U-fft for a NU-fft, that is even better, and will massively improve our computing efficiency (this software is often ran on OSU supercomputers)

This was written by Austin Gerlt, and glosses over a LOT of details, so feel free to contact me for any details, as this is a very quick-and-dirty explination of the full problem.

## Current To-Do List
---
[x] - Create To-Do list

[x] - Fork the repo

[x] - as a test and as practice, check this box, commit the change, push it to your repo, them make a pull request to push your change to the MESOOSU repo (Austin can walk you through all this over zoom if you want) 

[] - run "non_uniform_test_data.py", which should read in "spinodal.png" and create two text files from it: one with regular data and one with irregular data. the format is "x-coordinate, y-coordinate, data value" for both. 

[] - figure out how to read in the data into nfft and have it return the 2D fourier transform. you can do this however you want: python, java, directly in command line, your choice. also figure out how to write the output to a text file. again, your choice on how.

[] - read in the regular grid data into the uniform solver, and both the regular and irregular into the non-uniform, and save the results. since the results are going to be complex numbers, find the magnitude of the results, then plot them side by side so we can see how different they are (this should be 3 2-by-2 plots, each with several sharp peaks throughout.)


## Future To-Do List (after we can show that NUFFT works and withminimal slowdown)
---

[] -successfully compile M-THEMES as is. THE BEST WAY TO DO THIS IS TO CHANGE THE INSTALL.SH FILE. Edit it so instead of installing directories locally, it just works on your machine specifically (might need to still make local installs since it uses such outdated mpi and fftw libraries, but the install.sh deletes, recompiles, and relinks those everytime, which takes hours and is pointless.)

[] - once compiled as-is, attempt replacing fftw2 with nfft

[] - once working with nfft, try replacing nfft's fft solver with the nu-fft solver


## Meeting Notes
---

7/16/22:
- Successfully compiled nfft and fftw
- forked a copy of M-THEMES
- discussed in broad terms the Background and Objectives defined in this document
- Need to schedule regular meetings, but also let Steve and Austin know whenever you get M-THEMES compiled as-is 

8/4/22:
- NFFT is working, git repo is functioning
- need demo of I/O process for current NFFT (ie, how do I input data, and what do i get out)
- Austin gave Arvind some python code to generate regular and irregular data from the same dataset.
- next goal is proof of concept and validation for super basic NFFT example before going onto replacing FFTW2 with NFFT
- Need to schedule regular meetings, but also let Steve and Austin know whenever you get M-THEMES compiled as-is 

## Resources
---
Arvind, Add whatever links you want here, especially links to guides on how to compile NFFT and the examples you ran. this will be accessed by lots of people down the road who will also use this software.



-  [Link to nfft manual](https://www-user.tu-chemnitz.de/~potts/nfft/guide/nfft3.pdf)


COMPASS software package
========================

**COMPASS: Convex Optimization for Microbial Profiling by Aggregating Short Sequence reads**

The COMPASS repository contains programs for reconstruction of microbial species identities and frequencies from massively parallel sequencing short-reads data.

The repository includes source code, executables and example datasets, and can be used for microbial profiling of experimental and simulated reads, for simulating reads, and for evaluating profiling results.

Installation: 
-------------

Download the COMPASS directory from https://github.com/NoamShental/COMPASS
Use the "Download ZIP" button on the right hand side of the page. It's about 300M, which are mostly database. Once done, uncompress it, and follow the introductory files listed below.

You will have to change the path variable to your COMPASS location. For example, if COMPASS was saved under your home directory then enter the following in each introductory file:
userDir = '~/COMPASS/'; 

Requirements:
-------------
The package is written in Matlab and uses mex files which are currently only compiled for Linux 64bit. The software requires Matlab's Parallel Computing Toolbox. 
If a Mac/Windows version is needed, please contact shental@openu.ac.il


How to start? Introductory files

We included 3 introductory files. These are located in the directory "mFiles". Please follow these examples to better understand how to apply COMPASS.

1.	Example of a single simulation: An example of running a specific simulation - from creating the reads to applying COMPASS. See file:  "example_of_a_single_simulation.m". We recommend following this example first, section by section.

2.	Example of simulations similar to those performed to produce Figure 3 in the manuscript: namely varying the number of reads, read length and number of bacteria – See file: "example_simulations.m"

3.	Example of analysing the Drosophila larva sample L2 using COMPASS – See file: " example_experimental_reads_sample_L2.m "  



Directories in the COMPASS package:
-----------------------------------

**mFiles** – Matlab files and mex files.

**bin** - Compiled executibles (t.b.d.)

**results** – A directory used to save simulation results for the examples. Note that we have already saved several example results that appear in the introductory examples.

**database** – the Greengenes database used in a COMPASS format. Two databases are provided – the full 16S rRNA gene database and the database of 750bp long sequences covering variable regions V3-V6 (see manuscript)

**experimentalReads** –  example files for larva sample L2. Other experimental data presented in the paper are available at the MG-RAST website: http://metagenomics.anl.gov/linkin.cgi?project=5237

**cvx** – the optimization software used, downloaded from the CVX website (http://cvxr.com/cvx/)

**mothur** – the mothur software (http://www.mothur.org/) used to calculated distances between correct and reconstructed bacteria.




Other issues
------------

**Types of inputs to COMPASS**

COMPASS receives two inputs – a database and the set of reads. This package already contains the Greengenes database used in the paper. As for the reads – these can originate from two different sources – either simulated reads or from an experimental source, i.e. fasta files. These two options are described below.

Simulated reads - Reads can be simulated using the file  "createReads_package.m", that appears as part of the introductory examples.

Experimental reads – Reads should be prepared in a specific way for COMPASS. For an example of loading a fasta file and preparing the needed input for COMPASS, see file – " example_read_fasta.m "



**Computing weighted specificity and weighted sensitivity**

The file "example_evaluating_simulation_results.m" presents an example of calculating weighted specificity and sensitivity for a specific simulation, as used in Figure 3 and Figure 4 in the manuscript.


Acknowledgment
--------------

COMPASS was developed by: Noam Shental, Or Zuk, Amnon Amir and Amit Zeisel, as part of work on the papers,

[1] “High Resolution Microbial Community Reconstruction by Integrating Short Reads from Multiple 16S rRNA Regions“ A. Amir, A. Zeisel, O. Zuk, M. Elgart, S. Stern, O. Shamir, P.J. Turnbaugh, Y. Soen and N. Shental (submitted).

[2] ["Accurate Profiling of Microbial Communities from Massively Parallel Sequencing using Convex Optimization"](http://arxiv.org/abs/1309.6919), O. Zuk, A. Amir, A. Zeisel, O. Shamir and N. Shental ([SPIRE13](http://u.cs.biu.ac.il/~porately/spire2013/)).

Please cite the above papers if using the package.

For support, any questions or comments, please contact:

Noam Shental: shental@openu.ac.il

Or Zuk: or.zuk@mail.huji.ac.il



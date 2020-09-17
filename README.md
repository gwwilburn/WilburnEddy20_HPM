# HPM Proof of Principle Code

This is software for a proof of principle hidden Potts model implementation. We describe hidden Potts Models and the performance of this method in RNA remote homology search here: https://www.biorxiv.org/content/biorxiv/early/2020/06/23/2020.06.23.168153.full.pdf



## Contents and Dependencies:

Here is what is included.
* ./HPMHomology/ : code to build HPMs, align and score sequences using HPMs.
* ./CMIS/ : code to score and align sequences with Infernal CMs using importance sampling.
* ./data/ : benchmark datasets for tRNA, Twister ribozyme, and SAM riboswitch benchmarks.
* ./TwisterBenchmark_commands.txt: instructions for how to recreate the Twister ribozyme benchmark results.
* ./tRNABenchmark_commands.txt: instructions for how to recreate the tRNA benchmark results.
* ./SAMBenchmark_commands.txt: instructions for how to recreate the SAM riboswitch benchmark results.
* ./SimulationControl_commands.txt: instructions for how to create the synthetic sequence positive control results.
* ./clone.sh : bash script for cloning external dependencies from GitHub.
* ./compile.sh : bash script for compiling.

Here are the external dependencies, which are not included but needed to run the code.
* ./hmmer4/ : HMMER version 4 (in development) (must be cloned from GitHub, see below)
* ./infernal-1.1.3/ Infernal version 1.1.3 (must be cloned from GitHub, see below)


## Compiling instructions

### Programs required to compile:
 * git
 * make
 * gcc
 * autoconf

### Quick install

To clone Infernal and HMMER4, run
> ./clone.sh

Next, to compile, run:
> ./compile.sh

Or follow the step-by-step instructions below.

### Step-by-step install

#### Instructions for cloning Infernal from Git
> git clone https://github.com/EddyRivasLab/infernal.git
>
> mv infernal infernal-1.1.3
>
> cd infernal-1.1.3
>
> git checkout a358a49
>
> git clone https://github.com/EddyRivasLab/hmmer.git
>
> git clone https://github.com/EddyRivasLab/easel.git
>
> cd hmmer
>
> git checkout bcd05b0
>
> cd ../easel
>
> git checkout fa17f53
>
> cd ../..


#### Instructions for cloning HMMER4 from git
> git clone https://github.com/EddyRivasLab/hmmer.git
>
> mv hmmer hmmer4
>
> cd hmmer4
>
> git checkout 7994cc7
>
> git clone https://github.com/EddyRivasLab/easel.git
>
> cd easel
>
> git checkout 66745ea
>
> cd ../..

####  Instructions to compile Infernal
> cd infernal-1.1.3/
>
> autoconf
>
> ./configure
>
> make
>
> cd ..

#### Instructions to compile HMMER4
> cd hmmer4/easel
>
> make
>
> cd ..
>
> ln -s ${PWD}/easel lib/easel
>
> autoconf
>
> ./configure
>
> make
>
> cd nwo/
>
> make
>
> cd ../..

#### Instructions to compile HPMHomology code
> cd HPMHomology/src/
>
> make
>
> g++ -O3 -std=c++0x -o gremlin_h2_henikoff gremlin_h2_henikoff.cpp
>
> cd ../..

#### Instructions to compile CMIS code
> cd CMIS/src/
>
> make
>
> cd ../..

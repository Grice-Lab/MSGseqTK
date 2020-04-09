MSGseqTK introduction
=====================
MSGseqTK is a Metagenomics ShotGun Sequencing Toolkit for NGS read mapping and cleaning,
based on a novel incremental and parallel FMD-index.
The core algorithm is an improved bidirectional FM-index (FMD-index) in which traditionally not-allowed
null terminals ('\0') are used as sentinal between genomes/chromosomes so an imcremental, parallel and on-the-fly
FMD-index merge algorithm is possible.
This new algorithm allows MSGseqTK build, update and merge very large databases at whole metagenomic-scale with reasonable time and space.

The main program 'MSGseqTK' takes single or paired-end NGS FASTA/FASTQ files and map/align them to a pre-built metagenomics database.
It includes native supports for multi-threading (w/ OpenMP), .gz or .bz2 compressed inputs (w/ zlib), and direct write into binary BAM (or text SAM) output (w/ HTSLIB).

Implementation
--------------
MSGseqTK is written in pure C++11, and built with the GNU Autotools (autoconfig/automake), and can be easily installed under Linux, Windows and Mac OS X.

Download
--------
You can download the latest release from GitHub at: https://github.com/Grice-Lab/MSGseqTK/releases.
You can clone or fork and pull the source codes from GitHub at: https://github.com/Grice-Lab/MSGseqTK.

Dependencies
------------
MSGseqTK depends on the popular head-only C++ libraries Boost and Eigen3. They are available and often pre-installed on most Linux distributions, and can be easily installed on Windows and Mac OS X.
The HTSLIB C library (from SAMtools) is dependent for directly writing alignment into SAM/BAM files. Many Linux OS have HTSLIB in their software repositories. Try your `yum`, `dnf` or `apt-get` for details.
The ZLIB and Boost-IOSTREAMS libraries are optionally dependent for handling GZIP/GZIP2 compressed files, but are not required and can be disabled.

Installation
------------
1. Configure installation, by running the command
```bash
./configure
```
You may consider providing additional options, such as `--prefix`, `--exec-prefix`,
`--with-zlib`, `--with-boost`, `--with-htslib`, etc.

2. Compile and link, by running the command
```bash
make
```
Look for errors and try to resolve them yourself first (by Google).
Contact us only if you are sure it is a bug in our programs.

3. (Optionally) Test, by running the command
```bash
make check
```
It may take a while depending on your processor's speed.

4. Install
```bash
make install
```
You may need root privilege to do it, such as using `sudo`.

Pre-built databases
-------------------
You need to build an MSGseqTK database before using it tools for mapping and clean metagenomic NGS reads.
You can build your own database using `msgseqtk-build`, or alternatively [download the pre-built databases](https://www.med.upenn.edu/gricelab/msgseqtk.html#databases "Pre-built databases").

Documentations
--------------
Please check the help and documentations at MSGseqTK's [home page](https://www.med.upenn.edu/gricelab/msgseqtk.html "MSGseqTK home")


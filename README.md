Pandoo
============

To get help::

    pandoo -h

![Pandoo](http://bioinformatics.mdu.unimelb.edu.au/~schultzm/pandoo/pipeline.png)

Authors
-------

Mark B Schultz, Anders Gonçalves da Silva, Jason Kwong, Torsten Seemann

Introduction
------------

**Pandoo** is a command-line tool for for exploring and characterising bacterial whole genome DNA sequence data.  It is a computational pipeline written in python and is scalable via implementation using the ruffus pipeline library.  **Ruffus** handles task-scheduling and task-parallelisation during the run.  Pandoo is particularly useful for pathogenic species as it performs AMR and MLST profiling but in theory could be used for any bacterial species against any sequence database.  It was originally written to characterise and QC assemblies and reads of bacterial isolates at the Microbiological Diagnostic Unit Public Health Laboratory, Victoria, Australia.  Input is a tab-delimited text file that points the software to the assembly and paired-end read files for each isolate.  Specifically, the file has no header, one line per isolate, and four columns per line in the following column order::  

    +-------------+------------------------------------------------+----------------------------+---------------------------+
    |isolate_name | full_path_to_assembly (i.e., a 'contigs' file) | full_path_to_paired_reads1 | full_path_to_paired_reads2|
    +-------------+------------------------------------------------+----------------------------+---------------------------+

An example **isolates.tab** file looks like this::  

    isolate-1	/path/to/isolate-1/contigs.fa	/path/to/isolate-1/reads-1.fq.gz	/path/to/isolate-1/reads-2.fq.gz
    isolate-2	/path/to/isolate-2/contigs.fa	/path/to/isolate-2/reads-1.fq.gz	/path/to/isolate-2/reads-2.fq.gz
    isolate-3	/path/to/isolate-3/contigs.fa	/path/to/isolate-3/reads-1.fq.gz	/path/to/isolate-3/reads-2.fq.gz
    isolate-4	/path/to/isolate-4/contigs.fa	/path/to/isolate-4/reads-1.fq.gz	/path/to/isolate-4/reads-2.fq.gz
    isolate-5	/path/to/isolate-5/contigs.fa	/path/to/isolate-5/reads-1.fq.gz	/path/to/isolate-5/reads-2.fq.gz


Installing dependencies
-----------------------

The following packages need to be installed before pip3 installing Pandoo.  To install dependencies, do::  

    cpan -i Moo
    cpan -i List::MoreUtils
    cpan -i Bio::Perl
    brew tap homebrew/science
    brew tap tseemann/homebrew-bioinformatics-linux
    brew update
    brew install mlst --HEAD
    brew install abricate --HEAD
    brew install seqtk
    brew install mummer
    brew install bowtie2
    brew install cd-hit
    brew install ariba
    brew install kraken

Follow the instructions at https://ccb.jhu.edu/software/kraken to set up the databases.

Additionally, mashtree.pl needs to be installed. Follow the instructions at https://github.com/lskatz/mashtree
Add mashtree.pl to your path and ensure that mashtree.pl can be executed by typing on the command line:
    mashtree.pl


Installing Pandoo
-----------------------

To perform any of these install steps **for all users, remove '--user'**.  The final symlink step is not required if installing for all users.  Pandoo is written for **python3** and installation requires **pip3** and **setuptools**.  To install the latest 'stable' version of pandoo for the current user only, do::  

    pip3 install pandoo --user

To upgrade::  

    pip3 install pandoo --user --upgrade

To install the latest, potentially unstable, bleeding-edge version::  

    pip3 install --user https://github.com/schultzm/pandoo/zipball/master

**If installing via the '--user' option**
Check where the executable is::  

    which pandoo  # ~/.local/bin/pandoo

Check where the site-packages are::  

    python3 -m site --user-site  # ~/.local/lib/python3.6/site-packages

Now, symlink the packaged databases in site-packages above to the folder containing the executable shown above::  

    ln -s ~/.local/lib/python3.6/site-packages/Pandoo/CARD/ ~/.local/bin
    ln -s ~/.local/lib/python3.6/site-packages/Pandoo/VFDB/ ~/.local/bin
    ln -s ~/.local/lib/python3.6/site-packages/Pandoo/VFDB_core/ ~/.local/bin
    ln -s ~/.local/lib/python3.6/site-packages/Pandoo/plasmidfinder/ ~/.local/bin


**To uninstall**, do::  

    pip3 uninstall pandoo
    # Remove the symlinks
    rm -r ~/.local/bin/CARD
    rm -r ~/.local/bin/VFDB*
    rm -r ~/.local/bin/plasmidfinder

Now, check that the dependencies are in the path using::

    pandoo check

Quickstart tutorial
-------------------

**For large jobs, run in screen mode**, otherwise skip this step::  

    screen -SL screenname 

All screen-output within the screen started above, will be saved to::  

    screenlog.%n #  Where %n is the number of the screen.

After exiting the screen, the screenlog.0 can be viewed as the run progresses using standard command line actions::  

    tail screenlog.0
    less screenlog.0
    watch "tail screenlog.0"

**To get help for pandoo**, just do::  

    pandoo -h

Output looks like::  

    usage: pandoo <command> <options>

    This is a tool for exploring your bacterial genome data. Given some assemblies
    and/or paired-end read sets, run a pipeline of software tools to generate an
    NJ tree from assemblies and a complementary table of metadata/results (contig
    and read QC, mlst, species ID, resistance genes, virulence factors, plasmid
    replicon types).

    optional arguments:
      -h, --help     show this help message and exit
      -v, --version  Print version and quit.

    Commands:
  
        check        Check pipeline dependencies
        input        Generate input table
        run          Run the pandoo pipeline
        merge        Merge two metadata tables

Notice, above, four modules: one each for **check**, **input**, **run** and **merge**.  Each can be run independently. 

**check module**
This module will check if the required softwares are installed and executable as per the calls to these programs used by pandoo.

**input module**

This module is used to generate the isolates.tab file.  Final output from this command is sent to **standard out (stdout)**.  To capture the information from stdout, redirect it to a file using '> isolates.tab'::  

    pandoo input -h
    pandoo input -i isolates.txt > isolates.tab

**run module**

This module is used to run the analysis pipeline.  In the example below, output will be a folder called **results** and we have selected to use the **tree option with -t** using a **JC** model of evolution (default for the -t option, but user can also choose from **Raw** or **Kimura**)::  

    pandoo run -h
    pandoo run -i isolates.tab -o results -t

In the results folder there is a sub-folder for each isolate containing the results for each isolate.  Also within the results folder there are three files::  

    results/isolates_metadataAll.csv
    results/isolates_metadataAll_simplified.csv
    results/isolates_mashtree.tre

**merge module**

This module is used to join existing metadata tables (e.g., LIMS table) with one of the output tables from the run module.  In this example, an excel file (AMR\_ongoing\_20170307.xlsx, with a header on line 5 and first column containing isolate names that were expanded with a wildcard search during the **pandoo input** step above) is joined with the results/isolates\_metadataAll\_simplified.csv table.  Again, final output is sent to stdout, which in this example is redirected to file using '> results/AMR\_ongoing\_20170307\_join\_isolates\_metadataAll\_simplified.csv'::  

    pandoo merge -h
    pandoo merge -l AMR_ongoing_20170307.xlsx -r results/isolates_metadataAll_simplified.csv > results/AMR_ongoing_20170307_join_isolates_metadataAll_simplified.csv

Why use Pandoo?
---------------------

Pandoo provides a single table of metadata and a tree that the metadata can be plotted next to using **phandango**  

https://jameshadfield.github.io/phandango/

![alt tag](http://bioinformatics.mdu.unimelb.edu.au/~schultzm/pandoo/phandango.png)

The single NJ tree for the whole isolate set is inferred using the program **quicktree** from a distance matrix computed by **andi** using any of the evolutionary models JC, Raw or Kimura.  Andi infers the distance matrix from the assemblies (typically contigs.fa files).  The tree will not include isolates for which the assembly file is missing.

The summary table (csv) for all isolates combines the results from:

1. inferred species from running kraken on the reads  
2. inferred species call from running kraken on the contigs (assemblies)  
3. Inferred consensus species call from a consensus of the best hit from kraken on the reads and kraken on the contigs  
4. Based on the species call, runs mlst using the appropriate scheme (if available) or autodetects scheme  
5. Gene content profiles using unlimited number of user databases of resistance genes, plasmid rep genes, virulence genes, etc., using abricate (contigs, BLAST contigs against database, assembly based) and ariba (reads, MiniMapping to database, mapping based)  
6. QC metrics using seqtk for reads and contigs  
7. Reports software versions and paths to databases used in the analysis (for repeatability)  
8. The pipeline is modular in that the user can choose not to perform the tree inference step with andi plus quicktree and/or the can choose not to perform the read mapping step using ariba  
9. A flowchart is produced for the run (however, if the total path length of the results folders combined exceeds 16384 characters then the flowchart cannot be drawn)  
10. The user can supply reads and/or contigs for each file.  The final tree will only include taxa for which contigs have been supplied  


Assumptions
-----------

First you need to download assemblies or perform the assemblies yourself from readsets (using e.g., unicycler https://github.com/rrwick/Unicycler, which uses SPAdes (http://bioinf.spbau.ru/spades) or MegaHit (https://github.com/voutcn/megahit)).  It doesn't necessarily make sense to supply reads for an isolate but contigs that have been assembled using a readset other than the one supplied.  If you don't have reads, leave the columns blank for that isolate (for example, if you just want to characterise assemblies downloaded from NCBI GenBank). If you don't have contigs and only have reads, leave the column blank for contigs (but without contigs there will be no tree).  

Licence
-------

Pandoo is:

| Copyright (C) 2017 Mark B Schultz  
| https://github.com/schultzm/  
| email: dr.mark.schultz@gmail.com  

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.
You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

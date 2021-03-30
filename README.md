# CCD

Crystallisation construct designer

## Prerequisites

The following prerequisites are required by ccd and must be installed manually:

* iupred
* predator (source code can be found at http://ftp.ebi.ac.uk/pub/software/unix/predator/)
* ncoils (e.g. `apt-get install ncoils` on a Debian-like OS)
* muscle (e.g. `apt-get install muscle` on a Debian-like OS)
* rabbitmq (e.g. `apt-get install rabbitmq-server` on a Debian-like OS)
* ncbi-blast (e.g. apt-get install ncbi-blast+ on a Debian-like OS)

Development prerequisites:

* git
* virtualenv
* virtualenvwrapper

# Installation of a development environment

Clone the repository from gitlab and go to the ccd directory:

    git clone http://gitlab.xtal.nki.nl/ccd/ccd.git
    cd ccd

Create a python3.6 virtual environment containing the dependences found in requirements.txt:
We suggest using Anaconda or virtualenv

Create environment with Anaconda(tested with version >= 4.6.2):

    conda create -n ccd3 python=3.6
    conda activate ccd3

then install the dependencies, starting with numpy

    pip install -r requirements


```
To perform homology searches, CCD uses a local mirror of the uniprot database.
Download the swissprot database files to a folder of choice, then decompress it
and compile it for blast. About 350 Mb of disk space are required).  
    
    #download (if outside europe, use ftp://uniprot.org instead)
    wget ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    #decompress
    gunzip uniprot_sprot.fasta.gz
    #compile
    makeblastdb -dbtype prot -in uniprot_sprot.fasta -out local_swissprot
    
(Optional) Download the Trembl database (vastly bigger, non curated, ~70 Gb)
    
    wget ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
    gunzip uniprot_trembl.fasta.gz
    makeblastdb -dbtype prot -in uniprot_trembl.fasta -out local_trembl
    
(Optional) For PDB support in similarity searches, run the following script. Note that this script needs to 
download and process a file of ~120 Mb, so might require some time on slow connections.We would also suggest 
to setup a periodic database update to include all new PDB releases, via i.e. a cron job.


    python ccd/update_blast_db.py

by default, the script will install the databases in ./ccd/static/blastdb; use --destination_folder=path to specify an
installation folder.
    

Modify the settings in `config/settings.example` to your needs and rename the
file to `config/ccd_settings.cfg`. The secret key can be generated using
a python shell, for example:

	```python
	import os
	os.urandom(96)

Start the app!

    ./run.sh



Modify the settings in `config/settings.example` to your needs and rename the
file to `config/ccd_settings.cfg`. Make sure you generate a new secret key,
and fill in a valid email to allow Uniprot, NCBI and EBI to notify you in case of
problems.

Optionally, copy the settings file to `/etc/ccd.cfg`.

If needed, configure a reverse proxy. `config/apache.conf` is an example for
apache.

Optionally, install `supervisor` and use the `config/supervisor.example` to
make a system service ccd managed by supervisor (edit where necessary).

## Deployment

A fabric file `fabfile.py` has been provided to automate deployment steps. Assuming a production version is running on <server> by user <user> using `supervisor` without virtual environment, issue the following to release a new version:

    fab --user=<user> --hosts=<server> deploy

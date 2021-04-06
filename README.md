# CCD core

CCD stands for Crystallisation Construct Designer. 

CCD is a web application that facilitates protein crystallography experiments. I developed it during my employment at the Netherlands Cancer Institute. In this repository you can find the core of the application.

You can visit the full site https://ccd.rhpc.nki.nl

Basic web software used in this project:
* Flask
* JQuery
* AJAX
* Google Material Design
* Bootstrap

# Installation for development on your computer

## Prerequisites

The following prerequisites are required by CCD and must be installed manually if not present:

* predator (already uploaded in this repository's deps directory)
* uipred (already uploaded in this repository's deps directory)
* ncoils (e.g. `apt-get install ncoils` on a Debian-like OS)
* git (e.g. `apt-get install git` on a Debian-like OS)
* python 3.8

## Steps for setting up the development environment

Clone the repository from github and go to the ccd directory:

    git clone https://github.com/gdamaskos/flask_with_ajax_and_matrix_visualization.git
    
    cd flask_with_ajax_and_matrix_visualization

Create a python3.8 virtual environment containing the dependences found in requirements.txt:

    python3 -m venv ~/virtual_environments/ccd3
    source ~/virtual_environments/ccd3/bin/activate

then install python dependencies:

    pip install numpy
    pip install wheel
    pip install -r requirements.txt

Modify the settings in `config/settings.example` to your needs and rename the
file to `config/ccd_settings.cfg`. Make sure you generate a new secret key,
and fill in a valid email to allow Uniprot, NCBI and EBI to notify you in case of
problems.

The secret key can be generated using a python shell, for example:

	$ python
	>>> import os
	>>> os.urandom(96)

Start the app!

    ./run_autoreload.sh


If needed, configure a reverse proxy. `config/apache.conf` is an example for
apache.

Optionally, copy the settings file to `/etc/ccd.cfg`.

Optionally, install `supervisor` and use the `config/supervisor.example` to
make a system service ccd managed by supervisor (edit where necessary).

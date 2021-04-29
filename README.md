# CCD core

CCD stands for Crystallisation Construct Designer. 

CCD is a web application that facilitates protein crystallography experiments. In this repository you can find the core of the application.

You can visit the complete site https://ccd.rhpc.nki.nl

Basic web software used in this project:
* Flask
* JQuery
* AJAX
* Google Material Design
* Bootstrap

# Installation for development on your computer

## Prerequisites

First of all, you must be running a Debian Operating System like Ubuntu.

The following prerequisites are required by CCD and must be installed manually if not present:

* predator (already uploaded in this repository's deps directory)
* uipred (already uploaded in this repository's deps directory)
* ncoils (e.g. `apt-get install ncoils` on a Debian-like OS)
* muscle (e.g. `apt-get install muscle` on a Debian-like OS)
* git (e.g. `apt-get install git` on a Debian-like OS)
* python 3.8, follow this tutorial to prepare https://www.digitalocean.com/community/tutorials/how-to-install-python-3-and-set-up-a-programming-environment-on-ubuntu-20-04-quickstart

## Steps for setting up the development environment

Clone the repository from github and go to the ccd directory:

    git clone https://github.com/gdamaskos/flask_with_ajax_and_matrix_visualization.git
    
    cd flask_with_ajax_and_matrix_visualization

Create a python3.8 virtual environment containing the dependences found in requirements.txt:

    python3 -m venv ~/virtual_environments/ccd3
    source ~/virtual_environments/ccd3/bin/activate

Install python dependencies:

    pip install wheel
    pip install -r requirements.txt

Modify the settings in `config/settings.example` to your needs and rename the
file to `config/ccd_settings.cfg`. Make sure you generate a new secret key,
and fill in a valid email to allow Uniprot to notify you in case of
problems.

The secret key can be generated using a python shell, for example:

	$ python
	>>> import os
	>>> os.urandom(96)

Run the tests

    ./test.sh

Start the app!

    ./run_autoreload.sh

Visit the address:

    http://localhost:5000

# Further options for deployment

If you would like to create `systemd` services to run CCD automatically see the `config/celery.service` and `config.ccd3.service`.

As an alternative, you can use python `supervisor` and use the `config/supervisor.conf` to make a system service ccd managed by supervisor (edit where necessary).

If needed, you cab configure a reverse proxy using `config/apache.conf` as an example for apache web server.

The contact page (available only on developer versions for ccd.rhpc.nki.nl) requires a working Google maps API. 
(see https://developers.google.com/maps/api-key-best-practices for how to obtain and configure one)
In short:
1. Generate an API key
2. Restrict it by HTML referrer to the website you are running CCD from + /contact (e.g. ccd.rhpc.nki.nl/contact)
3. Enable "Map Javascript API".
4. Add your API key to contact.html template in the place of the current key.

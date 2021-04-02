"""Make deployment easier.
To execute all tasks at once (update source, update requirements,
restart gunicorn), issue

fab --user=<user> --hosts=host1,...
"""


from fabric.api import cd, prefix, run, task
from fabric.contrib.files import exists


def update_source():
    """Updates the source code on the host.
    If the git repository already exists on the host, `git pull` is used;
    otherwise, `git clone` is used.
    """
    if exists('/var/www/apacheroot/html/ccd'):
        with cd('/var/www/apacheroot/html/ccd'):
            run('git pull')
    else:
        with cd('/var/www/apacheroot/html'):
            run('git clone https://gitlab.rhpc.nki.nl/ccd/ccd.git')

def update_requirements(within_venv=False):
    """Updates the virtualenv on the host."""
    with cd('/var/www/apacheroot/html/ccd'):
        if within_venv:
            with prefix('source venv/bin/activate'):
                run('pip install -r requirements')
        else:
            run('pip install -r requirements')


@task
def restart_gunicorn():
    """Restarts gunicorn on the host using supervisorctl."""
    run('supervisorctl restart ccd')


@task
def deploy():
    update_source()
    update_requirements()
    restart_gunicorn()

import logging
import os
from logging import FileHandler
from flask import Flask
from flask_debugtoolbar import DebugToolbarExtension
from celery import Celery

# Create the top-level logger.
_log = logging.getLogger(__name__)
log_file = 'ccd.log'
_file = FileHandler(log_file, 'a')
_log.addHandler(_file)
_formatter = logging.Formatter(
        '%(asctime)s | %(levelname)-7s | %(message)s '
        '[in %(pathname)s:%(lineno)d]')
_log.setLevel(logging.INFO)
_file.setFormatter(_formatter)

_log.info('Creating app')
app = Flask(__name__)
app.config.from_envvar('CCD_SETTINGS')
_log.info('Imported settings from {}'.format(os.environ['CCD_SETTINGS']))

# Celery configuration to use RabbitMQ broker for messaging
app.config['CELERY_BROKER_URL'] = 'amqp://localhost'
app.config['CELERY_RESULT_BACKEND'] = 'rpc://'
celery = Celery(app.name, backend=app.config['CELERY_RESULT_BACKEND'], broker=app.config['CELERY_BROKER_URL'], include=['ccd.celery_tasks'])
celery.conf.update(app.config)

# Use ProxyFix to correct URL's when redirecting.
from ccd.middleware import ReverseProxied
app.wsgi_app = ReverseProxied(app.wsgi_app)

# Initialize extensions
toolbar = DebugToolbarExtension()

if app.debug:
    _log.setLevel(logging.DEBUG)
    _log.debug('Initializing toolbar')
    toolbar.init_app(app)

# Template reload on change, if the following two lines are removed we need to refresh manually the jinja templates,
# every time we change them by restarting the server.
app.jinja_env.auto_reload = True
app.config['TEMPLATES_AUTO_RELOAD'] = True

_log.debug('Created app')

# The views should be loaded
# after the app has been created this is not accepted in pep8
import ccd.views

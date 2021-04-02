#!/usr/bin/env bash
export CCD_SETTINGS=../config/ccd_settings.cfg
celery worker --app=ccd:celery --loglevel=info --purge & # celery for process queueing
gunicorn -b 0.0.0.0:5000 --log-level debug --log-file "-" -k gevent ccd:app

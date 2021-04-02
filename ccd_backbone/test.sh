#!/usr/bin/env bash
export CCD_SETTINGS=../config/ccd_settings.cfg
python -m unittest discover tests

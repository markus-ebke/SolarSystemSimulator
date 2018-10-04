#!/bin/bash
# pep8: no line length check
pep8 --ignore=E501 SolarSystemSimulator04.py > SSSim04_lint.txt

# pylint: disable "Invalid Name", "Missing Docstring", "Line to long" and "Class has no __init__ method"
pylint --disable=C0103,C0111,C0301,W0232 SolarSystemSimulator04.py >> SSSim04_lint.txt

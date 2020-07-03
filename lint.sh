#!/bin/bash
mkdir lint

# pycodestyle (pep8 check): no line length check
echo pycodestyle
pycodestyle solar_system_simulator/__init__.py > lint/init.txt
pycodestyle solar_system_simulator/calculation.py > lint/calculation.txt
pycodestyle solar_system_simulator/opengl.py > lint/opengl.txt
pycodestyle solar_system_simulator/operators.py > lint/operators.txt
pycodestyle solar_system_simulator/panels.py > lint/panels.txt
pycodestyle solar_system_simulator/properties.py > lint/properties.txt

# pylint: disable "Invalid Name", "Missing Docstring", "Line to long" and "Class has no __init__ method"
echo pylint
pylint --disable=C0103,C0111,C0301,W0232 solar_system_simulator/__init__.py >> lint/init.txt
pylint --disable=C0103,C0111,C0301,W0232 solar_system_simulator/calculation.py >> lint/calculation.txt
pylint --disable=C0103,C0111,C0301,W0232 solar_system_simulator/opengl.py >> lint/opengl.txt
pylint --disable=C0103,C0111,C0301,W0232 solar_system_simulator/operators.py >> lint/operators.txt
pylint --disable=C0103,C0111,C0301,W0232 solar_system_simulator/panels.py >> lint/panels.txt
pylint --disable=C0103,C0111,C0301,W0232 solar_system_simulator/properties.py >> lint/properties.txt

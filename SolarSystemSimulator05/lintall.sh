#!/bin/bash
mkdir lint
# pep8: no line length check
pep8 --ignore=E501 __init__.py > lint/init_lint.txt
pep8 --ignore=E501 calculation.py > lint/calculation_lint.txt
pep8 --ignore=E501 opengl.py > lint/opengl.txt
pep8 --ignore=E501 operators.py > lint/operators_lint.txt
pep8 --ignore=E501 panels.py > lint/panels_lint.txt
pep8 --ignore=E501 properties.py > lint/properties_lint.txt

# pylint: disable "Invalid Name", "Missing Docstring", "Line to long" and "Class has no __init__ method"
pylint --disable=C0103,C0111,C0301,W0232 __init__.py >> lint/init_lint.txt
pylint --disable=C0103,C0111,C0301,W0232 calculation.py >> lint/calculation_lint.txt
pylint --disable=C0103,C0111,C0301,W0232 opengl.py >> lint/opengl.txt
pylint --disable=C0103,C0111,C0301,W0232 operators.py >> lint/operators_lint.txt
pylint --disable=C0103,C0111,C0301,W0232 panels.py >> lint/panels_lint.txt
pylint --disable=C0103,C0111,C0301,W0232 properties.py >> lint/properties_lint.txt

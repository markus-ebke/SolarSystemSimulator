#!/bin/bash
dst_dir=/home/markus/.config/blender/2.68/scripts/addons_extern/SolarSystemSimulator05
mkdir ${dst_dir}
cp __init__.py ${dst_dir}
cp calculation.py ${dst_dir}
cp opengl.py ${dst_dir}
cp operators.py ${dst_dir}
cp panels.py ${dst_dir}
cp properties.py ${dst_dir}

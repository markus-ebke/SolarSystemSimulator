# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>

# Blenderartist Thread:
# http://blenderartists.org/forum/showthread.php?267761-Solar-System-Simulator-WIP

bl_info = {
    "name": "Solar System Simulator",
    "description": "Simulation of solar systems using Kepler's laws of planetary motion",
    "author": "Markus Ebke",
    "version": (0, 6),
    "blender": (2, 72, 0),
    "location": "Properties > Physics (for objects); Properties > Scene (global settings); 3DView Toolbar",
    "warning": "In developement, version from 2014-11-05",
    "wiki_url": "",
    "tracker_url": "",
    "category": "Object"}

if "bpy" in locals():
    import imp
    imp.reload(calculation)
    imp.reload(operators)
    imp.reload(panels)
    imp.reload(properties)
    imp.reload(opengl)
    print("Reloaded SSSim modules")
else:
    from . import calculation, operators, panels, properties, opengl
    print("Imported SSSim modules")

import bpy
from bpy.props import PointerProperty


### Registration ###

def register():
    bpy.utils.register_module(__name__)

    obj = bpy.types.Object
    obj.sssim_obj = PointerProperty(type=properties.SSSIMObject)
    obj.sssim_orbit = PointerProperty(type=properties.SSSIMOrbit)
    obj.sssim_rotation = PointerProperty(type=properties.SSSIMRotation)
    obj.sssim_calc = PointerProperty(type=properties.SSSIMCalculation)
    obj.sssim_surface = PointerProperty(type=properties.SSSIMSurface)
    bpy.types.Scene.sssim_scn = PointerProperty(type=properties.SSSIMScene)

    # We need to add custom drivers for the location and rotation
    drv = bpy.app.driver_namespace
    drv["eval_planet_orbit"] = calculation.eval_planet_orbit
    drv["eval_planet_rotation"] = calculation.eval_planet_rotation

    bpy.types.PHYSICS_PT_add.append(panels.physics_panel)


def unregister():
    bpy.types.PHYSICS_PT_add.remove(panels.physics_panel)

    # Remove the drivers if any
    drv = bpy.app.driver_namespace
    if "eval_planet_orbit" in drv:
        del drv["eval_planet_orbit"]
    if "eval_planet_rotation" in drv:
        del drv["eval_planet_rotation"]

    # Remove drawing callbacks in all scenes
    for scn in bpy.data.scenes:
        scn.sssim_scn.draw_orbit = False

    obj = bpy.types.Object
    del obj.sssim_obj
    del obj.sssim_orbit
    del obj.sssim_rotation
    del obj.sssim_calc
    del obj.sssim_surface
    del bpy.types.Scene.sssim_scn

    bpy.utils.unregister_module(__name__)


if __name__ == "__main__":
    register()

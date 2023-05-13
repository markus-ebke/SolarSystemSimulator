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

bl_info = {
    "name": "Solar System Simulator",
    "description": "Simulate solar systems using Kepler's laws",
    "author": "Markus Ebke",
    "version": (0, 8),
    "blender": (2, 80, 0),
    "location": "Properties > Physics; Properties > Scene; 3D-View Toolbar",
    "warning": "",
    "doc_url": "https://blenderartists.org/t/solar-system-simulator/553099",
    "tracker_url": "",
    "category": "Object"}

"""
# support reloading
if "bpy" in locals():
    import importlib
    if "calculation" in locals():
        importlib.reload(calculation)
    if "operators" in locals():
        importlib.reload(operators)
    if "panels" in locals():
        importlib.reload(panels)
    if "properties" in locals():
        importlib.reload(properties)
    if "opengl" in locals():
        importlib.reload(opengl)
"""

import bpy
from bpy.props import PointerProperty

from . import calculation, operators, panels, properties
print("Successfully imported Solar System Simulator modules")


# =============================================================================
# Registration
# =============================================================================
classes = (
    operators.SCENE_OT_add_sim_time_fcurve,
    operators.OBJECT_OT_sssim_to_fcurve,
    operators.OBJECT_OT_sssim_clear_fcurve,
    operators.SCENE_OT_update_sssim_drivers,
    operators.OBJECT_OT_sssim_create_center,
    operators.OBJECT_OT_sssim_create_planet,
    operators.OBJECT_OT_sssim_create_surface,
    operators.OBJECT_OT_sssim_bake_all,
    operators.OBJECT_OT_sssim_bake_clear,
    panels.SSSIM_PT_object,
    panels.SSSIM_PT_orbit,
    panels.SSSIM_PT_rotation,
    panels.SSSIM_PT_surface,
    panels.SSSIM_PT_calculation,
    panels.SSSIM_PT_scene,
    panels.SSSIM_PT_tools,
    properties.SSSIMObject,
    properties.SSSIMOrbit,
    properties.SSSIMRotation,
    properties.SSSIMCalculation,
    properties.SSSIMSurface,
    properties.SSSIMScene,
)


def register():
    for cla in classes:
        bpy.utils.register_class(cla)

    obj = bpy.types.Object
    obj.sssim_obj = PointerProperty(type=properties.SSSIMObject)
    obj.sssim_orbit = PointerProperty(type=properties.SSSIMOrbit)
    obj.sssim_rotation = PointerProperty(type=properties.SSSIMRotation)
    obj.sssim_calc = PointerProperty(type=properties.SSSIMCalculation)
    obj.sssim_surface = PointerProperty(type=properties.SSSIMSurface)
    bpy.types.Scene.sssim_scn = PointerProperty(type=properties.SSSIMScene)

    # add custom drivers for the location and rotation to driver namespace
    drv = bpy.app.driver_namespace
    drv["eval_planet_orbit"] = calculation.eval_planet_orbit
    drv["eval_planet_rotation"] = calculation.eval_planet_rotation

    bpy.types.PHYSICS_PT_add.append(panels.physics_panel)


def unregister():
    bpy.types.PHYSICS_PT_add.remove(panels.physics_panel)

    # remove the drivers if any
    drv = bpy.app.driver_namespace
    if "eval_planet_orbit" in drv:
        del drv["eval_planet_orbit"]
    if "eval_planet_rotation" in drv:
        del drv["eval_planet_rotation"]

    # remove drawing callbacks in all scenes
    for scn in bpy.data.scenes:
        scn.sssim_scn.draw_orbit = False

    obj = bpy.types.Object
    del obj.sssim_obj
    del obj.sssim_orbit
    del obj.sssim_rotation
    del obj.sssim_calc
    del obj.sssim_surface
    del bpy.types.Scene.sssim_scn

    for cla in classes:
        bpy.utils.unregister_class(cla)


if __name__ == "__main__":
    register()

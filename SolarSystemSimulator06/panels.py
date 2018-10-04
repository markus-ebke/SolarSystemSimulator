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

import bpy
from math import degrees

from .operators import has_location_fcurve, has_rotation_fcurve


### UI Panels ###

# Panel for general Object settings
class SSSIMPanelObject(bpy.types.Panel):
    bl_label = "Solar System Simulator: Object"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "physics"

    @classmethod
    def poll(cls, context):
        return context.object.sssim_obj.use_sssim

    def draw(self, context):
        layout = self.layout
        obj = context.object
        simobj = obj.sssim_obj
        info = simobj.show_info

        if bpy.ops.scene.add_eval_time_fcurve.poll():
            box = layout.box()
            box.label("You need to add an Evaluation Time F-Curve",
                      icon='ERROR')
            box.operator("scene.add_eval_time_fcurve", icon='FCURVE')
            box.label("You can find the Settings in the Scene Tab",
                      icon='SCENE_DATA')

        layout.active = simobj.use_sssim
        layout.prop(simobj, "show_info", icon='INFO')
        layout.prop(simobj, "object_type", expand=True)

        # planet only props
        if simobj.object_type == 'PLANET':
            simorbit = obj.sssim_orbit
            layout.prop_search(simorbit, "center_name", context.scene, "objects")
            if not simorbit.center_is_valid:
                layout.label("No valid center", icon='ERROR')

        # props for center and planet but not surface
        if simobj.object_type in ('CENTER', 'PLANET'):
            row = layout.row(align=True)
            row.prop(simobj, "mass_mantissa")
            row.prop(simobj, "mass_exp", slider=True)
            if info:
                mantissa = round(simobj.mass_mantissa, 3)
                exp = simobj.mass_exp
                layout.label("Mass: {} * 10^{} kg".format(mantissa, exp))
        elif simobj.object_type == 'SURFACE':
            simsur = obj.sssim_surface
            layout.prop_search(simsur, "parent_name", context.scene, "objects")
            if not simsur.parent_is_valid:
                layout.label("Parent object is not valid", icon="ERROR")


# Panel for Orbit settings
class SSSIMPanelOrbit(bpy.types.Panel):
    bl_label = "Solar System Simulator: Orbit"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "physics"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        obj = context.object
        simobj = obj.sssim_obj
        simorbit = obj.sssim_orbit
        if simobj.use_sssim and simobj.object_type == "PLANET":
            return simorbit.center_object is not None

    def draw(self, context):
        layout = self.layout
        obj = context.object
        simorbit = obj.sssim_orbit
        info = obj.sssim_obj.show_info

        if has_location_fcurve(obj):
            # Deactivate layout if location fcurve found
            layout.label(text="Found Orbit F-Curve", icon='INFO')
            layout.active = False
        elif not obj.sssim_calc.use_driver:
            # There is no orbit because driver not active
            layout.label(text="Orbit Driver not Active", icon='ERROR')

        # Shape and Size:
        simscn = context.scene.sssim_scn
        self.draw_shape(simorbit, simscn, layout.box(), info)

        layout.separator()

        split = layout.split()

        # Orientation:
        col = split.column()
        self.draw_orientation(simorbit, col.box())

        # Position:
        col = split.column()
        self.draw_position(simorbit, col.box(), info)

        layout.separator()

        # Orbital Period:
        self.draw_orbital_period(simorbit, layout.box(), info)

    def draw_shape(self, simorbit, simscn, layout, info):
        layout.label(text="Shape and Size of the Orbit:")

        split = layout.row().split(percentage=0.6)
        row = split.row()
        row.prop(simorbit, "eccentricity", slider=True)

        row = split.row()
        ecc = simorbit.eccentricity
        if info:
            if ecc == 0:
                row.label(text="Circular Orbit", icon='SPHERECURVE')
            else:
                row.label(text="Elliptical Orbit", icon='ROOTCURVE')

        row = layout.row(align=True)
        if ecc == 0:
            row.prop(simorbit, "distance_mantissa", text="Radius (km)")
        else:
            row.prop(simorbit, "distance_mantissa")
        row.prop(simorbit, "distance_exp", slider=True)
        if info:
            axis_mantissa = round(simorbit.distance_mantissa, 3)
            axis_exp = simorbit.distance_exp
            # lengths in km / length_mult = lengths in Blender Units
            axis_bu = round(simorbit.semi_major_axis / simscn.length_mult, 2)
            data = "{} * 10^{} km ( ={:,} BU)"
            data = data.format(axis_mantissa, axis_exp, axis_bu)
            if ecc == 0:
                msg = "Radius: {}".format(data)
                layout.label(text=msg, icon='CURVE_BEZCIRCLE')
            else:
                msg = "Semi-Major Axis: {}".format(data)
                layout.label(text=msg, icon='CURVE_BEZCIRCLE')

    def draw_orientation(self, simorbit, layout):
        layout.label(text="Orientation of Orbital Plane:")
        subcol = layout.column(align=True)
        subcol.prop(simorbit, "inclination")
        subcol.prop(simorbit, "asc_node")
        subcol.prop(simorbit, "arg_periapsis")

    def draw_position(self, simorbit, layout, info):
        layout.label(text="Position:")
        subcol = layout.column()
        subcol.prop(simorbit, "time_offset", slider=True)
        if info:
            # True Anomaly (angle travelled)
            angle = degrees(simorbit.true_anomaly)
            subcol.label(text="True Anomaly: {0}°".format(round(angle, 3)),
                         icon='FORCE_HARMONIC')

    def draw_orbital_period(self, simorbit, layout, info):
        layout.label(text="Orbital Period:", icon='SORTTIME')

        if info:
            orb_sec = round(simorbit.orbital_period)
            # orb_years = orbital period in julian years (365.25 days)
            orb_years = round(simorbit.orbital_period / 86400 / 365.25, 3)
            orb_frames = round(simorbit.orbital_period_frames, 2)

            if simorbit.use_user_period and simorbit.use_frames:
                # Don't show frames again when the user uses frames
                orb_text = "{0:,} s = {1} y".format(orb_sec, orb_years)
            else:
                orb_text = "{0:,} s = {1} y (= {2} frames)"
                orb_text = orb_text.format(orb_sec, orb_years, orb_frames)
            layout.label(text=orb_text)

        row = layout.row()
        row.prop(simorbit, "use_user_period")
        if simorbit.use_user_period:
            row.prop(simorbit, "use_frames")

            if simorbit.use_frames:
                layout.prop(simorbit, "user_period_frames")
            else:
                layout.prop(simorbit, "user_period_seconds")


# Panel for self rotation
class SSSIMPanelRotation(bpy.types.Panel):
    bl_label = "Solar System Simulator: Rotation"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "physics"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        sim_obj = context.object.sssim_obj
        return sim_obj.use_sssim and sim_obj.object_type != "SURFACE"

    def draw_header(self, context):
        sim_rot = context.object.sssim_rotation
        self.layout.prop(sim_rot, "use_rotation", text="")

    def draw(self, context):
        layout = self.layout
        obj = context.object
        sim_rot = obj.sssim_rotation
        info = obj.sssim_obj.show_info
        layout.active = sim_rot.use_rotation

        if has_rotation_fcurve(obj):
            # Deactivate layout if rotation fcurve found
            layout.label(text="Found Rotation F-Curve", icon='INFO')
            layout.active = False
        elif not obj.sssim_calc.use_driver:
            # No rotation because driver not active
            layout.label(text="Rotation Driver not Active", icon='ERROR')

        if info:
            layout.label(text="Rotation:", icon='PREVIEW_RANGE')

            rot_sec = round(sim_rot.rotation_period, 2)
            rot_hours = round(sim_rot.rotation_period / 3600, 2)
            rot_frames = round(sim_rot.rotation_period_frames, 2)

            # Don't show frames again when the user uses frames
            if sim_rot.use_frames:
                rot_text = "{0:,} s = {1} h".format(rot_sec, rot_hours)
            else:
                rot_text = "{0:,} s = {1} h (= {2} frames)"
                rot_text = rot_text.format(rot_sec, rot_hours, rot_frames)
            layout.label(text=rot_text)

        layout.prop(sim_rot, "use_frames")
        if sim_rot.use_frames:
            layout.prop(sim_rot, "user_period_frames")
        else:
            layout.prop(sim_rot, "user_period_seconds")

        row = layout.row()
        row.label(text="Rotation Axis Direction:")
        sub = row.column()
        sub.active = obj.sssim_orbit.center_object is not None
        sub.prop(sim_rot, "relative_to_orbit")

        row = layout.row(align=True)
        row.prop(sim_rot, "axis_tilt")
        row.prop(sim_rot, "axis_direction")


# Panel for Planet surfaces
class SSSIMPanelSurface(bpy.types.Panel):
    bl_label = "Solar System Simulator: Surface"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "physics"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        simobj = context.object.sssim_obj
        return simobj.use_sssim and simobj.object_type == 'SURFACE'

    def draw(self, context):
        layout = self.layout
        obj = context.object
        info = obj.sssim_obj.show_info
        simsur = obj.sssim_surface

        row = layout.row()
        row.prop(simsur, "calc_size")
        if info:
            row.label(text="Radius in BU: {:.5}".format(simsur.radius_bu))

        if simsur.calc_size:
            # show radius controls
            row = layout.row(align=True)
            row.prop(simsur, "radius_mantissa")
            row.prop(simsur, "radius_exp", slider=True)
        else:
            # show dimension controls of surface object
            layout.label(text="Dimensions of {}:".format(obj.name))
            layout.prop(obj, "dimensions", text="")


# Panel for Calculation settings
class SSSIMPanelCalculation(bpy.types.Panel):
    bl_label = "Solar System Simulator: Calculation"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "physics"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        sim_obj = context.object.sssim_obj
        return sim_obj.use_sssim and sim_obj.object_type != "SURFACE"

    def draw(self, context):
        layout = self.layout
        obj = context.object
        simcalc = obj.sssim_calc

        layout.prop(simcalc, "use_driver", icon="DRIVER")

        col = layout.column()
        col.enabled = not simcalc.use_driver
        col.label(text="Create Simulation F-Curve:")

        split = col.split()
        col = split.column()
        sub = col.column(align=True)
        sub.operator("object.sssim_to_fcurve")
        sub.operator("object.sssim_clear_fcurve")

        if simcalc.calc_orbit:
            col.label(text="Calc Orbit: Yes")
        else:
            col.label(text="Calc Orbit: No")

        if simcalc.calc_rotation:
            col.label(text="Calc Rotaton: Yes")
        else:
            col.label(text="Calc Rotaton: No")

        col = split.column()
        sub = col.row()
        sub.active = simcalc.calc_orbit
        sub.prop(simcalc, "frame_step")

        col.prop(simcalc, "cyclic")
        sub = col.column(align=True)
        if not simcalc.cyclic:
            sub.prop(simcalc, "frame_start")
            sub.prop(simcalc, "frame_end")
        else:
            sub.label(text="Start Frame: 1")
            sub.label(text="End Frame:   {}".format(simcalc.real_frame_end))

        keys = simcalc.calc_needed_keys
        col.label(text="Needed Keyframes: {}".format(keys))


# Panel for scene properties in scene tab
class SSSIMPanelScene(bpy.types.Panel):
    bl_label = "Solar System Simulator: Scene Data"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "scene"

    def draw(self, context):
        layout = self.layout
        simscn = context.scene.sssim_scn

        add_fcurve = bpy.ops.scene.add_eval_time_fcurve.poll()
        if add_fcurve:
            err_text = "You need to add an Evaluation Time F-Curve"
            layout.label(err_text, icon='ERROR')

        row = layout.row(align=True)
        if add_fcurve:
            row.operator("scene.add_eval_time_fcurve", icon='FCURVE')
        row.prop(simscn, "eval_time")

        time = round(simscn.time)
        time_text = "Time: {0:,} sec".format(time)
        layout.label(time_text, icon='TIME')

        layout.separator()

        time_scale_text = "Time Scaling: 1 sec =  {0:,} sec"
        time_mult = round(simscn.time_mult, 2)
        layout.label(time_scale_text.format(time_mult))
        layout.prop(simscn, "time_exp", slider=True)

        layout.separator()

        len_scale_text = "Length Scaling: 1 Blender Unit = {0:,} km"
        length_mult = round(simscn.length_mult, 2)
        layout.label(len_scale_text.format(length_mult))
        layout.prop(simscn, "length_exp", slider=True)
        layout.prop(simscn, "planet_size_mult")

        layout.separator()

        row = layout.row(align=True)
        row.operator("scene.update_sssim_drivers", icon='FILE_REFRESH')
        row.prop(simscn, "draw_orbit", toggle=True, icon='FORCE_CURVE')


# Panel for scene properties in scene tab
class SSSIMPanelTools(bpy.types.Panel):
    bl_label = "Solar System Tools:"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_context = "objectmode"
    bl_category = "Tools"

    def draw(self, context):
        layout = self.layout
        col = layout.column(align=True)

        if bpy.ops.scene.add_eval_time_fcurve.poll():
            col.operator("scene.add_eval_time_fcurve", icon='FCURVE')

        col.operator("object.sssim_create_center", icon='LAMP_SUN')
        col.operator("object.sssim_create_planet", icon='WORLD')
        col.operator("object.sssim_create_surface", icon='SURFACE_NSPHERE')

        col = layout.column(align=True)
        row = col.row(align=True)
        row.operator("scene.sssim_bake_all", text="Bake All", icon='IPO')
        row.operator("scene.sssim_bake_clear", text="Clear All", icon='DRIVER')
        col.operator("scene.update_sssim_drivers", icon='FILE_REFRESH')

        simscn = context.scene.sssim_scn
        layout.prop(simscn, "draw_orbit", toggle=True, icon='FORCE_CURVE')


# adjust the physics-tab panel
def physics_panel(self, context):
    self.layout.prop(context.object.sssim_obj, "use_sssim", icon='LAMP_SUN')

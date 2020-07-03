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

from math import ceil, degrees

import bpy

from .operators import has_location_fcurve, has_rotation_fcurve


# =============================================================================
# UI panel in 3D Viewport
# =============================================================================

# panel in 3D Viewport, tools tab
class SSSIM_PT_tools(bpy.types.Panel):
    bl_label = "Solar System Simulator: Tools"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_context = "objectmode"
    bl_category = "Tool"

    def draw(self, context):
        layout = self.layout

        if "eval_planet_orbit" not in bpy.app.driver_namespace:
            box = layout.box()
            box.label(text="Simulation driver not in driver namespace",
                      icon='ERROR')
            box.operator("scene.update_sssim_drivers", icon='FILE_REFRESH')

        add_fcurve = bpy.ops.scene.add_sim_time_fcurve.poll()
        if add_fcurve:
            box = layout.box()
            err_text = "You need to add a simulation time F-Curve"
            box.label(text=err_text, icon='ERROR')
            box.operator("scene.add_sim_time_fcurve", icon='FCURVE')

        col = layout.column(align=True)
        col.operator("object.sssim_create_center", icon='LIGHT_SUN')
        col.operator("object.sssim_create_planet", icon='WORLD')
        col.operator("object.sssim_create_surface", icon='SHADING_WIRE')

        row = layout.row(align=True)
        row.operator("scene.sssim_bake_all", text="Bake All", icon='GRAPH')
        row.operator("scene.sssim_bake_clear", text="Clear All", icon='DRIVER')

        simscn = context.scene.sssim_scn
        layout.prop(simscn, "draw_orbit", toggle=True, icon='FORCE_CURVE')

        layout.operator("scene.update_sssim_drivers", icon='FILE_REFRESH')


# =============================================================================
# UI panel in scene tab
# =============================================================================

# panel for scene properties in scene tab
class SSSIM_PT_scene(bpy.types.Panel):
    bl_label = "Solar System Simulator: Scene Data"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "scene"

    def draw(self, context):
        layout = self.layout
        simscn = context.scene.sssim_scn

        if "eval_planet_orbit" not in bpy.app.driver_namespace:
            box = layout.box()
            box.label(text="Simulation driver not in driver namespace",
                      icon='ERROR')
            box.operator("scene.update_sssim_drivers", icon='FILE_REFRESH')

        add_fcurve = bpy.ops.scene.add_sim_time_fcurve.poll()
        if add_fcurve:
            box = layout.box()
            err_text = "You need to add a simulation time F-Curve"
            box.label(text=err_text, icon='ERROR')
            box.operator("scene.add_sim_time_fcurve", icon='FCURVE')

        time_text = "Time: {:.2f} s"
        layout.label(text=time_text.format(simscn.time), icon='TIME')

        row = layout.row(align=True)
        row.prop(simscn, "sim_time")
        row.prop(simscn, "time_exp", slider=True)

        time_scale_text = "Time Conversion: 1 Blender Second = {0:.0f} s"
        layout.label(text=time_scale_text.format(simscn.time_mult))

        layout.separator()

        layout.prop(simscn, "length_exp", slider=True)
        len_scale_text = "Length Conversion: 1 Blender Unit = {0:.0f} km"
        layout.label(text=len_scale_text.format(simscn.length_mult))

        layout.prop(simscn, "planet_size_mult")

        layout.separator()

        col = layout.column(align=True)
        col.prop(simscn, "draw_orbit", toggle=True, icon='FORCE_CURVE')
        col.operator("scene.update_sssim_drivers", icon='FILE_REFRESH')


# =============================================================================
# UI panel in physics tab
# =============================================================================

# add solar system simulator to the physics-tab panel
def physics_panel(self, context):
    self.layout.prop(context.object.sssim_obj, "use_sssim", icon='LIGHT_SUN')


# panel for general object settings
class SSSIM_PT_object(bpy.types.Panel):
    bl_label = "Solar System Simulator"
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

        if bpy.ops.scene.add_sim_time_fcurve.poll():
            box = layout.box()
            box.label(text="You need to add an simulation time F-Curve",
                      icon='ERROR')
            box.operator("scene.add_sim_time_fcurve", icon='FCURVE')
            box.label(text="You can find more settings in the scene Tab",
                      icon='SCENE_DATA')

        if "eval_planet_orbit" not in bpy.app.driver_namespace:
            box = layout.box()
            box.label(text="Simulation driver not in driver namespace",
                      icon='ERROR')
            box.operator("scene.update_sssim_drivers", icon='FILE_REFRESH')

        layout.active = simobj.use_sssim
        layout.prop(simobj, "show_info", icon='INFO')
        layout.prop(simobj, "object_type", expand=True)

        # planet only props
        if simobj.object_type == 'PLANET':
            simorbit = obj.sssim_orbit
            layout.prop_search(simorbit, "center_name",
                               context.scene, "objects")
            if not simorbit.center_is_valid:
                layout.label(text="No valid center", icon='ERROR')

        # props for center and planet but not surface
        if simobj.object_type in ('CENTER', 'PLANET'):
            if info:
                mantissa = round(simobj.mass_mantissa, 3)
                exp = simobj.mass_exp
                layout.label(text="Mass: {} * 10^{} kg".format(mantissa, exp))

            row = layout.row(align=True)
            row.prop(simobj, "mass_mantissa")
            row.prop(simobj, "mass_exp", slider=True)
        elif simobj.object_type == 'SURFACE':
            simsur = obj.sssim_surface
            layout.prop_search(simsur, "parent_name", context.scene, "objects")
            if not simsur.parent_is_valid:
                layout.label(text="Parent object is not valid", icon="ERROR")


# panel for orbit settings
class SSSIM_PT_orbit(bpy.types.Panel):
    bl_label = "Orbit"
    bl_parent_id = "SSSIM_PT_object"
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
        return False

    def draw(self, context):
        layout = self.layout
        obj = context.object
        simorbit = obj.sssim_orbit
        info = obj.sssim_obj.show_info

        if has_location_fcurve(obj):
            # deactivate layout if location fcurve found
            layout.label(text="Found Orbit F-curve", icon='LOCKED')
            return

        if not obj.sssim_calc.use_driver:
            # there is no orbit because driver not active
            layout.label(text="Orbit Driver not active", icon='ERROR')

        # shape and size
        simscn = context.scene.sssim_scn
        self.draw_shape(simorbit, simscn, layout, info)

        layout.separator()

        # orientation
        self.draw_orientation(simorbit, layout)

        # position
        self.draw_position(simorbit, layout, info)

        layout.separator()

        # orbital period
        self.draw_orbital_period(simorbit, layout, info)

    def draw_shape(self, simorbit, simscn, layout, info):
        layout.label(text="Shape and Size of Orbit:")

        layout.prop(simorbit, "eccentricity", slider=True)

        if info:
            axis_mantissa = round(simorbit.distance_mantissa, 3)
            axis_exp = simorbit.distance_exp
            # lengths in km / length_mult = lengths in Blender Units
            axis_bu = round(simorbit.semi_major_axis / simscn.length_mult, 2)
            data = "{} * 10^{} km ( = {:,} units)"
            data = data.format(axis_mantissa, axis_exp, axis_bu)
            if simorbit.eccentricity == 0:
                msg = "Radius: {}".format(data)
            else:
                msg = "Semi-Major Axis: {}".format(data)
            layout.label(text=msg, icon="DRIVER_DISTANCE")

        row = layout.row(align=True)
        if simorbit.eccentricity == 0:
            row.prop(simorbit, "distance_mantissa", text="Radius (km)")
        else:
            row.prop(simorbit, "distance_mantissa")
        row.prop(simorbit, "distance_exp", slider=True)

    def draw_orientation(self, simorbit, layout):
        layout.label(text="Orientation of Orbital Plane:")
        subcol = layout.column(align=True)
        subcol.prop(simorbit, "inclination")
        subcol.prop(simorbit, "asc_node")
        subcol.prop(simorbit, "arg_periapsis")

    def draw_position(self, simorbit, layout, info):
        layout.label(text="Position:")

        if info:
            # True Anomaly (angle travelled)
            angle = round(degrees(simorbit.true_anomaly), 2)
            layout.label(text="True Anomaly: {}°".format(angle),
                         icon="HANDLE_AUTO")

        layout.prop(simorbit, "time_offset", slider=True)

    def draw_orbital_period(self, simorbit, layout, info):
        row = layout.row()
        if simorbit.override_period:
            row.label(text="Orbital Period:")
        else:
            row.label(text="Orbital Period: Gravity")

        row.prop(simorbit, "override_period")

        if info:
            orb_sec = round(simorbit.orbital_period)
            # orb_years = orbital period in julian years (365.25 days)
            orb_years = round(simorbit.orbital_period / 86400 / 365.25, 3)
            orb_frames = round(simorbit.orbital_period_frames, 2)

            if simorbit.override_period and simorbit.use_frames:
                # don't show frames again when the user uses frames
                orb_text = "{0:,} s = {1} y".format(orb_sec, orb_years)
            else:
                orb_text = "{0:,} s = {1} y (= {2} frames)"
                orb_text = orb_text.format(orb_sec, orb_years, orb_frames)
            layout.label(text=orb_text, icon='RECOVER_LAST')

        if simorbit.override_period:
            layout.prop(simorbit, "use_frames")

            if simorbit.use_frames:
                layout.prop(simorbit, "override_period_frames")
            else:
                layout.prop(simorbit, "override_period_seconds")


# panel for self rotation
class SSSIM_PT_rotation(bpy.types.Panel):
    bl_label = "Rotation"
    bl_parent_id = "SSSIM_PT_object"
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
            # deactivate layout if rotation fcurve found
            layout.label(text="Found Rotation F-Curve", icon='LOCKED')
            return

        if not obj.sssim_calc.use_driver:
            # no rotation because driver not active
            layout.label(text="Rotation Driver not Active", icon='ERROR')

        layout.label(text="Rotation Period:")
        if info:
            rot_sec = round(sim_rot.rotation_period, 2)
            rot_hours = round(sim_rot.rotation_period / 3600, 2)
            rot_frames = round(sim_rot.rotation_period_frames, 2)

            # don't show frames again when the user uses frames
            if sim_rot.use_frames:
                rot_text = "{0:,} s = {1} h".format(rot_sec, rot_hours)
            else:
                rot_text = "{0:,} s = {1} h (= {2} frames)"
                rot_text = rot_text.format(rot_sec, rot_hours, rot_frames)
            layout.label(text=rot_text, icon='ORIENTATION_GIMBAL')

        layout.prop(sim_rot, "use_frames")
        if sim_rot.use_frames:
            layout.prop(sim_rot, "period_frames")
        else:
            layout.prop(sim_rot, "period_seconds")

        row = layout.row()
        row.label(text="Rotation Axis Direction:")
        sub = row.column()
        sub.active = obj.sssim_orbit.center_object is not None
        sub.prop(sim_rot, "relative_to_orbit")

        row = layout.row(align=True)
        row.prop(sim_rot, "axis_tilt")
        row.prop(sim_rot, "axis_direction")


# panel for planet surfaces
class SSSIM_PT_surface(bpy.types.Panel):
    bl_label = "Surface"
    bl_parent_id = "SSSIM_PT_object"
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
            row.label(text="Radius: {:.5} units".format(simsur.radius_bu))

        if simsur.calc_size:
            # show radius controls
            row = layout.row(align=True)
            row.prop(simsur, "radius_mantissa")
            row.prop(simsur, "radius_exp", slider=True)
        else:
            # show dimension controls of surface object
            layout.label(text="Dimensions of {}:".format(obj.name))
            layout.prop(obj, "dimensions", text="")


# panel for calculation settings
class SSSIM_PT_calculation(bpy.types.Panel):
    bl_label = "Calculation"
    bl_parent_id = "SSSIM_PT_object"
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
        info = obj.sssim_obj.show_info

        layout.prop(simcalc, "use_driver", icon="DRIVER")

        col = layout.column()
        col.enabled = not simcalc.use_driver
        col.label(text="Create Simulation F-Curve:")

        row = col.row(align=True)
        row.operator("object.sssim_to_fcurve")
        row.operator("object.sssim_clear_fcurve")

        if info:
            row = col.row()
            calc_orbit = "Yes" if simcalc.calc_orbit else "No"
            row.label(text="Calculate Orbit: {}".format(calc_orbit))

            calc_rotation = "Yes" if simcalc.calc_rotation else "No"
            row.label(text="Calculate Rotaton: {}".format(calc_rotation))

        col.prop(simcalc, "cyclic")
        sub = col.column(align=True)
        if not simcalc.cyclic:
            sub.prop(simcalc, "frame_start")
            sub.prop(simcalc, "frame_end")
        else:
            sub.label(text="Start Frame: {}".format(simcalc.real_frame_start))
            end_frame = round(simcalc.real_frame_end, 1)
            sub.label(text="End Frame:   {}".format(end_frame))
        sub.prop(simcalc, "frame_step")

        if info:
            start = simcalc.real_frame_start
            end = simcalc.real_frame_end
            step = simcalc.frame_step
            num_keys = ceil((end - start) / step) + 1
            col.label(text="Needed Keyframes: {}".format(num_keys))

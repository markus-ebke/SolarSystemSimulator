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

from itertools import chain
from math import sqrt

import bpy
from bpy.props import (
    BoolProperty,
    IntProperty,
    FloatProperty,
    StringProperty,
    EnumProperty)

from .calculation import eval_planet_orbit, eval_planet_rotation


# =============================================================================
# Add/remove relations between objects
# =============================================================================

def add_driver_loc(obj, scn):
    """Add or adjust the location-drivers for orbit"""
    loc_list = obj.driver_add("location")  # X-, Y- and Z-Driver in list

    for driver_fcurve in loc_list:
        index = driver_fcurve.array_index
        driver = driver_fcurve.driver

        # clear existing variables, if any
        for var in driver.variables.values():
            driver.variables.remove(var)

        # new variable: scene data
        var = driver.variables.new()
        var.name = "scn_name"
        var.targets[0].id_type = 'SCENE'
        var.targets[0].id = scn
        var.targets[0].data_path = "name"

        # setup scripted expression
        driver.type = 'SCRIPTED'
        driver.use_self = True
        expr = "eval_planet_orbit(self, scn_name, {})".format(index)
        driver.expression = expr


def add_driver_rot(obj, scn):
    """Add or adjust the rotation-drivers"""
    # Have to use a different mode because of rotation around Z
    obj.rotation_mode = 'ZYX'
    # X-, Y- and Z-driver in list
    rot_list = obj.driver_add("rotation_euler")

    # rotation driver:
    for driver_fcurve in rot_list:
        index = driver_fcurve.array_index
        driver = driver_fcurve.driver

        # clear existing variables, if any
        for var in driver.variables.values():
            driver.variables.remove(var)

        # new variable: scene data
        var = driver.variables.new()
        var.name = "scn_name"
        var.targets[0].id_type = 'SCENE'
        var.targets[0].id = scn
        var.targets[0].data_path = "name"

        # setup scripted expression
        driver.type = 'SCRIPTED'
        driver.use_self = True
        expr = "eval_planet_rotation(self, scn_name, {})".format(index)
        driver.expression = expr


def add_orbit_constraint(child_obj, parent_obj):
    """Add/adjust child-parent contraint for orbit"""
    orbitcon = child_obj.constraints.get("ORBIT")

    if orbitcon is None:
        orbitcon = child_obj.constraints.new('CHILD_OF')
        orbitcon.name = "ORBIT"

        # don't rotate with the center
        orbitcon.use_rotation_x = False
        orbitcon.use_rotation_y = False
        orbitcon.use_rotation_z = False

        # don't scale with the center
        orbitcon.use_scale_x = False
        orbitcon.use_scale_y = False
        orbitcon.use_scale_z = False

        msg = "Constraint {} to {}".format(child_obj.name, parent_obj.name)
    else:
        msg = "Constraint of {} adjusted".format(child_obj.name)

    orbitcon.target = parent_obj
    return msg


def add_surface(child_obj, parent_obj):
    """Add/adjust parenting for surface child_obj and parent parent_obj"""
    if parent_obj == child_obj:
        msg = "Surface not valid: child: {}, parent: {}"
        return msg.format(child_obj.name, parent_obj.name)
    child_obj.parent = parent_obj
    msg = "Surface: {} surface of {}"
    return msg.format(child_obj.name, parent_obj.name)


def remove_driver_loc(obj):
    """Remove the location-driver (= orbit)"""
    return obj.driver_remove("location")


def remove_driver_rot(obj):
    """Remove the rotation-driver"""
    return obj.driver_remove("rotation_euler")


def remove_orbit_constraint(obj):
    """Remove the child-parent constraint"""
    orbitcon = obj.constraints.get("ORBIT")
    if orbitcon:
        obj.constraints.remove(orbitcon)
        msg = "Orbit Constraint removed from {}".format(obj.name)
    else:
        msg = "No Orbit Constraint to remove from {}".format(obj.name)

    return msg


def remove_surface(child_obj):
    """Remove the child-parent constraint for surface"""
    parent_obj = child_obj.parent
    child_obj.parent = None

    if parent_obj:
        pname = parent_obj.name
        msg = "Surface {} removed from {}".format(child_obj.name, pname)
    else:
        msg = "Surface {} had no parent".format(child_obj.name)

    return msg


# =============================================================================
# F-Curve stuff
# =============================================================================

def get_fcurve(fcurves, search_data_path, array_index=-1):
    """Find in list fcurves the F-Curve with given data_path (and index)"""
    for fc in fcurves:
        if fc.data_path == search_data_path:
            if array_index in (-1, fc.array_index):
                return fc
    return None


def has_location_fcurve(obj):
    """Find out if the object has a location (=orbit) F-Curve"""
    if obj.animation_data:
        if obj.animation_data.action:
            fcurves = obj.animation_data.action.fcurves
            if fcurves:
                if get_fcurve(fcurves, "location"):
                    return True
    return False


def has_rotation_fcurve(obj):
    """Find out if the object has a rotation F-Curve"""
    if obj.animation_data:
        if obj.animation_data.action:
            fcurves = obj.animation_data.action.fcurves
            if fcurves:
                if get_fcurve(fcurves, "rotation_euler"):
                    return True
    return False


def key_insert(fcurve, start, end, step, get_value):
    """Insert keys in given range into the F-Curve

    fcurve = F-Curve where the keys are added
    start, end, step = where to insert the keys (along x-axis)
    get_value(frame), function which returns the y-value at the given frame
    returns number of inserted keyframes
    """
    frames = chain(range(start, int(end), step), [end])

    loops = 0
    for frame in frames:
        val = get_value(frame)
        fcurve.keyframe_points.insert(frame, val)

        loops += 1

    return loops


# =============================================================================
# Operators
# =============================================================================

class SCENE_OT_add_sim_time_fcurve(bpy.types.Operator):
    """Add a linear F-Curve to simulation time"""
    bl_idname = "scene.add_sim_time_fcurve"
    bl_label = "Add time F-Curve"

    @classmethod
    def poll(cls, context):
        # if sim_time is controlled by fcurve -> return False
        anim_data = context.scene.animation_data
        if anim_data and anim_data.action:
            if get_fcurve(anim_data.action.fcurves, "sssim_scn.sim_time"):
                return False
        return True

    def execute(self, context):
        scn = context.scene

        if scn.animation_data is None:
            scn.animation_data_create()

        anim_data = scn.animation_data
        act = bpy.data.actions.new("Simulation_Time")
        anim_data.action = act
        time_curve = act.fcurves.new(data_path="sssim_scn.sim_time")

        fmod = time_curve.modifiers.new("GENERATOR")
        fmod.mode = "POLYNOMIAL"
        fmod.poly_order = 1
        slope = scn.render.fps_base / scn.render.fps
        fmod.coefficients = (-slope, slope)  # first frame: time = 0
        return {'FINISHED'}


# create and use F-Curves for the simualation instead of drivers
class OBJECT_OT_sssim_to_fcurve(bpy.types.Operator):
    """Bake the simulation to F-Curves for location and/or rotation"""
    bl_idname = "object.sssim_to_fcurve"
    bl_label = "SSSim to F-Curve"

    @classmethod
    def poll(cls, context):
        obj = context.object
        if obj:
            simcalc = obj.sssim_calc
            # can we bake something?
            return simcalc.calc_orbit or simcalc.calc_rotation
        return False

    def execute(self, context):
        scn = context.scene
        obj = context.object

        # there must be a sim_time F-Curve
        if bpy.ops.scene.add_sim_time_fcurve.poll():
            msg = "No simulation time F-Curve found, check the scene tab"
            self.report({'ERROR'}, msg)
            return {'CANCELLED'}

        scn_act = scn.animation_data.action
        sim_time_fcurve = get_fcurve(scn_act.fcurves, "sssim_scn.sim_time")
        time_multiply = scn.sssim_scn.time_mult

        # set up animation data if needed
        if obj.animation_data is None:
            obj.animation_data_create()
        anim_data = obj.animation_data

        # use existing action to keep everything clean
        previous = bpy.data.actions.get("SSSim_%s_Action" % obj.name)
        if previous:
            act = previous
        else:
            act = bpy.data.actions.new(name="SSSim_%s_Action" % obj.name)

        anim_data.action = act

        # time_func(frame) returns the time at this frame
        def time_func(frame):
            return sim_time_fcurve.evaluate(frame) * time_multiply

        simcalc = obj.sssim_calc
        if simcalc.calc_orbit:
            self.orbit_to_fcurve(act, scn, obj, time_func)
        if simcalc.calc_rotation:
            self.rotation_to_fcurve(act, scn, obj, time_func)

        return {'FINISHED'}

    def orbit_to_fcurve(self, act, scn, obj, get_time):
        simcalc = obj.sssim_calc
        start = simcalc.real_frame_start
        end = simcalc.real_frame_end
        step = simcalc.frame_step

        for index in range(0, 3):
            # value_func(frame) returns the value at this frame
            def value_func(frame):
                return eval_planet_orbit(obj, scn.name, index, get_time(frame))

            # remove existing F-Curve if any
            fc = get_fcurve(act.fcurves, "location", index)
            if fc:
                act.fcurves.remove(fc)
            fc = act.fcurves.new(data_path="location", index=index)

            n = key_insert(fc, start, end, step, value_func)

            if simcalc.cyclic:
                fc.modifiers.new(type='CYCLES')

        print("Created {:3} keyframes for location".format(n))

    def rotation_to_fcurve(self, act, scn, obj, get_time):
        def value_func(time):
            return eval_planet_rotation(obj, scn.name, 2, get_time(time))

        # calculate rotation around z-axis only, index = 2
        # rotation of x- and y-axis stay fixed
        rot = eval_planet_rotation(obj, scn.name, index=None, time=0)
        rot_x, rot_y = rot[0], rot[1]
        obj.rotation_euler.x = rot_x
        obj.rotation_euler.y = rot_y

        # remove existing z-axis curve if any
        fc = get_fcurve(act.fcurves, "rotation_euler", 2)
        if fc:
            act.fcurves.remove(fc)

        fc = act.fcurves.new(data_path="rotation_euler", index=2)
        fc.extrapolation = 'LINEAR'

        simcalc = obj.sssim_calc
        start = simcalc.real_frame_start
        end = simcalc.real_frame_end
        step = simcalc.frame_step

        if simcalc.cyclic:
            # we need just two keyframes
            end = 1 + step

        n = key_insert(fc, start, end, step, value_func)

        print("Created {:3} keyframes for rotation".format(n))


class OBJECT_OT_sssim_clear_fcurve(bpy.types.Operator):
    """Clear location and rotation F-Curves"""
    bl_idname = "object.sssim_clear_fcurve"
    bl_label = "Clear SSSim F-Curve"

    @classmethod
    def poll(cls, context):
        # search for location and rotation fcurve
        obj = context.object
        if obj:
            return has_location_fcurve(obj) or has_rotation_fcurve(obj)
        return False

    def execute(self, context):
        fcurves = context.object.animation_data.action.fcurves
        for i in range(3):
            loc_curve = get_fcurve(fcurves, "location", i)
            rot_curve = get_fcurve(fcurves, "rotation_euler", i)
            if loc_curve:
                fcurves.remove(loc_curve)
            if rot_curve:
                fcurves.remove(rot_curve)

        return {'FINISHED'}


class SCENE_OT_update_sssim_drivers(bpy.types.Operator):
    """Update the simulation drivers, neccessary after file reload"""
    bl_idname = "scene.update_sssim_drivers"
    bl_label = "Update Drivers"

    def execute(self, context):
        # update driver namespace
        bpy.app.driver_namespace["eval_planet_orbit"] = eval_planet_orbit
        bpy.app.driver_namespace["eval_planet_rotation"] = eval_planet_rotation

        # update all objects in current scene
        scn = context.scene
        number_updates = 0
        for obj in scn.objects:
            simobj = obj.sssim_obj
            simcalc = obj.sssim_calc
            if simobj.use_sssim and simcalc.use_driver:
                if obj.sssim_orbit.center_is_valid:
                    # update location driver
                    add_driver_loc(obj, scn)
                if obj.sssim_rotation.use_rotation:
                    # update rotation driver
                    add_driver_rot(obj, scn)
                number_updates += 1

        number_objects = len(scn.objects)
        msg = "Updated {} of {} objects".format(number_updates, number_objects)
        self.report({'INFO'}, msg)

        return {'FINISHED'}


class OBJECT_OT_add_orbit_curve(bpy.types.Operator):
    """Create Bezier curves that match the orbits of the selected planets"""
    bl_idname = "object.sssim_add_orbit_curve"
    bl_label = "Add Orbit Curves"
    bl_options = {'REGISTER', 'UNDO'}

    bevel_depth: FloatProperty(
        name="Bevel Depth",
        description="Radius of the Bevel Geometry",
        min=0.0)

    @classmethod
    def poll(cls, context):
        return bool(context.selected_objects)

    def execute(self, context):
        simscn = context.scene.sssim_scn
        for obj in context.selected_objects:
            if not obj.sssim_obj.use_sssim:
                continue

            # can only add orbits to planets
            if not obj.sssim_obj.object_type == 'PLANET':
                continue
            simorbit = obj.sssim_orbit

            # planet without a center has no orbit
            if not simorbit.center_object:
                continue

            # create a new Bezier circle curve
            bpy.ops.curve.primitive_bezier_circle_add()
            curve = context.active_object
            curve.name = "{}_orbit".format(obj.name)
            curve.data.bevel_depth = self.bevel_depth

            # switch to edit mode, set handle type so that scaling is proportional
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.curve.handle_type_set(type='ALIGNED')

            # resize circle to ellipse and translate focal point to origin
            ecc = simorbit.eccentricity
            scale = simorbit.semi_major_axis / simscn.length_mult
            bpy.ops.transform.resize(value=(scale, sqrt(1 - ecc**2) * scale, scale))
            bpy.ops.transform.translate(value=(-ecc * scale, 0, 0))

            # switch back to object mode
            bpy.ops.object.mode_set(mode='OBJECT')

            # perform rotations to match orientation of the orbit
            bpy.ops.transform.rotate(value=simorbit.inclination, orient_axis='X',
                                     orient_type='GLOBAL',
                                     constraint_axis=(True, False, False))
            bpy.ops.transform.rotate(value=simorbit.asc_node, orient_axis='Z',
                                     orient_type='GLOBAL',
                                     constraint_axis=(False, False, True))
            bpy.ops.transform.rotate(value=simorbit.arg_periapsis, orient_axis='Z',
                                     orient_type='LOCAL',
                                     constraint_axis=(False, False, True))

            # parent orbit curve to the center object (i.e. copy location)
            add_orbit_constraint(curve, simorbit.center_object)

        return {'FINISHED'}

    def invoke(self, context, event):
        wm = context.window_manager
        return wm.invoke_props_dialog(self, width=200)

    def draw(self, context):
        layout = self.layout

        layout.prop(self, "bevel_depth")


# =============================================================================
# Operators in 3D Viewport toolbar
# =============================================================================

class OBJECT_OT_sssim_create_center(bpy.types.Operator):
    bl_idname = "object.sssim_create_center"
    bl_label = "Create Center"
    bl_description = "Create a central object for the Solar System"

    center_name: StringProperty(
        name="Center Name",
        description="Name of the center",
        default="Center")

    mass_mantissa: FloatProperty(
        name="Mass (kg)",
        description="Mass of the center (in kilogram)",
        soft_min=0.001, soft_max=1000, default=1)

    mass_exp: IntProperty(
        name="Mass Exponent",
        description="Base-10 exponent of the mass",
        soft_min=0, soft_max=50, default=27)

    # Rotation
    use_rotation: BoolProperty(
        name="Use Rotation",
        description="Make the center rotate around its axis",
        default=True)

    rotation_period: FloatProperty(
        name="Rotation Period (in seconds)",
        description="Time for one to rotation in seconds",
        soft_min=1, default=86400)

    # Surfaces
    with_surface: BoolProperty(
        name="Create Surface")

    mesh_type: EnumProperty(
        name="Mesh Type",
        items=(('UV_SPHERE', "UV Sphere", "UV sphere surface"),
               ('ICOSPHERE', "Icosphere", "Icosphere surface"),
               ),
        description=("Type of mesh, UV spheres are better for unwrapping, "
                     "Icospheres have an even distribution of vertices"),
        default='ICOSPHERE')

    subdivisions: IntProperty(
        name="Subdivisions",
        description="1 = coarse mesh, 7 = dense but round mesh",
        min=1, max=7, default=3)

    radius_mantissa: FloatProperty(
        name="Radius (in km)",
        description="The radius of the surface",
        soft_min=0.001, soft_max=1000, default=0.5)

    radius_exp: IntProperty(
        name="Radius Exponent",
        description="Base-10 exponent of the radius",
        soft_min=0, soft_max=8, default=6)

    def execute(self, context):
        # create emtpy
        center = bpy.data.objects.new(self.center_name, None)
        context.collection.objects.link(center)

        # set up center properties
        simobj = center.sssim_obj
        simobj.use_sssim = True
        simobj.object_type = 'CENTER'
        simobj.mass_mantissa = self.mass_mantissa
        simobj.mass_exp = self.mass_exp

        # set up rotation of the center
        simrot = center.sssim_rotation
        simrot.use_rotation = self.use_rotation
        simrot.use_frames = False
        simrot.period_seconds = self.rotation_period

        if self.with_surface:
            # create a surface for the center
            bpy.ops.object.sssim_create_surface(
                parent_name=center.name,
                mesh_type=self.mesh_type,
                subdivisions=self.subdivisions,
                radius_mantissa=self.radius_mantissa,
                radius_exp=self.radius_exp)

        # make the center object the active object
        center.select_set(True)
        context.view_layer.objects.active = center

        return {'FINISHED'}

    def invoke(self, context, event):
        wm = context.window_manager
        return wm.invoke_props_dialog(self, width=300)

    def draw(self, context):
        layout = self.layout
        layout.prop(self, "center_name")

        row = layout.row(align=True)
        row.prop(self, "mass_mantissa")
        row.prop(self, "mass_exp", slider=True)

        layout.separator()

        layout.prop(self, "use_rotation")
        layout.prop(self, "rotation_period")

        layout.separator()

        layout.prop(self, "with_surface")
        layout.prop(self, "mesh_type")
        col = layout.column(align=True)
        col.prop(self, "subdivisions")
        sub = col.row(align=True)
        sub.prop(self, "radius_mantissa")
        sub.prop(self, "radius_exp", slider=True)


class OBJECT_OT_sssim_create_planet(bpy.types.Operator):
    """Create a planet orbiting the active planet or center"""
    bl_idname = "object.sssim_create_planet"
    bl_label = "Create Planet"

    center_name: StringProperty(
        name="Center",
        description="Name of the center object")

    planet_name: StringProperty(
        name="Planet Name",
        description="Name of the new planet",
        default="Planet")

    mass_mantissa: FloatProperty(
        name="Mass (kg)",
        description="Mass of the planet in kilogram",
        soft_min=0.001, soft_max=1000, default=1)

    mass_exp: IntProperty(
        name="Mass Exponent",
        description="Base-10 exponent of the mass",
        soft_min=0, soft_max=50, default=24)

    # Orbit
    distance_mantissa: FloatProperty(
        name="Distance (in km)",
        description="Distance to the center object in kilometers",
        soft_min=0, soft_max=1000, default=10)

    distance_exp: IntProperty(
        name="Length Exponent",
        description="Base-10 exponent of the distance",
        soft_min=0, soft_max=12, default=6)

    # Rotation
    use_rotation: BoolProperty(
        name="Use Rotation",
        description="If the planet rotates",
        default=True)

    rotation_period: FloatProperty(
        name="Rotation Period (in seconds)",
        description="Time for one to rotation in seconds",
        soft_min=1, default=86400)

    # Surface
    with_surface: BoolProperty(
        name="Create Surface")

    mesh_type: EnumProperty(
        name="Mesh Type",
        items=(('UV_SPHERE', "UV Sphere", "UV sphere surface"),
               ('ICOSPHERE', "Icosphere", "Icosphere surface"),
               ),
        description=("Type of mesh, UV spheres are better for unwrapping, "
                     "Icospheres have an even distribution of vertices"),
        default='ICOSPHERE')

    subdivisions: IntProperty(
        name="Subdivisions",
        description="1 = coarse mesh, 7 = dense but round mesh",
        min=1, max=7, default=3)

    radius_mantissa: FloatProperty(
        name="Radius (in km)",
        description="The radius of the surface",
        soft_min=0.001, soft_max=1000, default=0.5)

    radius_exp: IntProperty(
        name="Radius Exponent",
        description="Base-10 exponent of the radius",
        soft_min=0, soft_max=8, default=6)

    def execute(self, context):
        # create emtpy object for planet
        planet = bpy.data.objects.new(self.planet_name, None)
        context.collection.objects.link(planet)

        # set up orbit properties
        simobj = planet.sssim_obj
        simobj.use_sssim = True
        simobj.object_type = 'PLANET'
        simobj.mass_mantissa = self.mass_mantissa
        simobj.mass_exp = self.mass_exp

        simorbit = planet.sssim_orbit
        center = context.scene.objects.get(self.center_name)
        if center:
            simorbit.center_name = center.name
        simorbit.distance_mantissa = self.distance_mantissa
        simorbit.distance_exp = self.distance_exp

        # set up rotation
        simrot = planet.sssim_rotation
        simrot.use_rotation = self.use_rotation
        simrot.use_frames = False
        simrot.period_seconds = self.rotation_period

        if self.with_surface:
            # use our create_surface operator,
            bpy.ops.object.sssim_create_surface(
                parent_name=planet.name,
                mesh_type=self.mesh_type,
                subdivisions=self.subdivisions,
                radius_mantissa=self.radius_mantissa,
                radius_exp=self.radius_exp)

        # make the planet object the active object
        planet.select_set(True)
        context.view_layer.objects.active = planet

        return {'FINISHED'}

    def invoke(self, context, event):
        # Try to use the active object as center
        obj = context.object
        if obj:
            simobj = obj.sssim_obj
            if simobj.use_sssim:
                # is a valid simulation object
                if simobj.object_type != 'SURFACE':
                    # use this object as center
                    self.center_name = obj.name
                else:
                    # object is surface, use parent
                    parent_name = obj.sssim_surface.parent_name
                    self.center_name = parent_name

        wm = context.window_manager
        return wm.invoke_props_dialog(self, width=300)

    def draw(self, context):
        layout = self.layout
        layout.prop_search(self, "center_name", context.scene, "objects")

        layout.prop(self, "planet_name")

        row = layout.row(align=True)
        row.prop(self, "mass_mantissa")
        row.prop(self, "mass_exp", slider=True)

        row = layout.row(align=True)
        row.prop(self, "distance_mantissa")
        row.prop(self, "distance_exp", slider=True)

        layout.separator()

        layout.prop(self, "use_rotation")
        layout.prop(self, "rotation_period")

        layout.separator()

        layout.prop(self, "with_surface")
        layout.prop(self, "mesh_type")
        col = layout.column(align=True)
        col.prop(self, "subdivisions")
        sub = col.row(align=True)
        sub.prop(self, "radius_mantissa")
        sub.prop(self, "radius_exp", slider=True)


class OBJECT_OT_sssim_create_surface(bpy.types.Operator):
    """Create a surface for the active planet or center"""
    bl_idname = "object.sssim_create_surface"
    bl_label = "Create Surface"

    parent_name: StringProperty(
        name="Parent Name",
        description="Name of the parent center or planet object")

    mesh_type: EnumProperty(
        name="Mesh Type",
        items=(('UV_SPHERE', "UV Sphere", "UV sphere surface"),
               ('ICOSPHERE', "Icosphere", "Icosphere surface"),
               ),
        description=("Type of Mesh, UV spheres are better for unwrapping, "
                     "Icospheres have an even distribution of vertices"),
        default='ICOSPHERE')

    subdivisions: IntProperty(
        name="Subdivisions",
        description="1 = coarse mesh, 7 = dense but round mesh",
        min=1, max=7, default=3)

    radius_mantissa: FloatProperty(
        name="Radius (in km)",
        description="The radius of the surface",
        soft_min=0.001, soft_max=1000, default=0.5)

    radius_exp: IntProperty(
        name="Radius Exponent",
        description="Base-10 exponent of the radius",
        soft_min=0, soft_max=8, default=6)

    def execute(self, context):
        # Create the surface from a primitive
        if self.mesh_type == 'UV_SPHERE':
            u_res = 2 ** (self.subdivisions + 1)
            v_res = 2 ** self.subdivisions
            bpy.ops.mesh.primitive_uv_sphere_add(
                segments=u_res,
                ring_count=v_res,
                location=(0, 0, 0))
        elif self.mesh_type == 'ICOSPHERE':
            bpy.ops.mesh.primitive_ico_sphere_add(
                subdivisions=self.subdivisions,
                location=(0, 0, 0))

        # smooth looking surface
        bpy.ops.object.shade_smooth()

        # set up surface properties
        surobj = bpy.context.object
        surobj.sssim_obj.use_sssim = True
        surobj.sssim_obj.object_type = 'SURFACE'

        simsur = surobj.sssim_surface
        simsur.calc_size = True
        simsur.radius_mantissa = self.radius_mantissa
        simsur.radius_exp = self.radius_exp

        if self.parent_name:
            surobj.name = "{}_surface".format(self.parent_name)
            surobj.sssim_surface.parent_name = self.parent_name

        return {'FINISHED'}

    def invoke(self, context, event):
        # try to use the active object as parent
        obj = context.object
        if obj:
            simobj = obj.sssim_obj
            if simobj.use_sssim:
                # is a valid simulation object
                if simobj.object_type != 'SURFACE':
                    # use this object as center
                    self.parent_name = obj.name
                else:
                    # object is surface, use parent
                    parent_name = obj.sssim_surface.parent_name
                    self.parent_name = parent_name

        wm = context.window_manager
        return wm.invoke_props_dialog(self, width=300)

    def draw(self, context):
        layout = self.layout
        layout.prop_search(self, "parent_name", context.scene, "objects")

        layout.prop(self, "mesh_type")
        col = layout.column(align=True)
        col.prop(self, "subdivisions")
        sub = col.row(align=True)
        sub.prop(self, "radius_mantissa")
        sub.prop(self, "radius_exp", slider=True)


class OBJECT_OT_sssim_bake_all(bpy.types.Operator):
    """Bake all simulations in the current scene"""
    bl_idname = "scene.sssim_bake_all"
    bl_label = "SSSim Bake All"

    def execute(self, context):
        scn = context.scene
        # handy name for the operator
        sim_bake = bpy.ops.object.sssim_to_fcurve
        # the operator works on the active object, we will override the context

        n = 0
        for obj in scn.objects:
            if not obj.sssim_obj.use_sssim:
                continue

            # deactivate Driver if neccessary
            if obj.sssim_calc.use_driver:
                obj.sssim_calc.use_driver = False

            # now override the active object and bake
            with context.temp_override(object=obj):
                if sim_bake.poll():
                    print("Baking {}".format(obj.name))
                    result = sim_bake()
                    n += 1

                    # check if everything was OK
                    if result == {'CANCELLED'}:
                        print("Object {} failed to bake".format(obj.name))
                        return {'CANCELLED'}

        print("Baked {} objects".format(n))
        return {'FINISHED'}


class OBJECT_OT_sssim_bake_clear(bpy.types.Operator):
    """Clear all baked simulations and use drivers again"""
    bl_idname = "scene.sssim_bake_clear"
    bl_label = "SSSim Bake Clear"

    def execute(self, context):
        scn = context.scene
        sim_clear = bpy.ops.object.sssim_clear_fcurve

        num_obj = 0
        for obj in scn.objects:
            if not obj.sssim_obj.use_sssim:
                continue

            with context.temp_override(object=obj):
                if sim_clear.poll():
                    sim_clear()
                    num_obj += 1

            obj.sssim_calc.use_driver = True
        print("Cleared Bake of {} objects".format(num_obj))
        return {'FINISHED'}

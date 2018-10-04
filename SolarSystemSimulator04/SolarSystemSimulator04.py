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
    "description": "Simulation of solar systems using kepler's laws of planetary motion.",
    "author": "Markus Ebke",
    "version": (0, 4),
    "blender": (2, 66, 0),
    "location": "Properties > Physics (for objects), Properties > Scene (global settings)",
    "warning": "In developement, version from 5.04.2013, more info: http://blenderartists.org/forum/showthread.php?267761-Solar-System-Simulator-WIP",
    "wiki_url": "",
    "tracker_url": "",
    "category": "Object"}

import bpy
from bpy.props import (
    BoolProperty,
    IntProperty,
    FloatProperty,
    StringProperty)

from mathutils import Vector, Euler
from math import (
    pi,
    sqrt,
    sin,
    cos,
    acos,
    degrees,
    ceil,)


G = 6.673e-11  # Gravitational constant


def true_anomaly_get(time, T, e=0, time_offset=0):
    """Calculate the true anomaly (the angle since the nearest point of the orbit)

    time = time in seconds
    T = Orbital Period in seconds
    e = eccentricity of the orbit, e = 0: circular orbit, 0 < e < 1: elliptical
    time_offset = offset of the initial angle, in range 0 - 1,
                  offset * T = starting time in seconds if time=0
    """
    if T == 0:
        # Planet cannot have an orbit -> return static position dependent on offset
        return 2 * pi * time_offset

    # mean anomaly = "angle" travelled, in range 0 - 2*pi
    mean_anomaly = 2 * pi * (((time / T) + time_offset) % 1)

    if e == 0:
        # circle, no fancy stuff here
        return mean_anomaly
    else:
        # the true anomaly is the "angle" travelled (as seen from center)
        # since the nearest point (periapsis) of the orbit

        # mean_anomaly (in rad)-> eccentric anomaly (rad)-> true anomaly (rad)
        # calculate the eccentric anomaly according to a formular found at
        # http://www-spof.gsfc.nasa.gov/stargaze/Smotion.htm
        eccentric_anomaly = mean_anomaly
        e_old = 0
        while abs(e_old - eccentric_anomaly) > 1e-10:  # precision = 10 digits
            e_old = eccentric_anomaly
            eccentric_anomaly = mean_anomaly + e * sin(e_old)

        # calculate the true anomaly
        # formular from http://en.wikipedia.org/wiki/True_anomaly
        t_0 = cos(eccentric_anomaly) - e
        t_1 = 1 - e * cos(eccentric_anomaly)
        true_anomaly = acos(t_0 / t_1)  # true anomaly in range: 0 - pi

        if mean_anomaly > pi:
            # we are on the other side (angle > 180°) of the circle
            true_anomaly = 2 * pi - true_anomaly  # range: 0 - 2*pi
        return true_anomaly


### Driver functions ###

def eval_planet_orbit(scn_name, obj_name, index=None, time=None):
    """Evaluate the planets position, used by driver.
    scn_name = Name of a scene which contains the object
    obj_name = Name of the object to simulate
    index = index of the location channel
    time = time when to calculate, if not given use current scene time
    returns a 3-tuple with xyz-locations or, if index given, only one component
    """
    scn = bpy.data.scenes.get(scn_name)
    obj = bpy.data.objects.get(obj_name)
    if not obj or not scn:
        errmsg = "DRIVER ERROR: Invalid obj_name ({}) or scn_name ({}) variable"
        print(errmsg.format(obj_name, scn_name))
        return 0

    sssim_scn = scn.sssim_scn
    sssim_obj = obj.sssim_obj

    # time = time for which to calculate the planet position, in seconds
    if time is None:
        time = sssim_scn.time

    planet_loc = orbit_position(sssim_obj, sssim_scn, time)

    orbit_rot = orbit_rotate(planet_loc, sssim_obj)

    if index is None:
        return orbit_rot
    else:
        return orbit_rot[index]


def orbit_position(sssim_obj, sssim_scn, time):
    """Calculate the position of the planet in the xy-plane"""
    # a = (to BU adjusted) semi_major_axis
    a = sssim_obj.orbit_semi_major_axis / (10 ** sssim_scn.length_exp)

    # eccentricity = more circular (small e) or more elliptical (bigger e) orbit
    e = sssim_obj.orbit_eccentricity

    # orbital period = time to complete one orbit
    T = sssim_obj.orbital_period

    time_offset = sssim_obj.orbit_time_offset

    # calculate the position
    true_anomaly = true_anomaly_get(time, T, e, time_offset)
    if e == 0:
        # circular orbit, r = constant
        r = a
    else:
        # elliptical orbit, r = distance from sun (in BU)
        # from http://www-spof.gsfc.nasa.gov/stargaze/Smotion.htm
        r = a * (1 - e ** 2) / (1 + e * cos(true_anomaly))

    x = cos(true_anomaly) * r
    y = sin(true_anomaly) * r
    z = 0
    planet_loc = Vector((x, y, z))

    return planet_loc


def orbit_rotate(other, sssim_obj):
    """Rotates a mathutils value (Vector, Euler, ...) using orbital elements."""
    inc = sssim_obj.orbit_inclination
    asc_node = sssim_obj.orbit_asc_node
    arg_periapsis = sssim_obj.orbit_arg_periapsis

    # rotate mathutils value around global z axis (argument of periapsis)
    other.rotate(Euler((0, 0, arg_periapsis)))

    # rotate around the global x axis by (inclination)
    other.rotate(Euler((inc, 0, 0)))

    # rotate around global z again (rotate the ascending node)
    other.rotate(Euler((0, 0, asc_node)))

    return other


def eval_planet_rotation(scn_name, obj_name, index=None, time=None):
    """Evaluate the planets rotation, used by driver.
    scn_name = Name of a scene which contains the object
    obj_name = Name of the object to simulate
    index = index of the rotation channel, usually only z-axis (index=2) changes
    time = time when to calculate, if not given use current scene time
    returns an Euler in mode ZYX or, if index given, an angle in radians
    """
    scn = bpy.data.scenes.get(scn_name)
    obj = bpy.data.objects.get(obj_name)
    if not obj or not scn:
        errmsg = "DRIVER ERROR: Invalid obj_name ({}) or scn_name ({}) variable"
        print(errmsg.format(obj_name, scn_name))
        return 0

    sssim_scn = scn.sssim_scn
    sssim_obj = obj.sssim_obj
    tilt = sssim_obj.rotation_axis_tilt
    direction = sssim_obj.rotation_axis_direction

    # time = time in seconds, if None use current scene time
    if time is None:
        time = sssim_scn.time

    # rotation_period is also in seconds
    rotation_period = sssim_obj.rotation_period
    if rotation_period != 0:
        rot_z = 2 * pi * time / rotation_period
    else:
        # invalid input -> no rotation
        rot_z = 0

    planet_rot = Euler((tilt, 0.0, 0.0), 'ZYX')  # note that mode is 'ZYX'

    # rotate around global (not local) z-axis
    planet_rot.rotate(Euler((0.0, 0.0, direction), 'XYZ'))

    # rotate around local z-axis
    # NOTE: we won't use planet_rot.rotate_axis('Z', rot_z) because then
    # all rotation are between -180° and 180° and for the rotation around
    # z we need a continous motion with increasing z values
    planet_rot.z += rot_z

    if sssim_obj.rotation_relative_to_orbit:
        planet_rot = orbit_rotate(planet_rot, sssim_obj)

    if index is None:
        return planet_rot
    else:
        return planet_rot[index]


### Useful functions ###

def add_driver_loc(obj, scn):
    """Add or adjust the location-drivers for orbit"""
    loc_list = obj.driver_add("location")  # X-, Y- and Z-Driver in list

    # location driver:
    for driver_fcurve in loc_list:
        index = driver_fcurve.array_index
        driver = driver_fcurve.driver
        driver.show_debug_info = True

        driver.type = "SCRIPTED"
        expr = "eval_planet_orbit('{0}', '{1}', {2})"
        driver.expression = expr.format(scn.name, obj.name, index)

    return "location driver added/adjusted"


def add_driver_rot(obj, scn):
    """Add or adjust the rotation-driver"""
    obj.rotation_mode = 'ZYX'  # Have to use a different mode because of rotation around Z
    rot_list = obj.driver_add("rotation_euler")  # X-, Y- and Z-Driver in list

    # rotation driver:
    for driver_fcurve in rot_list:
        index = driver_fcurve.array_index
        driver = driver_fcurve.driver
        driver.show_debug_info = True

        driver.type = "SCRIPTED"
        expr = "eval_planet_rotation('{0}', '{1}', {2})"
        driver.expression = expr.format(scn.name, obj.name, index)

    return "rotation driver added/adjusted"


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


def add_surface_constraint(child_obj, parent_obj):
    """Add/adjust child-parent contraint for surface child_obj and parent parent_obj"""
    surfacecon = child_obj.constraints.get("SURFACE")

    if surfacecon is None:
        surfacecon = child_obj.constraints.new('CHILD_OF')
        surfacecon.name = "SURFACE"

        # don't use the size of the parent, we will set it over surface radius
        surfacecon.use_scale_x = False
        surfacecon.use_scale_y = False
        surfacecon.use_scale_z = False

    surfacecon.target = parent_obj
    return "Surface: Constraint {} to {}".format(child_obj.name, parent_obj.name)


def remove_driver_loc(obj):
    """Remove the location-driver (= Orbit)"""
    rem_loc = obj.driver_remove("location")
    if rem_loc:
        return "location driver removed"
    else:
        return "no location driver to remove"


def remove_driver_rot(obj):
    """Remove the rotation-driver"""
    rem_rot = obj.driver_remove("rotation_euler")
    if rem_rot:
        return "rotation driver removed"
    else:
        return "no rotation driver to remove"


def remove_orbit_constraint(obj):
    """Remove the child-parent constraint"""
    orbitcon = obj.constraints.get("ORBIT")
    if orbitcon:
        obj.constraints.remove(orbitcon)
        return "Orbit Constraint removed from {}".format(obj.name)
    else:
        return "No Orbit Constraint to remove from {}".format(obj.name)


def remove_surface_constraint(child_obj):
    """Remove the child-parent constraint for surface"""
    surfacecon = child_obj.constraints.get("SURFACE")
    if surfacecon:
        child_obj.constraints.remove(surfacecon)
        return "Surface Constraint removed from {}".format(child_obj.name)
    else:
        return "No Surface Constraint to remove from {}".format(child_obj.name)


def get_fcurve(fcurves, search_data_path, array_index=-1):
    """Find in list fcurves the fcurve with given data_path (and index)"""
    for f in fcurves:
        if f.data_path == search_data_path:
            if array_index == -1 or f.array_index == array_index:
                return f
    return None


def has_location_fcurve(obj):
    """Find out if the object has a location (=orbit) fcurve."""
    if obj.animation_data:
        if obj.animation_data.action:
            fcurves = obj.animation_data.action.fcurves
            if fcurves:
                if get_fcurve(fcurves, "location"):
                    return True
    return False


def has_rotation_fcurve(obj):
    """Find out if the object has a rotation fcurve."""
    if obj.animation_data:
        if obj.animation_data.action:
            fcurves = obj.animation_data.action.fcurves
            if fcurves:
                if get_fcurve(fcurves, "rotation_euler"):
                    return True
    return False


def key_insert(fcurve, start, end, step, get_value):
    """Insert keys and refining the curve (adjusts the handles).
    fcurve = FCurve where the keys are added
    start, end, step = where to insert the keys
    get_value(frame), function which returns the value at the given frame
    returns number of inserted keyframes
    """
    frame_list = list(range(start, end, step))
    frame_list.append(end)

    for frame in frame_list:
        # time difference for calculating the previous and next position
        delta = 0.4
        frame_pre = frame - delta * step
        frame_past = frame + delta * step

        # y-coordinates of left handle, key, right handle
        val_pre = get_value(frame_pre)
        val = get_value(frame)
        val_past = get_value(frame_past)

        key = fcurve.keyframe_points.insert(frame, val)
        key.handle_left_type = 'FREE'
        key.handle_left = (frame_pre, val_pre)
        key.handle_right_type = 'FREE'
        key.handle_right = (frame_past, val_past)

    return len(frame_list)


### Operators ###

class AddEvalTimeFCurve(bpy.types.Operator):
    bl_idname = "scene.add_eval_time_fcurve"
    bl_label = "Add Time FCurve"
    bl_description = "Add F-Curve to Evaluation Time"

    @classmethod
    def poll(cls, context):
        # if eval_time is controlled by fcurve -> return False
        anim_data = context.scene.animation_data
        if anim_data and anim_data.action:
            if get_fcurve(anim_data.action.fcurves, "sssim_scn.eval_time"):
                return False
        return True

    def execute(self, context):
        scn = context.scene

        if scn.animation_data is None:
            scn.animation_data_create()

        anim_data = scn.animation_data
        anim_data.action = bpy.data.actions.new("Evaluation_Time")
        time_curve = anim_data.action.fcurves.new(data_path="sssim_scn.eval_time")

        fmod = time_curve.modifiers.new("GENERATOR")
        fmod.mode = "POLYNOMIAL"
        fmod.poly_order = 1
        fmod.coefficients = (-1, 1)
        return {'FINISHED'}


# Create and use a fcurves for the simualation instead of drivers
class SSSimToFCurve(bpy.types.Operator):
    bl_idname = "object.sssim_to_fcurve"
    bl_label = "SSSim to F-Curve"
    bl_description = "Orbit (if any) and Rotation (if any) to F-Curve"

    calc_orbit = BoolProperty()
    calc_rotation = BoolProperty()
    cyclic = BoolProperty()
    frame_start = IntProperty()
    frame_end = IntProperty()
    frame_step = IntProperty()
    motion_path = BoolProperty()

    @classmethod
    def poll(cls, context):
        if context.object:
            # There is an active object
            sssim_obj = context.object.sssim_obj
            if sssim_obj.use_sssim:
                if sssim_obj.center_object or sssim_obj.use_rotation:
                    # sssim is activated and there is an orbit or rotation
                    return True
        return False

    def orbit_to_fcurve(self, act, scn, obj, get_time):
        for index in range(0, 3):
            # value_func(frame) returns the value at this frame
            def value_func(frame):
                return eval_planet_orbit(scn.name, obj.name, index, get_time(frame))

            # Remove the existing curve if any
            fc = get_fcurve(act.fcurves, "location", index)
            if fc:
                act.fcurves.remove(fc)
            fc = act.fcurves.new(data_path="location", index=index)

            n = key_insert(fc, self.frame_start, self.frame_end, self.frame_step, value_func)

            if self.cyclic:
                fc.modifiers.new(type='CYCLES')

        if self.motion_path:
            # set up motion paths for visualization
            mpath_calc = bpy.ops.object.paths_calculate
            if mpath_calc.poll():
                mpath_calc(start_frame=self.frame_start, end_frame=self.frame_end)

        print("Created {:3} keyframes for location".format(n))

    def rotation_to_fcurve(self, act, scn, obj, get_time):
        def value_func(time):
            return eval_planet_rotation(scn.name, obj.name, 2, get_time(time))

        # Calc rotation around z-axis only, index = 2
        # rotation of x- and y-axis stay fixed
        rot = eval_planet_rotation(scn.name, obj.name, index=None, time=0)
        rot_x, rot_y = rot[0], rot[1]
        obj.rotation_euler.x = rot_x
        obj.rotation_euler.y = rot_y

        # Remove the existing z-axis curve if any
        fc = get_fcurve(act.fcurves, "rotation_euler", 2)
        if fc:
            act.fcurves.remove(fc)
        fc = act.fcurves.new(data_path="rotation_euler", index=2)
        fc.extrapolation = 'LINEAR'

        if self.cyclic:
            start = 1
            end = 1 + self.frame_step
        else:
            start = self.frame_start
            end = self.frame_end

        n = key_insert(fc, start, end, self.frame_step, value_func)

        print("Created {:3} keyframes for rotation".format(n))

    def execute(self, context):
        obj = context.object
        scn = context.scene

        # There must be a eval_time f-curve
        if bpy.ops.scene.add_eval_time_fcurve.poll():
            msg = "No Evaluation Time F-Curve found, control the scene tab"
            self.report({'ERROR'}, msg)
            return {'CANCELLED'}

        eval_time_fcurve = get_fcurve(scn.animation_data.action.fcurves, "sssim_scn.eval_time")
        time_multiply = 10 ** scn.sssim_scn.time_exp / scn.render.fps

         # time_func(frame) returns the time at this frame
        time_func = lambda frame: eval_time_fcurve.evaluate(frame) * time_multiply

        # Set up animation data if needed
        if obj.animation_data is None:
            obj.animation_data_create()
        anim_data = obj.animation_data

        # Use existing action to keep it clean
        previous = bpy.data.actions.get("SSSim_%s_Action" % obj.name)
        if previous:
            act = previous
        else:
            act = bpy.data.actions.new(name="SSSim_%s_Action" % obj.name)

        anim_data.action = act

        if self.calc_orbit:
            #print("Calculating Orbit:......", end="  ")
            self.orbit_to_fcurve(act, scn, obj, time_func)
        if self.calc_rotation:
            #print("Calculating Rotation:...", end="  ")
            self.rotation_to_fcurve(act, scn, obj, time_func)

        return {'FINISHED'}

    def invoke(self, context, event):
        sssim_obj = context.object.sssim_obj
        self.calc_orbit = sssim_obj.center_object is not None
        self.calc_rotation = sssim_obj.use_rotation

        self.cyclic = sssim_obj.cyclic
        if self.cyclic:
            self.frame_start = 1
            # end = start + orbital_period
            self.frame_end = 1 + ceil(sssim_obj.orbital_period_frames)
        else:
            self.frame_start = sssim_obj.frame_start
            self.frame_end = sssim_obj.frame_end

        self.frame_step = sssim_obj.frame_step
        self.motion_path = sssim_obj.motion_path

        return self.execute(context)


class SSSimClearFCurve(bpy.types.Operator):
    bl_idname = "object.sssim_clear_fcurve"
    bl_label = "Clear SSSim F-Curve"
    bl_description = "Clear location and rotation F-Curves"

    @classmethod
    def poll(cls, context):
        # search for location and rotation fcurve
        obj = context.object
        if obj:
            if has_location_fcurve(obj) or has_rotation_fcurve(obj):
                return True
        return False

    def execute(self, context):
        if bpy.ops.object.paths_clear.poll():
            bpy.ops.object.paths_clear()

        fcurves = context.object.animation_data.action.fcurves
        for i in range(3):
            loc_curve = get_fcurve(fcurves, "location", i)
            rot_curve = get_fcurve(fcurves, "rotation_euler", i)
            if loc_curve:
                fcurves.remove(loc_curve)
            if rot_curve:
                fcurves.remove(rot_curve)

        return {'FINISHED'}


# Surface Operators Start -----------------------------------------------------
class SurfaceAdd(bpy.types.Operator):
    bl_idname = "object.sssim_surface_add"
    bl_label = "Add Planet Surface"
    bl_description = "Add a Surface Item to the Planet"

    @classmethod
    def poll(cls, context):
        if context.object.sssim_obj.use_sssim:
            return True
        return False

    def execute(self, context):
        # add list item
        sssim_obj = context.object.sssim_obj
        item = sssim_obj.surfaces.add()

        # set index to the last added item at the end of the list
        sssim_obj.active_surface_index = len(sssim_obj.surfaces) - 1

        return {'FINISHED'}


class SurfaceRemove(bpy.types.Operator):
    bl_idname = "object.sssim_surface_remove"
    bl_label = "Remove Planet Surface"
    bl_description = "Remove the Active Surface"

    @classmethod
    def poll(cls, context):
        obj = context.object
        return obj.sssim_obj.use_sssim and len(obj.sssim_obj.surfaces) > 0

    def execute(self, context):
        sssim_obj = context.object.sssim_obj
        surfaces = sssim_obj.surfaces
        active_index = sssim_obj.active_surface_index

        # clear item.name to reset surface object
        surfaces[active_index].name = ""

        # remove the list item and update index
        surfaces.remove(active_index)
        if active_index > 0:
            sssim_obj.active_surface_index = active_index - 1

        return {'FINISHED'}


class SurfaceMove(bpy.types.Operator):
    bl_idname = "object.sssim_surface_move"
    bl_label = "Move Planet Surfaces"
    bl_description = "Move the active surface up or down the list"

    @classmethod
    def poll(cls, context):
        obj = context.object
        if obj.sssim_obj.use_sssim and len(obj.sssim_obj.surfaces) > 0:
            return True
        return False

    direction = StringProperty(default='UP')

    def execute(self, context):
        sssim_obj = context.object.sssim_obj
        surfaces = sssim_obj.surfaces
        active_index = sssim_obj.active_surface_index

        if self.direction == 'UP' and active_index != 0:
            surfaces.move(active_index, active_index - 1)
            sssim_obj.active_surface_index = active_index - 1
        elif self.direction == 'DOWN' and (active_index < len(surfaces) - 1):
            surfaces.move(active_index, active_index + 1)
            sssim_obj.active_surface_index = active_index + 1

        return {'FINISHED'}

# Surface Operators End -------------------------------------------------------


class UpdateSSSimDrivers(bpy.types.Operator):
    bl_idname = "scene.update_sssim_drivers"
    bl_label = "Update Drivers"
    bl_description = "Update Solar System Simulation Drivers"

    def execute(self, context):
        # Update driver namespace
        bpy.app.driver_namespace["eval_planet_orbit"] = eval_planet_orbit
        bpy.app.driver_namespace["eval_planet_rotation"] = eval_planet_rotation

        # update all objects in current scene
        scn = context.scene
        number_updates = 0
        for obj in scn.objects:
            sssim_obj = obj.sssim_obj
            if sssim_obj.use_sssim and sssim_obj.use_driver:
                if sssim_obj.center_object:
                    # update location driver
                    add_driver_loc(obj, scn)
                if sssim_obj.use_rotation:
                    # update rotation driver
                    add_driver_rot(obj, scn)
                number_updates += 1

        number_objects = len(scn.objects)
        msg = "Updated {} of {} objects".format(number_updates, number_objects)
        self.report({'INFO'}, msg)

        return {'FINISHED'}


### Group Properties for Data ####

# Surface list items
class SSSIMSurfaceListItems(bpy.types.PropertyGroup):
    def update_object(self, context):
        # print("update surface object", self.name)
        # remove the constraint of the previous surface if any
        if self.constraint_object:
            con_obj = bpy.data.objects.get(self.constraint_object)
            if con_obj:
                msg = remove_surface_constraint(con_obj)
            else:
                msg = "WARNING: {} not found, could not delete surface constraint"
                msg = msg.format(self.constraint_object)
            self.constraint_object = ""
            print(msg)

        if self.object_is_valid:  # can add constraint
            surobj = self.surface_object
            parent = self.id_data
            msg = add_surface_constraint(child_obj=surobj, parent_obj=parent)
            print(msg)
            # update the name of the constraint object
            self.constraint_object = self.name

            # get the proportions
            if surobj.dimensions.length != 0:
                self.object_proportions = surobj.dimensions / max(surobj.dimensions)
            else:
                # We can't scale this object because all components are 0
                # but we have to update object_proportions
                self.object_proportions = surobj.dimensions

    def update_size(self, context):
        if self.object_is_valid and self.calc_radius:
            # dimensions = proportions * diameter, diameter = radius * 2
            dim = Vector(self.object_proportions) * self.radius_bu * 2
            self.surface_object.dimensions = dim

    name = StringProperty(
        name="Object Name",
        description="Name of the Surface Object",
        update=update_object)
    # Name of the constraint object, needed if item.name changes, not visible in UI
    constraint_object = StringProperty()

    # Properties for the objects scale
    calc_radius = BoolProperty(
        name="Calculate Radius",
        description="Calculate the Scale of the Object instead of using the Current Scale",
        default=False,
        update=update_size)
    # Used for scaling to maintain the original proportions -> none uniform objects are allowed
    object_proportions = bpy.props.FloatVectorProperty()
    radius_coef = FloatProperty(
        name="Radius (km)",
        description="The Radius of the Object",
        soft_min=0.001, soft_max=1000, default=1,
        update=update_size)
    radius_exp = IntProperty(
        name="Exponent",
        description="The Radius Exponent",
        soft_min=0, soft_max=8, default=3,
        update=update_size)

    def _get_object(self):
        # Returns the object from bpy.data.objects
        return bpy.data.objects.get(self.name)

    surface_object = property(_get_object)

    def _validate_object(self):
        surobj = self.surface_object
        if surobj and surobj != self.id_data:
            # center object can't be a surface object (else we get dependency cycle)
            if surobj.name != self.id_data.sssim_obj.center:
                # There is a valid object we can use as surface
                return True
        return False

    object_is_valid = BoolProperty(
        name="Object is valid",
        get=_validate_object)

    def _get_radius_km(self):
        if self.object_is_valid and not self.calc_radius:
            # use the current dimensions to get the radius in km
            sssim_scn = bpy.context.scene.sssim_scn
            mult = (10 ** sssim_scn.length_exp) / sssim_scn.planet_size_mult
            surobj = self.surface_object
            return max(surobj.dimensions) * mult / 2
        return self.radius_coef * 10 ** self.radius_exp

    radius_km = FloatProperty(
        name="Radius in km",
        get=_get_radius_km)

    def _get_radius_bu(self):
        sssim_scn = bpy.context.scene.sssim_scn
        mult = sssim_scn.planet_size_mult / (10 ** sssim_scn.length_exp)
        return self.radius_km * mult

    radius_bu = FloatProperty(
        name="Radius in Blender Units",
        get=_get_radius_bu)


# SSSIM Data for Objects
class SSSIMObjectData(bpy.types.PropertyGroup):
    # ------------------------------------------------
    # Update methods:
    # ------------------------------------------------
    def update_relations(self, context):
        # Orbit Constraint
        if self.center != "":
            msg_con = add_orbit_constraint(self.id_data, self.center_object)
        else:
            msg_con = remove_orbit_constraint(self.id_data)

        # Location driver
        if self.center != "" and self.use_driver:
            msg_loc = add_driver_loc(self.id_data, context.scene)
        else:
            msg_loc = remove_driver_loc(self.id_data)

        # Rotation driver
        if self.use_rotation and self.use_driver:
            msg_rot = add_driver_rot(self.id_data, context.scene)
        else:
            msg_rot = remove_driver_rot(self.id_data)

        if self.use_driver:
            # Delete Fcurves if any
            if bpy.ops.object.sssim_clear_fcurve.poll():
                bpy.ops.object.sssim_clear_fcurve()

        print("Orbit of {}:".format(self.id_data.name))
        print("    ", msg_con)
        print("    ", msg_loc)
        print("    ", msg_rot)

    # Test if center is valid
    def update_center(self, context):
        obj = context.object
        if self.center == "":
            self.update_relations(context)
        else:
            if self.center_recursive_control(obj):
                self.update_relations(context)
            else:
                print("ERROR: property 'center' has no valid value:", self.center)
                self.center = ""  # will update again and remove relations

    # recursive search for cyclic dependencies of center objects
    def center_recursive_control(self, ref_obj):
        if self.center == "":
            # no center object -> no cycle possible
            return True
        else:
            center = bpy.data.objects.get(self.center)
            if center is None or center.name == ref_obj.name:
                # center invalid or cycle detected
                return False
            else:
                # no link back to ref_obj -> no cycle
                return center.sssim_obj.center_recursive_control(ref_obj)

    def check_calc_frames(self, context):
        # end always later than start
        if self.frame_end <= self.frame_start:
            self.frame_end = self.frame_start + 1

    def update_surface_index(self, context):
        if len(self.surfaces) == 0:
            self.active_surface_index = 0
        elif self.active_surface_index >= len(self.surfaces):
            self.active_surface_index = len(self.surfaces) - 1

    # ------------------------------------------------
    # Properties:
    # ------------------------------------------------
    use_sssim = BoolProperty(
        name="Use Solar System Simulator",
        description="Use the Solar System Simulator",
        default=False,
        update=update_relations)

    ### Calculation Panel ###
    use_driver = BoolProperty(
        name="Use Driver",
        description="Use Driver for Calculation of Orbit and Rotation",
        default=True,
        update=update_relations)
    cyclic = BoolProperty(
        name="Cyclic",
        description="Add Cyclic F-Curve Modifier to Repeat the Motion",
        default=True)
    frame_start = IntProperty(
        name="Start Frame",
        description="First frame for the calculation",
        default=1,
        update=check_calc_frames)
    frame_end = IntProperty(
        name="End Frame",
        description="Last frame for the calculation",
        default=100,
        update=check_calc_frames)
    frame_step = IntProperty(
        name="Frame Step",
        description="Number of frames between calculated positions",
        min=1, default=5)
    motion_path = BoolProperty(
        name="Show Motion Path",
        description="Calculate and show the object's Motion Path",
        default=True)

    ### Object Panel ###
    show_info = BoolProperty(
        name="Show Info",
        description="Show additional info",
        default=True)
    center = StringProperty(
        name="Center Object",
        description="Name of the object as the orbit's center",
        update=update_center)
    mass_coef = FloatProperty(
        name="Mass (kg)",
        description="Mass of the object (in kg)",
        soft_min=0.001, soft_max=1000, default=1)
    mass_exp = IntProperty(
        name="Mass Exponent",
        description="Exponent for mass",
        soft_min=0, soft_max=50, default=24)

    ### Orbit Panel ###
    # Shape and Size
    orbit_eccentricity = FloatProperty(
        name="Eccentricity",
        description="The 'ellipticalness' of the orbit, 0 equals a circle",
        min=0, max=0.99, default=0)
    orbit_distance_coef = FloatProperty(
        name="Semi-Major Axis (km)",
        description="The Distance to the center, in a Circle: Radius, in an Ellipse: Semi Major Axis",
        soft_min=0.001, soft_max=1000, default=1)
    orbit_distance_exp = IntProperty(
        name="Length Exponent",
        description="Exponent of the Length to make it bigger",
        soft_min=0, soft_max=12, default=6)

    # Orientation of the orbital plane
    orbit_inclination = FloatProperty(
        name="Inclination",
        description="Inclination (= tilt against the xy-plane)",
        min=0, max=2 * pi,
        unit='ROTATION')
    orbit_asc_node = FloatProperty(
        name="Ascending Node",
        description="Longitude of the Ascending Node (= rotation of the orbit around the global z-axis)",
        min=0, max=2 * pi,
        unit='ROTATION')
    orbit_arg_periapsis = FloatProperty(
        name="Argument of Periapsis",
        description="Argument of Periapsis (= rotation of the nearest orbit point around the local z-axis)",
        min=0, max=2 * pi,
        unit='ROTATION')

    # Position
    orbit_time_offset = FloatProperty(
        name="Time Offset",
        description="Temporal Offset (in Percent) of the Orbital Period",
        min=0, max=1,
        subtype='PERCENTAGE')

    # Orbital Period specified by user, in frames
    use_user_orbital_period = BoolProperty(
        name="Use User Orbital Period",
        description="Use a Custom Orbital Period",
        default=False)
    use_orbital_frames = BoolProperty(
        name="Use Frames",
        description="Use Frames instead of Seconds for the Orbital Period",
        default=True)
    user_orbital_period_seconds = FloatProperty(
        name="Orbital Period (seconds)",
        description="Custom Orbital Period in Seconds",
        min=0, default=10000000)
    user_orbital_period_frames = FloatProperty(
        name="Orbital Period (frames)",
        description="Custom Orbital Period in Frames",
        min=0, default=100)

    ### Rotation properties ###
    use_rotation = BoolProperty(
        name="Use Rotation",
        description="Enable Rotation of the object",
        default=False,
        update=update_relations)
    rotation_use_frames = BoolProperty(
        name="Use Frames",
        description="Use Frames instead of Seconds for the Rotation Period",
        default=True,
        update=update_relations)
    user_rotation_period_frames = FloatProperty(
        name="Rotation Period (frames)",
        description="Time for one Rotation in Frames",
        min=0, default=100)
    user_rotation_period_seconds = FloatProperty(
        name="Rotation Period (seconds)",
        description="Time for one Rotation in Seconds",
        min=0, default=86400)
    rotation_axis_tilt = FloatProperty(
        name="Axis Tilt",
        description="Tilt of the rotation axis",
        min=0, max=pi, default=0,
        unit='ROTATION')
    rotation_axis_direction = FloatProperty(
        name="Tilt Direction",
        description="Direction in which to tilt",
        min=0, max=2 * pi, default=0,
        unit='ROTATION')
    rotation_relative_to_orbit = BoolProperty(
        name="Relative to Orbit",
        description="Axis direction relative to orientation of the orbital plane",
        default=True)

    ### Surface Properties ###
    surfaces = bpy.props.CollectionProperty(type=SSSIMSurfaceListItems)

    active_surface_index = IntProperty(
        name="Active Surface Index",
        min=0,
        update=update_surface_index)

    # ------------------------- Read-Only Properties ---------------------------

    def _get_center_object(self):
        if self.center != "":
            return bpy.data.objects.get(self.center)
        return None

    center_object = property(_get_center_object)

    def _get_surface(self):
        # Returns the active list item, type: SSSIMSurfaceListItems
        if len(self.surfaces) > 0:
            return self.surfaces[self.active_surface_index]
        return None

    active_surface = property(_get_surface)

    # The mass of the object
    def _mass_calc(self):
        return self.mass_coef * 10 ** self.mass_exp

    mass = FloatProperty(name="Mass", get=_mass_calc)

    # Semi Major Axis = half the longest axis of an ellipse
    def _semi_major_axis_calc(self):
        return self.orbit_distance_coef * 10 ** self.orbit_distance_exp

    orbit_semi_major_axis = FloatProperty(name="Semi Major Axis", get=_semi_major_axis_calc)

    # Semi-Minor Axis
    def _semi_minor_axis_calc(self):
        return self.orbit_semi_major_axis * sqrt(1 - self.orbit_eccentricity ** 2)

    orbit_semi_minor_axis = FloatProperty(name="Semi Minor Axis", get=_semi_minor_axis_calc)

    # Periapsis = Nearest Point in Orbit
    def _periapsis_calc(self):
        return (1 - self.orbit_eccentricity) * self.orbit_semi_major_axis

    orbit_periapsis = FloatProperty(name="Periapsis", get=_periapsis_calc)

    # Apoapsis = Furthest Point in Orbit
    def _apoapsis_calc(self):
        return (1 + self.orbit_eccentricity) * self.orbit_semi_major_axis

    orbit_apoapsis = FloatProperty(name="Apoapsis", get=_apoapsis_calc)

    # Orbital Period, time in seconds to complete one orbit, used to calc the orbit
    def _orbital_period_calc(self):
        if self.center_object:
            if self.use_user_orbital_period:
                if self.use_orbital_frames:
                    scn = bpy.context.scene
                    fps = scn.render.fps
                    return self.user_orbital_period_frames * 10 ** scn.sssim_scn.time_exp / fps
                else:
                    return self.user_orbital_period_seconds
            else:
                center = self.center_object
                a = self.orbit_semi_major_axis * 1000  # meter = kilometer * 1000
                m1 = self.mass
                m2 = center.sssim_obj.mass

                return sqrt((4 * pi ** 2 * a ** 3) / (G * (m1 + m2)))
        return 0  # No center, no orbit

    orbital_period = FloatProperty(name="Orbital Period in Seconds", get=_orbital_period_calc)

    # Orbital Period in frames, just for the UI
    def _orbital_period_frames_calc(self):
        scn = bpy.context.scene
        fps = scn.render.fps
        return self.orbital_period * fps / (10 ** scn.sssim_scn.time_exp)

    orbital_period_frames = FloatProperty(name="Orbital Period in Frames", get=_orbital_period_frames_calc)

    # True anomaly, angle of the current position since periapsis
    def _true_anomaly_calc(self):
        time = bpy.context.scene.sssim_scn.time
        T = self.orbital_period
        e = self.orbit_eccentricity
        time_offset = self.orbit_time_offset
        return true_anomaly_get(time, T, e, time_offset)

    orbit_true_anomaly = FloatProperty(name="True Anomaly", get=_true_anomaly_calc)

    # The Planet's rotation period in seconds
    def _rot_period_get(self):
        if self.use_rotation:
            if self.rotation_use_frames:
                scn = bpy.context.scene
                fps = scn.render.fps
                return self.user_rotation_period_frames * 10 ** scn.sssim_scn.time_exp / fps
            else:
                return self.user_rotation_period_seconds
        else:
            return 0  # No rotation

    rotation_period = FloatProperty(name="Rotation Period in seconds", get=_rot_period_get)

    # Rotation Period in frames, just for the UI
    def _rot_period_frames_calc(self):
        if self.use_rotation and self.rotation_use_frames:
            return self.user_rotation_period_frames
        scn = bpy.context.scene
        fps = scn.render.fps
        return self.rotation_period * fps / (10 ** scn.sssim_scn.time_exp)

    rotation_period_frames = FloatProperty(name="Rotation Period in Frames", get=_rot_period_frames_calc)


# SSSIM Scaling Data and Time for the Scene (access over scene.sssim_scn)
class SSSIMSceneData(bpy.types.PropertyGroup):
    def update_objects(self, context):
        for obj in context.scene.objects:
            if obj.sssim_obj.use_sssim:
                # update the location and rotation of the object
                obj.update_tag({'OBJECT'})

                # update the surfaces if any
                if len(obj.sssim_obj.surfaces) > 0:
                    for sur in obj.sssim_obj.surfaces:
                        sur.update_size(context)

        # redraw all visible areas
        for area in context.screen.areas:
            area.tag_redraw()

    eval_time = FloatProperty(
        name="Evaluation Time",
        description="Main control for the simulation time (in frames)",
        min=0, default=0)

    time_exp = IntProperty(
        name="Time Exponent",
        description="Exponent for converting simualtion time to animation time",
        min=0, max=10, default=6,
        update=update_objects)

    length_exp = IntProperty(
        name="Length Exponent",
        description="Exponent for scaling lengths",
        min=0, max=10, default=6,
        update=update_objects)

    planet_size_mult = IntProperty(
        name="Planet Size Multiplier",
        description="Correct the size of planet surfaces to make them better visible",
        min=1, soft_max=1000, default=1,
        update=update_objects)

    # simulation time in seconds:
    def _time_get(self):
        scn = self.id_data
        return self.eval_time * (10 ** self.time_exp) / scn.render.fps
    time = FloatProperty(
        name="Time",
        get=_time_get)


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
        sssim_obj = context.object.sssim_obj
        info = sssim_obj.show_info

        if bpy.ops.scene.add_eval_time_fcurve.poll():
            box = layout.box()
            box.label(text="You need to add an Evaluation Time F-Curve", icon='ERROR')
            box.operator("scene.add_eval_time_fcurve", icon='FCURVE')
            box.label(text="You can find the Settings in the Scene Tab", icon='SCENE_DATA')

        if context.object.constraints.get("SURFACE"):
            layout.label("This Object is a Surface or has a Surface Constraint", icon='ERROR')

        layout.active = sssim_obj.use_sssim
        layout.prop(sssim_obj, "show_info", icon='INFO')
        layout.prop_search(sssim_obj, "center", context.scene, "objects")

        row = layout.row(align=True).split(percentage=0.6)
        row.prop(sssim_obj, "mass_coef")
        row.prop(sssim_obj, "mass_exp", slider=True)
        if info:
            mass_coef = round(sssim_obj.mass_coef, 3)
            mass_exp = sssim_obj.mass_exp
            layout.label(text="Mass: {} * 10^{} kg".format(mass_coef, mass_exp))


# Panel for Orbit settings
class SSSIMPanelOrbit(bpy.types.Panel):
    bl_label = "Solar System Simulator: Orbit"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "physics"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        sssim_obj = context.object.sssim_obj
        return sssim_obj.use_sssim and sssim_obj.center_object

    def draw(self, context):
        layout = self.layout
        sssim_obj = context.object.sssim_obj
        sssim_scn = context.scene.sssim_scn
        info = sssim_obj.show_info

        if has_location_fcurve(context.object):
            # Deactivate layout if location fcurve found
            layout.label(text="Found Orbit F-Curve", icon='INFO')
            layout.active = False
        elif not sssim_obj.use_driver:
            # There is no orbit because driver not active
            layout.label(text="Orbit Driver not Active", icon='ERROR')

        # Shape and Size:
        box = layout.box()
        box.label(text="Shape and Size of the Orbit:")

        split = box.row().split(percentage=0.6)
        row = split.row()
        row.prop(sssim_obj, "orbit_eccentricity", slider=True)

        row = split.row()
        ecc = sssim_obj.orbit_eccentricity
        if info:
            if ecc == 0:
                row.label(text="Circular Orbit", icon='SPHERECURVE')
            else:
                row.label(text="Elliptical Orbit", icon='ROOTCURVE')

        row = box.row(align=True).split(percentage=0.6)
        if ecc == 0:
            row.prop(sssim_obj, "orbit_distance_coef", text="Radius (km)")
        else:
            row.prop(sssim_obj, "orbit_distance_coef")
        row.prop(sssim_obj, "orbit_distance_exp", slider=True)
        if info:
            a_coef = round(sssim_obj.orbit_distance_coef, 3)
            a_exp = sssim_obj.orbit_distance_exp
            # lengths in km / div_bu = lengths in Blender Units
            div_bu = 10 ** sssim_scn.length_exp
            a_bu = round(sssim_obj.orbit_semi_major_axis / div_bu, 2)
            data = "{} * 10^{} km ( ={:,} Blender Units)".format(a_coef, a_exp, a_bu)
            if ecc == 0:
                msg = "Radius: {}".format(data)
                box.label(text=msg, icon='CURVE_BEZCIRCLE')
            else:
                msg = "Semi-Major Axis: {}".format(data)
                box.label(text=msg, icon='CURVE_BEZCIRCLE')

        layout.separator()

        split = layout.split()
        col = split.column()

        # Orientation:
        box = col.box()
        box.label(text="Orientation of the Orbital Plane:")
        subcol = box.column(align=True)
        subcol.prop(sssim_obj, "orbit_inclination")
        subcol.prop(sssim_obj, "orbit_asc_node")
        subsub = subcol.row()
        subsub.active = sssim_obj.orbit_eccentricity != 0
        subsub.prop(sssim_obj, "orbit_arg_periapsis")

        col = split.column()

        # Position:
        box = col.box()
        box.label(text="Position:")
        subcol = box.column()
        subcol.prop(sssim_obj, "orbit_time_offset", slider=True)
        if info:
            # True Anomaly (angle travelled)
            angle = degrees(sssim_obj.orbit_true_anomaly)
            subcol.label(text="True Anomaly: {0}°".format(round(angle, 3)), icon='FORCE_HARMONIC')

        layout.separator()

        # Orbital Period:
        box = layout.box()
        box.label(text="Orbital Period:", icon='SORTTIME')

        if info:
            orb_sec = round(sssim_obj.orbital_period)
            orb_days = round(sssim_obj.orbital_period / 86400, 3)
            orb_frames = round(sssim_obj.orbital_period_frames, 2)

            if sssim_obj.use_user_orbital_period and sssim_obj.use_orbital_frames:
                # Don't show frames again when the user uses frames
                orb_text = "{0:,} s = {1} d".format(orb_sec, orb_days)
            else:
                orb_text = "{0:,} s = {1} d (= {2} frames)"
                orb_text = orb_text.format(orb_sec, orb_days, orb_frames)
            box.label(text=orb_text)

        row = box.row()
        row.prop(sssim_obj, "use_user_orbital_period")
        if sssim_obj.use_user_orbital_period:
            row.prop(sssim_obj, "use_orbital_frames")

            if sssim_obj.use_orbital_frames:
                box.prop(sssim_obj, "user_orbital_period_frames")
            else:
                box.prop(sssim_obj, "user_orbital_period_seconds")


# Panel for self rotation
class SSSIMPanelRotation(bpy.types.Panel):
    bl_label = "Solar System Simulator: Rotation"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "physics"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return context.object.sssim_obj.use_sssim

    def draw_header(self, context):
        sssim_obj = context.object.sssim_obj
        self.layout.prop(sssim_obj, "use_rotation", text="")

    def draw(self, context):
        layout = self.layout
        sssim_obj = context.object.sssim_obj
        info = sssim_obj.show_info
        layout.active = sssim_obj.use_rotation

        if has_rotation_fcurve(context.object):
            # Deactivate layout if rotation fcurve found
            layout.label(text="Found Rotation F-Curve", icon='INFO')
            layout.active = False
        elif not sssim_obj.use_driver:
            # No rotation because driver not active
            layout.label(text="Rotation Driver not Active", icon='ERROR')

        if info:
            layout.label(text="Rotation:", icon='PREVIEW_RANGE')

            rot_sec = round(sssim_obj.rotation_period, 2)
            rot_hours = round(sssim_obj.rotation_period / 3600, 2)
            rot_frames = round(sssim_obj.rotation_period_frames, 2)

            # Don't show frames again when the user uses frames
            if sssim_obj.rotation_use_frames:
                rot_text = "{0:,} s = {1} h".format(rot_sec, rot_hours)
            else:
                rot_text = "{0:,} s = {1} h (= {2} frames)"
                rot_text = rot_text.format(rot_sec, rot_hours, rot_frames)
            layout.label(text=rot_text)

        layout.prop(sssim_obj, "rotation_use_frames")
        if sssim_obj.rotation_use_frames:
            layout.prop(sssim_obj, "user_rotation_period_frames")
        else:
            layout.prop(sssim_obj, "user_rotation_period_seconds")

        row = layout.row()
        row.label(text="Rotation Axis Direction:")
        sub = row.column()
        sub.active = sssim_obj.center_object is not None
        sub.prop(sssim_obj, "rotation_relative_to_orbit")

        row = layout.row(align=True)
        row.prop(sssim_obj, "rotation_axis_tilt")
        row.prop(sssim_obj, "rotation_axis_direction")


# Panel for Planet surfaces:
class SSSIMPanelSurface(bpy.types.Panel):
    bl_label = "Solar System Simulator: Surfaces"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "physics"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return context.object.sssim_obj.use_sssim

    def draw(self, context):
        layout = self.layout

        sssim_obj = context.object.sssim_obj
        info = sssim_obj.show_info
        surface = sssim_obj.active_surface

        rows = 2
        if surface:
            rows = 4

        if info:
            text = "Items: {}, Active Index: {}"
            text = text.format(len(sssim_obj.surfaces), sssim_obj.active_surface_index)
            layout.label(text=text)

        row = layout.row()
        row.template_list("SSSIMSurfaceList", "", sssim_obj, "surfaces", sssim_obj, "active_surface_index", rows=rows)
        col = row.column(align=True)
        col.operator("object.sssim_surface_add", text="", icon='ZOOMIN')
        col.operator("object.sssim_surface_remove", text="", icon='ZOOMOUT')
        if surface:
            col.separator()
            col.operator("object.sssim_surface_move", icon='TRIA_UP', text="").direction = 'UP'
            col.operator("object.sssim_surface_move", icon='TRIA_DOWN', text="").direction = 'DOWN'

        if surface:
            self.draw_item(context, layout, surface, info)

    def draw_item(self, context, layout, item, info):
        layout.prop_search(item, "name", context.scene, "objects")

        if item.object_is_valid:
            row = layout.row()
            row.prop(item, "calc_radius")
            if info:
                row.label(text="Radius in BU: {:.5}".format(item.radius_bu))

            if item.calc_radius:
                # show radius controls
                row = layout.row(align=True)
                row.prop(item, "radius_coef")
                row.prop(item, "radius_exp", slider=True)
            else:
                # show dimension controls of surface object
                surobj = item.surface_object
                layout.label(text="Dimensions of {}:".format(surobj.name))
                layout.prop(surobj, "dimensions", text="")
        else:
            if info:
                msg = "'{}' is not a valid surface object".format(item.name)
                layout.label(text=msg, icon='ERROR')


# The Surface List as drawn in the UI
class SSSIMSurfaceList(bpy.types.UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            if item.object_is_valid:
                surobj = item.surface_object

                # if the object has an active material, use that icon
                if surobj.active_material:
                    icon = layout.icon(surobj.active_material)

                # use the first half for icon and name
                split = layout.split(0.5)
                row = split.row()
                row.label(text=item.name, translate=False, icon_value=icon)

                # second half is for radius and visibility in 3d-view
                row = split.row()
                radius_km = round(item.radius_km, 3)
                row.label(text="{:,} km".format(radius_km))
                row.prop(surobj, "hide", text="", emboss=False)
            else:
                # If the object is not valid show info text
                layout.label(text="<No Valid Object>")
        elif self.layout_type in {'GRID'}:
            layout.alignment = 'center'
            # the (valid) object has an active material, show that icon
            if item.object_is_valid and item.surface_object.active_material:
                icon = layout.icon(item.surface_object.active_material)
            layout.label(text="", icon_value=icon)


# Panel for Calculation settings
class SSSIMPanelCalculation(bpy.types.Panel):
    bl_label = "Solar System Simulator: Calculation"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "physics"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return context.object.sssim_obj.use_sssim

    def draw(self, context):
        layout = self.layout
        sssim_obj = context.object.sssim_obj
        calc_orbit = sssim_obj.center_object is not None
        calc_rot = sssim_obj.use_rotation
        layout.prop(sssim_obj, "use_driver", icon="DRIVER")

        col = layout.column()
        col.enabled = not sssim_obj.use_driver
        col.label(text="Create Simulation F-Curve:")

        split = col.split()
        col = split.column()
        sub = col.column(align=True)
        sub.operator("object.sssim_to_fcurve")
        sub.operator("object.sssim_clear_fcurve")

        if calc_orbit:
            col.label(text="Calc Orbit: Yes")
        else:
            col.label(text="Calc Orbit: No")

        if calc_rot:
            col.label(text="Calc Rotaton: Yes")
        else:
            col.label(text="Calc Rotaton: No")

        sub = col.row()
        sub.active = calc_orbit
        sub.prop(sssim_obj, "motion_path")

        col = split.column()
        sub = col.row()
        sub.active = calc_orbit
        sub.prop(sssim_obj, "frame_step")

        col.prop(sssim_obj, "cyclic")
        sub = col.column(align=True)
        if not sssim_obj.cyclic:
            sub.prop(sssim_obj, "frame_start")
            sub.prop(sssim_obj, "frame_end")
            start = sssim_obj.frame_start
            end = sssim_obj.frame_end
        else:
            start = 1
            end = 1 + ceil(sssim_obj.orbital_period_frames)
            sub.label(text="Start Frame: 1")
            sub.label(text="End Frame:   {}".format(end))

        keys_needed = ceil((end - start) / sssim_obj.frame_step) + 1
        col.label(text="Needed Keyframes: {}".format(keys_needed))


# Panel for scene properties in scene tab
class SSSIMPanelScene(bpy.types.Panel):
    bl_label = "Solar System Simulator: Scene Data"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "scene"

    def draw(self, context):
        layout = self.layout
        sssim_scn = context.scene.sssim_scn

        if bpy.ops.scene.add_eval_time_fcurve.poll():
            layout.label(text="You need to add an Evaluation Time F-Curve", icon='ERROR')

        row = layout.row(align=True)
        if bpy.ops.scene.add_eval_time_fcurve.poll():
            row.operator("scene.add_eval_time_fcurve", icon='FCURVE')
        row.prop(sssim_scn, "eval_time")

        layout.label(text="Time: {0:,} sec".format(round(sssim_scn.time)), icon='TIME')

        layout.separator()

        layout.label(text="Time Scaling: 1 sec =  {0:,} sec".format(10 ** sssim_scn.time_exp))
        layout.prop(sssim_scn, "time_exp", slider=True)

        layout.separator()

        layout.label(text="Length Scaling: 1 Blender Unit = {0:,} km".format(10 ** sssim_scn.length_exp))
        layout.prop(sssim_scn, "length_exp", slider=True)
        layout.prop(sssim_scn, "planet_size_mult")

        layout.separator()
        layout.operator("scene.update_sssim_drivers", icon='FILE_REFRESH')


# adjust the physics-tab panel
def physics_panel(self, context):
    self.layout.prop(context.object.sssim_obj, "use_sssim", icon="LAMP_SUN")


### Registration ###

def register():
    bpy.utils.register_module(__name__)

    bpy.types.Object.sssim_obj = bpy.props.PointerProperty(type=SSSIMObjectData)
    bpy.types.Scene.sssim_scn = bpy.props.PointerProperty(type=SSSIMSceneData)

    # We need to add custom drivers for the location and rotation
    bpy.app.driver_namespace["eval_planet_orbit"] = eval_planet_orbit
    bpy.app.driver_namespace["eval_planet_rotation"] = eval_planet_rotation

    bpy.types.PHYSICS_PT_add.append(physics_panel)


def unregister():
    bpy.types.PHYSICS_PT_add.remove(physics_panel)

    # Remove the drivers if any
    if "eval_planet_orbit" in bpy.app.driver_namespace:
        del bpy.app.driver_namespace["eval_planet_orbit"]
    if "eval_planet_rotation" in bpy.app.driver_namespace:
        del bpy.app.driver_namespace["eval_planet_rotation"]

    del bpy.types.Object.sssim_obj
    del bpy.types.Scene.sssim_scn

    bpy.utils.unregister_module(__name__)


if __name__ == "__main__":
    print("----------SOLARSYSTEMSIMULATOR-START----------")  # script starts
    register()

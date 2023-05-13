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

from math import pi, sqrt

import bpy
from bpy.props import (
    BoolProperty,
    IntProperty,
    FloatProperty,
    StringProperty,
    EnumProperty)
from mathutils import Vector

from . import operators, calculation
from .opengl import draw_3d_callback, draw_2d_callback


# handles for OpenGl drawing, used in SSSIMScene, must be a global variable
handle_3d = handle_2d = None

# printing to console
verbose = False


def update_relations(self, context):
    """Update all relations of this object with other objects.
    updates - orbit contraints
            - drivers for location and rotation (deletes fcurves if neccessary)
            - surface parenting
    """
    obj = self.id_data

    if verbose:
        print("Update relations of", self)
        print("Orbit of {}:".format(obj.name))

    simobj = obj.sssim_obj
    simorbit = obj.sssim_orbit
    simrot = obj.sssim_rotation
    simcalc = obj.sssim_calc
    simsur = obj.sssim_surface

    is_planet = simobj.object_type == 'PLANET' and simorbit.center_is_valid
    is_surface = simobj.object_type == 'SURFACE' and simsur.parent_is_valid

    # orbit constraint
    if is_planet:
        msg_con = operators.add_orbit_constraint(obj, simorbit.center_object)
    else:
        msg_con = operators.remove_orbit_constraint(obj)

    if verbose:
        print("    ", msg_con)

    # location driver
    if is_planet and simcalc.use_driver:
        operators.add_driver_loc(obj, context.scene)
        msg_loc = "location driver added/adjusted"
    else:
        rem_loc = operators.remove_driver_loc(obj)
        if rem_loc:
            msg_loc = "location driver removed"
        else:
            msg_loc = "no location driver to remove"

    if verbose:
        print("    ", msg_loc)

    # rotation driver
    if simrot.use_rotation and simcalc.use_driver and not is_surface:
        operators.add_driver_rot(obj, context.scene)
        msg_rot = "rotation driver added/adjusted"
    else:
        rem_rot = operators.remove_driver_rot(obj)
        if rem_rot:
            msg_rot = "rotation driver removed"
        else:
            msg_rot = "no rotation driver to remove"

    if verbose:
        print("    ", msg_rot)

    # delete F-Curves if any
    if simcalc.use_driver and bpy.ops.object.sssim_clear_fcurve.poll():
        bpy.ops.object.sssim_clear_fcurve()
        if verbose:
            print("    ", "cleared simulation fcurves")

    # surface
    if is_surface:
        msg_sur = operators.add_surface(obj, simsur.parent_object)
    else:
        msg_sur = operators.remove_surface(obj)
        # will reset size if necessary
        if obj.sssim_surface.calc_size:
            obj.sssim_surface.calc_size = False
    if verbose:
        print("    ", msg_sur)


# =============================================================================
# Scene data
# =============================================================================

# SSSIM scaling and time settings for the scene (access via scene.sssim_scn)
class SSSIMScene(bpy.types.PropertyGroup):
    # ------------------------------------------------
    # Update methods:
    # ------------------------------------------------
    def update_objects(self, context):
        """Update orbits of planets and size of surfaces"""
        for obj in context.scene.objects:
            simobj = obj.sssim_obj
            if simobj.use_sssim:
                # update the location and rotation of the object
                obj.update_tag(refresh={'OBJECT'})
                # update size of surfaces
                if simobj.object_type == 'SURFACE':
                    obj.sssim_surface.update_size(context)

        # redraw all visible areas
        for area in context.screen.areas:
            area.tag_redraw()

    def update_handle(self, context):
        """If draw_orbit then add handle (if neccessary), else remove it"""
        global handle_3d, handle_2d

        if self.draw_orbit and handle_3d is None:
            # the arguments for draw_orbit_callback
            args = (self, context)
            handle_3d = bpy.types.SpaceView3D.draw_handler_add(
                draw_3d_callback, args,
                'WINDOW', 'POST_VIEW')
            handle_2d = bpy.types.SpaceView3D.draw_handler_add(
                draw_2d_callback, args,
                'WINDOW', 'POST_PIXEL')
            # print("OpenGL drawing handle added")
        elif not self.draw_orbit and handle_3d is not None:
            bpy.types.SpaceView3D.draw_handler_remove(handle_3d, 'WINDOW')
            bpy.types.SpaceView3D.draw_handler_remove(handle_2d, 'WINDOW')
            handle_3d = handle_2d = None
            # print("OpenGL drawing handle removed")

    # ------------------------------------------------
    # Properties:
    # ------------------------------------------------

    sim_time: FloatProperty(
        name="Simulation Time",
        description="Main control of the simulation time (in seconds)",
        min=0, default=0)

    time_exp: IntProperty(
        name="Time Exponent",
        description=("Base-10 exponent for converting simulation time "
                     "to Blender time"),
        min=0, max=10, default=6,
        update=update_objects)

    length_exp: IntProperty(
        name="Length Exponent",
        description="Base-10 exponent for scaling lengths to Blender units",
        min=0, max=10, default=6,
        update=update_objects)

    planet_size_mult: IntProperty(
        name="Planet Size Multiplier",
        description=("Increase the size of planet surfaces "
                     "to make them better visible"),
        min=1, soft_max=1000, default=1,
        update=update_objects)

    draw_orbit: BoolProperty(
        name="Draw Orbit",
        description="Activate/deactivate drawing orbits in the viewport",
        default=False,
        update=update_handle)

    # ---------------- Read-only properties ----------------

    # time_mult:
    # multiplier for converting animation seconds to simulation seconds
    def _get_time_mult(self):
        return 10 ** self.time_exp

    time_mult: FloatProperty(
        get=_get_time_mult)

    # length_mult:
    # multiplier for converting Blender Units to km using length_exp
    def _get_length_mult(self):
        return 10 ** self.length_exp

    length_mult: FloatProperty(
        get=_get_length_mult)

    # time:
    # simulation time in seconds
    def _get_time(self):
        return self.sim_time * self.time_mult

    time: FloatProperty(
        get=_get_time)


# =============================================================================
# Object data
# =============================================================================

# SSSIM object settings
class SSSIMObject(bpy.types.PropertyGroup):
    # ------------------------------------------------
    # Update methods:
    # ------------------------------------------------

    # There are no specific methods for obj.sssim_obj,
    # but we use update_relations a lot

    # ------------------------------------------------
    # Properties:
    # ------------------------------------------------
    use_sssim: BoolProperty(
        name="Solar System Simulator",
        description="Use the Solar System Simulator",
        default=False,
        update=update_relations)

    # object panel
    object_type: EnumProperty(
        name="Type",
        items=(('CENTER', "Center", "Solar system center, has no orbit."),
               ('PLANET', "Planet", ("Planets orbit around a center or other "
                                     "planets.")),
               ('SURFACE', "Surface", "Surface of planets or centers."),
               ),
        default='CENTER',
        description="The type of the simulated object",
        update=update_relations)
    show_info: BoolProperty(
        name="Show Info",
        description="Show additional info",
        default=True)
    mass_mantissa: FloatProperty(
        name="Mass (kg)",
        description="Mass of the object (in kilogram)",
        soft_min=0.001, soft_max=1000, default=1)
    mass_exp: IntProperty(
        name="Mass Exponent",
        description="Base-10 exponent of the mass",
        soft_min=0, soft_max=50, default=24)

    # ---------------- Read-only properties ----------------

    # The mass of the object in kg
    def _mass_calc(self):
        return self.mass_mantissa * 10 ** self.mass_exp

    mass: FloatProperty(
        name="Mass",
        get=_mass_calc)

    # ----------- Helpful Functions -----------
    # nothing


# for objects with sssim_obj.object_type == 'PLANET'
class SSSIMOrbit(bpy.types.PropertyGroup):
    # ------------------------------------------------
    # Properties:
    # ------------------------------------------------
    center_name: StringProperty(
        name="Center Object",
        description="Name of the object at the orbit's focal point",
        update=update_relations)

    # shape and size
    eccentricity: FloatProperty(
        name="Eccentricity",
        description="The ellipticalness of the orbit, 0 equals a circle",
        min=0, max=0.999, default=0)
    distance_mantissa: FloatProperty(
        name="Semi-Major Axis (km)",
        description=("The distance to the center, in a circle: radius, "
                     "in an ellipse: half of the longest axis"),
        soft_min=0.001, soft_max=1000, default=1)
    distance_exp: IntProperty(
        name="Length Exponent",
        description="Base-10 exponent of the distance",
        soft_min=0, soft_max=12, default=6)

    # orientation of the orbital plane
    inclination: FloatProperty(
        name="Inclination",
        description="Tilt against the global xy-plane",
        min=-pi / 2, max=3 / 2 * pi,
        unit='ROTATION')
    asc_node: FloatProperty(
        name="Ascending Node",
        description=("Longitude of the ascending node "
                     "(= rotates of the orbit around the global z-axis)"),
        min=0, max=2 * pi,
        unit='ROTATION')
    arg_periapsis: FloatProperty(
        name="Argument of Periapsis",
        description=("Argument of periapsis (= rotates the nearest point in "
                     "the orbit around the local z-axis)"),
        min=0, max=2 * pi,
        unit='ROTATION')

    # position
    time_offset: FloatProperty(
        name="Time Offset",
        description=("Timeshift the starting point by the given fraction "
                     "of the orbital period"),
        min=0, max=1)

    # orbital period specified by user
    override_period: BoolProperty(
        name="Override Period",
        description=("Use a custom orbital period, otherwise gravity will "
                     "determine the period"),
        default=False)
    use_frames: BoolProperty(
        name="Use Frames",
        description="Input orbital period in frames",
        default=True)
    override_period_seconds: FloatProperty(
        name="Orbital Period (seconds)",
        description="Custom orbital period in seconds",
        min=0, default=10000000)
    override_period_frames: FloatProperty(
        name="Orbital Period (frames)",
        description="Custom orbital period in frames",
        min=0, default=50)

    # ---------------- Read-only properties ----------------

    # get the center object
    def _get_center_object(self):
        is_planet = self.id_data.sssim_obj.object_type == 'PLANET'
        if is_planet and self.center_name != "":
            return bpy.data.objects.get(self.center_name)
        return None

    center_object = property(_get_center_object)

    # recursive search for cyclic dependencies of center objects
    def center_recursive_control(self, ref_obj, ind=0):
        center = self.center_object
        self._print_verbose("Checking  {}".format(center), ind)
        if center is None:
            self._print_verbose("object is None => invalid", ind)
            return False

        if center.sssim_obj.object_type == 'CENTER':
            self._print_verbose("object is 'CENTER' => valid", ind)
            # no further center object -> no cycle possible
            return True

        if center.sssim_obj.object_type == 'SURFACE':
            # center cannot be surface object
            msg = "Object {} is type 'SURFACE' => invalid center"
            self._print_verbose(msg.format(center.name), ind)
            return False

        self._print_verbose("object is type 'PLANET'", ind)
        if center.name == ref_obj.name:
            self._print_verbose("Detected Cycle!", ind)
            # cycle detected
            return False

        self._print_verbose("This object is valid, continue", ind)
        # no link back to ref_obj -> no cycle -> go on
        center_control = center.sssim_orbit.center_recursive_control
        return center_control(ref_obj, ind + 1)

    # is the center valid
    def _center_valid(self):
        return self.center_recursive_control(self.id_data)

    center_is_valid: BoolProperty(get=_center_valid)

    # semi-major axis = half the longest axis of an ellipse
    def _semi_major_axis_calc(self):
        return self.distance_mantissa * 10 ** self.distance_exp

    semi_major_axis: FloatProperty(
        get=_semi_major_axis_calc)

    # semi-minor axis = helf the shortest distance
    def _semi_minor_axis_calc(self):
        return self.semi_major_axis * sqrt(1 - self.eccentricity ** 2)

    semi_minor_axis: FloatProperty(
        get=_semi_minor_axis_calc)

    # periapsis = nearest point of orbit
    def _periapsis_calc(self):
        return (1 - self.eccentricity) * self.semi_major_axis

    periapsis: FloatProperty(
        get=_periapsis_calc)

    # apoapsis = furthest point of orbit
    def _apoapsis_calc(self):
        return (1 + self.eccentricity) * self.semi_major_axis

    apoapsis: FloatProperty(
        get=_apoapsis_calc)

    # orbital period (in seconds) = time to complete one orbit
    def _orbital_period_calc(self):
        obj = self.id_data
        center_object = obj.sssim_orbit.center_object
        if not center_object:
            return 0  # no center, no orbit

        if self.override_period:
            if self.use_frames:
                scn = bpy.context.scene
                fps = scn.render.fps / scn.render.fps_base
                return self.override_period_frames * scn.sssim_scn.time_mult / fps
            else:
                return self.override_period_seconds

        axis = self.semi_major_axis
        mass = obj.sssim_obj.mass
        cmass = center_object.sssim_obj.mass
        return calculation.orbital_period(axis, mass, cmass)

    orbital_period: FloatProperty(
        get=_orbital_period_calc)

    # orbital period in frames, just for the UI
    def _orbital_period_frames_calc(self):
        scn = bpy.context.scene
        fps = scn.render.fps / scn.render.fps_base
        return self.orbital_period / scn.sssim_scn.time_mult * fps

    orbital_period_frames: FloatProperty(
        get=_orbital_period_frames_calc)

    # true anomaly, angle of the current position since periapsis
    def _true_anomaly_calc(self):
        time = bpy.context.scene.sssim_scn.time
        period = self.orbital_period
        ecc = self.eccentricity
        time_offset = self.time_offset
        return calculation.true_anomaly(time, period, ecc, time_offset)

    true_anomaly: FloatProperty(
        get=_true_anomaly_calc)

    # ----------- Helpful Functions -----------
    # for printing information, depending on verbosity level,
    # is used only in center_recursive_control
    verbose: BoolProperty(default=False)

    def _print_verbose(self, text, indent=0):
        if indent > 0:
            ind = "    " * indent
        else:
            ind = ""

        if self.verbose:
            print(ind, text)


class SSSIMRotation(bpy.types.PropertyGroup):
    # ------------------------------------------------
    # Properties:
    # ------------------------------------------------
    use_rotation: BoolProperty(
        name="Use Rotation",
        description="Enable rotation of the object around its axis",
        default=False,
        update=update_relations)
    use_frames: BoolProperty(
        name="Use Frames",
        description="Input rotation period in frames instead of seconds",
        default=False,
        update=update_relations)
    period_frames: FloatProperty(
        name="Rotation Period (frames)",
        description="Time in frames for one rotation",
        min=0, default=100)
    period_seconds: FloatProperty(
        name="Rotation Period (seconds)",
        description="Time in seconds for one rotation",
        min=0, default=86400)
    axis_tilt: FloatProperty(
        name="Axis Tilt",
        description="Tilt of the rotation axis",
        min=0, max=pi, default=0,
        unit='ROTATION')
    axis_direction: FloatProperty(
        name="Tilt Direction",
        description="Direction of tilting",
        min=0, max=2 * pi, default=0,
        unit='ROTATION')
    relative_to_orbit: BoolProperty(
        name="Relative to Orbit",
        description="Orient the rotation axis relative to the orbital plane",
        default=True)

    # ---------------- Read-only properties ----------------

    # the planet's rotation period in simulation seconds
    def _rot_period_get(self):
        if not self.use_rotation:
            return 0  # no rotation

        if self.use_frames:
            scn = bpy.context.scene
            fps = scn.render.fps / scn.render.fps_base
            return self.period_frames * scn.sssim_scn.time_mult / fps

        return self.period_seconds

    rotation_period: FloatProperty(
        get=_rot_period_get)

    # rotation period in frames, just for the UI
    def _rot_period_frames_calc(self):
        if self.use_rotation and self.use_frames:
            return self.period_frames

        scn = bpy.context.scene
        return self.rotation_period / scn.sssim_scn.time_mult

    rotation_period_frames: FloatProperty(
        get=_rot_period_frames_calc)


class SSSIMCalculation(bpy.types.PropertyGroup):
    # ------------------------------------------------
    # Update methods:
    # ------------------------------------------------
    def check_calc_frames(self, context):
        # end always later than start
        if self.frame_end <= self.frame_start:
            self.frame_end = self.frame_start + 1

    # ------------------------------------------------
    # Properties:
    # ------------------------------------------------
    use_driver: BoolProperty(
        name="Use Driver",
        description="Use driver for calculation of orbit and rotation",
        default=True,
        update=update_relations)
    cyclic: BoolProperty(
        name="Cyclic",
        description="Add cyclic F-Curve modifier to repeat the motion",
        default=True)
    frame_start: IntProperty(
        name="Start Frame",
        description="First frame for the calculation",
        default=1,
        update=check_calc_frames)
    frame_end: IntProperty(
        name="End Frame",
        description="Last frame for the calculation",
        default=100,
        update=check_calc_frames)
    frame_step: IntProperty(
        name="Frame Step",
        description="Size of timestep between calculated positions",
        min=1, default=20)

    # ---------------- Read-only properties ----------------

    # the real start frame for baking
    def _determine_frame_start(self):
        if self.cyclic:
            return 1
        return self.frame_start

    real_frame_start: IntProperty(
        get=_determine_frame_start)

    # the real end frame
    def _determine_frame_end(self):
        if self.cyclic:
            simorbit = self.id_data.sssim_orbit
            # end = start + orbital_period
            return 1 + simorbit.orbital_period_frames

        return self.frame_end

    real_frame_end: FloatProperty(
        get=_determine_frame_end)

    # Does the object have an orbit?
    def _calc_orbit(self):
        obj = self.id_data
        active = obj.sssim_obj.use_sssim
        valid_type = obj.sssim_obj.object_type == 'PLANET'
        can_orbit = obj.sssim_orbit.center_is_valid
        return active and valid_type and can_orbit

    calc_orbit: BoolProperty(
        get=_calc_orbit)

    # Can we calculate the rotation?
    def _calc_rotation(self):
        obj = self.id_data
        active = obj.sssim_obj.use_sssim
        valid_type = obj.sssim_obj.object_type in ('CENTER', 'PLANET')
        can_rotate = obj.sssim_rotation.use_rotation
        return active and valid_type and can_rotate

    calc_rotation: BoolProperty(
        get=_calc_rotation)


# for objects with sssim_obj.object_type == 'SURFACE'
class SSSIMSurface(bpy.types.PropertyGroup):
    # ------------------------------------------------
    # Update methods:
    # ------------------------------------------------
    def update_parent(self, context):
        obj = self.id_data
        if self.parent_is_valid:
            parent = obj.sssim_surface.parent_object
            msg_sur = operators.add_surface(obj, parent)
        else:
            msg_sur = operators.remove_surface(obj)

        if verbose:
            print("    ", msg_sur)

    def update_size(self, context):
        if self.calc_size:
            # dimensions = proportions * diameter, diameter = radius * 2
            dim = Vector(self.object_proportions) * self.radius_bu * 2
            self.id_data.dimensions = dim

    def reset_size(self, context):
        obj = self.id_data
        if self.calc_size and obj.sssim_obj.object_type == 'SURFACE':
            self.original_scale = obj.scale
            if obj.dimensions.length != 0:
                self.object_proportions = obj.dimensions / max(obj.dimensions)
            else:
                # we can't scale this object because all components are 0
                # but we have to update object_proportions
                self.object_proportions = obj.dimensions
            self.update_size(context)
        else:
            obj.scale = self.original_scale

    # ------------------------------------------------
    # Properties:
    # ------------------------------------------------
    parent_name: StringProperty(
        name="Parent Object",
        description="Name of the object as the surface's parent",
        update=update_parent)

    calc_size: BoolProperty(
        name="Calculate Size",
        description=("Calculate the size of the object instead of "
                     "using the current size"),
        default=False,
        update=reset_size)

    # used for scaling to maintain the original proportions
    # => none-uniform objects are allowed
    object_proportions: bpy.props.FloatVectorProperty()
    # to reset the scale
    original_scale: bpy.props.FloatVectorProperty()

    radius_mantissa: FloatProperty(
        name="Radius (in km)",
        description="The radius of the surface",
        soft_min=0.001, soft_max=1000, default=1,
        update=update_size)
    radius_exp: IntProperty(
        name="Exponent",
        description="Base-10 exponent of the radius",
        soft_min=0, soft_max=8, default=3,
        update=update_size)

    # ---------------- Read-only properties ----------------

    # return the parent object
    def _get_parent_object(self):
        is_surface = self.id_data.sssim_obj.object_type == 'SURFACE'
        if is_surface and self.parent_name != "":
            return bpy.data.objects.get(self.parent_name)
        return None

    parent_object = property(_get_parent_object)

    # is the parent object valid (= simulation object with right type)
    def _parent_valid(self):
        par = self.parent_object
        if par is not None:
            # so the object exists, has it the right type?
            simobj = par.sssim_obj
            if simobj.use_sssim and simobj.object_type != 'SURFACE':
                return True
        return False

    parent_is_valid: BoolProperty(
        get=_parent_valid)

    # radius in km
    def _get_radius_km(self):
        if self.calc_size:
            return self.radius_mantissa * 10 ** self.radius_exp

        # use the current dimensions to get the radius in km
        sssim_scn = bpy.context.scene.sssim_scn
        mult = sssim_scn.length_mult / sssim_scn.planet_size_mult
        return max(self.id_data.dimensions) * mult / 2

    radius_km: FloatProperty(
        get=_get_radius_km)

    # radius in Blender units from radius_km
    def _get_radius_bu(self):
        sssim_scn = bpy.context.scene.sssim_scn
        div = sssim_scn.length_mult / sssim_scn.planet_size_mult
        return self.radius_km / div

    radius_bu: FloatProperty(
        get=_get_radius_bu)

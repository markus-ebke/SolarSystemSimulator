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

from math import pi, sin, cos, acos, sqrt
import bpy
from mathutils import Vector, Euler

G = 6.673e-11  # Gravitational constant


def orbital_period(axis, mass1, mass2):
    """Calc the orbital period

    axis = length of the semi-major-axis in kilometers
    mass1, mass2 = Mass of the two bodies in kg
    return the orbital period in seconds.
    """
    axis *= 1000  # meter = kilometer * 1000

    return sqrt((4 * pi ** 2 * axis ** 3) / (G * (mass1 + mass2)))


def true_anomaly(time, period, ecc=0, time_offset=0):
    """Calc the true anomaly (the angle since the nearest point of the orbit)

    time = time in seconds
    period = Orbital Period in seconds
    ecc = eccentricity of the orbit, ecc = 0: circular, 0 < ecc < 1: elliptical
    time_offset = offset of the initial angle, in range 0 - 1,
                  offset * period = starting time in seconds if time=0
    return angle between 0 and 2pi.
    """
    if period == 0:
        # planet cannot have an orbit
        # => return static position dependent on offset
        return 2 * pi * time_offset

    # mean anomaly = "angle" travelled, in range 0 - 2*pi
    mean_anomaly = 2 * pi * (((time / period) + time_offset) % 1)

    if ecc == 0:
        # circle, no fancy stuff here
        return mean_anomaly

    # the true anomaly is the "angle" travelled (as seen from center)
    # since the nearest point (periapsis) of the orbit

    # mean_anomaly (in rad)-> eccentric anomaly (rad)-> true anomaly (rad)
    # calculate the eccentric anomaly according to a formular found at
    # http://www-spof.gsfc.nasa.gov/stargaze/Smotion.htm
    eccentric_anomaly = mean_anomaly
    ecc_old = 0
    while abs(ecc_old - eccentric_anomaly) > 1e-10:  # precision = 10 digits
        ecc_old = eccentric_anomaly
        eccentric_anomaly = mean_anomaly + ecc * sin(ecc_old)

    # calculate the true anomaly
    # formular from http://en.wikipedia.org/wiki/True_anomaly
    t_0 = cos(eccentric_anomaly) - ecc
    t_1 = 1 - ecc * cos(eccentric_anomaly)
    angle = acos(t_0 / t_1)  # true anomaly in range: 0 - pi

    if mean_anomaly > pi:
        # we are on the other side (angle > 180) of the circle
        angle = 2 * pi - angle  # range: 0 - 2*pi

    return angle


def orbit_position(simorbit, simscn, time):
    """Calculate the position of the planet in the xy-plane"""
    # sma = (to BU adjusted) semi_major_axis
    sma = simorbit.semi_major_axis / simscn.length_mult

    # eccentricity = more circular (small e) or more elliptical (bigger e)
    ecc = simorbit.eccentricity

    # orbital period = time to complete one orbit
    period = simorbit.orbital_period

    time_offset = simorbit.time_offset

    # calculate the position
    theta = true_anomaly(time, period, ecc, time_offset)
    if ecc == 0:
        # circular orbit, distance is constant
        dist = sma
    else:
        # elliptical orbit, dist = distance from sun (in BU)
        # from http://www-spof.gsfc.nasa.gov/stargaze/Smotion.htm
        dist = sma * (1 - ecc ** 2) / (1 + ecc * cos(theta))

    pos_x = cos(theta) * dist
    pos_y = sin(theta) * dist
    pos_z = 0
    position = Vector((pos_x, pos_y, pos_z))

    return position


def orbit_orientation(other, simorbit):
    """Rotates a mathutils value (Vector, Euler, ...) using orbital elements.

    simorbit contains the needed information,
    return a (rotated) Vector or Euler.
    """
    inc = simorbit.inclination
    asc_node = simorbit.asc_node
    arg_periapsis = simorbit.arg_periapsis

    # rotate mathutils value around global z axis (argument of periapsis)
    other.rotate(Euler((0, 0, arg_periapsis)))

    # rotate around the global x axis by (inclination)
    other.rotate(Euler((inc, 0, 0)))

    # rotate around global z again (rotate the ascending node)
    other.rotate(Euler((0, 0, asc_node)))

    return other


# =============================================================================
# Driver functions
# =============================================================================
def eval_planet_orbit(obj, scn_name, index=None, time=None):
    """Evaluate the planets position, used by driver.

    obj = object to be simulated
    sssim_scn = scene data
    index = index of the location channel
    time = time when to calculate, if not given use current scene time
            time is in seconds of the simulation
    returns a 3-tuple with xyz-locations or, if index given, only one component
    """
    simorbit = obj.sssim_orbit
    scn = bpy.data.scenes.get(scn_name)
    if not scn:
        errmsg = "DRIVER ERROR: Invalid scene name {}"
        print(errmsg.format(scn_name))
        return 0
    sssim_scn = scn.sssim_scn

    # time = time for which to calculate the planet position, in seconds
    if time is None:
        time = sssim_scn.time

    planet_loc = orbit_position(simorbit, sssim_scn, time)
    orbit_rot = orbit_orientation(planet_loc, simorbit)

    return orbit_rot if index is None else orbit_rot[index]


def eval_planet_rotation(obj, scn_name, index=None, time=None):
    """Evaluate the planets rotation, used by driver.

    scn_name = Name of a scene which contains the object
    obj_name = Name of the object to simulate
    index = index of the rotation channel,
            usually only z-axis (index=2) changes
    time = time when to calculate, if not given use current scene time
            time is in seconds of the simulation
    returns an Euler in mode ZYX or, if index given, an angle in radians
    """
    simrot = obj.sssim_rotation
    scn = bpy.data.scenes.get(scn_name)
    if not scn:
        errmsg = "DRIVER ERROR: Invalid scene name {}"
        print(errmsg.format(scn_name))
        return 0
    sssim_scn = scn.sssim_scn

    # time = time in seconds, if None use current scene time
    if time is None:
        time = sssim_scn.time

    # rotation_period is also in seconds
    rotation_period = simrot.rotation_period
    if rotation_period != 0:
        rot_z = 2 * pi * time / rotation_period
    else:
        # invalid input -> no rotation
        rot_z = 0

    tilt = simrot.axis_tilt
    planet_rot = Euler((tilt, 0.0, 0.0), 'ZYX')  # note that mode is 'ZYX'

    # rotate around global (not local) z-axis
    direction = simrot.axis_direction
    planet_rot.rotate(Euler((0.0, 0.0, direction), 'XYZ'))

    # rotate around local z-axis
    # NOTE: we won't use planet_rot.rotate_axis('Z', rot_z) because then
    # all rotations are between -180 and 180 and for the rotation around
    # z we need a continous motion with increasing z values
    planet_rot.z += rot_z

    if simrot.relative_to_orbit and obj.sssim_obj.object_type == 'PLANET':
        planet_rot = orbit_orientation(planet_rot, obj.sssim_orbit)

    return planet_rot if index is None else planet_rot[index]

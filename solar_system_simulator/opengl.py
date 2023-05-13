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

"""
OpenGL Drawing code for viewing the orbit and other stuff in the 3DView.
draw_3d_callback is the function used as SpaceView3D draw handler.
"""

from math import degrees, cos

import blf
import gpu
from gpu_extras.batch import batch_for_shader
from mathutils import Vector, Euler
from bpy_extras.view3d_utils import location_3d_to_region_2d

from .calculation import eval_planet_orbit, orbit_orientation

shader = None
scn = None  # scene
vl = None  # view layer

# number of points used to draw the orbit and an arc
NUM_ORBIT = 128
NUM_ARC = 32

# elements that should be drawn on the 2D screen (not the 3D view) go here
DRAW_2D = []

# Only one function (to_screen_coord) need these variables, I don't wanted
# to use them as parameter for a lot of these functions, so I made them global
# they are set when we call draw_3d_callback
region = None
rv3d = None


def draw_3d_callback(self, context):
    """Draw inside the 3D Viewport"""
    # set global variables
    global scn, region, rv3d, vl, shader
    scn = context.scene
    vl = context.view_layer
    region = context.region
    rv3d = context.space_data.region_3d

    # drawing shader
    shader = gpu.shader.from_builtin('3D_UNIFORM_COLOR')
    shader.bind()
    gpu.state.blend_set("ALPHA")

    # default settings
    default_point_size = 1  # there is not gpu.state.point_size_get()
    default_line_width = gpu.state.line_width_get()

    # draw info for every sssim object
    for obj in vl.objects:
        if obj.sssim_obj.use_sssim and obj.visible_get():
            if obj.sssim_obj.object_type == 'PLANET':
                show_orbit(obj)
            if obj.sssim_obj.object_type == 'CENTER':
                show_center_info(obj)

    # restore opengl defaults
    gpu.state.point_size_set(default_point_size)
    gpu.state.line_width_set(default_line_width)
    gpu.state.blend_set("NONE")


def draw_2d_callback(self, context):
    """Draw on the 2D screen of the 3D Viewport"""
    while DRAW_2D:
        fct, args = DRAW_2D.pop()
        fct(*args)


# =============================================================================
# These functions say what should be drawn
# =============================================================================

def show_orbit(obj):
    # BLENDER PYTHON AND THE HOLY OPENGL [PART I]
    #
    # GALAHAD: There it is!
    # ARTHUR: The show_orbit function!
    # ROBIN: Oh, great.

    simorbit = obj.sssim_orbit
    center = simorbit.center_object
    center_loc = real_location(center)

    # we want to draw more for the selected active object
    is_active = obj.select_get() and vl.objects.active == obj

    # mark the current position with a point and the object's name
    pos = real_location(obj)
    draw_point(pos)
    draw_text(obj.name, pos)

    # draw the orbit for every planet
    points = []
    period = simorbit.orbital_period
    for i in range(NUM_ORBIT):
        # get position relative to center
        pos = orbit_position(obj, time=period * i / NUM_ORBIT)
        points.append(center_loc + pos)

    # color of the orbit, active orbit is draw orange
    if is_active:
        # that is the color active objects are usually shown in
        col = (1.0, 0.66, 0.25, 1.0)
    else:
        col = (0.5, 0.5, 0.7, 1.0)

    # [PART II]
    # ARTHUR: Look! There's the old draw_line function from scene 24!
    # BEDEVERE: What is he doing here?
    # ARTHUR: He is the keeper of the Loop Drawing.
    #         He asks each traveller five parameters--
    # GALAHAD: Three parameters.
    # ARTHUR: Three parameters. He who answers the five parameters--
    # GALAHAD: Three parameters.
    # ARTHUR: Three parameters may draw in safety.
    # ROBIN: What if you get a parameter wrong?
    # ARTHUR: Then you are cast into the Gorge of Eternal Artifacts.
    # ROBIN: Oh, I won't go.
    # GALAHAD: Who's going to answer the parameters?
    draw_line(points, width=2, color=col, closed=True)

    if is_active:
        show_orbit_extra_info(obj)
        show_rotation_angles(obj)


def show_orbit_extra_info(obj):
    simorbit = obj.sssim_orbit
    center = simorbit.center_object
    center_loc = real_location(center)
    period = simorbit.orbital_period

    # compensate for the time offset
    periapsis_time = -simorbit.time_offset * period

    if period != 0:
        show_true_anomaly(obj, center, periapsis_time, period)

    if simorbit.eccentricity > 0:
        # show periapsis with distance
        peri_pos = center_loc + orbit_position(obj, periapsis_time)
        draw_point(peri_pos, size=10)
        draw_line((center_loc, peri_pos), color=(1.0, 1.0, 0.5, 1.0))
        text = "Periapsis: {:,} km".format(round(simorbit.periapsis))
        draw_text(text, peri_pos)

        # show apoapsis with distance
        time = periapsis_time + period / 2
        apo_pos = center_loc + orbit_position(obj, time)
        draw_point(apo_pos, size=10)
        draw_line((center_loc, apo_pos), color=(1.0, 0.8, 0.5, 1.0))
        text = "Apoapsis: {:,} km".format(round(simorbit.apoapsis))
        draw_text(text, apo_pos)
    else:
        # draw radius text between center and starting position
        pos = orbit_position(obj, time=0)
        draw_line((center_loc, center_loc + pos))
        text = "Radius: {:,} km".format(round(simorbit.semi_major_axis))
        draw_text(text, center_loc + 0.8 * pos)

    axis_bu = simorbit.semi_major_axis / scn.sssim_scn.length_mult
    angle_radius = axis_bu / 3
    show_orbit_angles(simorbit, angle_radius, center_loc)


def show_true_anomaly(obj, center, periapsis_time, period):
    """Draw the true anomaly = the angle from the periapsis

    Note: This is not perfect, we get some overlap when time_offset > 0
    """
    simorbit = obj.sssim_orbit
    center_loc = real_location(center)

    points = []
    # time of this period = current time - time for full periods
    time_since_periapsis = (scn.sssim_scn.time - periapsis_time) % period
    for i in range(NUM_ORBIT + 1):
        time = periapsis_time + time_since_periapsis * i / NUM_ORBIT
        position = orbit_position(obj, time)
        points.append(center_loc + position)

    draw_sector(center_loc, points, color=(0.0, 0.3, 1.0, 0.2))

    if simorbit.true_anomaly != 0:
        # draw text in the circular sector
        angle = degrees(simorbit.true_anomaly)
        text = "True anomaly: {}°".format(round(angle, 2))
        # the text position should be inside the arc
        # between the center and the middle point of the arc
        midarc = orbit_position(obj, periapsis_time + time_since_periapsis / 2)
        text_pos = center_loc + midarc / 2

        draw_text(text, text_pos, align="CENTER")


def show_orbit_angles(simorbit, radius, center_loc):
    """Draw the Inclination, Ascending Node, Argument of Periapsis"""
    inc = simorbit.inclination
    asc_node = simorbit.asc_node
    arg_peri = simorbit.arg_periapsis

    if inc != 0:
        col = (1.0, 1.0, 0.0, 0.5)
        text = "Inclination: {}°"

        def inc_rotate(angle):
            vec = Vector((0, 1, 0))
            vec.rotate(Euler((angle, 0, 0)))
            vec.rotate(Euler((0, 0, asc_node)))
            return vec

        show_angle(inc, radius, inc_rotate, center_loc, col, text)

    if asc_node != 0:
        col = (0.0, 1.0, 0.5, 0.5)
        text = "Ascending Node: {}°"

        def asc_rotate(angle):
            vec = Vector((1, 0, 0))
            vec.rotate(Euler((0, 0, angle)))
            return vec

        show_angle(asc_node, radius, asc_rotate, center_loc, col, text)

    if arg_peri != 0:
        col = (0.0, 0.7, 1.0, 0.5)
        text = "Argument of Periapsis: {}°"

        def peri_rotate(angle):
            vec = Vector((1, 0, 0))
            vec.rotate(Euler((0, 0, angle)))
            vec.rotate(Euler((inc, 0, 0)))
            vec.rotate(Euler((0, 0, asc_node)))
            return vec

        show_angle(arg_peri, radius, peri_rotate, center_loc, col, text)


def show_rotation_angles(obj):
    simorbit = obj.sssim_orbit
    simrot = obj.sssim_rotation

    if obj.sssim_obj.object_type == 'PLANET':
        # when using planets, better scale with distance to center
        axis_bu = simorbit.semi_major_axis / scn.sssim_scn.length_mult
        radius = axis_bu / 5
    else:
        radius = 1

    # the center of the angle is the object location
    obj_loc = real_location(obj)

    # orbit_normal = axis with no tilt, normal of the orbital plane
    # is used as helper vector
    orbit_normal = Vector((0, 0, 1))
    if simrot.relative_to_orbit:
        orbit_orientation(orbit_normal, simorbit)

    tilt = simrot.axis_tilt
    direction = simrot.axis_direction

    # draw the tilt of the rotation axis
    if tilt != 0:
        col = (1.0, 0.5, 0.0, 0.5)
        text = "Tilt: {}°"

        def tilt_rotate(angle):
            vec = Vector((0, 0, 1))
            vec.rotate(Euler((angle, 0, 0), 'ZYX'))
            vec.rotate(Euler((0, 0, direction)))
            if simrot.relative_to_orbit:
                vec = orbit_orientation(vec, simorbit)
            return vec

        show_angle(tilt, radius, tilt_rotate, obj_loc, col, text)

        # draw the orbit normal
        # the normal is drawn in blue, just like the z-axis
        col = (0.0, 0.0, 1.0, 1.0)
        up = obj_loc + orbit_normal * radius
        draw_line((obj_loc, up), color=col)

        # draw the Rotation Axis
        axis_vec = tilt_rotate(tilt)
        draw_line((obj_loc, obj_loc + axis_vec * radius))
        draw_text("Rotation Axis", obj_loc + axis_vec * radius)

        # draw the direction where the rotation axis is pointing
        if direction != 0:
            col = (1.0, 0.1, 0.0, 0.5)
            text = "Direction: {}°"

            def direction_rotate(angle):
                vec = Vector((0, 0, 1))
                vec.rotate(Euler((tilt, 0, 0), 'ZYX'))
                vec.rotate(Euler((0, 0, angle)))
                if simrot.relative_to_orbit:
                    vec = orbit_orientation(vec, simorbit)
                return vec - orbit_normal * cos(tilt)

            # we want a disk section, not a cone section so the center of
            # the angle must be translated in direction of the orbit_normal
            hub = obj_loc  # + orbit_normal * cos(tilt) * radius
            show_angle(direction, radius, direction_rotate, hub, col, text)


def show_center_info(obj):
    draw_point(obj.location, size=10)
    draw_text(obj.name, obj.location)

    # we want to draw more for the selected active object
    is_active = obj.select_get() and vl.objects.active == obj

    if is_active:
        show_rotation_angles(obj)


# =============================================================================
# Helpful functions
# =============================================================================

def real_location(obj, time=None):
    """Position of obj at time (in seconds) relative to origin."""
    if obj.sssim_obj.object_type != 'PLANET':  # that's easy
        return obj.location

    # So, the object is orbiting around another object,
    # what is the location of the center object?
    center = obj.sssim_orbit.center_object

    if center is None:
        center_loc = Vector((0, 0, 0))
    else:
        center_loc = real_location(center, time)

    return center_loc + orbit_position(obj, time)


def orbit_position(obj, time=None):
    """Position of an orbiting object relative to center.

    This wraps eval_planet_orbit in a more convenient function.
    """
    return eval_planet_orbit(obj, scn.name, index=None, time=time)


def show_angle(angle, radius, rotate_func, center_loc, color, text):
    """Draw an angle defined by rotate_func with text around center_loc.

    'angle' is an angle in radians.
    'radius' is the radius of the drawn angle.
    'rotate_func' is a function with one parameter, an angle, which returns a
    unit vector.
    Example:
        def rot_func(angle):
            vec = Vector((0, 0, 1))
            vec.rotate(Euler((angle, 0, 0)))
            return vec

    'center_loc' is a 3D-Vector, the location around which the angle is drawn.
    'color' is the color of the angle.
    'text' is a string written inside the angle, which can be formated with
    the angle (will be converted to degrees).
    For example text="Angle: {}°" will be shown as "Angle: 42°".
    """
    arc_points = []
    for i in range(NUM_ARC):
        part_angle = angle * i / (NUM_ARC - 1)
        vec = rotate_func(part_angle)
        arc_points.append(center_loc + vec * radius)

    # draw the sector
    draw_sector(center_loc, arc_points, color=color)

    # draw text
    angle_deg = degrees(angle)
    text = text.format(round(angle_deg, 2))
    midarc = arc_points[NUM_ARC // 2 - 1]
    text_pos = 0.3 * center_loc + 0.7 * midarc
    draw_text(text, text_pos, align='CENTER')


# =============================================================================
# These functions actually do the drawing
# =============================================================================

def draw_point(position, size=5, color=(1.0, 0.2, 0.2, 1.0)):
    """Draw single red point"""

    # bgl.glEnable(bgl.GL_POINT_SMOOTH)  # TODO fix
    gpu.state.point_size_set(size)
    shader.uniform_float("color", color)
    batch = batch_for_shader(shader, 'POINTS', {"pos": [position]})
    batch.draw(shader)

    # bgl.glDisable(bgl.GL_POINT_SMOOTH)  # TODO fix


def draw_line(points, width=2, color=(1.0, 1.0, 1.0, 1.0), closed=False):
    """Draw polyline, points = [(x1, y1, z1), (x2, y2, z2),...]"""
    # [PART V]
    # BRIDGEKEEPER: Stop! What... is your name?
    # GALAHAD: 'Sir Galahad of Camelot'.
    # BRIDGEKEEPER: What... is your line width?
    # GALAHAD: width.
    gpu.state.line_width_set(width)

    # BRIDGEKEEPER: What... is your favorite color?
    # GALAHAD: bgl.glColor4f(*color). No, shader.unif-- auuuuuuuugh!
    shader.uniform_float("color", color)

    mode = 'LINE_LOOP' if closed else 'LINE_STRIP'
    batch = batch_for_shader(shader, mode, {"pos": points})
    batch.draw(shader)


def draw_sector(center_point, arc_points, color=(1.0, 1.0, 1.0, 0.2)):
    """Draw a circular sector along arc_points around center_point"""
    shader.uniform_float("color", color)

    # The center_points is the hub of the triangle fan.
    # Don't forget it, it's not just a flesh wound.
    vertices = [center_point] + arc_points
    indices = [(0, idx, idx+1) for idx in range(len(arc_points))]
    batch = batch_for_shader(shader, 'TRIS', {"pos": vertices},
                             indices=indices)

    batch.draw(shader)


def draw_text(text, position, align="LEFT"):
    pos_2d = location_3d_to_region_2d(region, rv3d, position)
    if pos_2d is None:
        # we cannot draw because the point is not on the screen
        return

    # draw later when we are in POST_PIXEL mode
    # XXX this is ugly, but I don't know how to fix this 3D-2D mashing
    args = (text, pos_2d, align)
    DRAW_2D.append((_draw_text_2d, args))


def _draw_text_2d(text, position, align="LEFT"):
    """Draw text at given screen position aligned to the left"""

    # [PART VI]
    # BRIDGEKEEPER: Hee hee heh. Stop! What... is your name?
    # ARTHUR: It is 'Arthur', King of the Britons.
    # BRIDGEKEEPER: What... is your quest?
    # ARTHUR: To draw text on the screen.
    font_id = 0

    if align == "CENTER":
        # get length of text, place text .5 of length to the left of the coord.
        text_width, text_height = blf.dimensions(font_id, text)
        position.x -= text_width / 2

    blf.position(font_id, position.x, position.y, 0)
    blf.size(font_id, 12, 72)

    # BRIDGEKEEPER: What... is the OpenGL command for drawing text?
    # ARTHUR: What do you mean? Blender's bfl module can do that for me!
    # BRIDGEKEEPER: Huh? I-- I didn't know that. Auuuuuuuugh!

    blf.draw(font_id, text)

    # BEDEVERE: How do know so much about blender modules?
    # ARTHUR: Well, you have to know these things when you want
    #         to draw in Blender's 3DView, you know.

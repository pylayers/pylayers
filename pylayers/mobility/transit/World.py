from SimPy.SimulationRT import *
# from Tkinter import *
from math import pi, sin, cos
import os
from pylayers.mobility.transit.vec3 import vec3
import pdb
# This has to be 50% more than influence radius for the sub-tiling to work
# 1.71 from Pascal Stucki, Obstacles in Pedestrian Simulations
tile_size = 1.5
DIST_OBST = 2.0  # Distance from where wall have no influence on boids

_world = None


def world(**args):
    """ defining world

    Notes
    -----
   
    `Obstacles in pedestrian simulations <http://e-collection.library.ethz.ch/eserv/eth:27090/eth-27090-01.pdf>`_

    """

    global _world
    if _world is None:
        _world = World(**args)
    return _world


def near(boid, items, distance):
    tile_x, tile_y = boid.tile
    tile_distance = int(distance / tile_size)
    all_items = []
    for xx in range(tile_x - tile_distance, tile_x + tile_distance + 1):
        for yy in range(tile_y - tile_distance, tile_y + tile_distance + 1):
            all_items += items.get((xx, yy), [])
    return all_items


class World:
    """ Class World 
    
    Methods
    -------

    boids
    add_boid
    remove_boid
    update_boid
    obstacles
    add_wall

    """
    def __init__(self, **args):
        # self.tk = TkWorld(**args)
        self._boids = {}
        self._obstacles = {}
        self._zones = {}

    def boids(self, boid, distance=2):
        """
        Parameters
        ----------

        boid 
        distance : int 
            default 2

        Returns
        -------

        other_boids

        """

        other_boids = near(boid, self._boids, distance)
        other_boids.remove(boid)

        return other_boids

    def add_boid(self, boid):
        """

        Parameters
        ----------

        boid 

        """
        tile = (int(boid.position.x / tile_size), int(
            boid.position.y / tile_size))
        if tile not in self._boids:
            self._boids[tile] = [boid]
        else:
            self._boids[tile].append(boid)
        boid.tile = tile

    def remove_boid(self, boid):
        """
        Parameters
        ----------

        boid

        """
        self._boids[boid.tile].remove(boid)

    def update_boid(self, boid):
        tile = (int(boid.position.x / tile_size), int(
            boid.position.y / tile_size))
        if tile is not boid.tile:
            self._boids[boid.tile].remove(boid)
            self.add_boid(boid)

    def obstacles(self, boid):
        return near(boid, self._obstacles, DIST_OBST)

    def add_wall(self, *wall):
        the_obstacles = self._obstacles
        for ii in range(0, len(wall) - 1):
            line_start, line_end = wall[ii], wall[ii + 1]
            start_tile_x, start_tile_y = int(line_start[0] /
                                             tile_size), int(line_start[1] / tile_size)
            end_tile_x, end_tile_y = int(line_end[0] /
                                         tile_size), int(line_end[1] / tile_size)
            if end_tile_x < start_tile_x:
                end_tile_x, start_tile_x = start_tile_x, end_tile_x
            if end_tile_y < start_tile_y:
                end_tile_y, start_tile_y = start_tile_y, end_tile_y
            for xx in range(start_tile_x - 1, end_tile_x + 2):
                for yy in range(start_tile_y - 1, end_tile_y + 2):
                    tile = (xx, yy)
                    if tile in the_obstacles:
                        the_obstacles[tile].append((line_start, line_end))
                    else:
                        the_obstacles[tile] = [(line_start, line_end)]

    def zones(self, boid):
        """
        Parameters
        ----------

        boid 

        """
        tile = (int(boid.position.x / tile_size), int(
            boid.position.y / tile_size))
        return self._zones.get(tile, [])

    def add_zone(self, zone):
        """
        Parameters 
        ----------

        zone 

        """
        the_zones = self._zones
        start_tile_x, start_tile_y = int(zone.lower_left.x/ tile_size), int(zone.lower_left.y / tile_size)
        end_tile_x, end_tile_y = int(zone.upper_right.x / tile_size), int(zone.upper_right.y / tile_size)
        for xx in range(start_tile_x, end_tile_x + 1):
            for yy in range(start_tile_y, end_tile_y + 1):
                tile = (xx, yy)
                if tile in the_zones:
                    the_zones[tile].append(zone)
                else:
                    the_zones[tile] = [zone]

    def to_local(self, boid, point):
        xx, yy, unused_z = point - boid.position
        return vec3(boid.localy.y * xx - boid.localy.x * yy,
                    -boid.localx.y * xx + boid.localx.x * yy)

# 
# class TkWorld:
#     def __init__(self, **args):
#         defaults = {'width': 100, 'height': 200, 'geometry': '-20-90', 'scale': 2.5,
#                     'x_offset': 0, 'y_offset': 0, 'vectors': 0, 'collision_vectors': 0}
#         self.main = Tk()
#         for key, value in defaults.items():
#             if key.upper() in os.environ:
#                 setattr(self, key, float(os.environ[key.upper()]))
#             elif key in args:
#                 setattr(self, key, args[key])
#             else:
#                 setattr(self, key, value)
#         self.main.geometry(self.geometry)
#         self.canvas = Canvas(width=self.width * self.scale,
#                              height=self.height * self.scale)
#         self.canvas.pack(side='left')
#         # Tk has (0.0) in upper left, we need (0,0) in lower left
#         self.y_offset = self.height - self.y_offset
# 
#     def x_(self, xx):
#         return (xx + self.x_offset) * self.scale
# 
#     def y_(self, yy):
#         return (self.y_offset - yy) * self.scale
# 
#     def s_(self, ss):
#         return ss * self.scale

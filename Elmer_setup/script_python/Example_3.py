#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 14:02:02 2018

@author: Schmulius
"""


# Standard library imports
from os.path import join, dirname

# Enthought library imports
from mayavi import *
from mayavi.sources.vtk_xml_file_reader import VTKXMLFileReader
from mayavi.modules.outline import Outline
from mayavi.modules.glyph import Glyph
from mayavi.modules.vector_cut_plane import VectorCutPlane
from mayavi.modules.vectors import Vectors
from mayavi.filters.mask_points import MaskPoints


def glyph():
    """The script itself.  We needn't have defined a function but
    having a function makes this more reusable.
    """
    # 'mayavi' is always defined on the interpreter.
    # Create a new VTK scene.
    #mayavi.new_scene()

    # Read a VTK (old style) data file.
    r = VTKXMLFileReader()
    r.initialize(join(dirname(mayavi.__file__),
                      'examples', 'data', 'fire_ug.vtu'))
    mayavi.add_source(r)

    # Create an outline and a vector cut plane.
    mayavi.add_module(Outline())

    v = VectorCutPlane()
    mayavi.add_module(v)
    v.glyph.color_mode = 'color_by_scalar'

    # Now mask the points and show glyphs (we could also use
    # Vectors but glyphs are a bit more generic)
    m = MaskPoints()
    m.filter.set(on_ratio=10, random_mode=True)
    mayavi.add_filter(m)

    g = Glyph()
    mayavi.add_module(g)
    # Note that this adds the module to the filtered output.
    g.glyph.scale_mode = 'scale_by_vector'
    # Use arrows to view the scalars.
    g.glyph.glyph_source = g.glyph.glyph_list[1]


if __name__ == '__main__':
    glyph()
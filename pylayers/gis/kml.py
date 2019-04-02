import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from simplekml import (Kml, OverlayXY, ScreenXY, Units, RotationXY,
                       AltitudeMode, Camera)

def make_kml(extent,figs,colorbar=None, **kw):
    """
    Parameters
    ----------

    extent : tuple
        (lm,lM,Lm,LM)
        lower left corner longitude
        upper right corner longitude
        lower left corner Latitude
        upper right corner Latitude

    altitude : float
    altitudemode :
    roll : float
    tilt : float
    visibility : int

    """

    lm = extent[0]
    lM = extent[1]
    Lm = extent[2]
    LM = extent[3]

    kml = Kml()

    altitude = kw.pop('altitude', 2e6)
    roll = kw.pop('roll', 0)
    tilt = kw.pop('tilt', 0)
    altitudemode = kw.pop('altitudemode', AltitudeMode.relativetoground)

    camera = Camera(latitude  = np.mean([Lm, LM]),
                    longitude = np.mean([lm, lM]),
                    altitude = altitude,
                    roll=roll,
                    tilt=tilt,
                    altitudemode=altitudemode)

    kml.document.camera = camera

    draworder = 0

    for fig in figs:  # NOTE: Overlays are limited to the same bbox.
        draworder += 1
        ground = kml.newgroundoverlay(name='GroundOverlay')
        ground.draworder = draworder
        ground.visibility = kw.pop('visibility', 1)
        ground.name = kw.pop('name', 'overlay')
        ground.color = kw.pop('color', '9effffff')
        ground.atomauthor = kw.pop('author', 'author')
        ground.latlonbox.rotation = kw.pop('rotation', 0)
        ground.description = kw.pop('description', 'matplotlib figure')
        ground.gxaltitudemode = kw.pop('gxaltitudemode', 'clampToSeaFloor')
        ground.icon.href = fig
        ground.latlonbox.east = lM
        ground.latlonbox.south = Lm
        ground.latlonbox.north = LM
        ground.latlonbox.west = lm

    if colorbar:  # Options for colorbar are hard-coded (to avoid a big mess).
        screen = kml.newscreenoverlay(name='ScreenOverlay')
        screen.icon.href = colorbar
        screen.overlayxy = OverlayXY(x=0, y=0,
                                     xunits=Units.fraction,
                                     yunits=Units.fraction)
        screen.screenxy = ScreenXY(x=0.015, y=0.075,
                                   xunits=Units.fraction,
                                   yunits=Units.fraction)
        screen.rotationXY = RotationXY(x=0.5, y=0.5,
                                       xunits=Units.fraction,
                                       yunits=Units.fraction)
        screen.size.x = 0
        screen.size.y = 0
        screen.size.xunits = Units.fraction
        screen.size.yunits = Units.fraction
        screen.visibility = 1

    kmzfile = kw.pop('kmzfile', 'overlay.kmz')
    kml.savekmz(kmzfile)



def gearth_fig(extent,extent_c):
    """google earth figure

    Parameters
    ----------

    """
    Dx = extent_c[1]-extent_c[0]
    Dy = extent_c[3]-extent_c[2]
    aspect = Dy/Dx

    if aspect < 1.0:
        figsize = (10.0 / aspect, 10.0)
    else:
        figsize = (10.0, 10.0 * aspect)

    figsize=(10,10)

    fig = plt.figure(figsize=figsize)

    if aspect < 1.0:
        lax = [0, 0, 1.0/aspect,1]
    else:
        lax = [0, 0, 1.0,aspect]

    lax = [0,0,1,1] 
    ax = fig.add_axes(lax)
    print(figsize,aspect,lax)

    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])

    return fig, ax

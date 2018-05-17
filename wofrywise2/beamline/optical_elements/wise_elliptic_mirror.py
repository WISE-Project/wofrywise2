import numpy

from syned.beamline.shape import Ellipse

from wofrywise2.beamline.wise_optical_element import WiseOpticalElement

from wiselib2.Fundation import OpticalElement, PositioningDirectives
from wiselib2.Optics import MirrorElliptic


class WiseEllipticMirror(WiseOpticalElement):
    def __init__(self,
                 name="Undefined",
                 elliptic_mirror = MirrorElliptic(f1 = 98, f2 = 1.2, Alpha = numpy.deg2rad(2.5), L = 0.4),
                 position_directives=PositioningDirectives()):
        #:TODO boundary shape must be checked, is actually useless right now
        super(WiseEllipticMirror, self).__init__(name=name,
                                                 boundary_shape=Ellipse(min_ax_left=-0.5*elliptic_mirror.f1 + position_directives.XYCentre[1],
                                                                        min_ax_right=0.5*elliptic_mirror.f1 + position_directives.XYCentre[1],
                                                                        maj_ax_bottom=-0.5*elliptic_mirror.f2 + position_directives.XYCentre[0],
                                                                        maj_ax_top=-0.5*elliptic_mirror.f2 + position_directives.XYCentre[0]),
                                                 wise_optical_element = OpticalElement(Element=elliptic_mirror,
                                                                                       PositioningDirectives=position_directives,
                                                                                       Name=name,
                                                                                       IsSource=False))

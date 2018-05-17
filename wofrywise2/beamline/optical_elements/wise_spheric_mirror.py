
from syned.beamline.shape import Sphere

from wofrywise2.beamline.wise_optical_element import WiseOpticalElement

from wiselib2.Fundation import OpticalElement, PositioningDirectives
from wiselib2.Optics import MirrorSpheric


class WiseSphericMirror(WiseOpticalElement):
    def __init__(self, name="Undefined", spheric_mirror = MirrorSpheric(), position_directives=PositioningDirectives()):
        #:TODO boundary shape must be checked, is actually useless right now
        super(WiseSphericMirror, self).__init__(name=name,
                                                boundary_shape=Sphere(radius=spheric_mirror.f2),
                                                wise_optical_element = OpticalElement(Element=spheric_mirror,
                                                                                      PositioningDirectives=position_directives,
                                                                                      Name=name,
                                                                                      IsSource=False))

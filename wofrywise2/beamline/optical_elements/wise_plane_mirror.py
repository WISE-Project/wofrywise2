import numpy

from syned.beamline.shape import Plane

from wofrywise2.beamline.wise_optical_element import WiseOpticalElement

from wiselib2.Fundation import OpticalElement, PositioningDirectives
from wiselib2.Optics import MirrorPlane

class WisePlaneMirror(WiseOpticalElement):
    def __init__(self,
                 name="Undefined",
                 plane_mirror = MirrorPlane(L=0.4, AngleGrazing = numpy.deg2rad(2.5)),
                 position_directives=PositioningDirectives()):
        super(WisePlaneMirror, self).__init__(name=name,
                                              boundary_shape=Plane(),
                                              wise_optical_element = OpticalElement(Element=plane_mirror,
                                                                                    PositioningDirectives=position_directives,
                                                                                    Name=name,
                                                                                    IsSource=False))

import numpy

from wofrywise2.beamline.wise_optical_element import WiseOpticalElement

from wiselib2.Fundation import OpticalElement, PositioningDirectives
from wiselib2.Optics import Detector

class WiseDetector(WiseOpticalElement):
    def __init__(self,
                 name="Undefined",
                 detector = Detector(L=400e-6, AngleGrazing = numpy.deg2rad(90) ),
                 position_directives=PositioningDirectives()):
        super(WiseDetector, self).__init__(name=name,
                                           boundary_shape=None,
                                           wise_optical_element = OpticalElement(Element=detector,
                                                                                 PositioningDirectives=position_directives,
                                                                                 Name=name,
                                                                                 IsSource=False))

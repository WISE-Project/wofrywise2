
from wofrywise2.beamline.wise_optical_element import WiseOpticalElement

from wiselib2.Fundation import OpticalElement, PositioningDirectives
from wiselib2.Optics import SourceGaussian


class WiseGaussianSource(WiseOpticalElement):
    def __init__(self, name="Undefined", source_gaussian = SourceGaussian(Waist0=0.0, Lambda=10.0), position_directives=PositioningDirectives()):
        #:TODO boundary shape must be checked, is actually useless right now
        super(WiseGaussianSource, self).__init__(name=name,
                                                 boundary_shape=None,
                                                 wise_optical_element = OpticalElement(Element=source_gaussian,
                                                                                       PositioningDirectives=position_directives,
                                                                                       Name=name,
                                                                                       IsSource=True))

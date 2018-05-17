
from syned.beamline.beamline_element import BeamlineElement
from wofrywise2.beamline.wise_optical_element import WiseOpticalElement

class WiseBeamlineElement(BeamlineElement):

    def __init__(self, optical_element=WiseOpticalElement()):
        super(WiseBeamlineElement, self).__init__(optical_element=optical_element, coordinates=None)


    def get_coordinates(self):
        raise NotImplementedError("this method cannot be used in WISE 2")


from syned.beamline.optical_element import OpticalElement as SynedOpticaElement

from wiselib2.Fundation import OpticalElement


class WiseOpticalElement(SynedOpticaElement):

    def __init__(self, name="Undefined", boundary_shape=None, wise_optical_element=None):
        super(WiseOpticalElement, self).__init__(name=wise_optical_element.Name if name is None else name,
                                                 boundary_shape=boundary_shape)

        self.wise_optical_element = wise_optical_element




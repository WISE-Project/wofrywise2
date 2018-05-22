import numpy

import scipy.constants as codata
angstroms_to_eV = codata.h*codata.c/codata.e*1e10

from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
from wofry.propagator.propagator import Propagator1D, PropagationParameters, PropagationElements, PropagationManager, InteractiveMode

from wofrywise2.propagator.wavefront1D.wise_wavefront import WiseWavefront
from wofrywise2.beamline.wise_beamline_element import WiseBeamlineElement
from wofrywise2.beamline.wise_optical_element import WiseOpticalElement
from wofrywise2.beamline.optical_elements.wise_detector import WiseDetector

from wiselib2 import Fundation, Optics

class WisePropagationElements(PropagationElements):

    __wise_propagation_elements = None

    def __init__(self):
        super(WisePropagationElements, self).__init__()

        self.__wise_propagation_elements = Fundation.BeamlineElements()

    def add_beamline_element(self, beamline_element=WiseBeamlineElement()):
        super(WisePropagationElements, self).add_beamline_element(beamline_element)

        self.__wise_propagation_elements.Append(beamline_element.get_optical_element().wise_optical_element)

    def insert_beamline_element(self, index, new_element=WiseBeamlineElement(), mode=PropagationElements.INSERT_BEFORE):
        super(WisePropagationElements, self).insert_beamline_element(index, new_element, mode)

        self.__wise_propagation_elements.Insert(new_element.get_optical_element().wise_optical_element,
                                                ExistingName=self.get_wise_propagation_element(index).Name,
                                                Mode=mode+1)


    def add_beamline_elements(self, beamline_elements=[]):
        for beamline_element in beamline_elements:
            self.add_beamline_element(beamline_element)

    def get_wise_propagation_element(self, index):
        return self.get_propagation_element(index).get_optical_element().wise_optical_element

    def get_wise_propagation_elements(self):
        return self.__wise_propagation_elements

class WisePropagator(Propagator1D):

    HANDLER_NAME = "WISE2_PROPAGATOR"

    def get_handler_name(self):
        return self.HANDLER_NAME

    def do_propagation(self, parameters=PropagationParameters()):
        wavefront = parameters.get_wavefront()

        if not wavefront is None:
            is_generic_wavefront = isinstance(wavefront, GenericWavefront1D)
        else:
            is_generic_wavefront = False

        if not is_generic_wavefront and not wavefront is None:
            if not isinstance(wavefront, WiseWavefront): raise ValueError("Wavefront cannot be managed by this propagator")

        if is_generic_wavefront:
            wavefront = WiseWavefront.fromGenericWavefront(wavefront)

        wise_propagation_elements = parameters.get_PropagationElements()

        optical_element_end = wise_propagation_elements.get_propagation_element(-1).get_optical_element()
        oeEnd = optical_element_end.wise_optical_element

        if parameters.get_additional_parameter("single_propagation") == True:
            if oeEnd.IsSource: raise ValueError("Computation is impossibile: Optical Element is the Source")

            oeStart = wise_propagation_elements.get_wise_propagation_element(-2)
        else:
            oeStart = wise_propagation_elements.get_wise_propagation_element(0)

        beamline = wise_propagation_elements.get_wise_propagation_elements()
        if isinstance(optical_element_end, WiseDetector): beamline.RefreshPositions()
        beamline.ComputationSettings.NPools = int(parameters.get_additional_parameter("NPools"))

        if PropagationManager.Instance().get_interactive_mode() == InteractiveMode.ENABLED or parameters.get_additional_parameter("is_full_propagator"):
            beamline.ComputeFields(oeStart=oeStart, oeEnd=oeEnd, Verbose=False)

            result = WiseWavefront(wise_computation_results=oeEnd.ComputationResults)
        elif PropagationManager.Instance().get_interactive_mode() == InteractiveMode.DISABLED:
            result = wavefront
        else:
            result = None

        if is_generic_wavefront:
            return None if result is None else result.toGenericWavefront()
        else:
            return result



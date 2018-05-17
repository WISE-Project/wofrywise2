import scipy.constants as codata
angstroms_to_eV = codata.h*codata.c/codata.e*1e10

from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
from wofry.propagator.propagator import Propagator1D, PropagationParameters, PropagationElements

from wofrywise2.propagator.wavefront1D.wise_wavefront import WiseWavefront
from wiselib2 import Fundation

class WisePropagator(Propagator1D):

    HANDLER_NAME = "WISE2_PROPAGATOR"

    def get_handler_name(self):
        return self.HANDLER_NAME

    def do_propagation(self, parameters=PropagationParameters()):
        propagation_elements = parameters.get_PropagationElements()

        wavefront = parameters.get_wavefront()

        is_generic_wavefront = False

        if not wavefront is None:
            is_generic_wavefront = isinstance(wavefront, GenericWavefront1D)

            if is_generic_wavefront:
                wavefront = WiseWavefront.fromGenericWavefront(wavefront)
            else:
                if not isinstance(wavefront, WiseWavefront): raise ValueError("Wavefront cannot be managed by this propagator")

        computation_result = None if wavefront is None else parameters.get_wavefront().wise_computation_result

        NPools = int(parameters.get_additional_parameter("NPools"))

        beamline = Fundation.BeamlineElements()

        if propagation_elements.get_propagation_elements_number() == 1:
            wise_optical_element = propagation_elements.get_propagation_element(0).get_optical_element().wise_optical_element

            if wise_optical_element.IsSource:
                raise ValueError("Computation is impossibile: Optical Element is the Source")

            if wise_optical_element.Parent is None and not wise_optical_element.Parent.IsSource:
                wise_optical_element.Parent = Fundation.OpticalElement(Name="Dummy")

            if wise_optical_element.Parent.ComputationResults.Field is None and not wise_optical_element.Parent.IsSource:
                wise_optical_element.Parent.ComputationResults = computation_result

                raise ValueError("Computation is impossibile: Parent Optical Element has no computed Field and Optical Element is not the Source")

            beamline.Append(wise_optical_element)
            beamline.RefreshPositions()
            beamline.ComputeFields(oeStart = wise_optical_element.Parent, oeEnd = wise_optical_element, NPools=NPools)

        elif propagation_elements.get_propagation_elements_number() > 1:

            first_wise_optical_element = propagation_elements.get_propagation_element(0).get_optical_element().wise_optical_element

            if first_wise_optical_element.ComputationResult is None and first_wise_optical_element.IsSource:
                first_wise_optical_element.ComputationResult = computation_result

            if first_wise_optical_element.Parent is None and not first_wise_optical_element.Parent.IsSource:
                first_wise_optical_element.Parent = Fundation.OpticalElement(Name="Dummy")

            if first_wise_optical_element.Parent.ComputationResults.Field is None and not first_wise_optical_element.Parent.IsSource:
                first_wise_optical_element.Parent.ComputationResults = computation_result

            for beamline_element in propagation_elements.get_propagation_elements():
                beamline.Append(beamline_element.get_optical_element().wise_optical_element)

            beamline.RefreshPositions()

            beamline.ComputeFields(oeStart=parameters.get_PropagationElements[0].wise_optical_element,
                                   oeEnd=parameters.get_PropagationElements()[-1].wise_optical_element,
                                   NPools = NPools)

        else:
            raise ValueError("Computation is impossibile: PropagationParameters object contains no elements")

        result = WiseWavefront(wise_computation_result=parameters.get_PropagationElements()[-1].wise_optical_element.ComputationResult)

        if is_generic_wavefront:
            return result.toGenericWavefront()
        else:
            return result


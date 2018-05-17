import scipy.constants as codata
angstroms_to_eV = codata.h*codata.c/codata.e*1e10

from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
from wofry.propagator.propagator import Propagator1D, PropagationParameters

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

        beamline = Fundation.BeamlineElements()
        beamline.ComputationSettings.NPools = int(parameters.get_additional_parameter("NPools"))

        oeEnd = parameters.get_PropagationElements().get_propagation_element(-1).get_optical_element().wise_optical_element

        if propagation_elements.get_propagation_elements_number() == 1:
            wise_optical_element = propagation_elements.get_propagation_element(0).get_optical_element().wise_optical_element

            if wise_optical_element.IsSource:
                raise ValueError("Computation is impossibile: Optical Element is the Source")

            if wise_optical_element.Parent is None:
                wise_optical_element.Parent = Fundation.OpticalElement(Name="Dummy", PositioningDirectives=Fundation.PositioningDirectives())

            if wise_optical_element.Parent.ComputationResults.Field is None:
                if computation_result is None:  raise ValueError("Computation is impossibile: Parent Optical Element has no computed Field and Optical Element is not the Source")

                wise_optical_element.Parent.ComputationResults = computation_result

            beamline.Append(wise_optical_element)

            oeStart = wise_optical_element.Parent

        elif propagation_elements.get_propagation_elements_number() > 1:
            first_wise_optical_element = propagation_elements.get_propagation_element(0).get_optical_element().wise_optical_element

            if not first_wise_optical_element.IsSource:
                if first_wise_optical_element.Parent is None:
                    first_wise_optical_element.Parent = Fundation.OpticalElement(Name="Dummy", PositioningDirectives=Fundation.PositioningDirectives())

                if first_wise_optical_element.Parent.ComputationResults.Field is None:
                    if computation_result is None:  raise ValueError("Computation is impossibile: Parent Optical Element has no computed Field and Optical Element is not the Source")

                    first_wise_optical_element.Parent.ComputationResults = computation_result
            else:
                if first_wise_optical_element.ComputationResults.Field is None: first_wise_optical_element.ComputationResults = computation_result

            for beamline_element in propagation_elements.get_propagation_elements():
                beamline.Append(beamline_element.get_optical_element().wise_optical_element)

            oeStart=parameters.get_PropagationElements().get_propagation_element(0).get_optical_element().wise_optical_element
        else:
            raise ValueError("Computation is impossibile: PropagationParameters object contains no elements")
        
        beamline.RefreshPositions()
        beamline.ComputeFields(oeStart=oeStart, oeEnd=oeEnd)

        result = WiseWavefront(wise_computation_results=oeEnd.ComputationResults)
        
        if is_generic_wavefront:
            return result.toGenericWavefront()
        else:
            return result


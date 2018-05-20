import numpy

import scipy.constants as codata
angstroms_to_eV = codata.h*codata.c/codata.e*1e10

from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
from wofry.propagator.propagator import Propagator1D, PropagationParameters, PropagationElements

from wofrywise2.propagator.wavefront1D.wise_wavefront import WiseWavefront
from wofrywise2.beamline.wise_beamline_element import WiseBeamlineElement

from wiselib2 import Fundation, Optics

class WisePropagationElements(PropagationElements):
    __wise_propagation_elements = None

    def __init__(self):
        super(WisePropagationElements, self).__init__()

        self.__wise_propagation_elements = Fundation.BeamlineElements()

    def add_beamline_element(self, beamline_element=WiseBeamlineElement()):
        super(WisePropagationElements, self).add_beamline_element(beamline_element)

        self.__wise_propagation_elements.Append(beamline_element.get_optical_element().wise_optical_element)

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

        wise_propagation_elements = parameters.get_PropagationElements()

        oeEnd = wise_propagation_elements.get_wise_propagation_element(-1)

        if parameters.get_additional_parameter("single_propagation") == True:
            if oeEnd.IsSource:
                raise ValueError("Computation is impossibile: Optical Element is the Source")

            if oeEnd.Parent is None:
                oeEnd.Parent = Fundation.OpticalElement(Name="Dummy",
                                                        PositioningDirectives= Fundation.PositioningDirectives(ReferTo = Fundation.PositioningDirectives.ReferTo.AbsoluteReference,
                                                                                                               XYCentre = [0,0],
                                                                                                               Angle = numpy.deg2rad(0)))
                oeEnd.Parent.ComputationResults = wavefront.wise_computation_result

            if not oeEnd.Parent.IsSource and oeEnd.Parent.ComputationResults.Field is None:
                if wavefront.wise_computation_result is None:  raise ValueError("Computation is impossibile: Parent Optical Element has no computed Field and Optical Element is not the Source")

                oeEnd.Parent.ComputationResults = wavefront.wise_computation_result

            oeStart = oeEnd.Parent
        else:
            oeStart = wise_propagation_elements.get_wise_propagation_element(0)

            if not oeStart.IsSource:
                if oeStart.Parent is None:
                    oeStart.Parent = Fundation.OpticalElement(Name="Dummy",
                                                            PositioningDirectives= Fundation.PositioningDirectives(ReferTo = Fundation.PositioningDirectives.ReferTo.AbsoluteReference,
                                                                                                                   XYCentre = [0,0],
                                                                                                                   Angle = numpy.deg2rad(0)))
                    oeStart.Parent.ComputationResults = wavefront.wise_computation_result

                if oeStart.Parent.ComputationResults.Field is None:
                    if wavefront.wise_computation_result is None:  raise ValueError("Computation is impossibile: Parent Optical Element has no computed Field and Optical Element is not the Source")

                    oeStart.Parent.ComputationResults = wavefront.wise_computation_result
            else:
                if oeStart.ComputationResults.Field is None: oeStart.ComputationResults = wavefront.wise_computation_result

            oeStart = wise_propagation_elements.get_wise_propagation_element(0)

        beamline = wise_propagation_elements.get_wise_propagation_elements()

        if isinstance(oeEnd, Optics.MirrorElliptic): beamline.RefreshPositions()
        beamline.ComputationSettings.NPools = int(parameters.get_additional_parameter("NPools"))
        beamline.ComputeFields(oeStart=oeStart, oeEnd=oeEnd)

        result = WiseWavefront(wise_computation_results=oeEnd.ComputationResults)

        if is_generic_wavefront:
            return result.toGenericWavefront()
        else:
            return result

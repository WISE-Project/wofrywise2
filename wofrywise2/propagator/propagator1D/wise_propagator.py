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


class DummyElement(Optics.TransmissionMask):
    def __init__(self, electric_field):
        super(DummyElement, self).__init__(L=1)

        self.electric_field = electric_field

    def EvalField(self, x1, y1, Lambda, E0, NPools = 3,  Options = ['HF']):
        return self.electric_field

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
        beamline = wise_propagation_elements.get_wise_propagation_elements()

        oeEnd = wise_propagation_elements.get_wise_propagation_element(-1)

        if parameters.get_additional_parameter("single_propagation") == True:
            if oeEnd.IsSource:
                raise ValueError("Computation is impossibile: Optical Element is the Source")

            if oeEnd.Parent is None:
                dummy_optical_element = Fundation.OpticalElement(Name="Dummy",
                                                                 Element=DummyElement(wavefront.wise_computation_result.Field),
                                                                 PositioningDirectives= Fundation.PositioningDirectives(ReferTo = Fundation.PositioningDirectives.ReferTo.AbsoluteReference,
                                                                                                                        XYCentre = [0,0],
                                                                                                                        Angle = numpy.deg2rad(0)))

                #dummy_optical_element.ComputationSettings.Ignore = True
                dummy_optical_element.ComputationResults = wavefront.wise_computation_result
                dummy_optical_element.ComputationSettings.UseCustomSampling = oeEnd.ComputationSettings.UseCustomSampling
                dummy_optical_element.ComputationSettings.NSamples = oeEnd.ComputationSettings.NSamples

                beamline.Insert(dummy_optical_element, ExistingName=oeEnd.Name, Mode=Fundation.INSERT_MODE.Before)
            if not oeEnd.Parent.IsSource and oeEnd.Parent.ComputationResults.Field is None:
                if wavefront.wise_computation_result is None:  raise ValueError("Computation is impossibile: Parent Optical Element has no computed Field and Optical Element is not the Source")

                oeEnd.Parent.ComputationResults = wavefront.wise_computation_result

            oeStart = oeEnd.Parent
        else:
            oeStart = wise_propagation_elements.get_wise_propagation_element(0)

            if not oeStart.IsSource:
                if oeStart.Parent is None:
                    oeStart.Parent = Fundation.OpticalElement(Name="Dummy",
                                                              IsSource = True,
                                                              PositioningDirectives= Fundation.PositioningDirectives(ReferTo = Fundation.PositioningDirectives.ReferTo.AbsoluteReference,
                                                                                                                     XYCentre = [0,0],
                                                                                                                     Angle = numpy.deg2rad(0)))
                    oeStart.Parent.ComputationResults = wavefront.wise_computation_result
                    oeStart.Parent.ComputationSettings.UseCustomSampling = oeStart.ComputationSettings.UseCustomSampling
                    oeStart.Parent.ComputationSettings.NSamples = oeStart.ComputationSettings.NSamples

                if oeStart.Parent.ComputationResults.Field is None:
                    if wavefront.wise_computation_result is None:  raise ValueError("Computation is impossibile: Parent Optical Element has no computed Field and Optical Element is not the Source")

                    oeStart.Parent.ComputationResults = wavefront.wise_computation_result
            else:
                if oeStart.ComputationResults.Field is None: oeStart.ComputationResults = wavefront.wise_computation_result

            oeStart = wise_propagation_elements.get_wise_propagation_element(0)

        if isinstance(oeEnd, Optics.MirrorElliptic): beamline.RefreshPositions()
        beamline.ComputationSettings.NPools = int(parameters.get_additional_parameter("NPools"))
        beamline.ComputeFields(oeStart=oeStart, oeEnd=oeEnd)

        result = WiseWavefront(wise_computation_results=oeEnd.ComputationResults)

        if is_generic_wavefront:
            return result.toGenericWavefront()
        else:
            return result

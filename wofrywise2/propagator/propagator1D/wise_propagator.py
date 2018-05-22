import numpy

import scipy.constants as codata
angstroms_to_eV = codata.h*codata.c/codata.e*1e10

from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
from wofry.propagator.propagator import Propagator1D, PropagationParameters, PropagationElements

from wofrywise2.propagator.wavefront1D.wise_wavefront import WiseWavefront
from wofrywise2.beamline.wise_beamline_element import WiseBeamlineElement
from wofrywise2.beamline.wise_optical_element import WiseOpticalElement

from wiselib2 import Fundation, Optics

class WisePropagationElements(PropagationElements):

    __wise_propagation_elements = None

    def __init__(self):
        super(WisePropagationElements, self).__init__()

        self.__wise_propagation_elements = Fundation.BeamlineElements()

    def add_beamline_element(self, beamline_element=WiseBeamlineElement()):
        super(WisePropagationElements, self).add_beamline_element(beamline_element)

        self.__wise_propagation_elements.Append(beamline_element.get_optical_element().wise_optical_element)

    def insert_beamline_element(self, index, new_element=WiseBeamlineElement(), mode=Fundation.INSERT_MODE.Before):
        if not (mode == Fundation.INSERT_MODE.After or mode == Fundation.INSERT_MODE.Before): raise ValueError("Fork mode is not supported")

        self.__wise_propagation_elements.Insert(new_element.get_optical_element().wise_optical_element,
                                                ExistingName=self.get_wise_propagation_element(index).Name,
                                                Mode=mode)

        if mode == Fundation.INSERT_MODE.Before:
            if index == 0:
                super(WisePropagationElements, self).__propagation_elements = [new_element] + super(WisePropagationElements, self).__propagation_elements
            else:
                super(WisePropagationElements, self).__propagation_elements.insert(index, new_element)
        elif mode == Fundation.INSERT_MODE.After:
            if index == len(super(WisePropagationElements, self).__propagation_elements) - 1:
                super(WisePropagationElements, self).__propagation_elements = super(WisePropagationElements, self).__propagation_elements + [new_element]
            else:
                super(WisePropagationElements, self).__propagation_elements.insert(index+1, new_element)

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

        oeEnd = wise_propagation_elements.get_wise_propagation_element(-1)

        if parameters.get_additional_parameter("single_propagation") == True:
            if oeEnd.IsSource:
                raise ValueError("Computation is impossibile: Optical Element is the Source")

            if oeEnd.Parent is None:
                wise_propagation_elements.insert_beamline_element(index = -1,
                                                                  new_element=WiseBeamlineElement(optical_element=WiseOpticalElement(wise_optical_element=get_dummy_element(wavefront, oeEnd))),
                                                                  mode=Fundation.INSERT_MODE.Before)

            oeStart = wise_propagation_elements.get_wise_propagation_element(-2)
        else:
            oeStart = wise_propagation_elements.get_wise_propagation_element(0)

            if not oeStart.IsSource and oeStart.Parent is None:
                wise_propagation_elements.insert_beamline_element(index = 0,
                                                                  new_element=WiseBeamlineElement(optical_element=WiseOpticalElement(wise_optical_element=get_dummy_element(wavefront, oeStart))),
                                                                  mode=Fundation.INSERT_MODE.Before)

                oeStart = wise_propagation_elements.get_wise_propagation_element(0)

        beamline = wise_propagation_elements.get_wise_propagation_elements()
        if isinstance(oeEnd.CoreOptics, Optics.Detector): beamline.RefreshPositions()
        beamline.ComputationSettings.NPools = int(parameters.get_additional_parameter("NPools"))
        beamline.ComputeFields(oeStart=oeStart, oeEnd=oeEnd)

        result = WiseWavefront(wise_computation_results=oeEnd.ComputationResults)

        if is_generic_wavefront:
            return result.toGenericWavefront()
        else:
            return result


class DummyElement(Optics.SourceGaussian):
    def __init__(self, Lambda, electric_field):
        super(DummyElement, self).__init__(Lambda, 1e-6)

        self.electric_field = electric_field

    def EvalField(self, x1, y1, Lambda, NPools=3, **kwargs):
        return self.electric_field

def get_dummy_element(wavefront, oe):
    dummy_optical_element = Fundation.OpticalElement(Name="Dummy",
                                                     IsSource=True,
                                                     Element=DummyElement(wavefront.wise_computation_result.Lambda,
                                                                          wavefront.wise_computation_result.Field),
                                                     PositioningDirectives= Fundation.PositioningDirectives(ReferTo = Fundation.PositioningDirectives.ReferTo.AbsoluteReference,
                                                                                                            XYCentre = [0,0],
                                                                                                            Angle = numpy.deg2rad(0)))

    dummy_optical_element.ComputationResults = wavefront.wise_computation_result
    dummy_optical_element.ComputationSettings.UseCustomSampling = oe.ComputationSettings.UseCustomSampling
    dummy_optical_element.ComputationSettings.NSamples = oe.ComputationSettings.NSamples

    return dummy_optical_element

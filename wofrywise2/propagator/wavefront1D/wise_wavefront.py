import numpy

from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D, WavefrontDimension
from wofry.propagator.decorators import WavefrontDecorator

from wiselib2.Fundation import ComputationResults

class WiseWavefront(WavefrontDecorator):

    def __init__(self,
                 wise_computation_results=ComputationResults()):

        self.wise_computation_result = wise_computation_results

    def get_dimension(self):
        return WavefrontDimension.ONE

    def toGenericWavefront(self):
        wavelength = self.wise_computation_result.Lambda
        position = self.wise_computation_result.S
        electric_field = numpy.real(self.wise_computation_result.Field) + 1j*numpy.imag(self.wise_computation_result.Field)

        return GenericWavefront1D.initialize_wavefront_from_arrays(x_array=position, y_array=electric_field, wavelength=wavelength)

    @classmethod
    def fromGenericWavefront(cls, wavefront):
        wise_computation_result = ComputationResults()
        wise_computation_result.Lambda = wavefront.get_wavelength()
        wise_computation_result.S = wavefront.get_abscissas()
        wise_computation_result.Field = wavefront.get_complex_amplitude()

        return WiseWavefront(wise_computation_results=wise_computation_result)

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 11:56:47 2017

@author: Mic

- Funziona
- la lascio che c'è da correggere in EllipticalMirror pTan e quindi VersosNorm
- la lascio che c'è ancora l'inifto bisticco tra XYCentre e XYCentre, XYF1 e XYF1
"""

import importlib

import wiselib2.Rayman as rm
import wiselib2.Fundation as Fundation
import wiselib2.Optics as Optics
import wiselib2.ToolLib as tl
import wiselib2.FermiSource as Fermi

importlib.reload(Fundation)
importlib.reload(Optics)
importlib.reload(tl)
importlib.reload(rm)
importlib.reload(Fermi)

from matplotlib import pyplot as plt

from numpy import *

from wofry.propagator.propagator import PropagationParameters, PropagationManager

from wofrywise2.beamline.wise_beamline_element import WiseBeamlineElement
from wofrywise2.beamline.light_sources.wise_gaussian_source import WiseGaussianSource
from wofrywise2.beamline.optical_elements.wise_plane_mirror import WisePlaneMirror
from wofrywise2.beamline.optical_elements.wise_elliptic_mirror import WiseEllipticMirror
from wofrywise2.beamline.optical_elements.wise_detector import WiseDetector

from wofrywise2.propagator.wavefront1D.wise_wavefront import WiseWavefront
from wofrywise2.propagator.propagator1D.wise_propagator import WisePropagator, WisePropagationElements

def plot(oe, id):
    S = oe.ComputationResults.S
    E = oe.ComputationResults.Field

    plt.figure(id)
    plt.plot(S*1e3, abs(E)**2/max(abs(E)**2))
    plt.xlabel('mm')
    plt.title('|E| (' + oe.Name + ')')

print(__name__)
if __name__ == '__main__':

    PropagationManager.Instance().add_propagator(WisePropagator())


    tl.Debug.On = True
    N = 7000
    UseCustomSampling = True
    # SOURCE
    #==========================================================================
    Lambda = 32e-9
    Waist0 = Fermi.Waist0E(Lambda)

    s_k = Optics.SourceGaussian(Lambda, Waist0)      # Kernel delle ottiche
    s_pd = Fundation.PositioningDirectives(            # Direttive di posizionamento
                        ReferTo = Fundation.PositioningDirectives.ReferTo.AbsoluteReference,
                        XYCentre = [0.0, 0.0],
                        Angle = 0.0)


    s = WiseGaussianSource(name='source',
                           source_gaussian=s_k,
                           position_directives=s_pd)


    # PM1A (h)
    #==========================================================================
    pm1a_k = Optics.MirrorPlane(L=0.4, AngleGrazing = deg2rad(2.5) )
    pm1a_pd = Fundation.PositioningDirectives(
                                    ReferTo = 'upstream',
                                    PlaceWhat = 'centre',
                                    PlaceWhere = 'centre',
                                    Distance = 48.0901)

    pm1a = WisePlaneMirror(name= 'pm1a',
                           plane_mirror=pm1a_k,
                           position_directives=pm1a_pd)

    pm1a.wise_optical_element.ComputationSettings.Ignore = False          # Lo user decide di non simulare lo specchio ()
    pm1a.wise_optical_element.ComputationSettings.UseCustomSampling = UseCustomSampling # l'utente decide di impostare a mano il campionamento
    pm1a.wise_optical_element.ComputationSettings.NSamples = N

    # PM1B (h)
    #==========================================================================
    pm1b_k = Optics.MirrorPlane(L=0.4, AngleGrazing = deg2rad(5.0) )
    pm1b_pd = Fundation.PositioningDirectives(
                                    ReferTo = 'upstream',
                                    PlaceWhat = 'centre',
                                    PlaceWhere = 'centre',
                                    Distance = 6.2388)
    pm1b = WisePlaneMirror(name= 'pm1b',
                           plane_mirror=pm1b_k,
                           position_directives=pm1b_pd)

    pm1b.wise_optical_element.ComputationSettings.Ignore = False          # Lo user decide di non simulare lo specchio ()
    pm1b.wise_optical_element.ComputationSettings.UseCustomSampling = UseCustomSampling # l'utente decide di impostare a mano il campionamento
    pm1b.wise_optical_element.ComputationSettings.NSamples = N

    # KB(h)
    #==========================================================================
    f1 = 98
    f2 = 1.2
    GrazingAngle = deg2rad(2.5)
    L = 0.4

    kb_k = Optics.MirrorElliptic(f1 = f1, f2 = f2 , L= L, Alpha = GrazingAngle)
    kb_pd = Fundation.PositioningDirectives(
                        ReferTo = 'upstream',
                        PlaceWhat = 'centre',
                        PlaceWhere = 'centre',
                        Distance=44.9751)#49.9099)

    kb = WiseEllipticMirror(name='kb',
                            elliptic_mirror=kb_k,
                            position_directives=kb_pd)

    #----- Impostazioni KB
    kb.wise_optical_element.CoreOptics.ComputationSettings.UseFigureError = False
    kb.wise_optical_element.CoreOptics.ComputationSettings.UseRoughness = False
    kb.wise_optical_element.CoreOptics.ComputationSettings.UseSmallDisplacements = False # serve per traslare/ruotare l'EO
    kb.wise_optical_element.CoreOptics.SmallDisplacements.Rotation = deg2rad(0)
    kb.wise_optical_element.CoreOptics.SmallDisplacements.Trans = 0 # Transverse displacement (rispetto al raggio uscente, magari faremo scegliere)
    kb.wise_optical_element.CoreOptics.SmallDisplacements.Long = 0 # Longitudinal displacement (idem)
    # aggiungo figure error
    kb.wise_optical_element.CoreOptics.FigureErrorLoad(File = "/Users/admin/Documents/workspace/OASYS-Develop/1.0/wiselib2/Examples/DATI/kbv.txt",
                  Step = 2e-3, # passo del file
                  AmplitudeScaling = 1*1e-3 # fattore di scala
                  )
    kb.wise_optical_element.ComputationSettings.UseCustomSampling = UseCustomSampling # l'utente decide di impostare a mano il campionamento
    kb.wise_optical_element.ComputationSettings.NSamples = N

    # detector (h)
    #==========================================================================
    d_k = Optics.Detector(
                        L=400e-6,
                        AngleGrazing = deg2rad(90))
    d_pd = Fundation.PositioningDirectives(
                        ReferTo = 'upstream',
                        PlaceWhat = 'centre',
                        PlaceWhere = 'downstream focus',
                        Distance = 0)
    d = WiseDetector(name = 'detector',
                     detector=d_k,
                     position_directives=d_pd)

    d.wise_optical_element.ComputationSettings.UseCustomSampling = UseCustomSampling
    d.wise_optical_element.ComputationSettings.NSamples = N                          # come sopra. In teoria il campionamento può essere specificato elemento per elmeento



    beamline = WisePropagationElements()
    wavefront = WiseWavefront(wise_computation_results=None)

    beamline.add_beamline_element(WiseBeamlineElement(optical_element=s))
    beamline.add_beamline_element(WiseBeamlineElement(optical_element=pm1a))

    parameters = PropagationParameters(wavefront=wavefront, propagation_elements=beamline)
    parameters.set_additional_parameters("single_propagation", True)
    parameters.set_additional_parameters("NPools", 5)

    wavefront = PropagationManager.Instance().do_propagation(propagation_parameters=parameters, handler_name=WisePropagator.HANDLER_NAME)

    if not pm1a.wise_optical_element.ComputationSettings.Ignore: plot(pm1a.wise_optical_element, 11)

    beamline.add_beamline_element(WiseBeamlineElement(optical_element=pm1b))

    parameters = PropagationParameters(wavefront=wavefront, propagation_elements=beamline)
    parameters.set_additional_parameters("single_propagation", True)
    parameters.set_additional_parameters("NPools", 5)

    wavefront = PropagationManager.Instance().do_propagation(propagation_parameters=parameters, handler_name=WisePropagator.HANDLER_NAME)

    if not pm1b.wise_optical_element.ComputationSettings.Ignore: plot(pm1b.wise_optical_element, 12)

    beamline.add_beamline_element(WiseBeamlineElement(optical_element=kb))

    parameters = PropagationParameters(wavefront=wavefront, propagation_elements=beamline)
    parameters.set_additional_parameters("single_propagation", True)
    parameters.set_additional_parameters("NPools", 5)

    wavefront = PropagationManager.Instance().do_propagation(propagation_parameters=parameters, handler_name=WisePropagator.HANDLER_NAME)

    plot(kb.wise_optical_element, 22)

    # DETECTOR
    beamline.add_beamline_element(WiseBeamlineElement(optical_element=d))

    parameters = PropagationParameters(wavefront=wavefront, propagation_elements=beamline)
    parameters.set_additional_parameters("single_propagation", True)
    parameters.set_additional_parameters("NPools", 5)

    wavefront = PropagationManager.Instance().do_propagation(propagation_parameters=parameters, handler_name=WisePropagator.HANDLER_NAME)

    plot(d.wise_optical_element, 33)

    print(beamline.get_wise_propagation_elements()) # comodo per controllare la rappresentazione interna di Beamline Element

    parameters = PropagationParameters(wavefront=WiseWavefront(wise_computation_results=None), propagation_elements=beamline)
    parameters.set_additional_parameters("single_propagation", False)
    parameters.set_additional_parameters("NPools", 5)

    wavefront = PropagationManager.Instance().do_propagation(propagation_parameters=parameters, handler_name=WisePropagator.HANDLER_NAME)

    plot(d.wise_optical_element, 44)

    plt.show()




#%%








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

from wofry.propagator.propagator import PropagationElements, PropagationParameters, PropagationManager

from wofrywise2.beamline.wise_beamline_element import WiseBeamlineElement
from wofrywise2.beamline.light_sources.wise_gaussian_source import WiseGaussianSource
from wofrywise2.beamline.optical_elements.wise_plane_mirror import WisePlaneMirror
from wofrywise2.beamline.optical_elements.wise_elliptic_mirror import WiseEllipticMirror
from wofrywise2.beamline.optical_elements.wise_detector import WiseDetector

from wofrywise2.propagator.wavefront1D.wise_wavefront import WiseWavefront
from wofrywise2.propagator.propagator1D.wise_propagator import WisePropagator

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
                        XYCentre = [0,0],
                        Angle = deg2rad(0))


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
                                    Distance = 48090.1)

    pm1a = WisePlaneMirror(name= 'pm1a',
                           plane_mirror=pm1a_k,
                           position_directives=pm1a_pd)

    pm1a.wise_optical_element.ComputationSettings.Ignore = True          # Lo user decide di non simulare lo specchio ()

    # KB(h)
    #==========================================================================
    f1 = 98
    f2 = 1.2
    GrazingAngle = deg2rad(2.5)
    L = 0.4

    #ob = Optics.Obstruction()


    kb_k = Optics.MirrorElliptic(f1 = f1, f2 = f2 , L= L, Alpha = GrazingAngle)
    kb_pd = Fundation.PositioningDirectives(
                        ReferTo = 'source',
                        PlaceWhat = 'upstream focus',
                        PlaceWhere = 'centre')

    kb = WiseEllipticMirror(name='kb',
                            elliptic_mirror=kb_k,
                            position_directives=kb_pd)

    #----- Impostazioni KB
    kb.wise_optical_element.CoreOptics.ComputationSettings.UseFigureError = True
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
                        AngleGrazing = deg2rad(90) )
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



    beamline = PropagationElements()

    beamline.add_beamline_element(WiseBeamlineElement(optical_element=s))
    beamline.add_beamline_element(WiseBeamlineElement(optical_element=pm1a))
    beamline.add_beamline_element(WiseBeamlineElement(optical_element=kb))
    beamline.add_beamline_element(WiseBeamlineElement(optical_element=d))

    parameters = PropagationParameters(wavefront=WiseWavefront(wise_computation_results=None),
                                       propagation_elements=beamline)

    parameters.set_additional_parameters("NPools", 5)

    wavefront = PropagationManager.Instance().do_propagation(propagation_parameters=parameters, handler_name=WisePropagator.HANDLER_NAME)

    assert (wavefront.wise_computation_result == d.wise_optical_element.ComputationResults)

    """
    s = elemento di inizio; d = elemento finale. Se non specificati, li trova lui
    quinti 
    t.ComputeFields() fa tutta la beamline
     Utile sarebbero i comandi Oasys (adeguali alla notazione che immagino esista già)
     Propagate To Next => t.ComputeFields(oe,oe.Children[0])
     Propagate From Previous => t.ComputeFields(oe.Parent,oe)
    
    .Children è una lista e contiene la possibilità di diramare il path ottico, 
    cosa che di fatto ora non ho mai usato.
    Va da sé che il "Child" è quindi Children[0]
    
    .ComputeFields riempie l'oggetto ComputationResults
    """

    # to fasten....

    kb = kb.wise_optical_element
    d = d.wise_optical_element

    #%%
    if 1==1:
        #------------------Intensità normalizzata su specchio
        S = kb.ComputationResults.S
        E = kb.ComputationResults.Field
        I = abs(E)**2
        I = I/max(I)

        plt.figure(11)
        plt.plot(S*1e3, abs(E)**2/max(abs(E)**2))
        plt.xlabel('mm')
        plt.title('|E| (mirror)')

        #------------------Intensità normalizzata nel fuoco
        Sd = d.ComputationResults.S
        Ed = d.ComputationResults.Field
        Id = abs(Ed)**2
        Id = Id/max(Id)

        plt.figure(13)
        plt.plot(Sd*1e3, abs(Ed))
        plt.xlabel('mm')
        plt.title('|E|^2 (detector)')



#%% Caustica

    DefocusList = linspace(-6e-3, 1e-3,   21)
    DefocusList_mm = DefocusList * 1e3

    import copy
    kb_copy = copy.deepcopy(kb)

    ResultList, HewList,SigmaList, More = Fundation.FocusSweep(kb_copy, DefocusList,
                                                            DetectorSize = 200e-6,
                                                            NPools = 4)

    N = len(ResultList)

    # Plotta il campo sui detector a varie distanze
    if 1==1:
        plt.figure(23)
        for Res in ResultList:
            plt.plot(Res.S *1e6, abs(Res.Field))
            plt.title('Campo')
            plt.xlabel('um')



#%% Plot della HEW

    plt.figure(32)
    plt.plot(DefocusList_mm, HewList,'.')
    plt.plot(DefocusList_mm, 2*0.68* SigmaList,'x')

    plt.xlabel('defocus (mm)')
    plt.ylabel('Hew')
    plt.legend(['Hew', '0.68 * 2 Sigma'])


    plt.show()
#%%








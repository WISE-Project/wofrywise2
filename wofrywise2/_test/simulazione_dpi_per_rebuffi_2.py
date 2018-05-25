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

from wiselib2.must import *
from wiselib2.Fundation import OpticalElement

from matplotlib import pyplot as plt

from numpy import *

def plot(oe, id):

    S = oe.ComputationResults.S
    E = oe.ComputationResults.Field

    plt.figure(id)
    plt.plot(S*1e3, abs(E)**2/max(abs(E)**2))
    plt.xlabel('mm')
    plt.title('|E| (' + oe.Name + ')')

print(__name__)
if __name__ == '__main__':

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
    s = OpticalElement(s_k,
                       PositioningDirectives = s_pd,
                       Name = 'source',
                       IsSource = True)

    # PM1A (h)
    #==========================================================================
    pm1a_k = Optics.MirrorPlane(L=0.4, AngleGrazing = deg2rad(2.5) )
    pm1a_pd = Fundation.PositioningDirectives(
                                    ReferTo = 'upstream',
                                    PlaceWhat = 'centre',
                                    PlaceWhere = 'centre',
                                    Distance = 48.0901)
    pm1a = OpticalElement(pm1a_k, 
                          PositioningDirectives = pm1a_pd,
                          Name = 'pm1a')
    pm1a.ComputationSettings.Ignore = False          # Lo user decide di non simulare lo specchio ()
    pm1a.ComputationSettings.UseCustomSampling = UseCustomSampling # l'utente decide di impostare a mano il campionamento
    pm1a.ComputationSettings.NSamples = N

    pm1a.CoreOptics.ComputationSettings.UseFigureError = True
    pm1a.CoreOptics.ComputationSettings.UseRoughness = False
    # aggiungo figure error
    pm1a.CoreOptics.FigureErrorLoad(File = "/Users/admin/Oasys/Mirror_figure_error.dat",
                  Step = 1e-3, # passo del file
                  AmplitudeScaling = 1.0 # fattore di scala
                  )


    # PM1A (h)
    #==========================================================================
    pm1b_k = Optics.MirrorPlane(L=0.4, AngleGrazing = deg2rad(5.0) )
    pm1b_pd = Fundation.PositioningDirectives(
                                    ReferTo = 'upstream',
                                    PlaceWhat = 'centre',
                                    PlaceWhere = 'centre',
                                    Distance = 6.2388)
    pm1b = OpticalElement(pm1b_k,
                          PositioningDirectives = pm1b_pd,
                          Name = 'pm1b')
    pm1b.ComputationSettings.Ignore = False          # Lo user decide di non simulare lo specchio ()
    pm1b.ComputationSettings.UseCustomSampling = UseCustomSampling # l'utente decide di impostare a mano il campionamento
    pm1b.ComputationSettings.NSamples = N
    pm1b.CoreOptics.ComputationSettings.UseFigureError = True
    pm1b.CoreOptics.ComputationSettings.UseRoughness = False
    # aggiungo figure error
    pm1b.CoreOptics.FigureErrorLoad(File = "/Users/admin/Oasys/Mirror_figure_error.dat",
                  Step = 1e-3, # passo del file
                  AmplitudeScaling = 1.0 # fattore di scala
                  )

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
    kb = OpticalElement(                                
                        kb_k, 
                        PositioningDirectives = kb_pd, 
                        Name = 'kb')

    #----- Impostazioni KB
    kb.CoreOptics.ComputationSettings.UseFigureError = True
    kb.CoreOptics.ComputationSettings.UseRoughness = False 
    kb.CoreOptics.ComputationSettings.UseSmallDisplacements = False # serve per traslare/ruotare l'EO
    kb.CoreOptics.SmallDisplacements.Rotation = deg2rad(0)
    kb.CoreOptics.SmallDisplacements.Trans = 0 # Transverse displacement (rispetto al raggio uscente, magari faremo scegliere)
    kb.CoreOptics.SmallDisplacements.Long = 0 # Longitudinal displacement (idem)
    # aggiungo figure error
    kb.CoreOptics.FigureErrorLoad(File = "/Users/admin/Desktop/Lavoro/Private/Progetti/OASYS/per luca/DATI/kbv.txt",
                  Step = 2e-3, # passo del file
                  AmplitudeScaling = 1*1e-3 # fattore di scala 
                  )
    kb.ComputationSettings.UseCustomSampling = UseCustomSampling # l'utente decide di impostare a mano il campionamento
    kb.ComputationSettings.NSamples = N

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
    d = OpticalElement(
                        d_k, 
                        PositioningDirectives = d_pd, 
                        Name = 'detector')
    d.ComputationSettings.UseCustomSampling = UseCustomSampling 
    d.ComputationSettings.NSamples = N                          # come sopra. In teoria il campionamento può essere specificato elemento per elmeento
    

    # Assemblamento beamline
    #==========================================================================
    t = Fundation.BeamlineElements()
    t.Append(s)
    t.Append(pm1a)         # per ora lo lasciamo commentato, devo aggiustare una cosa che si è rotta 2 gg fa
    t.RefreshPositions()

    t.ComputationSettings.NPools = 5
    t.ComputeFields(oeStart=s, oeEnd=pm1a, Verbose = False)

    if not pm1a.ComputationSettings.Ignore: plot(pm1a, 11)

    t.Append(pm1b)         # per ora lo lasciamo commentato, devo aggiustare una cosa che si è rotta 2 gg fa
    t.RefreshPositions()

    t.ComputationSettings.NPools = 5
    t.ComputeFields(oeStart=pm1a, oeEnd=pm1b, Verbose = False)

    if not pm1b.ComputationSettings.Ignore: plot(pm1b, 12)

    t.Append(kb)
    t.RefreshPositions()

    t.ComputeFields(oeStart=pm1b, oeEnd=kb, Verbose = False)

    plot(kb, 22)

    t.Append(d)
    t.RefreshPositions()

    t.ComputeFields(oeStart=kb, oeEnd=d, Verbose = False)

    print(d.ComputationResults.Action, d.ComputationResults.Lambda)

    plot(d, 33)

    print(t) # comodo per controllare la rappresentazione interna di Beamline Element

    t.ComputeFields(oeStart=s, oeEnd=d, Verbose = False)

    plot(d, 44)

    plt.show()
#%%
    

     
    
    
    
    

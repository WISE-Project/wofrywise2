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
                        XYCentre = [0,0],
                        Angle = deg2rad(0))
    s = OpticalElement(                                    # Optical Element (la cosa più vicina al pupolo Oasys)
                        s_k, 
                        PositioningDirectives = s_pd, 
                        Name = 'source', IsSource = True)


    # PM1A (h)
    #==========================================================================
    pm1a_k = Optics.MirrorPlane(L=0.4, AngleGrazing = deg2rad(2.5) )
    pm1a_pd = Fundation.PositioningDirectives(
                                    ReferTo = 'upstream',
                                    PlaceWhat = 'centre',
                                    PlaceWhere = 'centre',
                                    Distance = 48090.1)
    pm1a = OpticalElement(pm1a_k, 
                            PositioningDirectives = pm1a_pd, 
                            Name = 'pm1a')
    pm1a.ComputationSettings.Ignore = False          # Lo user decide di non simulare lo specchio ()
    pm1a.ComputationSettings.UseCustomSampling = UseCustomSampling # l'utente decide di impostare a mano il campionamento
    pm1a.ComputationSettings.NSamples = N

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
    t = None
    t = Fundation.BeamlineElements()
    t.Append(s)
    t.Append(pm1a)         # per ora lo lasciamo commentato, devo aggiustare una cosa che si è rotta 2 gg fa
    t.Append(kb)
    t.Append(d)
    t.RefreshPositions()

    print(t) # comodo per controllare la rappresentazione interna di Beamline Element

    
    #%%      Calcolo il campo: sorgente -> specchio

    t.ComputationSettings.NPools = 5
    t.ComputeFields(oeStart=s, oeEnd=d, Verbose = False)
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
    if 1==0:
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
    

     
    
    
    
    

#!/usr/bin/python3
# -*- coding: utf-8 -*-

from MSPyBentley import *
from MSPyBentleyGeom import *
from MSPyECObjects import *
from MSPyDgnPlatform import *
from MSPyDgnView import *
from MSPyMstnPlatform import *


'''
falls triangle noch nicht installiert ist, MicroStation 2024 oder 2025 beenden und triangle installieren
- Eingabeaufforderung öffnen
- Verzeichnis wechseln:
    cd C:/ProgramData/Bentley/PowerPlatformPython/python
- richtiges Verzeichnis kontrollieren ! , danach eingeben:
    python -m pip install triangle
- triangle wird installiert
    siehe auch     https://pypi.org/project/triangle/
    documentation: https://rufat.be/triangle/
                   https://www.cs.cmu.edu/~quake/triangle.html
- Microstation neu starten

## sonstiges: https://people.sc.fsu.edu/~jburkardt/c_src/triangle_shewchuk/triangle_shewchuk.html

## utm reduktion: https://www.bfrvermessung.de/materialien-1/etrs89/utm-projektmassstab-und-planungskoordinatensystem

## Python Beispiele: https://github.com/BentleySystems/MicroStationPython/tree/main/MSPythonSamples
'''


import sys

try:
  import triangle
except:
  print("triangle konnte nicht geladen werden. Bitte installieren")
  sys.exit()

import re
import numpy as np
import collections
import itertools
import os

from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import time
import datetime as dt
from tkinter import filedialog
from tkinter import scrolledtext
from tkinter import colorchooser
from tkinter import Frame
from datetime import datetime

import tkinter as tk
from tkinter import ttk


def open_file_dialog():
    file_ = filedialog.askopenfilename(initialdir = dgn_pfad, title="1 von 4: Waehle Umringkanten - line segments", filetypes=[("Text files", "*.tmp"), ("All files", "*.*")])
    return file_

def open_file_dialog_2():
    file_ = filedialog.askopenfilename(initialdir = dgn_pfad, title="2. von 4. Waehle Innenkanten - line segments", filetypes=[("Text files", "*.tmp"), ("All files", "*.*")])
    return file_

def open_file_dialog_3():
    file_ = filedialog.askopenfilename(initialdir = dgn_pfad, title="4. von 4. Waehle Lochkoordinaten - line in holes", filetypes=[("Text files", "*.tmp"), ("All files", "*.*")])
    return file_

def daten_einlesen(datei):
    liste_local = []
    with open(datei) as f:
        for zeile in f:
            liste_local.append(zeile.strip())      #  am Ende wird \n entfernt sowie Leerzeichen am Enfang und Ende
        return liste_local

def open_files_dialog():
    files_ = filedialog.askopenfilenames(initialdir = dgn_pfad, title="3. von 4. Waehle eine oder mehrere Punktdateien - coordinates", filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
    return files_

def open_files_dialog_from_ex():
    files_ = filedialog.askopenfilenames(initialdir = dgn_pfad, title="Waehle eine oder mehrere Koordinatendateien fuer den Import", filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
    return files_

def file_save(koordinaten_):
    f = filedialog.asksaveasfile(mode='w', initialdir = dgn_pfad, initialfile = label_ex.cget('text'), defaultextension=".tmp")
    if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
        return
    for zeile_ in koordinaten_:
        f.write(zeile_ + '\n')
    f.close() 

def file_save_2(koordinaten_, nummer1_):
    f = filedialog.asksaveasfile(mode='w', initialdir = dgn_pfad, initialfile = label_ex_datei.cget('text'), defaultextension=".txt")
    if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
        return
    nummer = nummer1_
    for zeile_ in koordinaten_:
        punkt_zeile = ' ' + str(nummer) + ' , '
        f.write(punkt_zeile + zeile_ + '\n')
        nummer = nummer + 1
    f.close() 
    
def file_save_3(text_):
    f = filedialog.asksaveasfile(mode='w', initialdir = dgn_pfad, initialfile = 'KGE_DGM_Volumenberechnung.txt', defaultextension=".txt")
    if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
        return
    for zeile_ in text_:
        f.write(zeile_)
    f.close()    

def file_save_4(text_):
    f = filedialog.asksaveasfile(mode='w', initialdir = dgn_pfad, initialfile = 'KGE_DGM_to_LandXML.xml', defaultextension=".xml")
    if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
        return
    for zeile_ in text_:
        f.write(zeile_+'\n')
    f.close()    

def dateien_einlesen(dateien):
    liste_local = []
    for datei in dateien:
        with open(datei) as f:
            for zeile in f:
                liste_local.append(zeile.strip())      #  am Ende wird \n entfernt sowie Leerzeichen am Enfang und Ende
        f.close()
    return liste_local

def triangle_det_2d(dreieck):
    #der absolutwert der determinante ist der doppelte flächeninhalt des dreiecks (=parallelogramm)
    #das vorzeichen gibt die anordnung der punkte an: plus, wenn im gegenuhrzeigersinn (GUZ CCW)
    #
    #print(dreieck)
    #p1 = dreieck[0]
    #p2 = dreieck[1]
    #p3 = dreieck[2]
    #print(p1, p2, p3)

    x1 = dreieck[0][0]
    y1 = dreieck[0][1]
    #z1 = dreieck[0][2]

    x2 = dreieck[1][0]
    y2 = dreieck[1][1]
    #z2 = dreieck[1][2]

    x3 = dreieck[2][0]
    y3 = dreieck[2][1]
    #z3 = dreieck[2][2]        

    #print(x1, y1)
    #print(x2, y2)
    #print(x3, y3)

    flaechehor = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
    # print(flaechehor)
    return flaechehor

def flaeche_masche(dreieck):
 
    x1 = dreieck[0][0]
    y1 = dreieck[0][1]
    z1 = dreieck[0][2]

    x2 = dreieck[1][0]
    y2 = dreieck[1][1]
    z2 = dreieck[1][2]

    x3 = dreieck[2][0]
    y3 = dreieck[2][1]
    z3 = dreieck[2][2]        

    #print(x1, y1, z1)
    #print(x2, y2, z2)
    #print(x3, y3, z3)

    a1 = x2 - x1
    a2 = y2 - y1
    a3 = z2 - z1

    b1 = x3 - x1
    b2 = y3 - y1
    b3 = z3 - z1

    flaeche = 0.5 * ((a2 * b3 - a3 * b2)**2 + (a3 * b1 -  a1 * b3)**2 + (a1 * b2 - a2 * b1)**2)**0.5
    #print(flaeche)
    return flaeche


def daten_aufbereiten(liste_, kontrolle_):
    punkte_xyz = []
    punkte_xy = []
    x_ = []
    y_ = []
    z_ = []
    mark_ = []

    for j,zeile in enumerate(liste_):
        #finde alle dezimalwerte vom format " 40.5 +79.9 123.12345  1.0 -5.67 "
        #überlese " a45 52 1. bc"
        werte = re.findall(r"[-+]?\d+\.\d+",zeile)
        if len(werte) != 3:
            #keine 3 dezimalwerte für x, y, z
            #print('Die Zeile : ', zeile, '  wird überlesen.', ' Anzahl Dezimalzahlen : ', len(werte))
            continue
        else:
            #print(werte)
            x = float(werte[0])
            y = float(werte[1])
            z = float(werte[2])
            #print(x, y, z)

            #überlese identische punkte ?
            # ist bei scipy.spatial.Delaunay und triangle nicht notwendig - werden nicht berücksichtigt
            if ([x,y] in kontrolle_):continue

            punkte_xyz.append([x, y, z])
            punkte_xy.append([x, y])
            x_.append(x)
            y_.append(y)
            z_.append(z)
            mark_.append(0)
    return punkte_xyz, punkte_xy, x_, y_, z_, mark_

def trianguliere(punkte_xy_,  mark_, segmente_, holes_):
    #
    #print('Anzahl Punkte in Koordinatenliste : ', len(punkte_xy_))
    file_text.insert(tk.END, 'Anzahl Punkte in Koordinatenliste : ' + str(len(punkte_xy_)) +'\n')
    # Triangle
    #print("Die Triangulierung wird gestartet")
    file_text.insert(tk.END, 'Die Triangulierung wird gestartet ' + '\n')
    root.update_idletasks()
    #
    startzeit_3 = time.perf_counter()
    if len(segmente_) > 0 and len(holes_) > 0 :
        tri_ = dict(vertices=punkte_xy_ , vertex_markers=mark_ , segments=segmente_, holes = holes_)
        tri2 = triangle.triangulate(tri_, 'p')
    else:
        if len(segmente_) > 0:
            tri_ = dict(vertices=punkte_xy_ , vertex_markers=mark_ ,  segments=segmente_)
            tri2 = triangle.triangulate(tri_, 'p')
        else:
            tri_ = dict(vertices=punkte_xy_)
            tri2 = triangle.triangulate(tri_)

    #print('Die Triangulierung erfolgte in : ', (time.perf_counter() - startzeit_3) * 1000, ' Milli-Sekunden')
    file_text.insert(tk.END, 'Die Triangulierung erfolgte in    :  {0:.3f}'.format((time.perf_counter() - startzeit_3) * 1000) + ' Milli-Sekunden' + '\n')

    try:
        tri = tri2['triangles'].tolist()
        #print('Anzahl Dreicksmaschen in der Triangulation : ', len(tri2['triangles'].tolist()))
        file_text.insert(tk.END, 'Anzahl Dreicksmaschen in der Triangulation : ' + str(len(tri)) +'\n')
    except:
        tri = []
        #print("Keine Dreiecksmaschen vorhanden nach der Triangulierung")
        file_text.insert(tk.END, 'Keine Dreiecksmaschen vorhanden nach der Triangulierung' +'\n')
    #print('Die ersten 5 Maschen sind : ', tri[0:5] )
    root.update_idletasks()
    return tri


def maschen_bilden(tri_, punkte_xyz_):
    summe_grund_fl = 0
    summe_ober_fl = 0
    maschen_koordinaten = []
    for i,masche in enumerate(tri_[:]):
        #print(f'die ',i+1 ,'. Masche in tri : ',masche)
        koordinaten_masche=[]
        for j,punkt_id in enumerate(masche):
            #print('Punkt_id ist : ',punkt_id)

            punkt_xyz = punkte_xyz_[punkt_id]
            #print('Die Koordinate in punkte_xyz ist : ',punkt_xyz)

            koordinaten_masche.append(punkt_xyz)

            #x_p = punkt_xyz[0]
            #y_p = punkt_xyz[1]
            #z_p = punkt_xyz[2]
            #print(f' {x_p = }  {y_p = }  {z_p = }')

            #koordinaten_p = tri.points[tri.simplices[i, j], :] # ist die 2d-Koordinate im np_array
            #print(koordinaten_p)
            
        #print('Die Koordinaten der ', i+1, '. Masche sind : ',koordinaten_masche)
        
        # der nachfolgende test kann entfallen, weil:
        #  in scipy und triangle gilt: for 2-D, the points are oriented counterclockwise.  
        #test_im_guz = triangle_det_2d(koordinaten_masche)
        #if test_im_guz < 0.0: print('Punktanordnung im UZ : ', test_im_guz)

        # flaechenberechnungen
        fl_grund_ = abs(triangle_det_2d(koordinaten_masche))

        fl_masche_ = flaeche_masche(koordinaten_masche)

        summe_grund_fl = summe_grund_fl + fl_grund_

        summe_ober_fl = summe_ober_fl + fl_masche_

        maschen_koordinaten.append(koordinaten_masche)
    #print(maschen_koordinaten[0:5])
    #print('Anzahl Maschen in maschen_koordinaten      : ', len(maschen_koordinaten))

    #print('Grundfläche des DGM : {0:.3f}'.format(summe_grund_fl))

    #print('Oberfläche  des DGM : {0:.3f}'.format(summe_ober_fl))
    file_text.insert(tk.END, 'Anzahl Maschen in maschen_koordinaten      : ' + str(len(maschen_koordinaten)) + '\n')
    file_text.insert(tk.END, 'Grundflaeche des DGM : {0:.3f}'.format(summe_grund_fl)  + '\n')
    file_text.insert(tk.END, 'Oberflaeche  des DGM : {0:.3f}'.format(summe_ober_fl)  + '\n')
    root.update_idletasks()
    return maschen_koordinaten


def bereich(x_, y_, z_):
    x_min = min(x_)
    y_min = min(y_)
    z_min = min(z_)

    x_max = max(x_)
    y_max = max(y_)
    z_max = max(z_)

    diff = []
    diff.append(abs(x_max - x_min))
    diff.append(abs(y_max - y_min))
    diff.append(abs(z_max - z_min))
    diff_max = max(diff)

    #print('Koordinatenbereich des DGM :')
    #print(f' {x_min = } {y_min = } {z_min = }')
    #print(f' {x_max = } {y_max = } {z_max = }')
    #print(f' { diff_max = }')
    file_text.insert(tk.END, 'Koordinatenbereich des Roh-DGM :'  + '\n')
    file_text.insert(tk.END, f' {x_min = } {y_min = } {z_min = }'  + '\n')
    file_text.insert(tk.END, f' {x_max = } {y_max = } {z_max = }'  + '\n')
    file_text.insert(tk.END, '\n')    
    return x_min, y_min, z_min, x_max, y_max, z_max, diff_max

def plot_2d(tri_):
    #if len(maschen_koordinaten) < 20000:
    #_ = delaunay_plot_2d(tri_)
    #_ = convex_hull_plot_2d(hull)
    #plt.show()
    pass

def plot_3d(maschen_koordinaten, x_min, y_min, z_min, x_max, y_max, z_max, diff_max):
    #if len(maschen_koordinaten) < 20000:
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    #from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    #from matplotlib import cm

    ax.add_collection(Poly3DCollection(maschen_koordinaten))

    ax.set_xlim([int((x_max + x_min)/2.0 - diff_max/2.0 - 1.0),int((x_max + x_min)/2.0 + diff_max/2.0 + 1.0)])
    ax.set_ylim([int((y_max + y_min)/2.0 - diff_max/2.0 - 1.0),int((y_max + y_min)/2.0 + diff_max/2.0 + 1.0)])
    ax.set_zlim([int((z_max + z_min)/2.0 - diff_max/2.0 - 1.0),int((z_max + z_min)/2.0 + diff_max/2.0 + 1.0)])

    plt.show()

def plot_3d_ms2024(maschen_koordinaten_, mu_):
    global ACTIVEMODEL
    global selected_level
    global dgm_color
    global g_1mu
    global g_go    
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu
    
    #dgnModel = ACTIVEMODEL.GetDgnModel()
    #dgnfile = dgnModel.GetDgnFile()

    #rgb_color = hex_to_rgb(dgm_color)
    #rgb_str = rgb_to_string (rgb_color)
    #print("rgb_str = ", rgb_str)

    level_name = str(selected_level.get())
    levelList,levelListIDs = GetLevelList()
    levelID = 0
    if level_name in levelList:
        j = levelList.index(level_name)
        levelId = levelListIDs[j]
        
    dgm_color = w_color.get()

    z_min_limit = -9000.0

    PyCadInputQueue.SendKeyin("mark")
    button_undo_mark.config(state='active')

    zaehler_ = 0
    summe_grund_fl = 0
    summe_ober_fl = 0
    volumen_ueber_0 = 0  
    x_cad_ = [] 
    y_cad_ = [] 
    z_cad_ = [] 
    for masche in maschen_koordinaten_[:]:
        if ((masche[0][2] < z_min_limit) or (masche[1][2] < z_min_limit) or (masche[2][2] < z_min_limit)): continue
        test = createShapeElement(masche, mu_, levelId, dgm_color)

        fl_grund_ = abs(triangle_det_2d(masche))
        fl_masche_ = flaeche_masche(masche)
        summe_grund_fl = summe_grund_fl + fl_grund_
        summe_ober_fl = summe_ober_fl + fl_masche_
        volumen_ueber_0 = volumen_ueber_0 + fl_grund_ * ((masche[0][2] + masche[1][2] + masche[2][2]) / 3.0)
        x_cad_.append(masche[0][0])
        x_cad_.append(masche[1][0])
        x_cad_.append(masche[2][0])

        y_cad_.append(masche[0][1])
        y_cad_.append(masche[1][1])
        y_cad_.append(masche[2][1])

        z_cad_.append(masche[0][2])
        z_cad_.append(masche[1][2])
        z_cad_.append(masche[2][2])        
        zaehler_ = zaehler_ + 1

    #print('Anzahl Maschen in CAD uebernommen : ', zaehler_)
    file_text.insert(tk.END, 'Anzahl Maschen in CAD uebernommen : ' + str(zaehler_) + '\n')
    file_text.insert(tk.END, 'folgende Werte sind ohne Massstabskorrekturen (z.B. utm): ' + '\n')
    file_text.insert(tk.END, 'Grundflaeche des DGM              : {0:_.3f}'.format(summe_grund_fl)  + '\n')
    file_text.insert(tk.END, 'Oberflaeche  des DGM              : {0:_.3f}'.format(summe_ober_fl)  + '\n')
    file_text.insert(tk.END, 'Volumen ausgeglichen, Hoehe  0.00 : {0:_.3f}'.format(volumen_ueber_0)  + '\n')
    file_text.insert(tk.END, 'mittlere Hoehe                    : {0:_.6f}'.format(volumen_ueber_0/summe_grund_fl)  + '\n')
    file_text.insert(tk.END, 'Koordinatenbereich des DGM in CAD :'  + '\n')
    file_text.insert(tk.END, f' {min(x_cad_) = } {min(y_cad_) = } {min(z_cad_) = }'  + '\n')
    file_text.insert(tk.END, f' {max(x_cad_) = } {max(y_cad_) = } {max(z_cad_) = }'  + '\n')
    file_text.insert(tk.END, '\n')
    root.update_idletasks()     
    return True

# Create Shape Element function
def createShapeElement(masche, mu_, level_id_, dgm_color_ ):
    global ACTIVEMODEL
    global selected_level
    global dgm_color
    global g_1mu
    global g_go
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu

    shape_eeh = EditElementHandle()
    #global origin
    #print(g_go, g_go.x*-1/g_1mu, g_go.y*-1/g_1mu, g_go.z*-1/g_1mu)

    # Set Shape coordinates
    points = DPoint3dArray()
        
    #points.append (DPoint3d (masche[0][0] * mu_ + g_go.x, masche[0][1] * mu_ + g_go.y, masche[0][2] * mu_ + g_go.z))
    #points.append (DPoint3d (masche[1][0] * mu_ + g_go.x, masche[1][1] * mu_ + g_go.y, masche[1][2] * mu_ + g_go.z))
    #points.append (DPoint3d (masche[2][0] * mu_ + g_go.x, masche[2][1] * mu_ + g_go.y, masche[2][2] * mu_ + g_go.z))
    #points.append (DPoint3d (masche[0][0] * mu_ + g_go.x, masche[0][1] * mu_ + g_go.y, masche[0][2] * mu_ + g_go.z))

    points.append (DPoint3d (masche[0][0] * g_1mu + g_go.x, masche[0][1] * g_1mu + g_go.y, masche[0][2] * g_1mu + g_go.z))
    points.append (DPoint3d (masche[1][0] * g_1mu + g_go.x, masche[1][1] * g_1mu + g_go.y, masche[1][2] * g_1mu + g_go.z))
    points.append (DPoint3d (masche[2][0] * g_1mu + g_go.x, masche[2][1] * g_1mu + g_go.y, masche[2][2] * g_1mu + g_go.z))
    points.append (DPoint3d (masche[0][0] * g_1mu + g_go.x, masche[0][1] * g_1mu + g_go.y, masche[0][2] * g_1mu + g_go.z))

    
    # Create Shape element
    status = ShapeHandler.CreateShapeElement (shape_eeh, None, points, 
                                               ACTIVEMODEL.Is3d(), ACTIVEMODEL)

    if BentleyStatus.eSUCCESS != status:
        return False

    color = dgm_color_
    #rgb_color = hex_to_rgb(dgm_color)
    
    lineStyle = 0
    lineWeight = 0
    propertiesSetter = ElementPropertiesSetter()
    propertiesSetter.SetColor(color)
    propertiesSetter.SetLinestyle(lineStyle, None)
    propertiesSetter.SetWeight(lineWeight)
    propertiesSetter.SetLevel(level_id_)
    
    propertiesSetter.Apply(shape_eeh)        

    # Add the Shape element to model
    if BentleyStatus.eSUCCESS != shape_eeh.AddToModel():
        return False

    return True

# hoehenlinie in CAD
def createLineElement(punkt_von, punkt_nach, hl_ebene_id, color, style, weight):
    global ACTIVEMODEL
    global g_1mu
    global g_go
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu
    
    segment = DSegment3d(DPoint3d(punkt_von[0] * g_1mu + g_go.x, punkt_von[1] * g_1mu + g_go.y, punkt_von[2] * g_1mu + g_go.z),
                         DPoint3d(punkt_nach[0] * g_1mu + g_go.x, punkt_nach[1] * g_1mu + g_go.y, punkt_nach[2] * g_1mu + g_go.z))
    line_eeh = EditElementHandle()
    
    status = LineHandler.CreateLineElement(line_eeh, None, segment, dgnModel.Is3d(), dgnModel)

    if BentleyStatus.eSUCCESS != status:
        return False
    
    propertiesSetter = ElementPropertiesSetter()
    propertiesSetter.SetColor(color)
    propertiesSetter.SetLinestyle(style, None)
    propertiesSetter.SetWeight(weight)
    propertiesSetter.SetLevel(hl_ebene_id)
    
    propertiesSetter.Apply(line_eeh)        

    # Add the line element to model
    if BentleyStatus.eSUCCESS != line_eeh.AddToModel():
        return False
    
    return True

def createLineStringElement(punktliste, hl_ebene_id, color, style, weight):
    global ACTIVEMODEL
    global g_1mu
    global g_go
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu
    
    points = DPoint3dArray()
    
    for punkt in punktliste:
        #print(punkt)
        points.append(DPoint3d(punkt[0] * g_1mu + g_go.x, punkt[1] * g_1mu + g_go.y, punkt[2] * g_1mu + g_go.z))

    linestring_eeh = EditElementHandle()
    
    #status = LineHandler.CreateLineElement(line_eeh, None, segment, dgnModel.Is3d(), dgnModel)
    status = LineStringHandler.CreateLineStringElement(linestring_eeh, None, points, dgnModel.Is3d(), dgnModel)

    if BentleyStatus.eSUCCESS != status:
        return False
    
    propertiesSetter = ElementPropertiesSetter()
    propertiesSetter.SetColor(color)
    propertiesSetter.SetLinestyle(style, None)
    propertiesSetter.SetWeight(weight)
    propertiesSetter.SetLevel(hl_ebene_id)
    
    propertiesSetter.Apply(linestring_eeh)        

    # Add the linestring element to model
    if BentleyStatus.eSUCCESS != linestring_eeh.AddToModel():
        return False
    
    return True

def createShapeElement_2(punktliste, hl_ebene_id, color, style, weight):
    global ACTIVEMODEL
    global g_1mu
    global g_go
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu
    
    points = DPoint3dArray()
    
    for punkt in punktliste:
        #print(punkt)
        points.append(DPoint3d(punkt[0] * g_1mu + g_go.x, punkt[1] * g_1mu + g_go.y, punkt[2] * g_1mu + g_go.z))

    shape_eeh = EditElementHandle()
    
    #status = LineHandler.CreateLineElement(line_eeh, None, segment, dgnModel.Is3d(), dgnModel)
    #status = LineStringHandler.CreateLineStringElement(linestring_eeh, None, points, dgnModel.Is3d(), dgnModel)
    status = ShapeHandler.CreateShapeElement (shape_eeh, None, points, ACTIVEMODEL.Is3d(), ACTIVEMODEL)

    if BentleyStatus.eSUCCESS != status:
        return False
    
    propertiesSetter = ElementPropertiesSetter()
    propertiesSetter.SetColor(color)
    propertiesSetter.SetLinestyle(style, None)
    propertiesSetter.SetWeight(weight)
    propertiesSetter.SetLevel(hl_ebene_id)
    
    propertiesSetter.Apply(shape_eeh)        

    # Add the shape element to model
    if BentleyStatus.eSUCCESS != shape_eeh.AddToModel():
        return False
    
    return True

def createBsplineElement(punktliste, hl_ebene_id, color, style, weight):
    global ACTIVEMODEL
    global g_1mu
    global g_go
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu
    
    points = DPoint3dArray()
    
    for punkt in punktliste:
        #print(punkt)
        points.append(DPoint3d(punkt[0] * g_1mu + g_go.x, punkt[1] * g_1mu + g_go.y, punkt[2] * g_1mu + g_go.z))
    
    # Set the points to draw curve.
    curve = MSBsplineCurve.CreateFromPolesAndOrder(points, None, None, 3, False, False)
    # Draw the curve based of the points.
    ICurvePrimitive.CreateBsplineCurve(curve)
    
    bspLine_eeh = EditElementHandle()
    
    status = BSplineCurveHandler.CreateBSplineCurveElement(bspLine_eeh, None, curve, dgnModel.Is3d(), dgnModel)

    if BentleyStatus.eSUCCESS != status:
        return False

    propertiesSetter = ElementPropertiesSetter()
    propertiesSetter.SetColor(color)
    propertiesSetter.SetLinestyle(style, None)
    propertiesSetter.SetWeight(weight)
    propertiesSetter.SetLevel(hl_ebene_id)
    
    propertiesSetter.Apply(bspLine_eeh)        

    # Add the linestring element to model
    if BentleyStatus.eSUCCESS != bspLine_eeh.AddToModel():
        return False

    return True

def main():
    
    file_text.delete("1.0", tk.END)
    root.update_idletasks()

    liste_umring = []
    segments_o = []
    punkte_xyz_o = []
    punkte_xy_o = []
    x_o = []
    y_o = []
    z_o = []    
    mark_o=[]
    Datei_umring = open_file_dialog()
    if Datei_umring != '':
        liste_umring = daten_einlesen(Datei_umring)
        #print('Anzahl Zeilen im Umring : ', len(liste_umring)) 
        file_text.insert(tk.END, 'Anzahl Zeilen im Umring           : ' + str(len(liste_umring)) + '\n')
        root.update_idletasks()
        #
        punkte_xyz_o, punkte_xy_o, x_o, y_o, z_o, mark_o = daten_aufbereiten(liste_umring,[])
        # alle doppelten punkte entfernen:  list -> set -> list
        punkte_set=set()
        for punkt in punkte_xyz_o:
            punkte_set.add((punkt[0], punkt[1], punkt[2]))
            #print(punkt[0], punkt[1], punkt[2])
        punkte_list_o = []
        for item in punkte_set:punkte_list_o.append([item[0], item[1], item[2]])
        #print(len(punkte_xyz_o), len(punkte_set), len(punkte_list_o))
        #print(punkte_list_o[0:5], punkte_xyz_o[0:5])

        if len(punkte_xyz_o) > 0:
            lp_=len(punkte_xyz_o)
            for i in range(0, lp_, 2):
                punkt_1 = punkte_xyz_o[i]
                j=0
                if punkt_1 in punkte_list_o:
                    j = punkte_list_o.index(punkt_1)
                punkt_2 = punkte_xyz_o[i+1]
                k=0
                if punkt_2 in punkte_list_o:
                    k = punkte_list_o.index(punkt_2)
               
                segments_o.append((j, k))
                #print(punkt_1, punkt_2)
            #print(segments_o[:5], ' ... ', segments_o[-5:])
            punkte_xyz_o = []
            punkte_xy_o = []
            x_o = []
            y_o = []
            z_o = []
            mark_o = []
            for m in range(len(punkte_list_o)):
                punkt_ = punkte_list_o[m]
                punkte_xyz_o.append([punkt_[0], punkt_[1], punkt_[2]])
                punkte_xy_o.append([punkt_[0], punkt_[1]])
                x_o.append(punkt_[0])
                y_o.append(punkt_[1])
                z_o.append(punkt_[2])
                mark_o.append(1)
    # innere segmente (lines) einlesen
    segments_s = []
    punkte_xyz_o_s = punkte_xyz_o + []
    punkte_xyz_s = []
    punkte_xy_s = []
    x_s = []
    y_s = []
    z_s = []    
    mark_s=[]    
    Datei_innen = open_file_dialog_2()
    if Datei_innen != '':
        liste_innen = daten_einlesen(Datei_innen)
        #print('Anzahl Zeilen Segmente innen : ', len(liste_innen)) 
        file_text.insert(tk.END, 'Anzahl Zeilen Segmente innen      : ' + str(len(liste_innen)) + '\n')
        root.update_idletasks()
        #           
        punkte_xyz_s, punkte_xy_s, x_s, y_s, z_s, mark_s = daten_aufbereiten(liste_innen,[])
        # alle doppelten punkte entfernen:  list -> set -> list 
        punkte_set_s=set()
        for punkt in punkte_xyz_s:
            punkte_set_s.add((punkt[0], punkt[1], punkt[2]))
            #print(punkt[0], punkt[1], punkt[2])
        punkte_xyz_o_s = punkte_xyz_o + []                    # die umringskoordinaten werden am ende ergänzt 
        punkte_xy_s = []
        x_s = []
        y_s = []
        z_s = []    
        mark_s=[]
        for item in punkte_set_s:
            punkt_temp = [item[0], item[1], item[2]]
            if punkt_temp in punkte_xyz_o:continue             # überlese punkte im umring
            punkte_xyz_o_s.append([item[0], item[1], item[2]])
            punkte_xy_s.append([item[0], item[1]])
            x_s.append(item[0])
            y_s.append(item[1])
            z_s.append(item[2])
            mark_s.append(0)

        #print(len(punkte_xyz_s), len(punkte_set_s), len(punkte_list_o_s))
        #print(punkte_list_o_s[0:5], punkte_xyz_o_s[0:5])
        if len(punkte_xyz_s) > 0:
            lp_=len(punkte_xyz_s)
            for i in range(0, lp_, 2):               # jedes segment hat zwei koodinaten triple
                punkt_1 = punkte_xyz_s[i]
                j=0
                if punkt_1 in punkte_xyz_o_s:
                    j = punkte_xyz_o_s.index(punkt_1)
                punkt_2 = punkte_xyz_s[i+1]
                k=0
                if punkt_2 in punkte_xyz_o_s:
                    k = punkte_xyz_o_s.index(punkt_2)
               
                segments_s.append((j, k))
                #print(punkt_1, punkt_2)
            #print(segments_s[:5], ' ... ', segments_s[-5:])

    # koordinatenpunkte einlesen
    Dateien = open_files_dialog()

    if len(Dateien) == 0 and Datei_umring == '':
        #print('Programmabbruch: keine Datei gewählt')
        file_text.insert(tk.END, 'Programmabbruch: keine Koordinaten-Datei gewählt' + '\n')
        root.update_idletasks()
        return
        #sys.exit()
    
    #liste = daten_einlesen(Datei)
    liste_p = []
    liste_p = dateien_einlesen(Dateien)
    #print('Anzahl Zeilen in der Punkte-Datei:',len(liste_p))
    file_text.insert(tk.END, 'Anzahl Zeilen in der Punkte-Datei : ' + str(len(liste_p)) + '\n')
    root.update_idletasks()
    #
    liste = liste_p
    #           print('Die Zeilen 0:5 sind : ', liste[0:5])      
    punkte_xyz_i, punkte_xy_i, x_i, y_i, z_i, mark_i = daten_aufbereiten(liste, (punkte_xy_o + punkte_xy_s))

    punkte_xyz = punkte_xyz_o_s + punkte_xyz_i
    punkte_xy = punkte_xy_o + punkte_xy_s + punkte_xy_i
    x_ = x_o + x_s + x_i
    y_ = y_o + y_s + y_i
    z_ = z_o + z_s + z_i
    mark_ = mark_o + mark_s + mark_i
    segments_ = segments_o + segments_s

    #print(punkte_xy[:5], punkte_xy_o[:5])
    if len(punkte_xy) < 3:
        
        #print('Programmabbruch: Die Dateien sind keine Koordinatendateien')
        file_text.insert(tk.END, 'Programmabbruch: Die Dateien sind keine Koordinatendateien' + '\n')
        root.update_idletasks()
        return
        #sys.exit()
    #           print('Anzahl Punkte in Koordinatenliste : ', len(punkte_xyz))

    hole_points=[]    
    Datei_hole = open_file_dialog_3()
    if Datei_hole != '':
        liste_holes = daten_einlesen(Datei_hole)
        #print('Anzahl Zeilen in holes: ', len(liste_holes))
        file_text.insert(tk.END, 'Anzahl Zeilen in holes            : ' + str(len(liste_holes)) + '\n')
        root.update_idletasks()    
        punkte_xyz_h, punkte_xy_h, x_h, y_h, z_h, mark_h = daten_aufbereiten(liste_holes,[])
        hole_points = punkte_xy_h
        #print(hole_points)

    startzeit = time.perf_counter()
    tri = trianguliere(punkte_xy, mark_, segments_, hole_points)

    if len(tri) > 0:

        maschen_koordinaten = maschen_bilden(tri, punkte_xyz)

        x_min, y_min, z_min, x_max, y_max, z_max, diff_max = bereich(x_, y_, z_)
        #
        #print('Die Aufgabe erledigt in: ', (time.perf_counter() - startzeit) * 1000, ' Milli-Sekunden')
        file_text.insert(tk.END, 'Die Aufgabe erledigt in:    :  {0:.3f}'.format((time.perf_counter() - startzeit) * 1000) + ' Milli-Sekunden' + '\n')
        file_text.insert(tk.END, '\n')
        root.update_idletasks()
        
        #if len(maschen_koordinaten) < 10000: plot_2d(tri)
        #if len(maschen_koordinaten) < 10000: plot_3d(maschen_koordinaten, x_min, y_min, z_min, x_max, y_max, z_max, diff_max)
        
        # Global variable declaration
        global ACTIVEMODEL
        global dgnModel
        global g_1mu
        

        test = plot_3d_ms2024(maschen_koordinaten, g_1mu)
        PyCadInputQueue.SendReset()
        PyCadInputQueue.SendKeyin("FIT VIEW EXTENDED")
        PyCommandState.StartDefaultCommand()

        lift_window(root)
        
        #print('Die gesamte Aufgabe erledigt in: ', (time.perf_counter() - startzeit) * 1000, ' Milli-Sekunden')
        file_text.insert(tk.END, 'Die gesamte Aufgabe erledigt in:   :  {0:.3f}'.format((time.perf_counter() - startzeit) * 1000) + ' Milli-Sekunden' + '\n')
        file_text.insert(tk.END, '\n')
        root.update_idletasks()
        

    return 

def pick_color():
    global dgm_color
    dgm_color = colorchooser.askcolor()[1] 
    dgm_color = dgm_color.upper()
    if dgm_color:
        color_label.config(text = dgm_color)
        color_label.config(bg = dgm_color)

# RGB to String function
def rgb_to_string(rgb_tuple):
    return ', '.join(map(str, rgb_tuple))


# Hex to RGB function
def hex_to_rgb(hex_color):
    # Remove the hash symbol if present
    hex_color = hex_color.lstrip('#')
    
    # Convert hex to RGB
    rgb = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
    return rgb

# Get Level Name List
def GetLevelList ():
    level_List = []
    level_Id_List = []

    ACTIVEMODEL = ISessionMgr.ActiveDgnModelRef
    levelCache = ACTIVEMODEL.GetDgnFile().GetLevelCache()

    it = iter(levelCache)
    for level in it:
        levelName = level.GetName()
        levelID = level.GetLevelId()
        level_List.append(levelName)
        level_Id_List.append(levelID)

    return level_List, level_Id_List

# tkinter fenster ganz nach vorne bringen
def lift_window(window):
    window.attributes('-topmost', True)
    window.update_idletasks()  # get window on top
    window.attributes('-topmost', False)  # prevent permanent focus 
    window.focus_force()  # focus to the window

#tkinter sliderwert aendert sich
def farbe_aendern(i):
    dgm_color = w_color.get()
    #color_button.config(text = farben_[dgm_color], bg=farben_[dgm_color])
    color_button.config(text = dgm_color, bg=farben_[dgm_color])
    farbe.config(bg=farben_[dgm_color], text = dgm_color)
    
def farb_pick(k):
    global farbe_fuer
    if farbe_fuer == 0:
        dgm_color = k
        farbe.config(bg=farben_[dgm_color], text = dgm_color)
        w_color.set(k)
    elif farbe_fuer == 1:
        color_button_hl_0.config(bg=farben_[k], text = k)
        farbe.config(bg=farben_[k], text = k)
    elif farbe_fuer == 2:
        color_button_hl_2.config(bg=farben_[k], text = k)
        farbe.config(bg=farben_[k], text = k)
    elif farbe_fuer == 3:
        color_button_hl_3.config(bg=farben_[k], text = k)
        farbe.config(bg=farben_[k], text = k)
    elif farbe_fuer == 4:
        color_button_hl_4.config(bg=farben_[k], text = k)
        farbe.config(bg=farben_[k], text = k) 
    elif farbe_fuer == 5:
        color_button_um_ob.config(bg=farben_[k], text = k)
        farbe.config(bg=farben_[k], text = k) 
    elif farbe_fuer == 6:
        color_button_um_un.config(bg=farben_[k], text = k)
        farbe.config(bg=farben_[k], text = k) 
    elif farbe_fuer == 7:
        color_button_um_sei.config(bg=farben_[k], text = k)
        farbe.config(bg=farben_[k], text = k)                                              
    
def select_tab():
    global farbe_fuer
    farbe_fuer = 0
    notebook.select(frameFarben)
    
def select_tab_h0():
    global farbe_fuer
    farbe_fuer = 1
    #hl_0_Color = int(color_button_hl_0['text'])
    hl_0_Color = int(color_button_hl_0.cget('text'))
    farbe.config(bg=farben_[hl_0_Color], text = hl_0_Color)
    notebook.select(frameFarben)
    
def select_tab_h2():
    global farbe_fuer
    farbe_fuer = 2
    #hl_2_Color = int(color_button_hl_0['text'])
    hl_2_Color = int(color_button_hl_2.cget('text'))
    farbe.config(bg=farben_[hl_2_Color], text = hl_2_Color)
    notebook.select(frameFarben) 
    
def select_tab_h3():
    global farbe_fuer
    farbe_fuer = 3
    #hl_2_Color = int(color_button_hl_0['text'])
    hl_3_Color = int(color_button_hl_3.cget('text'))
    farbe.config(bg=farben_[hl_3_Color], text = hl_3_Color)
    notebook.select(frameFarben) 
    
def select_tab_h4():
    global farbe_fuer
    farbe_fuer = 4
    #hl_2_Color = int(color_button_hl_0['text'])
    hl_4_Color = int(color_button_hl_4.cget('text'))
    farbe.config(bg=farben_[hl_4_Color], text = hl_4_Color)
    notebook.select(frameFarben) 
    
def select_tab_um_ob():
    global farbe_fuer
    farbe_fuer = 5
    um_ob_Color = int(color_button_um_ob.cget('text'))
    farbe.config(bg=farben_[um_ob_Color], text = um_ob_Color)
    notebook.select(frameFarben) 
    
def select_tab_um_un():
    global farbe_fuer
    farbe_fuer = 6
    um_un_Color = int(color_button_um_un.cget('text'))
    farbe.config(bg=farben_[um_un_Color], text = um_un_Color)
    notebook.select(frameFarben)  

def select_tab_um_sei():
    global farbe_fuer
    farbe_fuer = 7
    um_sei_Color = int(color_button_um_sei.cget('text'))
    farbe.config(bg=farben_[um_sei_Color], text = um_sei_Color)
    notebook.select(frameFarben)    
                         
    
def select_tab_H():
    global farbe_fuer
    if farbe_fuer == 0:
        notebook.select(frameHaupt)   # DGM
    elif farbe_fuer == 1 or farbe_fuer == 2 or farbe_fuer == 3 or farbe_fuer == 4:
        notebook.select(frameHoehenlinien)  # hoehenlinien
    elif farbe_fuer == 5 or farbe_fuer == 6 or farbe_fuer == 7:
        notebook.select(frameUmring)  # umring        
    return
        
def undo_mark():
    PyCadInputQueue.SendKeyin("undo mark")
    button_undo_mark.config(state='disabled')
    lift_window(root)
    
def trim_ellipse_one_point():
    PyCadInputQueue.SendKeyin("trim break bypoint")
    PyCadInputQueue.SendReset()
    PyCommandState.StartDefaultCommand()
    #PyCadInputQueue.SendKeyin("null")
    lift_window(root)

def item_selected():
    for selected_item in tree_1.selection():
        item = tree_1.item(selected_item)
        record = item['values']
        # show a message
        #message=','.join(record)
        # showinfo(title='Information', message=','.join(record))
        #print(record)

# punkte als tuple punkt.x, punkt.y, punkt.z
# faktor fuer utm
def strecke_2d_f(punkt_a, punkt_b, faktor ):
    global faktor_projektion
    global hoehen_werte
    #faktor_projektion = ( (1.0 + ((ost_mittel_km - koord_meridian_km)**2)/(2.0 * radius_km**2)) * utm_faktor )
    dx_ = punkt_b[0] - punkt_a[0]
    dy_ = punkt_b[1] - punkt_a[1]
    #dz_ = punkt_b[2] - punkt_a[2]
    #print(dx_, dy_, dz_)
    hoehen_korrektur = str(options_utm_hoehe.get())
    
    if utm_var.get() == True and hoehen_korrektur == hoehen_werte[0]:
        utm_faktor = 0.9996
        undulation_km = float(geoid_un_eingabe.get())/1000.0
        radius_km = float(radius_eingabe.get())
        if radius_km < 5000 or radius_km > 7000:radius_km = 6382
        h_strecke_ellip_km = (punkt_a[2] + punkt_b[2])/2.0/1000.0 + undulation_km
        faktor_hoehe_strecke = 1.0 - (h_strecke_ellip_km / radius_km)
        faktor = 1.0 / (faktor_projektion - (h_strecke_ellip_km / radius_km) * utm_faktor)
        #print(faktor)
        
    s_2d = (dx_ * dx_ + dy_ * dy_ )**0.5
    s_2d_ = s_2d * faktor
    #print(s_2d_)
    return s_2d_

# punkte als tuple punkt.x, punkt.y, punkt.z
# faktor fuer 2d utm
def strecke_3d_f(punkt_a, punkt_b, faktor):
    global faktor_projektion
    global hoehen_werte
    dx_ = punkt_b[0] - punkt_a[0]
    dy_ = punkt_b[1] - punkt_a[1]
    dz_ = punkt_b[2] - punkt_a[2]
    #print(dx_, dy_, dz_)
    hoehen_korrektur = str(options_utm_hoehe.get())
    
    if utm_var.get() == True and hoehen_korrektur == hoehen_werte[0]:
        utm_faktor = 0.9996
        undulation_km = float(geoid_un_eingabe.get())/1000.0
        radius_km = float(radius_eingabe.get())
        if radius_km < 5000 or radius_km > 7000:radius_km = 6382
        h_strecke_ellip_km = (punkt_a[2] + punkt_b[2])/2.0/1000.0 + undulation_km
        faktor_hoehe_strecke = 1.0 - (h_strecke_ellip_km / radius_km)
        faktor = 1.0 / (faktor_projektion - (h_strecke_ellip_km / radius_km) * utm_faktor)
        #print(faktor) 
        
    s_2d = (dx_ * dx_ + dy_ * dy_)**0.5
    s_2d_ = s_2d * faktor
    # z-werte ohne faktor
    s_3d_ = (s_2d_ * s_2d_ + dz_ * dz_)**0.5
    #print(s_3d_)
    return s_3d_

def flaeche_heron(seite_a, seite_b, seite_c):
    # dreiecksflaeche nach satz des heron
    s = (seite_a + seite_b + seite_c) * 0.5
    fl_ = (s * (s - seite_a) * (s - seite_b) * (s - seite_c))**0.5
    #print(fl_)
    return fl_

    
# punkte als tuple punkt.x, punkt.y, punkt.z
def strecke_3d(punkt_a, punkt_b):
    dx_ = punkt_b[0] - punkt_a[0]
    dy_ = punkt_b[1] - punkt_a[1]
    dz_ = punkt_b[2] - punkt_a[2]
    #print(dx_, dy_, dz_)
    s_3d = (dx_ * dx_ + dy_ * dy_ + dz_ * dz_)**0.5
    #print(s_3d)
    return s_3d, (dx_/s_3d, dy_/s_3d, dz_/s_3d)

def zwischenpunkte(punkt_start_xyz_, punkt_ende_xyz_, intervall_):
    punkte_zwischen_ = []

    s_ges, deltas = strecke_3d(punkt_start_xyz_, punkt_ende_xyz_)

    #print(s_ges)
    #print(deltas)

    anz_punkte_zwi = int(s_ges / intervall_) + 1
    #print('3d Strecke : ', s_ges , 'Anzahl Punkte : ', anz_punkte_zwi)
    teil_strecke = s_ges/anz_punkte_zwi
    #print('statt Intervall : ', intervall, ' --> ', teil_strecke)

    for i in range(anz_punkte_zwi + 1):
        s_zwischen = teil_strecke * i
        #print(s_zwischen)
        dx_i = deltas[0] * s_zwischen
        dy_i = deltas[1] * s_zwischen
        dz_i = deltas[2] * s_zwischen
        #print(dx_i, dy_i, dz_i)
        x_i = punkt_start_xyz_[0] + dx_i
        y_i = punkt_start_xyz_[1] + dy_i
        z_i = punkt_start_xyz_[2] + dz_i
        #print(x_i, y_i, z_i)
        punkt_koordinate_s = " {0:.4f} , {1:.4f} , {2:.4f} ".format(x_i, y_i, z_i)
        #print(punkt_koordinate_s)
        punkte_zwischen_.append(punkt_koordinate_s)
    
    return punkte_zwischen_

#Function to select elements by its type 3 (line) or 4 (linestring)
def selectElementsbyType_3_4():
    levelAnzahlElemente=[]
    levelAnzahlLines=[]
    levelAnzahlLinestrings=[]
    for i in range(max(levelListIDs)+1):
        levelAnzahlElemente.append(0)
        levelAnzahlLines.append(0)
        levelAnzahlLinestrings.append(0)
    #print(levelAnzahlElemente)
    
    #Get active model
    ACTIVEMODEL = ISessionMgr.ActiveDgnModelRef
    dgnModel = ACTIVEMODEL.GetDgnModel()
    #name =  model.GetModelName()
    #Get all graphical elements from the model
    graphicalElements = dgnModel.GetGraphicElements()

    for perElementRef in graphicalElements:
        elementId = perElementRef.GetElementId()
        eeh = EditElementHandle(perElementRef, dgnModel)
        eh = ElementHandle(perElementRef)

        msElement = MSElement()
        msElement = eeh.GetElement ()

        isGraphics = msElement.ehdr.isGraphics
        isInvisible = msElement.hdr.dhdr.props.b.invisible

        if (isGraphics and not(isInvisible)):
            eleType = eh.GetElementType()
            levelId = msElement.ehdr.level
            if(eleType == 3) or (eleType == 4):
                #It will select highlight all elements added in selection set
                #selSetManager.AddElement(perElementRef,dgnModel)
                if levelId in levelListIDs:
                    levelAnzahlElemente[levelId] = levelAnzahlElemente[levelId] + 1
                    if eleType == 3:
                        levelAnzahlLines[levelId] = levelAnzahlLines[levelId] + 1
                    if eleType == 4:
                        levelAnzahlLinestrings[levelId] = levelAnzahlLinestrings[levelId] + 1
    # loesche alle items in treeview tree_1
    
    all_root_items = tree_1.get_children()
    tree_1.delete(*all_root_items)
    root.update_idletasks()


    # daten aktualisieren
    ebenen = []
    for i in range(len(levelListIDs)):
        levelId_ = levelListIDs[i]
        if levelAnzahlElemente[levelId_] > 0:
            ebenen.append((levelList[i], levelListIDs[i], levelAnzahlElemente[levelId_], levelAnzahlLines[levelId_], levelAnzahlLinestrings[levelId_]))

    # add data to the treeview
    for ebene in ebenen:
        tree_1.insert('', tk.END, values=ebene) 

    if len(ebenen) > 0: button_2.config(state='active')

    root.update_idletasks()   

    #print(levelAnzahlElemente)


#Function to select elements by its type 3 (line) or 4 (linestring) or 6 (shapes) or 16 (arcs) or 15 (circle ellipse)
def select_2_ElementsbyType_3_4():
    levelAnzahlElemente=[]
    levelAnzahlLines=[]
    levelAnzahlLinestrings=[]
    levelAnzahlShapes=[]
    levelAnzahlArcs=[]
    levelAnzahlEllipses=[]
    levelAnzahlBsplines=[]
    levelAnzahlComplex=[]
    for i in range(max(levelListIDs)+1):
        levelAnzahlElemente.append(0)
        levelAnzahlLines.append(0)
        levelAnzahlLinestrings.append(0)
        levelAnzahlShapes.append(0)
        levelAnzahlArcs.append(0)
        levelAnzahlEllipses.append(0)
        levelAnzahlBsplines.append(0)
        levelAnzahlComplex.append(0)
    #print(levelAnzahlElemente)
    
    #Get active model
    ACTIVEMODEL = ISessionMgr.ActiveDgnModelRef
    dgnModel = ACTIVEMODEL.GetDgnModel()
    #name =  model.GetModelName()
    #Get all graphical elements from the model
    graphicalElements = dgnModel.GetGraphicElements()

    for perElementRef in graphicalElements:
        elementId = perElementRef.GetElementId()
        eeh = EditElementHandle(perElementRef, dgnModel)
        eh = ElementHandle(perElementRef)

        msElement = MSElement()
        msElement = eeh.GetElement ()

        isGraphics = msElement.ehdr.isGraphics
        isInvisible = msElement.hdr.dhdr.props.b.invisible

        if (isGraphics and not(isInvisible)):
            eleType = eh.GetElementType()
            levelId = msElement.ehdr.level
            if(eleType == 3) or (eleType == 4) or (eleType == 6) or (eleType == 16) or (eleType == 15) or (eleType == 27) or (eleType == 12) or (eleType == 14):
                #It will select highlight all elements added in selection set
                #selSetManager.AddElement(perElementRef,dgnModel)
                if levelId in levelListIDs:
                    levelAnzahlElemente[levelId] = levelAnzahlElemente[levelId] + 1
                    if eleType == 3:
                        levelAnzahlLines[levelId] = levelAnzahlLines[levelId] + 1
                    if eleType == 4:
                        levelAnzahlLinestrings[levelId] = levelAnzahlLinestrings[levelId] + 1
                    if eleType == 6:
                        levelAnzahlShapes[levelId] = levelAnzahlShapes[levelId] + 1 
                    if eleType == 15:
                        levelAnzahlEllipses[levelId] = levelAnzahlEllipses[levelId] + 1 
                    if eleType == 16:
                        levelAnzahlArcs[levelId] = levelAnzahlArcs[levelId] + 1 
                    if eleType == 27:
                        levelAnzahlBsplines[levelId] = levelAnzahlBsplines[levelId] + 1
                    if eleType == 12 or  eleType == 14:
                        levelAnzahlComplex[levelId] = levelAnzahlComplex[levelId] + 1
                                                                                          
    # loesche alle items in treeview tree_2
    all_root_items = tree_2.get_children()
    tree_2.delete(*all_root_items)
    root.update_idletasks()


    # daten aktualisieren
    ebenen = []
    for i in range(len(levelListIDs)):
        levelId_ = levelListIDs[i]
        if levelAnzahlElemente[levelId_] > 0:
            ebenen.append((levelList[i], levelListIDs[i], levelAnzahlElemente[levelId_], levelAnzahlLines[levelId_],
             levelAnzahlLinestrings[levelId_], levelAnzahlShapes[levelId_], levelAnzahlArcs[levelId_], levelAnzahlEllipses[levelId_], levelAnzahlBsplines[levelId_], levelAnzahlComplex[levelId_]))

    # add data to the treeview
    for ebene in ebenen:
        tree_2.insert('', tk.END, values=ebene) 

    if len(ebenen) > 0:
        button_6.config(state='active')
        button_8.config(state='active')

    root.update_idletasks()   
   

def exportElementsbyType_3_4():
    global g_1mu
    global g_go    
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu    
    koordinaten_liste_=[]
    ebenenId_selected = []

    if len(tree_1.selection()) == 0: return

    for selected_item in tree_1.selection():
        item = tree_1.item(selected_item)
        record = item['values']
        id_=record[1]
        ebenenId_selected.append(id_)
    
    #print(ebenenId_selected)

    start = DPoint3d()
    end = DPoint3d()

    
    #Get active model
    ACTIVEMODEL = ISessionMgr.ActiveDgnModelRef
    dgnModel = ACTIVEMODEL.GetDgnModel()
    #Get all graphical elements from the model
    graphicalElements = dgnModel.GetGraphicElements()

    for perElementRef in graphicalElements:
        elementId = perElementRef.GetElementId()
        eeh = EditElementHandle(perElementRef, dgnModel)
        eh = ElementHandle(perElementRef)

        msElement = MSElement()
        msElement = eeh.GetElement ()

        isGraphics = msElement.ehdr.isGraphics
        isInvisible = msElement.hdr.dhdr.props.b.invisible

        if (isGraphics and not(isInvisible)):
            eleType = eh.GetElementType()
            levelId = msElement.ehdr.level
            if levelId in ebenenId_selected:
                if (eleType == 3) or (eleType == 4):
                    curve = ICurvePathQuery.ElementToCurveVector(eh)
                    curve.GetStartEnd(start, end)

                    if eleType == 3:
                        #print(" {0:.4f} , {1:.4f} , {2:.4f} ".format((start.x/g_1mu),(start.y/g_1mu),(start.z/g_1mu)))
                        #print(" {0:.4f} , {1:.4f} , {2:.4f} ".format((end.x/g_1mu), (end.y/g_1mu), (end.z/g_1mu)))
                        #print(g_go, g_go.x*-1/g_1mu, g_go.y*-1/g_1mu, g_go.z*-1/g_1mu)
                        koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((start.x/g_1mu + g_go.x*-1/g_1mu),(start.y/g_1mu + g_go.y*-1/g_1mu),(start.z/g_1mu + g_go.z*-1/g_1mu)))
                        koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((end.x/g_1mu + g_go.x*-1/g_1mu), (end.y/g_1mu + g_go.y*-1/g_1mu), (end.z/g_1mu + g_go.z*-1/g_1mu)))

                    if eleType == 4:
                        #print('linestring')
                        for element in curve:
                            points = element.GetLineString()
                            
                            for i in range(len(points)-1):
                                point_a = points[i]
                                point_b = points[i+1]
                                #print(" {0:.4f} , {1:.4f} , {2:.4f} ".format((point_a.x/g_1mu),(point_a.y/g_1mu),(point_a.z/g_1mu)))
                                #print(" {0:.4f} , {1:.4f} , {2:.4f} ".format((point_b.x/g_1mu), (point_b.y/g_1mu), (point_b.z/g_1mu)))
                                koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((point_a.x/g_1mu + g_go.x*-1/g_1mu),(point_a.y/g_1mu + g_go.y*-1/g_1mu),(point_a.z/g_1mu + g_go.z*-1/g_1mu)))
                                koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((point_b.x/g_1mu + g_go.x*-1/g_1mu),(point_b.y/g_1mu + g_go.y*-1/g_1mu),(point_b.z/g_1mu + g_go.z*-1/g_1mu)))

    #print(koordinaten_liste_[:])
    if len(koordinaten_liste_) > 0:
        file_save(koordinaten_liste_)



def export_2_ElementsbyType_3_4():
    flag_intervall = False
    global punkte_set
    global g_1mu
    global g_go    
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu    
    koordinaten_liste_=[]
    ebenenId_selected_ = []

    if len(tree_2.selection()) == 0: return

    interval_2 = get_intervall()
    #print(interval_2)

    if weitere_punkte_var.get():
        flag_intervall = True
 
    for selected_item in tree_2.selection():
        item = tree_2.item(selected_item)
        record = item['values']
        id_=record[1]
        ebenenId_selected_.append(id_)
    
    #print(ebenenId_selected)

    start = DPoint3d()
    end = DPoint3d()

    #Get active model
    ACTIVEMODEL = ISessionMgr.ActiveDgnModelRef
    dgnModel = ACTIVEMODEL.GetDgnModel()
    #Get all graphical elements from the model
    graphicalElements = dgnModel.GetGraphicElements()

    for perElementRef in graphicalElements:
        elementId = perElementRef.GetElementId()
        eeh = EditElementHandle(perElementRef, dgnModel)
        eh = ElementHandle(perElementRef)

        msElement = MSElement()
        msElement = eeh.GetElement ()

        isGraphics = msElement.ehdr.isGraphics
        isInvisible = msElement.hdr.dhdr.props.b.invisible

        if (isGraphics and not(isInvisible)):
            eleType = eh.GetElementType()
            levelId = msElement.ehdr.level
            if levelId in ebenenId_selected_:
                if (eleType == 3) or (eleType == 4) or (eleType == 6) or (eleType == 16) or (eleType == 27) or (eleType == 15) or (eleType == 12) or (eleType ==14):
                    curve = ICurvePathQuery.ElementToCurveVector(eh)
                    curve.GetStartEnd(start, end)

                    interval_0 = interval_2 * g_1mu
                    lineLength = curve.Length()
                    
                    anzPunkte = int(lineLength/interval_0) + 1
                    
                    #if eleType==12 : print(anzPunkte)
               
                    laengen = DoubleArray()
                    if anzPunkte > 0:
                        interval_1 = lineLength/anzPunkte
                        for i in range(anzPunkte + 1):
                            laengen.append(i * interval_1)
                        #print(laengen)
                    
                    neue_punkte = CurveLocationDetailArray()
                    ele_ = CurveLocationDetail()
                    
                    #if lineLength < interval_0: print('Laenge = ', lineLength/g_1mu, ' Intervall = ', interval_0)
                    if eleType == 3 and export_lines_var.get() == True:
                        #print('line')
                        ergebnis = curve.AddSpacedPoints(laengen, neue_punkte)
                        #print('ergebnis = ', ergebnis)
                        for i, ele_ in enumerate(neue_punkte):
                            #print(ele_)
                            #teilung_ = ele_.fraction
                            punkt_ = ele_.point
                            #print(i, teilung_, punkt_)
                            #print(" {0:.4f} , {1:.4f} , {2:.4f} ".format((punkt_.x/g_1mu),(punkt_.y/g_1mu),(punkt_.z/g_1mu)))
                            #if flag_intervall == True:
                                #koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((punkt_.x/g_1mu),(punkt_.y/g_1mu),(punkt_.z/g_1mu)))
                        #print(g_go, g_go.x*-1/g_1mu, g_go.y*-1/g_1mu, g_go.z*-1/g_1mu) 
                        # g_go_x_mu, g_go_y_mu,g_go_z_mu    
                        koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((start.x/g_1mu + g_go_x_mu),(start.y/g_1mu + g_go_y_mu),(start.z/g_1mu + g_go_z_mu)))
                        koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((end.x/g_1mu + g_go_x_mu), (end.y/g_1mu + g_go_y_mu), (end.z/g_1mu + g_go_z_mu)))
                        if flag_intervall == True:
                            zwischen_punkte = zwischenpunkte((start.x/g_1mu + g_go_x_mu, start.y/g_1mu + g_go_y_mu, start.z/g_1mu + g_go_z_mu), (end.x/g_1mu + g_go_x_mu, end.y/g_1mu + g_go_y_mu, end.z/g_1mu + g_go_z_mu), interval_2)
                            for punkt in zwischen_punkte:
                                punkte_set.add(punkt)                       
                    if eleType == 4 and export_linestrings_var.get() == True:
                        #print('linestring')
                        for element in curve:
                            points = element.GetLineString()
                            
                            for i in range(len(points)-1):
                                point_a = points[i]
                                point_b = points[i+1]
                                #print(" {0:.4f} , {1:.4f} , {2:.4f} ".format((point_a.x/g_1mu),(point_a.y/g_1mu),(point_a.z/g_1mu)))
                                #print(" {0:.4f} , {1:.4f} , {2:.4f} ".format((point_b.x/g_1mu), (point_b.y/g_1mu), (point_b.z/g_1mu)))
                                koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((point_a.x/g_1mu + g_go_x_mu),(point_a.y/g_1mu + g_go_y_mu),(point_a.z/g_1mu + g_go_z_mu)))
                                koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((point_b.x/g_1mu + g_go_x_mu),(point_b.y/g_1mu + g_go_y_mu),(point_b.z/g_1mu + g_go_z_mu)))
                                if flag_intervall == True:
                                    zwischen_punkte = zwischenpunkte((point_a.x/g_1mu + g_go_x_mu, point_a.y/g_1mu + g_go_y_mu, point_a.z/g_1mu + g_go_z_mu), (point_b.x/g_1mu + g_go_x_mu, point_b.y/g_1mu + g_go_y_mu, point_b.z/g_1mu + g_go_z_mu), interval_2)
                                    for punkt in zwischen_punkte:
                                        punkte_set.add(punkt)
                                        
                    if eleType == 6 and export_shapes_var.get() == True:
                        #print('shape')
                        for element in curve:
                            points = element.GetLineString()
                            
                            for i in range(len(points)-1):
                                point_a = points[i]
                                point_b = points[i+1]
                                #print(" {0:.4f} , {1:.4f} , {2:.4f} ".format((point_a.x/g_1mu),(point_a.y/g_1mu),(point_a.z/g_1mu)))
                                #print(" {0:.4f} , {1:.4f} , {2:.4f} ".format((point_b.x/g_1mu), (point_b.y/g_1mu), (point_b.z/g_1mu)))
                                koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((point_a.x/g_1mu + g_go_x_mu),(point_a.y/g_1mu + g_go_y_mu),(point_a.z/g_1mu + g_go_z_mu)))
                                koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((point_b.x/g_1mu + g_go_x_mu),(point_b.y/g_1mu + g_go_y_mu),(point_b.z/g_1mu + g_go_z_mu)))
                                if flag_intervall == True:
                                    zwischen_punkte = zwischenpunkte((point_a.x/g_1mu + g_go_x_mu, point_a.y/g_1mu + g_go_y_mu, point_a.z/g_1mu + g_go_z_mu), (point_b.x/g_1mu + g_go_x_mu, point_b.y/g_1mu + g_go_y_mu, point_b.z/g_1mu + g_go_z_mu), interval_2)
                                    for punkt in zwischen_punkte:
                                        punkte_set.add(punkt) 
                                        
                    if eleType == 16 and export_arcs_var.get() == True:
                        #print('arc')
                        ergebnis = curve.AddSpacedPoints(laengen, neue_punkte)
                        #print('ergebnis = ', ergebnis)
                        for i, ele_ in enumerate(neue_punkte):
                            #print(ele_)
                            #teilung_ = ele_.fraction
                            punkt_ = ele_.point
                            #print(i, teilung_, punkt_)
                            #print(" {0:.4f} , {1:.4f} , {2:.4f} ".format((punkt_.x/g_1mu),(punkt_.y/g_1mu),(punkt_.z/g_1mu)))
                            if flag_intervall == True:
                                koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((punkt_.x/g_1mu + g_go_x_mu),(punkt_.y/g_1mu + g_go_y_mu),(punkt_.z/g_1mu + g_go_z_mu)))
                            
                        koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((start.x/g_1mu + g_go_x_mu),(start.y/g_1mu + g_go_y_mu),(start.z/g_1mu + g_go_z_mu)))
                        koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((end.x/g_1mu + g_go_x_mu), (end.y/g_1mu + g_go_y_mu), (end.z/g_1mu + g_go_z_mu)))
                        
                    if eleType == 27 and export_bsplines_var.get() == True:
                        #print('bspline')
                        #print(" {0:.4f} , {1:.4f} , {2:.4f} ".format((start.x/g_1mu),(start.y/g_1mu),(start.z/g_1mu)))
                        #print(" {0:.4f} , {1:.4f} , {2:.4f} ".format((end.x/g_1mu), (end.y/g_1mu), (end.z/g_1mu)))
                        ergebnis = curve.AddSpacedPoints(laengen, neue_punkte)
                        #print('ergebnis = ', ergebnis)
                        for i, ele_ in enumerate(neue_punkte):
                            #print(ele_)
                            #teilung_ = ele_.fraction
                            punkt_ = ele_.point
                            #print(i, teilung_, punkt_)
                            #print(" {0:.4f} , {1:.4f} , {2:.4f} ".format((punkt_.x/g_1mu),(punkt_.y/g_1mu),(punkt_.z/g_1mu)))
                            if flag_intervall == True:
                                koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((punkt_.x/g_1mu + g_go_x_mu),(punkt_.y/g_1mu + g_go_y_mu),(punkt_.z/g_1mu + g_go_z_mu)))
                            
                        koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((start.x/g_1mu + g_go_x_mu),(start.y/g_1mu + g_go_y_mu),(start.z/g_1mu + g_go_z_mu)))
                        koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((end.x/g_1mu + g_go_x_mu), (end.y/g_1mu + g_go_y_mu), (end.z/g_1mu + g_go_z_mu)))
                        
                    if (eleType == 12 or eleType == 14) and export_complex_var.get() == True: 
                        #print('complex chain')
                        ergebnis = curve.AddSpacedPoints(laengen, neue_punkte)
                        #print('ergebnis = ', ergebnis)
                        for i, ele_ in enumerate(neue_punkte):
                            punkt_ = ele_.point
                            if flag_intervall == True:
                                koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((punkt_.x/g_1mu + g_go_x_mu),(punkt_.y/g_1mu + g_go_y_mu),(punkt_.z/g_1mu + g_go_z_mu)))
                                
                        koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((start.x/g_1mu + g_go_x_mu),(start.y/g_1mu + g_go_y_mu),(start.z/g_1mu + g_go_z_mu)))
                        koordinaten_liste_.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format((end.x/g_1mu + g_go_x_mu), (end.y/g_1mu + g_go_y_mu), (end.z/g_1mu + g_go_z_mu)))
                            
                        
                        
                        

    #print(koordinaten_liste_[:])
    for punkt_ in koordinaten_liste_:
        punkte_set.add(punkt_)
    
    #für den export vorbereiten
    punkte_liste_ = []

    for item in punkte_set:
        punkte_liste_.append(item)    

    
    if len(punkte_liste_) > 0:
        label_koor_datei.config(text=f"im Koordinatenspeicher sind : {len(punkte_set)} Punktkoordinaten")
        #print(len(punkte_liste_))
        #print(punkte_liste_[:10])

def exportKoordinaten():
    global punkte_set
    anzahl_p = len(punkte_set)
    if anzahl_p == 0: return
    nr = anzahl_p
    k=0
    while nr > 1:
        k = k + 1
        nr = nr / 10
    nr_1 = int(1 * 10**k + 1)
    
    punkte_liste_ = []
    for item in punkte_set:
        punkte_liste_.append(item)

    file_save_2(punkte_liste_, nr_1)

def importKoordinaten():
    global punkte_set
    koordinaten_liste_2 = []
    # koordinatenpunkte einlesen
    Dateien = open_files_dialog_from_ex()

    if len(Dateien) == 0:
        return  

    liste_p_ = dateien_einlesen(Dateien)

    punkte_xyz_ex, punkte_xy_ex, x_ex, y_ex, z_ex, mark_ex = daten_aufbereiten(liste_p_,[])

    for punkt in punkte_xyz_ex:
        #print(punkt[0],punkt[1],punkt[2])
        #print(" {0:.4f} , {1:.4f} , {2:.4f} ".format(punkt[0],punkt[1],punkt[2]))
        koordinaten_liste_2.append(" {0:.4f} , {1:.4f} , {2:.4f} ".format(punkt[0],punkt[1],punkt[2]))
        
    for punkt_ in koordinaten_liste_2[:]:
        #print(punkt_)
        punkte_set.add(punkt_)

    label_koor_datei.config(text=f"im Koordinatenspeicher sind : {len(punkte_set)} Punktkoordinaten")


# zur Volumenberechnung
def berechnung_loeschen():
    # volumenberechnung loeschen        
    file_text_vol.delete("1.0", tk.END)
    root.update_idletasks() 
    
def berechnung_kopf():
    # kopfzeilen eintragen
    now = dt.datetime.now()
    datum = now.strftime('%d. %B %Y')
    zeit = now.strftime('%H:%M:%S')
    file_text_vol.insert(tk.END, 'Volumenberechnung eines DGM' + '\n' + '\n')
    file_text_vol.insert(tk.END, 'Diese Berechnung erfolgte am ' + str(datum) + '  um  ' + str(zeit) + '\n')
    file_text_vol.insert(tk.END, 'durch :' + '\n' + '\n')
    dgnFile = ISessionMgr.GetActiveDgnFile()
    #beFileName = BeFileName(str(dgnFile.GetFileName()))
    #print(str(beFileName))
    dateiname = str(dgnFile.GetFileName())
    file_text_vol.insert(tk.END, 'Zeichnungsdatei : ' + dateiname + '\n')
    
    if len(tree_3.selection()) != 0:
        for selected_item in tree_3.selection():
            item = tree_3.item(selected_item)
            record = item['values']
            name_ = record[0]
            file_text_vol.insert(tk.END, 'DGM auf Ebene   : ' + str(name_) + '\n' + '\n')

    
    root.update_idletasks()
  
def berechnung_speichern():
    #volumenberechnung protokoll speichern
    berechnung = file_text_vol.get('1.0', 'end-1c')
    file_save_3(berechnung)
      

def select_3_ElementsbyType():
    #dgn file scannen nach dgm maschen
    global g_1mu
    global g_go
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu
    global maschen_eines_levels
    maschen_eines_levels=[]
    levelAnzahlMaschen=[]
    levelGrundFlaeche=[]
    levelVolumen=[]
    hoeheMittel=[]
    hoeheMin=[]
    hoeheMax=[]
    x_mittel=[]
    y_mittel=[]
    z_mittel=[]
    rechtsMin=[]
    rechtsMax=[]
    hochMin=[]
    hochMax=[]
    for i in range(max(levelListIDs)+1):
        maschen_eines_levels.append([])
        levelAnzahlMaschen.append(0)
        levelGrundFlaeche.append(0)
        levelVolumen.append(0)
        hoeheMittel.append(0)
        hoeheMin.append(0)
        hoeheMax.append(0)
        x_mittel.append(0)
        y_mittel.append(0)
        z_mittel.append(0)
        rechtsMin.append(0)
        rechtsMax.append(0)
        hochMin.append(0)
        hochMax.append(0)

    schritt = 2
    progressbar_1.step(schritt)
    root.update_idletasks()
    anzahl = 0    
    #Get active model
    ACTIVEMODEL = ISessionMgr.ActiveDgnModelRef
    dgnModel = ACTIVEMODEL.GetDgnModel()
    #name =  model.GetModelName()
    #Get all graphical elements from the model
    graphicalElements = dgnModel.GetGraphicElements()
    

    for perElementRef in graphicalElements:
        elementId = perElementRef.GetElementId()
        eeh = EditElementHandle(perElementRef, dgnModel)
        eh = ElementHandle(perElementRef)

        msElement = MSElement()
        msElement = eeh.GetElement ()

        isGraphics = msElement.ehdr.isGraphics
        isInvisible = msElement.hdr.dhdr.props.b.invisible

        if (isGraphics and not(isInvisible)):
            eleType = eh.GetElementType()
            levelId = msElement.ehdr.level
            if (eleType == 6) :
                curve = ICurvePathQuery.ElementToCurveVector(eh)
                for element in curve:
                    points = element.GetLineString()
                    if len(points) == 4:
                        if levelId in levelListIDs:
                            levelAnzahlMaschen[levelId] = levelAnzahlMaschen[levelId] + 1
                            masche = []
                            punkt_a = []
                            punkt_b = []
                            punkt_c = []
                            point_a = points[0]
                            point_b = points[1]
                            point_c = points[2]
                            #print(g_go, g_go.x*-1/g_1mu, g_go.y*-1/g_1mu, g_go.z*-1/g_1mu)
                            #global origin muss beruecksichtigt werden
                            #g_go_x_mu, g_go_y_mu, g_go_z_mu
                            #punkt_a = [point_a.x/g_1mu + g_go.x*-1/g_1mu, point_a.y/g_1mu + g_go.y*-1/g_1mu, point_a.z/g_1mu + g_go.z*-1/g_1mu]
                            #punkt_b = [point_b.x/g_1mu + g_go.x*-1/g_1mu, point_b.y/g_1mu + g_go.y*-1/g_1mu, point_b.z/g_1mu + g_go.z*-1/g_1mu]
                            #punkt_c = [point_c.x/g_1mu + g_go.x*-1/g_1mu, point_c.y/g_1mu + g_go.y*-1/g_1mu, point_c.z/g_1mu + g_go.z*-1/g_1mu]

                            punkt_a = [point_a.x/g_1mu + g_go_x_mu, point_a.y/g_1mu + g_go_y_mu, point_a.z/g_1mu + g_go_z_mu]
                            punkt_b = [point_b.x/g_1mu + g_go_x_mu, point_b.y/g_1mu + g_go_y_mu, point_b.z/g_1mu + g_go_z_mu]
                            punkt_c = [point_c.x/g_1mu + g_go_x_mu, point_c.y/g_1mu + g_go_y_mu, point_c.z/g_1mu + g_go_z_mu]
                                                        
                            h_min = punkt_a[2]
                            h_max = punkt_a[2]
                            if punkt_b[2] < h_min: h_min = punkt_b[2]
                            if punkt_c[2] < h_min: h_min = punkt_c[2]
                            if punkt_b[2] > h_max: h_max = punkt_b[2]
                            if punkt_c[2] > h_max: h_max = punkt_c[2]
                            
                            x_min = punkt_a[0]
                            x_max = punkt_a[0]
                            if punkt_b[0] < x_min: x_min = punkt_b[0]
                            if punkt_c[0] < x_min: x_min = punkt_c[0]
                            if punkt_b[0] > x_max: x_max = punkt_b[0]
                            if punkt_c[0] > x_max: x_max = punkt_c[0] 
                            
                            y_min = punkt_a[1]
                            y_max = punkt_a[1]
                            if punkt_b[1] < y_min: y_min = punkt_b[1]
                            if punkt_c[1] < y_min: y_min = punkt_c[1]
                            if punkt_b[1] > y_max: y_max = punkt_b[1]
                            if punkt_c[1] > y_max: y_max = punkt_c[1]                                                       
                            
                            masche.append(punkt_a)
                            masche.append(punkt_b)
                            masche.append(punkt_c)
                            
                            maschen_eines_levels[levelId].append(masche)

                            fl_grund = abs(triangle_det_2d(masche))
                            volumen_ueber_0 = fl_grund * ((masche[0][2] + masche[1][2] + masche[2][2]) / 3.0)
                            
                            levelGrundFlaeche[levelId] = levelGrundFlaeche[levelId] + fl_grund
                            levelVolumen[levelId] = levelVolumen[levelId] + volumen_ueber_0
                            
                            x_mittel[levelId]  = x_mittel[levelId] + (punkt_a[0] + punkt_b[0] + punkt_c[0]) / 3.0
                            y_mittel[levelId]  = y_mittel[levelId] + (punkt_a[1] + punkt_b[1] + punkt_c[1]) / 3.0
                            z_mittel[levelId]  = z_mittel[levelId] + (punkt_a[2] + punkt_b[2] + punkt_c[2]) / 3.0
                            
                            if levelAnzahlMaschen[levelId] == 1:
                                hoeheMin[levelId] = h_min
                                hoeheMax[levelId] = h_max
                                rechtsMin[levelId] = x_min
                                rechtsMax[levelId] = x_max
                                hochMin[levelId] = y_min
                                hochMax[levelId] = y_max
                                
                            if h_min < hoeheMin[levelId]: hoeheMin[levelId] = h_min
                            if h_max > hoeheMax[levelId]: hoeheMax[levelId] = h_max
                            
                            if x_min < rechtsMin[levelId]: rechtsMin[levelId] = x_min
                            if x_max > rechtsMax[levelId]: rechtsMax[levelId] = x_max
                            
                            if y_min < hochMin[levelId]: hochMin[levelId] = y_min
                            if y_max > hochMax[levelId]: hochMax[levelId] = y_max
                            
                            anzahl = anzahl + 1
                            if anzahl % 1000 == 0:
                                progressbar_1.step(schritt)
                                root.update_idletasks()
                            if anzahl % 45000 == 0:
                                schritt = schritt * -1
                        
    # loesche alle items in treeview tree_3
    all_root_items = tree_3.get_children()
    tree_3.delete(*all_root_items)
    root.update_idletasks()
                              
    # daten aktualisieren
    ebenen = []
    for i in range(len(levelListIDs)):
        levelId_ = levelListIDs[i]
        if levelAnzahlMaschen[levelId_] > 0:
            hoeheMittel = levelVolumen[levelId_] / levelGrundFlaeche[levelId_]
            ebenen.append((levelList[i], levelListIDs[i], levelAnzahlMaschen[levelId_],
             '{0:_.3f}'.format(levelGrundFlaeche[levelId_]), '{0:_.3f}'.format(levelVolumen[levelId_]),
              '{0:_.8f}'.format(hoeheMittel), '{0:_.3f}'.format(hoeheMin[levelId_]),
               '{0:_.3f}'.format(hoeheMax[levelId_]), '{0:_.3f}'.format(rechtsMin[levelId_]),
                '{0:_.3f}'.format(rechtsMax[levelId_]), '{0:_.3f}'.format(hochMin[levelId_]),
                 '{0:_.3f}'.format(hochMax[levelId_]), '{0:_.3f}'.format(x_mittel[levelId_]/levelAnzahlMaschen[levelId_]),
                  '{0:_.3f}'.format(y_mittel[levelId_]/levelAnzahlMaschen[levelId_]), '{0:_.3f}'.format(z_mittel[levelId_]/levelAnzahlMaschen[levelId_]) ))

    # add data to the treeview
    for ebene in ebenen:
        tree_3.insert('', tk.END, values=ebene)
    
    # resize columns 
    tree_3.column('anzahl_maschen_1',anchor='center', stretch=False, width=175, minwidth=100)   
    tree_3.column('grundflaeche_1',anchor='center', stretch=False, width=175, minwidth=100)
    tree_3.column('volumen_1',anchor='center', stretch=False, width=175, minwidth=100)
    tree_3.column('z_aus_1',anchor='center', stretch=False, width=175, minwidth=50)
    tree_3.column('z_min_1',anchor='center', stretch=False, width=120, minwidth=50)
    tree_3.column('z_max_1',anchor='center', stretch=False, width=120, minwidth=50)
    tree_3.column('x_min_1',anchor='center', stretch=False, width=120, minwidth=50)
    tree_3.column('x_max_1',anchor='center', stretch=False, width=120, minwidth=50)
    tree_3.column('y_min_1',anchor='center', stretch=False, width=130, minwidth=50)
    tree_3.column('y_max_1',anchor='center', stretch=False, width=130, minwidth=50)
    tree_3.column('x_mittel_1',anchor='center', stretch=False, width=120, minwidth=50)
    tree_3.column('y_mittel_1',anchor='center', stretch=False, width=130, minwidth=50)
    tree_3.column('z_mittel_1',anchor='center', stretch=False, width=120, minwidth=50)
    
    progressbar_1.stop()
    
    button_11.config(state='active')
        
    root.update_idletasks() 



def dgm_volumen_berechnen_2():
    #volumenberechnung fuer ausgewaehlte ebene
    global g_1mu
    global g_go
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu
    global faktor_projektion
    global hoehen_werte
    global maschen_eines_levels
    
    # keine ebene gewaehlt, abbruch
    if len(tree_3.selection()) == 0: return
    
    ebenenId_selected_ = []
    ebeneDatensatz_selected_ = []
    
    maschen = []
    
    for selected_item in tree_3.selection():
        item = tree_3.item(selected_item)
        record = item['values']
        #print(record)
        id_=record[1]
        ebenenId_selected_.append(id_)
        ebeneDatensatz_selected_.append(record)
        ost_mittel_km = float(record[12])/1000.0
        utm_zone = int(ost_mittel_km // 1000)
        koord_meridian_km = utm_zone * 1000 + 500
        if ost_mittel_km < 100:
            ost_mittel_km = 500
        z_mittel_km = float(record[14])/1000.0
    
    faktor = 1.0000
    faktor_projektion = 1.000
    if utm_var.get() == True:
        # utm koordinaten beruecksichtigen
        utm_faktor = 0.9996
        faktor = 1.0 / utm_faktor
        radius_km = float(radius_eingabe.get())
        if radius_km < 5000 or radius_km > 7000:
            radius_km = 6382
        #print(radius_km, ost_mittel_km)
        faktor_projektion = ( (1.0 + ((ost_mittel_km - koord_meridian_km)**2)/(2.0 * radius_km**2)) * utm_faktor )
        undulation_km = float(geoid_un_eingabe.get())/1000.0
        #mittlere Hoehe Gebiet
        h_ellip_km = z_mittel_km + undulation_km
        faktor_hoehe = 1.0 - (h_ellip_km / radius_km)
        faktor = 1.0 / ( (1.0 + ((ost_mittel_km - koord_meridian_km)**2)/(2.0 * radius_km**2) - (h_ellip_km / radius_km)) * utm_faktor )
            
    
    #print(faktor)
    abr_hoehe = float(abr_hoehe_eingabe.get())
    #print(abr_hoehe)
    
    schritt = 2
    progressbar_1.step(schritt)
    root.update_idletasks()
    anzahl = 0
    
    for masche in maschen_eines_levels[id_]:
        punkt_a = masche[0]
        punkt_b = masche[1]
        punkt_c = masche[2]
                                                                              
        # sortiere nach z_wert
            
        if punkt_b[2] < punkt_a[2]:
            punkt_temp = punkt_a
            punkt_a = punkt_b
            punkt_b = punkt_temp
            
        if punkt_c[2] < punkt_a[2]:
            punkt_temp = punkt_a
            punkt_a = punkt_c
            punkt_c = punkt_temp                                

        if punkt_c[2] < punkt_b[2]:
            punkt_temp = punkt_b
            punkt_b = punkt_c
            punkt_c = punkt_temp
                            
            
        #masche.append(punkt_a)
        #masche.append(punkt_b)
        #masche.append(punkt_c)
        
        #maschen.append(masche)

        anzahl = anzahl + 1
        if anzahl % 1000 == 0:
            progressbar_1.step(schritt)
            root.update_idletasks()
        if anzahl % 45000 == 0:
            schritt = schritt * -1

    grundflaeche_aus_ges = 0.0
    oberflaeche_aus_ges = 0.0
    volumen_aus_ges = 0.0
    
    grundflaeche_unter = 0.0
    oberflaeche_unter = 0.0
    volumen_unter = 0.0
    
    grundflaeche_ueber = 0.0
    oberflaeche_ueber = 0.0
    volumen_ueber = 0.0
    
    epsilon = 1e-8
    
    if utm_var.get() == True:
        faktor_hoehe_masche_min = faktor_hoehe
        faktor_hoehe_masche_max = faktor_hoehe
        
    hoehen_korrektur = str(options_utm_hoehe.get())
    
    for masche_ in maschen_eines_levels[id_]:       #maschen:
        
        punkt_1 = masche_[0]
        punkt_2 = masche_[1]
        punkt_3 = masche_[2]
        
        # Hoehenkorrektur Maschen Mittel
        if utm_var.get() == True and (hoehen_korrektur == hoehen_werte[1] or hoehen_korrektur == hoehen_werte[0]) :
            h_masche_ellip_km = (punkt_1[2] + punkt_2[2] + punkt_3[2])/3.0/1000 + undulation_km
            faktor_hoehe_masche = 1.0 - (h_masche_ellip_km / radius_km)
            #
            faktor = 1.0 / ( (1.0 + ((ost_mittel_km - koord_meridian_km)**2)/(2.0 * radius_km**2) - (h_masche_ellip_km / radius_km)) * utm_faktor )
            if faktor_hoehe_masche < faktor_hoehe_masche_min:
                faktor_hoehe_masche_min = faktor_hoehe_masche
            if faktor_hoehe_masche > faktor_hoehe_masche_max:
                faktor_hoehe_masche_max = faktor_hoehe_masche

        seite_a1_ = strecke_2d_f(punkt_1, punkt_2, faktor)
        seite_b1_ = strecke_2d_f(punkt_2, punkt_3, faktor)
        seite_c1_ = strecke_2d_f(punkt_1, punkt_3, faktor)
        
        flaeche_grund_ = flaeche_heron(seite_a1_, seite_b1_, seite_c1_)
        grundflaeche_aus_ges = grundflaeche_aus_ges + flaeche_grund_
 
        seite_a2_ = strecke_3d_f(punkt_1, punkt_2, faktor)
        seite_b2_ = strecke_3d_f(punkt_2, punkt_3, faktor)
        seite_c2_ = strecke_3d_f(punkt_1, punkt_3, faktor)

        flaeche_ober_ = flaeche_heron(seite_a2_, seite_b2_, seite_c2_)
        oberflaeche_aus_ges = oberflaeche_aus_ges + flaeche_ober_
        
        volumen_ = flaeche_grund_ * (punkt_1[2] + punkt_2[2] + punkt_3[2]) / 3.0
        volumen_aus_ges = volumen_aus_ges + volumen_
        volumen_aus_abrechnung = volumen_aus_ges - grundflaeche_aus_ges * abr_hoehe
        
        # volumen ober- oder unterhalb abr_hoehe
        if abr_hoehe <= (punkt_1[2] + epsilon):
            # 3 Punkte ueber Horizont
            #print('Fall 1: 3 Punkte ueber Horizont')
            grundflaeche_unter = grundflaeche_unter + 0.0
            oberflaeche_unter = oberflaeche_unter + 0.0
            volumen_unter = volumen_unter + 0.0
            
            grundflaeche_ueber = grundflaeche_ueber + flaeche_grund_ 
            oberflaeche_ueber = oberflaeche_ueber + flaeche_ober_
            volumen_ueber = volumen_ueber + (volumen_ - flaeche_grund_ * abr_hoehe)

        elif abr_hoehe >= (punkt_3[2] - epsilon):
            # 3 Punkte unter Horizont
            #print('Fall 2: 3 Punkte unter Horizont')
            grundflaeche_unter = grundflaeche_unter + flaeche_grund_ 
            oberflaeche_unter = oberflaeche_unter + flaeche_ober_
            volumen_unter = volumen_unter + (volumen_ - flaeche_grund_ * abr_hoehe)
            
            grundflaeche_ueber = grundflaeche_ueber + 0.0
            oberflaeche_ueber = oberflaeche_ueber + 0.0
            volumen_ueber = volumen_ueber + 0.0
            

        elif abs(abr_hoehe - punkt_2[2]) < epsilon:
            # Sonderfall Punkt 3 ueber Horizont - aber Horizont geht durch Punkt 2 !!
            #print('Fall 5: Punkt3  ueber Horizont , Horizont geht durch Punkt2')
            #file_text_vol.insert(tk.END, 'Fall 5: Punkt3  ueber Horizont , Horizont geht durch Punkt2'  + '\n')
            #root.update_idletasks()            
            hilfspunkt_1 = [0.0, 0.0, 0.0] 
            faktor_z = (abr_hoehe - punkt_1[2]) / (punkt_3[2] - punkt_1[2])

            hilfspunkt_1[0] = punkt_1[0] + ((punkt_3[0] - punkt_1[0]) * faktor_z)
            hilfspunkt_1[1] = punkt_1[1] + ((punkt_3[1] - punkt_1[1]) * faktor_z)
            hilfspunkt_1[2] = abr_hoehe
            
            #teil1 masche unterhalb abrechnungshoehe
            
            seite_a1_ = strecke_2d_f(punkt_1, punkt_2, faktor)
            seite_b1_ = strecke_2d_f(punkt_2, hilfspunkt_1, faktor)
            seite_c1_ = strecke_2d_f(punkt_1, hilfspunkt_1, faktor)
        
            flaeche_grund_ = flaeche_heron(seite_a1_, seite_b1_, seite_c1_)
            grundflaeche_unter = grundflaeche_unter + flaeche_grund_
 
            seite_a2_ = strecke_3d_f(punkt_1, punkt_2, faktor)
            seite_b2_ = strecke_3d_f(punkt_2, hilfspunkt_1, faktor)
            seite_c2_ = strecke_3d_f(punkt_1, hilfspunkt_1, faktor)

            flaeche_ober_ = flaeche_heron(seite_a2_, seite_b2_, seite_c2_)
            oberflaeche_unter = oberflaeche_unter + flaeche_ober_
        
            volumen_ = flaeche_grund_ * (punkt_1[2] + punkt_2[2] + hilfspunkt_1[2]) / 3.0
            volumen_unter = volumen_unter + (volumen_- flaeche_grund_ * abr_hoehe)

            #teil2 masche oberhalb abrechnungshoehe
            
            seite_a1_ = strecke_2d_f(hilfspunkt_1, punkt_2, faktor)
            seite_b1_ = strecke_2d_f(punkt_2, punkt_3, faktor)
            seite_c1_ = strecke_2d_f(hilfspunkt_1, punkt_3, faktor)
        
            flaeche_grund_ = flaeche_heron(seite_a1_, seite_b1_, seite_c1_)
            grundflaeche_ueber = grundflaeche_ueber + flaeche_grund_
 
            seite_a2_ = strecke_3d_f(hilfspunkt_1, punkt_2, faktor)
            seite_b2_ = strecke_3d_f(punkt_2, punkt_3, faktor)
            seite_c2_ = strecke_3d_f(hilfspunkt_1, punkt_3, faktor)

            flaeche_ober_ = flaeche_heron(seite_a2_, seite_b2_, seite_c2_)
            oberflaeche_ueber = oberflaeche_ueber + flaeche_ober_
        
            volumen_ = flaeche_grund_ * (punkt_3[2] + punkt_2[2] + hilfspunkt_1[2]) / 3.0
            volumen_ueber = volumen_ueber + (volumen_- flaeche_grund_ * abr_hoehe)
            
           
        elif abr_hoehe > punkt_2[2]:
            #print('Fall 4: Punkt3  ueber Horizont , Punkt2 und Punkt1 unter Horizont')
            #file_text_vol.insert(tk.END, 'Fall 4: Punkt3  ueber Horizont , Punkt2 und Punkt1 unter Horizont'  + '\n')
            #root.update_idletasks() 
            hilfspunkt_1 = [0.0, 0.0, 0.0] 
            faktor_z = (abr_hoehe - punkt_1[2]) / (punkt_3[2] - punkt_1[2])
            
            hilfspunkt_1[0] = punkt_1[0] + ((punkt_3[0] - punkt_1[0]) * faktor_z)
            hilfspunkt_1[1] = punkt_1[1] + ((punkt_3[1] - punkt_1[1]) * faktor_z)
            hilfspunkt_1[2] = abr_hoehe
                        
            hilfspunkt_2 = [0.0, 0.0, 0.0] 
            faktor_z = (abr_hoehe - punkt_2[2]) / (punkt_3[2] - punkt_2[2])
            
            hilfspunkt_2[0] = punkt_2[0] + ((punkt_3[0] - punkt_2[0]) * faktor_z)
            hilfspunkt_2[1] = punkt_2[1] + ((punkt_3[1] - punkt_2[1]) * faktor_z)
            hilfspunkt_2[2] = abr_hoehe
                        
            #teil1 masche 1 oberhalb abrechnungshoehe 
            seite_a1_ = strecke_2d_f(hilfspunkt_1, hilfspunkt_2, faktor)
            seite_b1_ = strecke_2d_f(punkt_3, hilfspunkt_2, faktor)
            seite_c1_ = strecke_2d_f(punkt_3, hilfspunkt_1, faktor)
                    
            flaeche_grund_ = flaeche_heron(seite_a1_, seite_b1_, seite_c1_)
            grundflaeche_ueber = grundflaeche_ueber + flaeche_grund_
             
            seite_a2_ = strecke_3d_f(hilfspunkt_1, hilfspunkt_2, faktor)
            seite_b2_ = strecke_3d_f(punkt_3, hilfspunkt_2, faktor)
            seite_c2_ = strecke_3d_f(punkt_3, hilfspunkt_1, faktor)
            
            flaeche_ober_ = flaeche_heron(seite_a2_, seite_b2_, seite_c2_)
            oberflaeche_ueber = oberflaeche_ueber + flaeche_ober_
                    
            volumen_ = flaeche_grund_ * (punkt_3[2] + hilfspunkt_2[2] + hilfspunkt_1[2]) / 3.0
            volumen_ueber = volumen_ueber + (volumen_- flaeche_grund_ * abr_hoehe)
                                                           
            #teil2 masche 2 unterhalb abrechnungshoehe 
            seite_a1_ = strecke_2d_f(punkt_1, punkt_2, faktor)
            seite_b1_ = strecke_2d_f(punkt_2, hilfspunkt_2, faktor)
            seite_c1_ = strecke_2d_f(punkt_1, hilfspunkt_2, faktor)
                    
            flaeche_grund_ = flaeche_heron(seite_a1_, seite_b1_, seite_c1_)
            grundflaeche_unter = grundflaeche_unter + flaeche_grund_
             
            seite_a2_ = strecke_3d_f(punkt_1, punkt_2, faktor)
            seite_b2_ = strecke_3d_f(punkt_2, hilfspunkt_2, faktor)
            seite_c2_ = strecke_3d_f(punkt_1, hilfspunkt_2, faktor)
            
            flaeche_ober_ = flaeche_heron(seite_a2_, seite_b2_, seite_c2_)
            oberflaeche_unter = oberflaeche_unter + flaeche_ober_
                    
            volumen_ = flaeche_grund_ * (punkt_1[2] + punkt_2[2] + hilfspunkt_2[2]) / 3.0
            volumen_unter = volumen_unter + (volumen_- flaeche_grund_ * abr_hoehe)
                       
            #teil3 masche 3 unterhalb abrechnungshoehe
            seite_a1_ = strecke_2d_f(punkt_1, hilfspunkt_1, faktor)
            seite_b1_ = strecke_2d_f(hilfspunkt_1, hilfspunkt_2, faktor)
            seite_c1_ = strecke_2d_f(punkt_1, hilfspunkt_2, faktor)

            flaeche_grund_ = flaeche_heron(seite_a1_, seite_b1_, seite_c1_)
            grundflaeche_unter = grundflaeche_unter + flaeche_grund_
 
            seite_a2_ = strecke_3d_f(punkt_1, hilfspunkt_1, faktor)
            seite_b2_ = strecke_3d_f(hilfspunkt_1, hilfspunkt_2, faktor)
            seite_c2_ = strecke_3d_f(punkt_1, hilfspunkt_2, faktor)
            
            flaeche_ober_ = flaeche_heron(seite_a2_, seite_b2_, seite_c2_)
            oberflaeche_unter = oberflaeche_unter + flaeche_ober_
            
            volumen_ = flaeche_grund_ * (punkt_1[2] + hilfspunkt_1[2] + hilfspunkt_2[2]) / 3.0
            volumen_unter = volumen_unter + (volumen_- flaeche_grund_ * abr_hoehe)
            
        elif abr_hoehe < punkt_2[2]:
            #print('Fall 3: Punkt2 und Punkt3  ueber Horizont ,  Punkt1 unter Horizont')
            #file_text_vol.insert(tk.END, 'Fall 3: Punkt2 und Punkt3  ueber Horizont ,  Punkt1 unter Horizont'  + '\n')
            #root.update_idletasks()
            hilfspunkt_1 = [0.0, 0.0, 0.0] 
            faktor_z = (abr_hoehe - punkt_1[2]) / (punkt_3[2] - punkt_1[2])
            
            hilfspunkt_1[0] = punkt_1[0] + ((punkt_3[0] - punkt_1[0]) * faktor_z)
            hilfspunkt_1[1] = punkt_1[1] + ((punkt_3[1] - punkt_1[1]) * faktor_z)
            hilfspunkt_1[2] = abr_hoehe
                        
            hilfspunkt_2 = [0.0, 0.0, 0.0] 
            faktor_z = (abr_hoehe - punkt_1[2]) / (punkt_2[2] - punkt_1[2])
            
            hilfspunkt_2[0] = punkt_1[0] + ((punkt_2[0] - punkt_1[0]) * faktor_z)
            hilfspunkt_2[1] = punkt_1[1] + ((punkt_2[1] - punkt_1[1]) * faktor_z)
            hilfspunkt_2[2] = abr_hoehe
                        
            #teil1 masche 1 unterhalb abrechnungshoehe 
            seite_a1_ = strecke_2d_f(hilfspunkt_1, hilfspunkt_2, faktor)
            seite_b1_ = strecke_2d_f(punkt_1, hilfspunkt_2, faktor)
            seite_c1_ = strecke_2d_f(punkt_1, hilfspunkt_1, faktor)
                    
            flaeche_grund_ = flaeche_heron(seite_a1_, seite_b1_, seite_c1_)
            grundflaeche_unter = grundflaeche_unter + flaeche_grund_
             
            seite_a2_ = strecke_3d_f(hilfspunkt_1, hilfspunkt_2, faktor)
            seite_b2_ = strecke_3d_f(punkt_1, hilfspunkt_2, faktor)
            seite_c2_ = strecke_3d_f(punkt_1, hilfspunkt_1, faktor)
            
            flaeche_ober_ = flaeche_heron(seite_a2_, seite_b2_, seite_c2_)
            oberflaeche_unter = oberflaeche_unter + flaeche_ober_
                    
            volumen_ = flaeche_grund_ * (punkt_1[2] + hilfspunkt_2[2] + hilfspunkt_1[2]) / 3.0
            volumen_unter = volumen_unter + (volumen_- flaeche_grund_ * abr_hoehe)
                                                           
            #teil2 masche 2 oberhalb abrechnungshoehe 
            seite_a1_ = strecke_2d_f(punkt_2, punkt_3, faktor)
            seite_b1_ = strecke_2d_f(punkt_2, hilfspunkt_1, faktor)
            seite_c1_ = strecke_2d_f(punkt_3, hilfspunkt_1, faktor)
                    
            flaeche_grund_ = flaeche_heron(seite_a1_, seite_b1_, seite_c1_)
            grundflaeche_ueber = grundflaeche_ueber + flaeche_grund_
             
            seite_a2_ = strecke_3d_f(punkt_2, punkt_3, faktor)
            seite_b2_ = strecke_3d_f(punkt_2, hilfspunkt_1, faktor)
            seite_c2_ = strecke_3d_f(punkt_3, hilfspunkt_1, faktor)
            
            flaeche_ober_ = flaeche_heron(seite_a2_, seite_b2_, seite_c2_)
            oberflaeche_ueber = oberflaeche_ueber + flaeche_ober_
                    
            volumen_ = flaeche_grund_ * (punkt_2[2] + punkt_3[2] + hilfspunkt_1[2]) / 3.0
            volumen_ueber = volumen_ueber + (volumen_- flaeche_grund_ * abr_hoehe)
                       
            #teil3 masche 3 oberhalb abrechnungshoehe
            seite_a1_ = strecke_2d_f(punkt_2, hilfspunkt_1, faktor)
            seite_b1_ = strecke_2d_f(hilfspunkt_1, hilfspunkt_2, faktor)
            seite_c1_ = strecke_2d_f(punkt_2, hilfspunkt_2, faktor)

            flaeche_grund_ = flaeche_heron(seite_a1_, seite_b1_, seite_c1_)
            grundflaeche_ueber = grundflaeche_ueber + flaeche_grund_
 
            seite_a2_ = strecke_3d_f(punkt_2, hilfspunkt_1, faktor)
            seite_b2_ = strecke_3d_f(hilfspunkt_1, hilfspunkt_2, faktor)
            seite_c2_ = strecke_3d_f(punkt_2, hilfspunkt_2, faktor)
            
            flaeche_ober_ = flaeche_heron(seite_a2_, seite_b2_, seite_c2_)
            oberflaeche_ueber = oberflaeche_ueber + flaeche_ober_
            
            volumen_ = flaeche_grund_ * (punkt_2[2] + hilfspunkt_1[2] + hilfspunkt_2[2]) / 3.0
            volumen_ueber = volumen_ueber + (volumen_- flaeche_grund_ * abr_hoehe)            

        anzahl = anzahl + 1
        if anzahl % 1000 == 0:
            progressbar_1.step(schritt)
            root.update_idletasks()
        if anzahl % 45000 == 0:
            schritt = schritt * -1        
        
    
    file_text_vol.insert(tk.END, 'Anzahl Maschen                    : ' + str(len(maschen_eines_levels[id_])) + '\n')
    if utm_var.get() == True:
        file_text_vol.insert(tk.END, 'utm Korrektur wird beruecksichtigt'  + '\n')
        file_text_vol.insert(tk.END, '       m_0 ist                    : {0:_.4f}'.format(utm_faktor) + '\n')
        file_text_vol.insert(tk.END, 'mittlerer Kruemmungsradius in km  : {0:_.1f}'.format(radius_km) + '\n')
        file_text_vol.insert(tk.END, '       x_mittel   in km           : {0:_.3f}'.format(float(record[12])/1000) + '\n')
        file_text_vol.insert(tk.END, '       utm zone                   : ' + str(utm_zone) + '\n')
        file_text_vol.insert(tk.END, 'mittlerer Ostwert in km           : {0:_.1f}'.format(ost_mittel_km) + '\n')
        file_text_vol.insert(tk.END, '       m_proj ist                 : {0:_.6f}'.format(faktor_projektion) + '\n')
        file_text_vol.insert(tk.END, '       z_mittel   in km           : {0:_.4f}'.format(z_mittel_km) + '\n')
        file_text_vol.insert(tk.END, '       undulation in km           : {0:_.4f}'.format(undulation_km) + '\n')
        file_text_vol.insert(tk.END, '       m_hoehe gebiet ist         : {0:_.6f}'.format(faktor_hoehe) + '\n')
        if hoehen_korrektur == hoehen_werte[0]:
            file_text_vol.insert(tk.END, 'die Umkehrung Hoehenreduktion erfolgt fuer ' + str(hoehen_werte[0])+ '\n')
            file_text_vol.insert(tk.END, '       m_hoehe_maschen_min ist    : {0:_.6f}'.format(faktor_hoehe_masche_min) + '\n')
            file_text_vol.insert(tk.END, '       m_hoehe_maschen_max ist    : {0:_.6f}'.format(faktor_hoehe_masche_max) + '\n')
            massstab_min = 1.0/(0.0 + faktor_projektion - (1.0 - faktor_hoehe_masche_min) * utm_faktor)
            massstab_max = 1.0/(0.0 + faktor_projektion - (1.0 - faktor_hoehe_masche_max) * utm_faktor)
            file_text_vol.insert(tk.END, ' Massstabsfaktor (Masche) max ist : {0:_.6f}'.format(massstab_min)  + '\n')
            file_text_vol.insert(tk.END, ' Massstabsfaktor (Masche) min ist : {0:_.6f}'.format(massstab_max)  + '\n') 
        elif hoehen_korrektur == hoehen_werte[1]:
            file_text_vol.insert(tk.END, 'die Umkehrung Hoehenreduktion erfolgt fuer ' + str(hoehen_werte[1])+ '\n')
            file_text_vol.insert(tk.END, '       m_hoehe_maschen_min ist    : {0:_.6f}'.format(faktor_hoehe_masche_min) + '\n')
            file_text_vol.insert(tk.END, '       m_hoehe_maschen_max ist    : {0:_.6f}'.format(faktor_hoehe_masche_max) + '\n')
            massstab_min = 1.0/(0.0 + faktor_projektion - (1.0 - faktor_hoehe_masche_min) * utm_faktor)
            massstab_max = 1.0/(0.0 + faktor_projektion - (1.0 - faktor_hoehe_masche_max) * utm_faktor)
            file_text_vol.insert(tk.END, ' Massstabsfaktor (Masche) max ist : {0:_.6f}'.format(massstab_min)  + '\n')
            file_text_vol.insert(tk.END, ' Massstabsfaktor (Masche) min ist : {0:_.6f}'.format(massstab_max)  + '\n')            
        else:
            file_text_vol.insert(tk.END, 'die Umkehrung Hoehenreduktion erfolgt fuer ' + str(hoehen_werte[2])+ '\n')
            file_text_vol.insert(tk.END, ' Massstabsfaktor (Gebiet) ist     : {0:_.6f}'.format(faktor)  + '\n')
    if utm_var.get() != True:
        file_text_vol.insert(tk.END, 'Massstabsfaktor ist               : {0:_.6f}'.format(faktor)  + '\n')
    file_text_vol.insert(tk.END, 'Abrechnungshoehe ist              : {0:_.8f}'.format(abr_hoehe)  + '\n')
    file_text_vol.insert(tk.END, 'Grundflaeche des DGM              : {0:_.3f}'.format(grundflaeche_aus_ges)  + '\n')
    file_text_vol.insert(tk.END, 'Oberflaeche  des DGM              : {0:_.3f}'.format(oberflaeche_aus_ges)  + '\n')
    file_text_vol.insert(tk.END, 'Volumen ausgegli.Abrechnungshoehe : {0:_.3f}'.format(volumen_aus_ges - grundflaeche_aus_ges * abr_hoehe)  + '\n')
    file_text_vol.insert(tk.END, 'mittlere Hoehe                    : {0:_.6f}'.format(volumen_aus_ges/grundflaeche_aus_ges)  + '\n' + '\n')
    file_text_vol.insert(tk.END, 'Grundflaeche unter Abr.Hoehe      : {0:_.3f}'.format(grundflaeche_unter)  + '\n')
    file_text_vol.insert(tk.END, 'Oberflaeche  unter Abr.Hoehe      : {0:_.3f}'.format(oberflaeche_unter)  + '\n')
    file_text_vol.insert(tk.END, 'Volumen      unter Abr.Hoehe      : {0:_.3f}'.format(volumen_unter)  + '\n' + '\n')
    file_text_vol.insert(tk.END, 'Grundflaeche ueber Abr.Hoehe      : {0:_.3f}'.format(grundflaeche_ueber)  + '\n')
    file_text_vol.insert(tk.END, 'Oberflaeche  ueber Abr.Hoehe      : {0:_.3f}'.format(oberflaeche_ueber)  + '\n')
    file_text_vol.insert(tk.END, 'Volumen      ueber Abr.Hoehe      : {0:_.3f}'.format(volumen_ueber)  + '\n')
        
    file_text_vol.insert(tk.END, '\n')
    root.update_idletasks()    
    progressbar_1.stop()
    
def select_4_ElementsbyType():
    #dgn file scannen nach dgm maschen
    global g_1mu
    global g_go
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu
    global maschen_eines_levels
    maschen_eines_levels=[]
    levelAnzahlMaschen=[]
    #levelGrundFlaeche=[]
    #levelVolumen=[]
    #hoeheMittel=[]
    hoeheMin=[]
    hoeheMax=[]
    #x_mittel=[]
    #y_mittel=[]
    #z_mittel=[]
    #rechtsMin=[]
    #rechtsMax=[]
    #hochMin=[]
    #hochMax=[]
    for i in range(max(levelListIDs)+1):
        maschen_eines_levels.append([])
        levelAnzahlMaschen.append(0)
        #levelGrundFlaeche.append(0)
        #levelVolumen.append(0)
        #hoeheMittel.append(0)
        hoeheMin.append(0)
        hoeheMax.append(0)
        #x_mittel.append(0)
        #y_mittel.append(0)
        #z_mittel.append(0)
        #rechtsMin.append(0)
        #rechtsMax.append(0)
        #hochMin.append(0)
        #hochMax.append(0)

    schritt = 2
    progressbar_4.step(schritt)
    root.update_idletasks()
    anzahl = 0    
    #Get active model
    ACTIVEMODEL = ISessionMgr.ActiveDgnModelRef
    dgnModel = ACTIVEMODEL.GetDgnModel()
    #name =  model.GetModelName()
    #Get all graphical elements from the model
    graphicalElements = dgnModel.GetGraphicElements()
    

    for perElementRef in graphicalElements:
        elementId = perElementRef.GetElementId()
        eeh = EditElementHandle(perElementRef, dgnModel)
        eh = ElementHandle(perElementRef)

        msElement = MSElement()
        msElement = eeh.GetElement ()

        isGraphics = msElement.ehdr.isGraphics
        isInvisible = msElement.hdr.dhdr.props.b.invisible

        if (isGraphics and not(isInvisible)):
            eleType = eh.GetElementType()
            levelId = msElement.ehdr.level
            if (eleType == 6) :
                curve = ICurvePathQuery.ElementToCurveVector(eh)
                for element in curve:
                    points = element.GetLineString()
                    if len(points) == 4:
                        if levelId in levelListIDs:
                            levelAnzahlMaschen[levelId] = levelAnzahlMaschen[levelId] + 1
                            masche = []
                            punkt_a = []
                            punkt_b = []
                            punkt_c = []
                            point_a = points[0]
                            point_b = points[1]
                            point_c = points[2]

                            punkt_a = [point_a.x/g_1mu + g_go_x_mu, point_a.y/g_1mu + g_go_y_mu, point_a.z/g_1mu + g_go_z_mu]
                            punkt_b = [point_b.x/g_1mu + g_go_x_mu, point_b.y/g_1mu + g_go_y_mu, point_b.z/g_1mu + g_go_z_mu]
                            punkt_c = [point_c.x/g_1mu + g_go_x_mu, point_c.y/g_1mu + g_go_y_mu, point_c.z/g_1mu + g_go_z_mu]
                                                        
                            h_min = punkt_a[2]
                            h_max = punkt_a[2]
                            if punkt_b[2] < h_min: h_min = punkt_b[2]
                            if punkt_c[2] < h_min: h_min = punkt_c[2]
                            if punkt_b[2] > h_max: h_max = punkt_b[2]
                            if punkt_c[2] > h_max: h_max = punkt_c[2]
                            

                            masche.append(punkt_a)
                            masche.append(punkt_b)
                            masche.append(punkt_c)
                            
                            maschen_eines_levels[levelId].append(masche)

                            if levelAnzahlMaschen[levelId] == 1:
                                hoeheMin[levelId] = h_min
                                hoeheMax[levelId] = h_max

                                
                            if h_min < hoeheMin[levelId]: hoeheMin[levelId] = h_min
                            if h_max > hoeheMax[levelId]: hoeheMax[levelId] = h_max

                            anzahl = anzahl + 1
                            if anzahl % 1000 == 0:
                                progressbar_4.step(schritt)
                                root.update_idletasks()
                            if anzahl % 45000 == 0:
                                schritt = schritt * -1
                        
    # loesche alle items in treeview tree_3
    all_root_items = tree_4.get_children()
    tree_4.delete(*all_root_items)
    root.update_idletasks()
                              
    # daten aktualisieren
    ebenen = []
    for i in range(len(levelListIDs)):
        levelId_ = levelListIDs[i]
        #print(len(maschen_eines_levels[levelId_]))
        #if len(maschen_eines_levels[levelId_]) > 0:print(maschen_eines_levels[levelId_][0:2])
        if levelAnzahlMaschen[levelId_] > 0:
            #hoeheMittel = levelVolumen[levelId_] / levelGrundFlaeche[levelId_]
            ebenen.append((levelList[i], levelListIDs[i], levelAnzahlMaschen[levelId_],
              '{0:_.3f}'.format(hoeheMin[levelId_]),
               '{0:_.3f}'.format(hoeheMax[levelId_]), '{0:_.3f}'.format(hoeheMax[levelId_] - hoeheMin[levelId_]) ))

    # add data to the treeview
    for ebene in ebenen:
        tree_4.insert('', tk.END, values=ebene)
    
    # resize columns 
    tree_4.column('anzahl_maschen_4',anchor='center', stretch=False, width=175, minwidth=100)   
    tree_4.column('z_min_4',anchor='center', stretch=False, width=120, minwidth=50)
    tree_4.column('z_max_4',anchor='center', stretch=False, width=120, minwidth=50)
    tree_4.column('diff_z_4',anchor='center', stretch=False, width=120, minwidth=50)

    progressbar_4.stop()
    
    #button_11.config(state='active')
        
    root.update_idletasks()


def select_5_ElementsbyType():
    #dgn file scannen nach dgm maschen
    global g_1mu
    global g_go
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu
    global maschen_eines_levels
    maschen_eines_levels=[]
    levelAnzahlMaschen=[]
    #levelGrundFlaeche=[]
    #levelVolumen=[]
    #hoeheMittel=[]
    hoeheMin=[]
    hoeheMax=[]
    #x_mittel=[]
    #y_mittel=[]
    #z_mittel=[]
    #rechtsMin=[]
    #rechtsMax=[]
    #hochMin=[]
    #hochMax=[]
    for i in range(max(levelListIDs)+1):
        maschen_eines_levels.append([])
        levelAnzahlMaschen.append(0)
        #levelGrundFlaeche.append(0)
        #levelVolumen.append(0)
        #hoeheMittel.append(0)
        hoeheMin.append(0)
        hoeheMax.append(0)
        #x_mittel.append(0)
        #y_mittel.append(0)
        #z_mittel.append(0)
        #rechtsMin.append(0)
        #rechtsMax.append(0)
        #hochMin.append(0)
        #hochMax.append(0)

    schritt = 2
    progressbar_5.step(schritt)
    root.update_idletasks()
    anzahl = 0    
    #Get active model
    ACTIVEMODEL = ISessionMgr.ActiveDgnModelRef
    dgnModel = ACTIVEMODEL.GetDgnModel()
    #name =  model.GetModelName()
    #Get all graphical elements from the model
    graphicalElements = dgnModel.GetGraphicElements()
    

    for perElementRef in graphicalElements:
        elementId = perElementRef.GetElementId()
        eeh = EditElementHandle(perElementRef, dgnModel)
        eh = ElementHandle(perElementRef)

        msElement = MSElement()
        msElement = eeh.GetElement ()

        isGraphics = msElement.ehdr.isGraphics
        isInvisible = msElement.hdr.dhdr.props.b.invisible

        if (isGraphics and not(isInvisible)):
            eleType = eh.GetElementType()
            levelId = msElement.ehdr.level
            if (eleType == 6) :
                curve = ICurvePathQuery.ElementToCurveVector(eh)
                for element in curve:
                    points = element.GetLineString()
                    if len(points) == 4:
                        if levelId in levelListIDs:
                            levelAnzahlMaschen[levelId] = levelAnzahlMaschen[levelId] + 1
                            masche = []
                            punkt_a = []
                            punkt_b = []
                            punkt_c = []
                            point_a = points[0]
                            point_b = points[1]
                            point_c = points[2]

                            punkt_a = [point_a.x/g_1mu + g_go_x_mu, point_a.y/g_1mu + g_go_y_mu, point_a.z/g_1mu + g_go_z_mu]
                            punkt_b = [point_b.x/g_1mu + g_go_x_mu, point_b.y/g_1mu + g_go_y_mu, point_b.z/g_1mu + g_go_z_mu]
                            punkt_c = [point_c.x/g_1mu + g_go_x_mu, point_c.y/g_1mu + g_go_y_mu, point_c.z/g_1mu + g_go_z_mu]
                                                        
                            h_min = punkt_a[2]
                            h_max = punkt_a[2]
                            if punkt_b[2] < h_min: h_min = punkt_b[2]
                            if punkt_c[2] < h_min: h_min = punkt_c[2]
                            if punkt_b[2] > h_max: h_max = punkt_b[2]
                            if punkt_c[2] > h_max: h_max = punkt_c[2]
                            

                            masche.append(punkt_a)
                            masche.append(punkt_b)
                            masche.append(punkt_c)
                            
                            maschen_eines_levels[levelId].append(masche)

                            if levelAnzahlMaschen[levelId] == 1:
                                hoeheMin[levelId] = h_min
                                hoeheMax[levelId] = h_max

                                
                            if h_min < hoeheMin[levelId]: hoeheMin[levelId] = h_min
                            if h_max > hoeheMax[levelId]: hoeheMax[levelId] = h_max

                            anzahl = anzahl + 1
                            if anzahl % 1000 == 0:
                                progressbar_5.step(schritt)
                                root.update_idletasks()
                            if anzahl % 45000 == 0:
                                schritt = schritt * -1
                        
    # loesche alle items in treeview tree_5
    all_root_items = tree_5.get_children()
    tree_5.delete(*all_root_items)
    root.update_idletasks()
                              
    # daten aktualisieren
    ebenen = []
    for i in range(len(levelListIDs)):
        levelId_ = levelListIDs[i]
        #print(len(maschen_eines_levels[levelId_]))
        #if len(maschen_eines_levels[levelId_]) > 0:print(maschen_eines_levels[levelId_][0:2])
        if levelAnzahlMaschen[levelId_] > 0:
            #hoeheMittel = levelVolumen[levelId_] / levelGrundFlaeche[levelId_]
            ebenen.append((levelList[i], levelListIDs[i], levelAnzahlMaschen[levelId_],
              '{0:_.3f}'.format(hoeheMin[levelId_]),
               '{0:_.3f}'.format(hoeheMax[levelId_]), '{0:_.3f}'.format(hoeheMax[levelId_] - hoeheMin[levelId_]) ))

    # add data to the treeview
    for ebene in ebenen:
        tree_5.insert('', tk.END, values=ebene)
    
    # resize columns 
    tree_5.column('anzahl_maschen_5',anchor='center', stretch=False, width=175, minwidth=100)   
    tree_5.column('z_min_5',anchor='center', stretch=False, width=120, minwidth=50)
    tree_5.column('z_max_5',anchor='center', stretch=False, width=120, minwidth=50)
    tree_5.column('diff_z_5',anchor='center', stretch=False, width=120, minwidth=50)

    progressbar_5.stop()
    
    #button_11.config(state='active')
        
    root.update_idletasks()
   
def are_ident_points(punkt_a_, punkt_b_):
    if abs(punkt_a_[0] - punkt_b_[0]) > 0.0001: return False
    if abs(punkt_a_[1] - punkt_b_[1]) > 0.0001: return False
    if abs(punkt_a_[2] - punkt_b_[2]) > 0.0001: return False
    return True

def hoehenlinien_berechnen_2():
    global g_1mu
    global g_go
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu
    global faktor_projektion
    global hoehenlinien
    global polylines
    global maschen_eines_levels
    
    # keine ebene gewaehlt, abbruch
    if len(tree_4.selection()) == 0: return
    
    ebenenId_selected_ = []
    ebeneDatensatz_selected_ = []
    
    maschen = []
    hoehenlinien = []
    polylines = []
    
    hoehen_intervall_1 = float(options_h_intervall.get())
    
    for selected_item in tree_4.selection():
        item = tree_4.item(selected_item)
        record = item['values']
        id_=record[1]
        ebenenId_selected_.append(id_)
        ebeneDatensatz_selected_.append(record)
        #print(record)
        #print(id_)
        #print(len(maschen_eines_levels[id_]))
        

    schritt = 2
    progressbar_4.step(schritt)
    root.update_idletasks()
    anzahl = 0
    
    for masche in maschen_eines_levels[id_]:
        punkt_a = masche[0]
        punkt_b = masche[1]
        punkt_c = masche[2]
        punkt_a[2] = punkt_a[2] + 0.005
        punkt_b[2] = punkt_b[2] + 0.005
        punkt_c[2] = punkt_c[2] + 0.005
    
        # von unten nach oben: 1(a) unten, 3(c) oben
                                      
        if punkt_b[2] < punkt_a[2]:
            punkt_temp = punkt_a
            punkt_a = punkt_b
            punkt_b = punkt_temp
                            
        if punkt_c[2] < punkt_a[2]:
            punkt_temp = punkt_a
            punkt_a = punkt_c
            punkt_c = punkt_temp 
                                    
        if punkt_c[2] < punkt_b[2]:
            punkt_temp = punkt_b
            punkt_b = punkt_c
            punkt_c = punkt_temp
        
        #hoehenlinien berechnen
        delta_x21 = punkt_b[0] - punkt_a[0]
        delta_y21 = punkt_b[1] - punkt_a[1]
        delta_z21 = punkt_b[2] - punkt_a[2]
        
        delta_x31 = punkt_c[0] - punkt_a[0]
        delta_y31 = punkt_c[1] - punkt_a[1]
        delta_z31 = punkt_c[2] - punkt_a[2]
        
        delta_x32 = punkt_c[0] - punkt_b[0]
        delta_y32 = punkt_c[1] - punkt_b[1]
        delta_z32 = punkt_c[2] - punkt_b[2]
        
        hoehe_start = int(punkt_a[2] / hoehen_intervall_1) * hoehen_intervall_1
        hoehe_start = hoehe_start - hoehen_intervall_1
        hoehe = hoehe_start
        
        # starthoehe festlegen
        while hoehe < punkt_a[2]:
            hoehe = hoehe + hoehen_intervall_1
            
        hoehe = round(hoehe, 3)
        
        # jetzt hoehenlinie bestimmen: [[hoehe, 1] ,[punkt_von] ,[punkt_nach]]
        # zuerst von unten von punkt_a bis punkt_b
        while hoehe <= punkt_b[2]:
            hoehenlinie = []
            if abs(delta_z31) >= 0.00001:
                x_von = punkt_a[0] + delta_x31 * (hoehe - punkt_a[2]) / delta_z31
                y_von = punkt_a[1] + delta_y31 * (hoehe - punkt_a[2]) / delta_z31
                z_von = hoehe
                punkt_von = [x_von, y_von, z_von]
                hoehenlinie.append([hoehe, 1])
                hoehenlinie.append(punkt_von)
            else:
                punkt_von = punkt_a
                punkt_von[2] = hoehe
                hoehenlinie.append([hoehe, 1])
                hoehenlinie.append(punkt_von)
            
            if abs(delta_z21) >= 0.00001:
                x_nach = punkt_a[0] + delta_x21 * (hoehe - punkt_a[2]) / delta_z21
                y_nach = punkt_a[1] + delta_y21 * (hoehe - punkt_a[2]) / delta_z21
                z_nach = hoehe
                punkt_nach = [x_nach, y_nach, z_nach]
                hoehenlinie.append(punkt_nach)
            else:
                punkt_nach = punkt_b
                punkt_nach[2] = hoehe
                hoehenlinie.append(punkt_nach)
                
            hoehenlinien.append(hoehenlinie)
            
            hoehe = hoehe + hoehen_intervall_1
            hoehe = round(hoehe, 3)
            
        while hoehe <= punkt_c[2]:
            hoehenlinie = []
            if abs(delta_z31) >= 0.00001:
                x_von = punkt_a[0] + delta_x31 * (hoehe - punkt_a[2]) / delta_z31
                y_von = punkt_a[1] + delta_y31 * (hoehe - punkt_a[2]) / delta_z31
                z_von = hoehe
                punkt_von = [x_von, y_von, z_von]
                hoehenlinie.append([hoehe, 1])
                hoehenlinie.append(punkt_von)
            else:                                
                punkt_von = punkt_c
                punkt_von[2] = hoehe
                hoehenlinie.append([hoehe, 1])
                hoehenlinie.append(punkt_von)
            
            if abs(delta_z32) >= 0.00001:
                x_nach = punkt_b[0] + delta_x32 * (hoehe - punkt_b[2]) / delta_z32
                y_nach = punkt_b[1] + delta_y32 * (hoehe - punkt_b[2]) / delta_z32
                z_nach = hoehe
                punkt_nach = [x_nach, y_nach, z_nach]
                hoehenlinie.append(punkt_nach)
            else:
                punkt_nach = punkt_b
                punkt_nach[2] = hoehe
                hoehenlinie.append(punkt_nach)
                
            hoehenlinien.append(hoehenlinie)
            
            hoehe = hoehe + hoehen_intervall_1
            hoehe = round(hoehe, 3)
            
        #ende hoehenlinien berechnen

        anzahl = anzahl + 1
        if anzahl % 1000 == 0:
            progressbar_4.step(schritt)
            root.update_idletasks()
        if anzahl % 45000 == 0:
            schritt = schritt * -1  

    hoehenlinien.sort()
    #die hoehenlinien sind nun nach hoehenwert [hoehe, 1] sortiert
    #und lagemaessig von links unten nach rechts oben
    
    label_anz_hoehenlinien.config(text=' Anzahl Hoehenlinien (segmente): ' + str(len(hoehenlinien)))
    
    #gruppieren der hoehenlinien nach hoehenwert [hoehe, 1]
    #jede hoehengruppe kann mehr als eine Polylinien enthalten
    hoehen_gruppen = itertools.groupby(hoehenlinien, key=lambda x: x[0])
    for key, gruppe in hoehen_gruppen:
        puffer = []
        lines = list(gruppe)
        #print(key, '--->',lines)
        anzahl_1en = len(lines)
        #print(f'{anzahl_1en = }')
        #kennzeichen für neue polylinie in der gruppe suchen:
        #falls nach einem durchlauf noch linien mit kennung 1 vorhanden sind
        neu = True
        
        while anzahl_1en > 0:
            for i, linie in enumerate(lines[:]):
                if linie[0][1] == 0: continue   #die linie ist bereits verarbeitet, also ueberspringen
                if neu:
                    punkte_dq = collections.deque([])
                    punkt_vorne = lines[i][1]
                    punkt_hinten = lines[i][2]
                    punkte_dq.append(punkt_vorne)
                    punkte_dq.append(punkt_hinten)
                    lines[i][0][1] = 0  # kennung 0 = linie als verarbeitet markieren
                    neu = False
                    #print()
                    #print(lines[i])
                    #print()
                    continue
                idx = i + 0
                #print(idx, linie)
                #if linie[0][1] == 0: continue   #ist bereits verarbeitet
                punkta = linie[1]
                punktb = linie[2]    
                #test fuer vorne anhaengen
                test_1 = are_ident_points(punkt_vorne, punkta)
                if test_1:
                    punkte_dq.appendleft(punktb)
                    punkt_vorne = punktb
                    #print(f'{test_1 = }')
                    lines[idx][0][1] = 0
                    #print(lines[idx])
                    continue        
                test_2 = are_ident_points(punkt_vorne, punktb)
                if test_2:
                    punkte_dq.appendleft(punkta)
                    punkt_vorne = punkta
                    #print(f'{test_2 = }')
                    lines[idx][0][1] = 0
                    #print(lines[idx])
                    continue        
                #test fuer hinten anhaengen
                test_3 = are_ident_points(punkt_hinten, punkta)
                if test_3:
                    punkte_dq.append(punktb)
                    punkt_hinten = punktb
                    #print(f'{test_3 = }')
                    lines[idx][0][1] = 0
                    #print(lines[idx])
                    continue
                test_4 = are_ident_points(punkt_hinten, punktb)
                if test_4:
                    punkte_dq.append(punkta)
                    punkt_hinten = punkta
                    #print(f'{test_4 = }')
                    lines[idx][0][1] = 0
                    #print(lines[idx])
                    continue
                #print(test_1, test_2, test_3, test_4)
                
            #print(punkte_dq)
            #im puffer werden die polygonlinien eines hoehenwertes zwischengespeichert
            #die gerade gefundene polygonlinie koennte ein teilabschnitt einer laengeren
            #polygonlinie sein. im puffer wird nun nachgeschaut, ob ein dort gespeichertes
            #teilstueck am anfang oder ende angehaengt werden kann.
            #ein passendes teilstueck wird als verarbeitet gekennzeichnet, weitere teilstuecke
            #werden im puffer gesucht und verarbeitet.
            #zum schluss wird die (eventuell) ergaenzte polygonlinie im puffer abgelegt
 
            if len(puffer) > 0:
                for j,zeile in enumerate(puffer):
                    #print('pufferzeile :',j, zeile)
                    test=zeile[0][1]
                    if test == 0: continue  # bereits verarbeitet
                    liste_punkte = zeile[1]
                    #print('puffer punkte in zeile :',j, liste_punkte)
                    #print(type(liste_punkte))
                    punkt_a = liste_punkte[0]
                    punkt_b = liste_punkte[-1]
                    #print('a=',punkt_a, 'b=',punkt_b, 'vo=',punkt_vorne,'hi=',punkt_hinten)
                    #test fuer vorne anhaengen
                    test_12 = are_ident_points(punkt_vorne, punkt_a)
                    if test_12:
                        punkte_dq.popleft()
                        punkte_dq.extendleft(liste_punkte)
                        punkt_vorne = punkt_b
                        #print(f'{test_12 = }')
                        puffer[j][0][1] = 0
                        #print(puffer[j])
                        continue
                    test_22 = are_ident_points(punkt_vorne, punkt_b)
                    if test_22:
                        punkte_dq.popleft()
                        #liste_punkte in umgekehrter reihenfolge anhangen
                        punkte_dq.extendleft(liste_punkte[::-1])
                        punkt_vorne = punkt_a
                        #print(f'{test_22 = }')
                        puffer[j][0][1] = 0
                        #print(puffer[j])
                        continue
                    #test fuer hinten anhaengen
                    test_32 = are_ident_points(punkt_hinten, punkt_a)
                    if test_32:
                        punkte_dq.pop()
                        punkte_dq.extend(liste_punkte)
                        punkt_hinten = punkt_b
                        #print(f'{test_32 = }')
                        puffer[j][0][1] = 0
                        #print(puffer[j])
                        continue
                    test_42 = are_ident_points(punkt_hinten, punkt_b)
                    if test_42:
                        punkte_dq.pop()
                        #liste_punkte in umgekehrter reihenfolge anhangen
                        punkte_dq.extend(liste_punkte[::-1])
                        punkt_hinten = punkt_a
                        #print(f'{test_42 = }')
                        puffer[j][0][1] = 0
                        #print(puffer[j])
                        continue                    

            #polylines.append([key, list(punkte_dq)])
            neu = True
            
            puffer.append([[key[0],2], list(punkte_dq)])
            
            #print('puffer = ',puffer)
            #for k,zeile_ in enumerate(puffer):
                #print('puffer zeile : ', k, zeile_, '\n')
                
            zaehler_kennung_1 = 0
            for line in lines:
                if line[0][1] == 1:
                    zaehler_kennung_1 = zaehler_kennung_1 + 1
            #print(f' { zaehler_kennung_1 = }')
            anzahl_1en = zaehler_kennung_1
            # -- hier ende der while schleife anzahl_1en

        for k,zeile_ in enumerate(puffer):
            test=zeile_[0][1]
            if test == 0: continue  # bereits verarbeitet
            polylines.append(zeile_)

    #print(polylines)
    #print()
    #for polyline in polylines:
    #    print(polyline)    
    
    #print(hoehenlinien[0])
    #print(hoehenlinien[1])
    #print(hoehenlinien[2])  
    #print(hoehenlinien[-1])
    
    label_anz_hoehenlinien.config(text=' Anzahl Hoehenlinien (segmente): ' + str(len(hoehenlinien)) + ' in Anzahl polylinien : ' + str(len(polylines)) )    
  
    progressbar_4.stop()
    
    if len(hoehenlinien) > 0:
        button_cad_hl.config(state='active')
    else:
        button_cad_hl.config(state='disabled')
        
    return    

def hoehenlinien_zeichnen():
    global ACTIVEMODEL
    global dgnModel
    global g_1mu
    global g_go
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu
    global levelList
    global levelListIDs        
    global hoehenlinien
    global polylines
    global hl_werte
    
    schritt = 2
    progressbar_4.step(schritt)
    root.update_idletasks()
    anzahl = 0    
    
    for selected_item in tree_4.selection():
        item = tree_4.item(selected_item)
        record = item['values']
        z_min_4 = float(record[3]) - float(options_h_intervall.get())
        z_max_4 = float(record[4]) + float(options_h_intervall.get())
        #ebenenId_selected_.append(id_)
        #ebeneDatensatz_selected_.append(record)
        
    anzahl_linien = 0
    anzahl_polylinien = 0
    
    level_name = str(selected_level_hl.get())
    levelID = 0
    if level_name in levelList:
        j = levelList.index(level_name)
        levelId = levelListIDs[j]
    hl_ebene_id = levelId
    
    hl_0_color = int(color_button_hl_0.cget('text'))
    hl_2_color = int(color_button_hl_2.cget('text'))
    hl_3_color = int(color_button_hl_3.cget('text'))
    hl_4_color = int(color_button_hl_4.cget('text'))
    
    hl_0_style = int(options_strichart_hl_0.get())
    hl_2_style = int(options_strichart_hl_2.get())
    hl_3_style = int(options_strichart_hl_3.get())
    hl_4_style = int(options_strichart_hl_4.get())
    
    hl_0_weight = int(options_strichbreite_hl_0.get())
    hl_2_weight = int(options_strichbreite_hl_2.get())
    hl_3_weight = int(options_strichbreite_hl_3.get())
    hl_4_weight = int(options_strichbreite_hl_4.get())
    
    hl_2_intervall = int(options_intervall_hl_2.get()) * float(options_h_intervall.get())
    hl_3_intervall = int(options_intervall_hl_3.get()) * float(options_h_intervall.get())
    hl_4_intervall = int(options_intervall_hl_4.get()) * float(options_h_intervall.get())
    
    hoehen_hl_2 = []
    hoehen_hl_3 = []
    hoehen_hl_4 = []
    
    hoehe_start_2 = int(z_min_4 / hl_2_intervall) * hl_2_intervall
    hoehe_start_3 = int(z_min_4 / hl_3_intervall) * hl_3_intervall
    hoehe_start_4 = int(z_min_4 / hl_4_intervall) * hl_4_intervall
    
    hoehe = hoehe_start_2
    hoehen_hl_2.append(hoehe)
    while hoehe < z_max_4:
        hoehe = hoehe + hl_2_intervall
        hoehen_hl_2.append(hoehe)
        
    hoehe = hoehe_start_3
    hoehen_hl_3.append(hoehe)
    while hoehe < z_max_4:
        hoehe = hoehe + hl_3_intervall
        hoehen_hl_3.append(hoehe)
        
    hoehe = hoehe_start_4
    hoehen_hl_4.append(hoehe)
    while hoehe < z_max_4:
        hoehe = hoehe + hl_4_intervall
        hoehen_hl_4.append(hoehe)
    
    hl_typ = 2   #linestrings default 
    hl_wert_s = str(options_hl_art.get())
    
    if hl_wert_s == hl_werte[0]:
        hl_typ = 1   #lines
    if hl_wert_s == hl_werte[1]:
        hl_typ = 2   #linestrings
    if hl_wert_s == hl_werte[2]:
        hl_typ = 3   #bsplines               
    
    #hl_typ = 1   #linien
    #hl_typ = 2   #linestrings
    #hl_typ = 3   #bsplines
    if len(hoehenlinien) > 0:
        PyCadInputQueue.SendKeyin("mark")
        button2_undo_mark.config(state='active')
        if hl_typ == 1:
            for hoehenlinie in hoehenlinien:
                anzahl_linien = anzahl_linien + 1
                hoehe = hoehenlinie[0][0]
                kennung =  hoehenlinie[0][1]
                punkt_von = hoehenlinie[1]
                punkt_nach = hoehenlinie[2]
                color = hl_0_color
                style = hl_0_style
                weight = hl_0_weight
                #if anzahl_linien == 1:
                    #print(hoehe)
                    #print(kennung)
                    #print(punkt_von)
                    #print(punkt_nach)
                if zeichne_intervall_hl_2_var.get() == True: 
                    if hoehe in hoehen_hl_2:
                        color = hl_2_color
                        style = hl_2_style
                        weight = hl_2_weight
                if zeichne_intervall_hl_3_var.get() == True: 
                    if hoehe in hoehen_hl_3:
                        color = hl_3_color
                        style = hl_3_style
                        weight = hl_3_weight
                if zeichne_intervall_hl_4_var.get() == True: 
                    if hoehe in hoehen_hl_4:
                        color = hl_4_color
                        style = hl_4_style
                        weight = hl_4_weight
                
                test = createLineElement(punkt_von, punkt_nach, hl_ebene_id, color, style, weight)
                
                anzahl = anzahl + 1
                if anzahl % 1000 == 0:
                    progressbar_4.step(schritt)
                    root.update_idletasks()
                if anzahl % 45000 == 0:
                    schritt = schritt * -1                
                
            #print(anzahl_linien)
            #print(hl_ebene_id)
            #print(hl_0_color, hl_2_color, hl_3_color, hl_4_color)
            #print(hl_0_style, hl_2_style, hl_3_style, hl_4_style)
            #print(hl_0_weight, hl_2_weight, hl_3_weight, hl_4_weight)
            #print(z_min_4, z_max_4)
            #print(hl_2_intervall, hl_3_intervall, hl_4_intervall)
            #print(hoehen_hl_2, hoehen_hl_3, hoehen_hl_4)
            #print(color, style, weight)
            progressbar_4.stop()
        elif hl_typ == 2:
            for polylinie in polylines:
                anzahl_polylinien = anzahl_polylinien + 1
                hoehe = polylinie[0][0]
                kennung =  polylinie[0][1]
                punktliste = polylinie[1]
                #if anzahl_polylinien < 3:
                    #print(hoehe, kennung)
                    #print(punktliste[0:3], '\n')
                    #punkt=punktliste[0]
                    #print(punkt)
                color = hl_0_color
                style = hl_0_style
                weight = hl_0_weight

                if zeichne_intervall_hl_2_var.get() == True: 
                    if hoehe in hoehen_hl_2:
                        color = hl_2_color
                        style = hl_2_style
                        weight = hl_2_weight
                if zeichne_intervall_hl_3_var.get() == True: 
                    if hoehe in hoehen_hl_3:
                        color = hl_3_color
                        style = hl_3_style
                        weight = hl_3_weight
                if zeichne_intervall_hl_4_var.get() == True: 
                    if hoehe in hoehen_hl_4:
                        color = hl_4_color
                        style = hl_4_style
                        weight = hl_4_weight
                
                test = createLineStringElement(punktliste, hl_ebene_id, color, style, weight)
                
                anzahl = anzahl + 1
                if anzahl % 1000 == 0:
                    progressbar_4.step(schritt)
                    root.update_idletasks()
                if anzahl % 45000 == 0:
                    schritt = schritt * -1            
            progressbar_4.stop()
            
        elif hl_typ == 3:
            for polylinie in polylines:
                anzahl_polylinien = anzahl_polylinien + 1
                hoehe = polylinie[0][0]
                kennung =  polylinie[0][1]
                punktliste = polylinie[1]
                #if anzahl_polylinien < 3:
                    #print(hoehe, kennung)
                    #print(punktliste[0:3], '\n')
                    #punkt=punktliste[0]
                    #print(punkt)
                color = hl_0_color
                style = hl_0_style
                weight = hl_0_weight

                if zeichne_intervall_hl_2_var.get() == True: 
                    if hoehe in hoehen_hl_2:
                        color = hl_2_color
                        style = hl_2_style
                        weight = hl_2_weight
                if zeichne_intervall_hl_3_var.get() == True: 
                    if hoehe in hoehen_hl_3:
                        color = hl_3_color
                        style = hl_3_style
                        weight = hl_3_weight
                if zeichne_intervall_hl_4_var.get() == True: 
                    if hoehe in hoehen_hl_4:
                        color = hl_4_color
                        style = hl_4_style
                        weight = hl_4_weight
                
                #test = createLineStringElement(punktliste, hl_ebene_id, color, style, weight)
                if len(punktliste) > 3:
                    test = createBsplineElement(punktliste, hl_ebene_id, color, style, weight)
                else:
                    test = createLineStringElement(punktliste, hl_ebene_id, color, style, weight)
                
                anzahl = anzahl + 1
                if anzahl % 1000 == 0:
                    progressbar_4.step(schritt)
                    root.update_idletasks()
                if anzahl % 45000 == 0:
                    schritt = schritt * -1            
            progressbar_4.stop()
                        

        PyCadInputQueue.SendReset()
        PyCadInputQueue.SendKeyin("FIT VIEW EXTENDED")
        PyCommandState.StartDefaultCommand()

        lift_window(root)    
    return    

def umring_berechnen():
    global g_1mu
    global g_go
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu
    global polylines
    global maschen_eines_levels
    
    # keine ebene gewaehlt, abbruch
    if len(tree_5.selection()) == 0: return
    
    ebenenId_selected_ = []
    ebeneDatensatz_selected_ = []
    
    maschen = []
    polylines = []
    
    seite = []
    seiten = []
    
    hoehe_grund = float(um_un_hoehe_eingabe.get())
    
    for selected_item in tree_5.selection():
        item = tree_5.item(selected_item)
        record = item['values']
        id_=record[1]
        ebenenId_selected_.append(id_)
        ebeneDatensatz_selected_.append(record)
        #print(record)
        #print(id_)
        #print(len(maschen_eines_levels[id_]))
        

    schritt = 2
    progressbar_5.step(schritt)
    root.update_idletasks()
    anzahl = 0
    
    for masche in maschen_eines_levels[id_]:
        #print(masche)
        seite_01=[masche[0], masche[1]]
        seite_02=[masche[0], masche[2]]
        seite_12=[masche[1], masche[2]]
        s_mitte_01=[round((masche[0][0] + masche[1][0])/2.0, 3),
                    round((masche[0][1] + masche[1][1])/2.0, 3),
                    round((masche[0][2] + masche[1][2])/2.0, 3)]
        s_mitte_02=[round((masche[0][0] + masche[2][0])/2.0, 3),
                    round((masche[0][1] + masche[2][1])/2.0, 3),
                    round((masche[0][2] + masche[2][2])/2.0, 3)]
        s_mitte_12=[round((masche[1][0] + masche[2][0])/2.0, 3),
                    round((masche[1][1] + masche[2][1])/2.0, 3),
                    round((masche[1][2] + masche[2][2])/2.0, 3)]    
        seiten.append([s_mitte_01, seite_01])
        seiten.append([s_mitte_02, seite_02])
        seiten.append([s_mitte_12, seite_12])
        
    seiten.sort()
        
    seiten_gruppen = itertools.groupby(seiten, key=lambda x: x[0])

    zwischenpuffer = []

    for key, gruppe in seiten_gruppen:
        zeilen = list(gruppe)
        #print(key, len(zeilen))
        if len(zeilen) == 1:
            zwischenpuffer.append([[2, 1],zeilen[0][1]])
            
    #print('zwischenpuffer ', len(zwischenpuffer))
    
    puffer = []
    lines = zwischenpuffer
    anzahl_1en = len(lines)
    #print(f'{anzahl_1en = }')
    #kennzeichen für neue polylinie in der gruppe suchen falls nach einem durchlauf noch linien mit kennung 1 vorhanden sind
    neu = True
    
    while anzahl_1en > 0:
        for i, linie in enumerate(lines[:]):
            key = linie[0]
            if linie[0][1] == 0: continue   #die linie ist bereits verarbeitet, also ueberspringen
            if neu:
                punkte_dq = collections.deque([])
                punkt_vorne = lines[i][1][0]
                punkt_hinten = lines[i][1][1]
                punkte_dq.append(punkt_vorne)
                punkte_dq.append(punkt_hinten)
                lines[i][0][1] = 0  # kennung 0 = linie als verarbeitet markieren
                neu = False
                #print()
                #print('++++',lines[i])
                #print()
                continue
            idx = i + 0
            #print(idx, linie)
            #if linie[0][1] == 0: continue   #ist bereits verarbeitet
            punkta = linie[1][0]
            punktb = linie[1][1]   
            #test fuer vorne anhaengen
            test_1 = are_ident_points(punkt_vorne, punkta)
            if test_1:
                punkte_dq.appendleft(punktb)
                punkt_vorne = punktb
                #print(f'{test_1 = }')
                lines[idx][0][1] = 0
                #print(lines[idx])
                continue        
            test_2 = are_ident_points(punkt_vorne, punktb)
            if test_2:
                punkte_dq.appendleft(punkta)
                punkt_vorne = punkta
                #print(f'{test_2 = }')
                lines[idx][0][1] = 0
                #print(lines[idx])
                continue        
            #test fuer hinten anhaengen
            test_3 = are_ident_points(punkt_hinten, punkta)
            if test_3:
                punkte_dq.append(punktb)
                punkt_hinten = punktb
                #print(f'{test_3 = }')
                lines[idx][0][1] = 0
                #print(lines[idx])
                continue
            test_4 = are_ident_points(punkt_hinten, punktb)
            if test_4:
                punkte_dq.append(punkta)
                punkt_hinten = punkta
                #print(f'{test_4 = }')
                lines[idx][0][1] = 0
                #print(lines[idx])
                continue
            #print(test_1, test_2, test_3, test_4)
            
        #print(punkte_dq)
            
        if len(puffer) > 0:
            for j,zeile in enumerate(puffer):
                #print('pufferzeile :',j, zeile)
                test=zeile[0][1]
                if test == 0: continue  # bereits verarbeitet
                liste_punkte = zeile[1]
                #print('puffer punkte in zeile :',j, liste_punkte)
                #print(type(liste_punkte))
                punkt_a = liste_punkte[0]
                punkt_b = liste_punkte[-1]
                #print('a=',punkt_a, 'b=',punkt_b, 'vo=',punkt_vorne,'hi=',punkt_hinten)
                #test fuer vorne anhaengen
                test_12 = are_ident_points(punkt_vorne, punkt_a)
                if test_12:
                    punkte_dq.popleft()
                    punkte_dq.extendleft(liste_punkte)
                    punkt_vorne = punkt_b
                    #print(f'{test_12 = }')
                    puffer[j][0][1] = 0
                    #print(puffer[j])
                    continue
                test_22 = are_ident_points(punkt_vorne, punkt_b)
                if test_22:
                    punkte_dq.popleft()
                    #punkte_dq.extendleft(liste_punkte.reverse())
                    punkte_dq.extendleft(liste_punkte[::-1])
                    punkt_vorne = punkt_a
                    #print(f'{test_22 = }')
                    puffer[j][0][1] = 0
                    #print(puffer[j])
                    continue
                #test fuer hinten anhaengen
                test_32 = are_ident_points(punkt_hinten, punkt_a)
                if test_32:
                    punkte_dq.pop()
                    punkte_dq.extend(liste_punkte)
                    punkt_hinten = punkt_b
                    #print(f'{test_32 = }')
                    puffer[j][0][1] = 0
                    #print(puffer[j])
                    continue
                test_42 = are_ident_points(punkt_hinten, punkt_b)
                if test_42:
                    punkte_dq.pop()
                    punkte_dq.extend(liste_punkte[::-1])
                    punkt_hinten = punkt_a
                    #print(f'{test_42 = }')
                    puffer[j][0][1] = 0
                    #print(puffer[j])
                    continue                    

        #polylines.append([key, list(punkte_dq)])
        neu = True
        
        puffer.append([[key[0],2], list(punkte_dq)])
        
        #print('puffer = ',puffer)
        #for k,zeile_ in enumerate(puffer):
            #print('puffer zeile : ', k, zeile_, '\n')
            
        
        zaehler_kennung_1 = 0
        for line in lines:
            if line[0][1] == 1:
                zaehler_kennung_1 = zaehler_kennung_1 + 1
        #print(f' { zaehler_kennung_1 = }')
        anzahl_1en = zaehler_kennung_1

    for k,zeile_ in enumerate(puffer):
        test=zeile_[0][1]
        if test == 0: continue  # bereits verarbeitet
        polylines.append(zeile_)
        
    anz_um_punkte = 0
    for polyline in polylines:
        #print('111',polyline, '\n')
        #print('222',polyline[1],'\n')
        polyline_oben = polyline[1]
        anz_um_punkte = anz_um_punkte + len(polyline_oben) - 1
        #polyline_unten = []
        #for punkt in polyline_oben:
            #punkt_x = punkt[0]
            #punkt_y = punkt[1]
            #punkt_z = hoehe_grund
            #print(punkt)
            #polyline_unten.append([punkt_x, punkt_y, punkt_z])
        #print('polylinie oben  : Anzahl Punkte ', len(polyline_oben), polyline_oben[0:2], polyline_oben[-1])
        #print('polylinie unten : Anzahl Punkte ', len(polyline_unten)) 
        
    label_anz_umring.config(text=' Anzahl Elemente Umring (polylines): ' + str(len(polylines)) + '   mit Seiten :  ' + str(anz_um_punkte) )    
    
    progressbar_5.stop()
    
    if len(polylines) > 0:
        button_cad_um.config(state='active')
    else:
        button_cad_um.config(state='disabled')   

    return 

def umring_zeichnen():
    global ACTIVEMODEL
    global dgnModel
    global g_1mu
    global g_go
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu
    global levelList
    global levelListIDs        
    global polylines
    global hl_werte
    
    schritt = 2
    progressbar_5.step(schritt)
    root.update_idletasks()
    anzahl = 0
    
    if len(polylines) == 0:
        return
    
    PyCadInputQueue.SendKeyin("mark")
    button3_undo_mark.config(state='active') 
        
    hoehe_grund = float(um_un_hoehe_eingabe.get())
    
    level_name = str(selected_level_um_ob.get())
    levelID = 0
    if level_name in levelList:
        j = levelList.index(level_name)
        levelId = levelListIDs[j]
    um_ob_ebene_id = levelId
    
    level_name = str(selected_level_um_un.get())
    levelID = 0
    if level_name in levelList:
        j = levelList.index(level_name)
        levelId = levelListIDs[j]
    um_un_ebene_id = levelId  
    
    level_name = str(selected_level_um_sei.get())
    levelID = 0
    if level_name in levelList:
        j = levelList.index(level_name)
        levelId = levelListIDs[j]
    um_sei_ebene_id = levelId        
    
    um_ob_color = int(color_button_um_ob.cget('text'))
    um_un_color = int(color_button_um_un.cget('text'))
    um_sei_color = int(color_button_um_sei.cget('text'))
    
        
    for polyline in polylines:
        #print('111',polyline, '\n')
        #print('222',polyline[1],'\n')
        polyline_oben = polyline[1]
        polyline_unten = []
        for punkt in polyline_oben:
            punkt_x = punkt[0]
            punkt_y = punkt[1]
            punkt_z = hoehe_grund
            #print(punkt)
            polyline_unten.append([punkt_x, punkt_y, punkt_z])
        #print('polylinie oben  : Anzahl Punkte ', len(polyline_oben), polyline_oben[0:2], polyline_oben[-1])
        #print('polylinie unten : Anzahl Punkte ', len(polyline_unten))
        shapes = []
        shapes_tri = []
        for i in range(len(polyline_oben) - 1):
            shape = []
            tri1 = []
            tri2 = []
            shape.append([polyline_unten[i], polyline_oben[i], polyline_oben[i+1],
                          polyline_unten[i+1], polyline_unten[i]])
            shapes.append(shape)
            tri1.append([polyline_unten[i], polyline_oben[i], polyline_unten[i+1], polyline_unten[i]])
            tri2.append([polyline_unten[i+1], polyline_oben[i], polyline_oben[i+1], polyline_unten[i+1]])
            shapes_tri.append(tri1)
            shapes_tri.append(tri2) 
                
        #for i in range(len(shapes)):
            #print('shape ',i, shapes[i])
            #j = 2 * i
            #print('tri1  ', j, shapes_tri[j])
            #print('tri2  ', j, shapes_tri[j+1])
        
        if zeichne_um_ob_var.get() == True:
            color = um_ob_color
            style = 0
            weight = 0
            test = createLineStringElement(polyline_oben, um_ob_ebene_id, color, style, weight)
            
        if zeichne_um_un_var.get() == True:
            color = um_un_color
            style = 0
            weight = 0
            punkt_1 = polyline_unten[0]
            punkt_2 = polyline_unten[-1]
            test_sh = are_ident_points(punkt_1, punkt_2) 
            if test_sh:
                test = createShapeElement_2(polyline_unten, um_un_ebene_id, color, style, weight)
            else:
                test = createLineStringElement(polyline_unten, um_un_ebene_id, color, style, weight)
                
        if zeichne_um_sei_var.get() == True:
            color = um_sei_color
            style = 0
            weight = 0
            
            sei_wert_s = str(options_sei_art.get())
            if sei_wert_s == um_sei_werte[0]:  #triangle
                for element in shapes_tri:
                    shape_tri = element[0]
                    #print(shape_tri)
                    test = createShapeElement_2(shape_tri, um_sei_ebene_id, color, style, weight)
            else:
                for element in shapes:
                    shape_quad = element[0]
                    #print(shape_quad)
                    test = createShapeElement_2(shape_quad, um_sei_ebene_id, color, style, weight)
            
    PyCadInputQueue.SendReset()
    PyCadInputQueue.SendKeyin("FIT VIEW EXTENDED")
    PyCommandState.StartDefaultCommand()

    lift_window(root)                       

    return

def export_dgm_to_lamdxml():
    global maschen_eines_levels
    
    # keine ebene gewaehlt, abbruch
    if len(tree_5.selection()) == 0: return
    
    ebenenId_selected_ = []
    ebeneDatensatz_selected_ = []

    for selected_item in tree_5.selection():
        item = tree_5.item(selected_item)
        record = item['values']
        ebenen_name_ = record[0]
        id_=record[1]
        ebenenId_selected_.append(id_)
        ebeneDatensatz_selected_.append(record)
        #print(record)
        #print(id_)
        #print(len(maschen_eines_levels[id_]))
        
    schritt = 2
    progressbar_5.step(schritt)
    root.update_idletasks()
    anzahl = 0
    
    current_dateTime = datetime.now()
    datum = current_dateTime.strftime("%Y-%m-%d")
    zeit = current_dateTime.strftime("%H:%M:%S")
    protokoll = []

    punkt_texte = []
    faces = []
    koordinaten_liste_ = []
    koordinaten_set_ = set()
    koordinaten_dict_ = {}

    for masche in maschen_eines_levels[id_]:
        #print(masche)
        #face = ''
        for punkt_ in masche[0:3]:
            koordinate_ = "{1:.4f} {0:.4f} {2:.4f}".format(punkt_[0], punkt_[1], punkt_[2])
            koordinaten_set_.add(koordinate_)

        anzahl = anzahl + 1
        if anzahl % 1000 == 0:
            progressbar_5.step(schritt)
            root.update_idletasks()
        if anzahl % 45000 == 0:
            schritt = schritt * -1 
            
    nummer = 0       
    for element in koordinaten_set_:
        nummer = nummer + 1
        koordinaten_dict_[element] = str(nummer)
        
    for masche in maschen_eines_levels[id_]:
        #print(masche)
        face = ''
        for punkt_ in masche[0:3]:
            koordinate_ = "{1:.4f} {0:.4f} {2:.4f}".format(punkt_[0], punkt_[1], punkt_[2])
            nummer = koordinaten_dict_[koordinate_]
            face = face + ' ' + nummer
            
        faces.append(face[1:]) #das erste leerzeichen ueberlesen
        
        anzahl = anzahl + 1
        if anzahl % 1000 == 0:
            progressbar_5.step(schritt)
            root.update_idletasks()
        if anzahl % 45000 == 0:
            schritt = schritt * -1                
    punkt_texte = []
    for koordinate, nummer in koordinaten_dict_.items():
        punkt_text = ' '*5*3 + '<P id="' + nummer + '">' + koordinate + '</P>'
        punkt_texte.append(punkt_text)
 
    zeile = ' '*0*3 + '<?xml version="1.0"?>'
    #print(zeile)
    protokoll.append(zeile)
    
    zeile = ' '*0*3 + '<LandXML xmlns="http://www.landxml.org/schema/LandXML-1.2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.landxml.org/schema/LandXML-1.2 http://www.landxml.org/schema/LandXML-1.2/LandXML-1.2.xsd" version="1.2" readOnly="false" language="English" ' #date="2018-08-27" time="17:56:22">'
    zeile = zeile + 'date="' + datum + '" time="' + zeit + '">'
    #print(zeile)
    protokoll.append(zeile)
    
    zeile = ' '*1*3 + '<Units>'
    #print(zeile)
    protokoll.append(zeile)
    
    zeile = ' '*2*3 + '<Metric linearUnit="meter" areaUnit="squareMeter" volumeUnit="cubicMeter"/>'
    #print(zeile)
    protokoll.append(zeile)
    
    zeile = ' '*1*3 + '</Units>'
    #print(zeile)
    protokoll.append(zeile)
    
    zeile = ' '*1*3 + '<Application name="KGE_DGM_Triangle_delau_pslg_lines_tk.py" URL="https://github.com/KGE-Brem/MicroStation_2024_2025_python_pgms"></Application>'
    #print(zeile)
    protokoll.append(zeile)
    
    zeile = ' '*1*3 + '<Surfaces>'
    #print(zeile)
    protokoll.append(zeile)
        
    #zeile = ' '*2*3 + '<Surface name="level_name">'
    zeile = ' '*2*3 + '<Surface name="' + ebenen_name_ + '">'
    #print(zeile)
    protokoll.append(zeile)
    
    zeile = ' '*3*3 + '<Definition surfType="TIN">'
    #print(zeile)
    protokoll.append(zeile)
    
    zeile = ' '*4*3 + '<Pnts>'
    #print(zeile)
    protokoll.append(zeile)
        
    for zeile in punkt_texte:
        #print(zeile)
        protokoll.append(zeile)
        
    zeile = ' '*4*3 + '</Pnts>'
    #print(zeile)
    protokoll.append(zeile)
    
    zeile = ' '*4*3 + '<Faces>'
    #print(zeile)
    protokoll.append(zeile)
    
    for face in faces:
        face_text = ' '*5*3 + '<F>' + face + '</F>'
        #print(face_text)
        protokoll.append(face_text)
        
    zeile = ' '*4*3 + '</Faces>'
    #print(zeile)
    protokoll.append(zeile)
    
    zeile = ' '*3*3 + '</Definition>'
    #print(zeile)
    protokoll.append(zeile)
    
    zeile = ' '*2*3 + '</Surface>'
    #print(zeile)
    protokoll.append(zeile)

    zeile = ' '*1*3 + '</Surfaces>'
    #print(zeile)
    protokoll.append(zeile)

    zeile = ' '*0*3 + '</LandXML>'
    #print(zeile)
    protokoll.append(zeile)

    #for zeile in protokoll:
        #print(zeile)
    file_save_4(protokoll)

    progressbar_5.stop()

    return

# 2 funktionen fuer zwischenpunkte - punktabstand
def get_intervall():
    _, wert = label_intervall.cget('text').split(':')
    return float(wert)

def callback(selection):
    label_intervall.config(text=f"Punktabstand ist ca: {selection}")
    
def callback_hoehe_hl1(selection):
    label_intervall_h2.config(text=" = " + str(int(options_intervall_hl_2.get()) * float(options_h_intervall.get())) + "  ")
    label_intervall_h3.config(text=" = " + str(int(options_intervall_hl_3.get()) * float(options_h_intervall.get())) + "  ")    
    label_intervall_h4.config(text=" = " + str(int(options_intervall_hl_4.get()) * float(options_h_intervall.get())) + "  ")
    
def callback_hoehe_hl2(selection):
    label_intervall_h2.config(text=" = " + str(int(options_intervall_hl_2.get()) * float(options_h_intervall.get())) + "  ")

def callback_hoehe_hl3(selection):
    label_intervall_h3.config(text=" = " + str(int(options_intervall_hl_3.get()) * float(options_h_intervall.get())) + "  ")

def callback_hoehe_hl4(selection):
    label_intervall_h4.config(text=" = " + str(int(options_intervall_hl_4.get()) * float(options_h_intervall.get())) + "  ")

def koordinaten_leeren():
    global punkte_set
    punkte_set=set()
    label_koor_datei.config(text=f"im Koordinatenspeicher sind : {len(punkte_set)} Punktkoordinaten")

if __name__ == '__main__':

    global ACTIVEMODEL
    global dgnModel
    global g_1mu
    global g_go
    global dgm_color
    global selected_level
    global levelList
    global levelListIDs
    global farben_
    global farbe_fuer
    global punkte_set
    global faktor_projektion
    global g_go_x_mu
    global g_go_y_mu
    global g_go_z_mu
    global hoehenlinien
    global polylines
    global maschen_eines_levels
    global dgn_pfad
    
    punkte_set=set()    # koordinaten zum export - keine doppelten
    
    farbe_fuer = 0  # fuer DGM
  
    # Get the active DGN model reference
    ACTIVEMODEL = ISessionMgr.ActiveDgnModelRef
    if ACTIVEMODEL != None:
            dgnModel = ACTIVEMODEL.GetDgnModel()
            modelInfo = dgnModel.GetModelInfo()
            g_1mu = modelInfo.GetUorPerStorage()
            g_go = modelInfo.GetGlobalOrigin()
            
            #print(g_go, g_go.x*-1/g_1mu, g_go.y*-1/g_1mu, g_go.z*-1/g_1mu)
            
            g_go_x_mu = g_go.x*-1/g_1mu
            g_go_y_mu = g_go.y*-1/g_1mu
            g_go_z_mu = g_go.z*-1/g_1mu

            dgnfile = dgnModel.GetDgnFile()
            
                      
            dgnFile_ = ISessionMgr.GetActiveDgnFile()
            dgn_dateiname = str(dgnFile_.GetFileName())
            dgn_name = os.path.basename(dgn_dateiname)
            dgn_pfad = os.path.dirname(dgn_dateiname)
            dgn_sep = os.path.sep
            #print(dgn_dateiname, dgn_name, dgn_pfad, dgn_sep)
            
            farben_ = []

            for color_ in range(0,255):
                eleinfoColor = DgnColorMap.ExtractElementColorInfo(color_ , dgnfile)
                rgbdef       = RgbColorDef
                rgbdef       = eleinfoColor[1].m_rgb #RGB in IntColorDef
                red   = rgbdef.red
                green = rgbdef.green
                blue  = rgbdef.blue
                farben_.append(f'#{red << 16 | green << 8 | blue:06x}')
            #print("color ", color_ ,"red,green,blue-->",red, green, blue) 
            #print(farben_[0:10])

            root = tk.Tk()
            root.title("KGE-DGM basierend auf Triangle (J.Shewchuk) und triangle (D.Rufat): hier PSLG: segmente (lines): umring und innen - punkte innen - lines in holes ")

            notebook = ttk.Notebook(root)
            notebook.pack(pady=10, expand=True)

            frameHaupt = ttk.Frame(notebook, width=600, height=480)
            frameFarben = ttk.Frame(notebook, width=600, height=480)
            frameExportLines = ttk.Frame(notebook, width=600, height=480)
            frameExportKoordinaten = ttk.Frame(notebook, width=600, height=480)
            frameVolumen = ttk.Frame(notebook, width=600, height=480)
            frameHoehenlinien = ttk.Frame(notebook, width=600, height=480)
            frameUmring = ttk.Frame(notebook, width=600, height=480)


            frameHaupt.pack(fill='both', expand=True)
            frameFarben.pack(fill='both', expand=True)
            frameExportLines.pack(fill='both', expand=True)
            frameExportKoordinaten.pack(fill='both', expand=True)
            frameVolumen.pack(fill='both', expand=True)
            frameHoehenlinien.pack(fill='both', expand=True)
            frameUmring.pack(fill='both', expand=True)

            notebook.add(frameHaupt, text='  Hauptprogramm  ')
            notebook.add(frameFarben, text='  Farbpicker aus Colortable  ')
            notebook.add(frameExportLines, text='  Segmente exportieren - aussen, innen, loecher  ')
            notebook.add(frameExportKoordinaten, text='  Koordinaten exportieren  ')
            notebook.add(frameVolumen, text=' Volumenberechnung ')
            notebook.add(frameHoehenlinien, text=' Hoehenlinien ')
            notebook.add(frameUmring, text=' Umring ')

            #------frameHaupt-------
            
            text_ = "zuerst Farb-Nr fuer DGM und DGM auf Ebene einstellen \n"
            text_ = text_ + " dann  -Im Programm fortfahren- klicken \n"
            text_ = text_ + " nacheinander auswaehlen bzw. abbrechen \n"
            text_ = text_ + "1. Koordinatendatei Segmente Umring (lines) \n"
            text_ = text_ + "2. Koordinatendatei Segmente Innen (lines) \n"
            text_ = text_ + "3. eine oder mehrere Koordinaten Dateien Punkte innen \n"
            text_ = text_ + "4. Koordinatendatei Linien innerhalb Loecher"


            hinweise_label = tk.Label(frameHaupt, text=text_)
            hinweise_label.pack(padx=20, pady=20)


            file_text = scrolledtext.ScrolledText(frameHaupt, wrap=tk.WORD, height=20, width=100)
            file_text.pack(padx=20, pady=20)

            
            color_label = tk.Label(frameHaupt, text = "Farb-Nr waehlen fuer das DGM :")
            color_label.pack()

            dgmColor = 9
            w_color = tk.Scale(frameHaupt, from_=0, to=254, length=1200, tickinterval=15, orient=tk.HORIZONTAL, command=farbe_aendern)
            w_color.set(dgmColor)
            w_color.pack()
            
            #color_button = tk.Button(frameHaupt, height=1, width=10, text = farben_[dgmColor], bg=farben_[dgmColor], command=select_tab)
            color_button = tk.Button(frameHaupt, height=1, width=10, text = dgmColor, bg=farben_[dgmColor], command=select_tab)
            color_button.pack(padx=20, pady=10)
            
            #--- FrameFarben - FrameFarben2  - colorpicker
            
            frameFarben2 = Frame(master=frameFarben)
            frameFarben2.pack()

            n=255 # anzahl farb-buttons
            i=0 # zeile 
            j=0 # spalte 

            buttons = []
            for k in range(n):
                but = tk.Button(frameFarben2, text=k,height=1,width=3, bg=farben_[k], command=lambda k=k: farb_pick(k)) 
                but.grid(row=i, column=j, padx=1, pady=1)
                buttons.append(but)
                j=j+1
                if(j%16==0):
                    i=i+1
                    j=0
            
            farbe = tk.Button(frameFarben, text=dgmColor, bg=farben_[dgmColor], height=1, width=6, padx=2, pady=2, command=select_tab_H )
            farbe.pack()

            #---- FrameExportLines
            
            frame3_1 = ttk.Frame(master=frameExportLines)
            frame3_2 = ttk.Frame(master=frameExportLines)
            frame3_1.pack(pady=20)
            frame3_2.pack(pady=20)


            # define columns
            columns = ('level_name', 'level_id', 'anzahl_elemente', 'lines', 'linestrings')

            tree_1 = ttk.Treeview(frame3_1, columns=columns, show='headings')

            # define headings
            tree_1.column('level_name',anchor='center', stretch=False, width=400)
            tree_1.heading('level_name', text='Level Name')
            tree_1.column('level_id',anchor='center', stretch=False, width=150)
            tree_1.heading('level_id', text='Level ID')
            tree_1.column('anzahl_elemente',anchor='center', stretch=False, width=150)
            tree_1.heading('anzahl_elemente', text='Anzahl Elemente')
            tree_1.column('lines',anchor='center', stretch=False, width=150)
            tree_1.heading('lines', text='Anzahl lines')
            tree_1.column('linestrings',anchor='center', stretch=False, width=150)
            tree_1.heading('linestrings', text='Anzahl linestrings')
 
            # Level
            levelList,levelListIDs = GetLevelList()

            # erste daten erzeugen
            ebenen = []
            for i in range(len(levelListIDs)):
                ebenen.append((levelList[i], levelListIDs[i], -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1))

            # und eintragen
            for ebene in ebenen:
                tree_1.insert('', tk.END, values=ebene)

            tree_1.grid(row=0, column=0, sticky='nsew')

            # add a scrollbar
            scrollbar = ttk.Scrollbar(frame3_1, orient=tk.VERTICAL, command=tree_1.yview)
            tree_1.configure(yscroll=scrollbar.set)
            scrollbar.grid(row=0, column=1, sticky='ns')
            

            button_1 = tk.Button(frame3_1, text='1. Zeichnung scannen (lines  linestrings)', command=selectElementsbyType_3_4) 
            button_1.grid(row=1, column=0)

            label_ex2 = tk.Label(frame3_1,text='2. Ebene(n) markieren          ')
            label_ex2.grid(row=3, column=0,pady = 20)

            label_ex3 = tk.Label(frame3_2,text='3. Art der Segmente waehlen :   ')
            label_ex3.pack(anchor="w", padx=10, pady=10)

            exports = ["segmente_umring", "segmente_innen", "lines_in_holes"]
            label_ex = tk.Label(frame3_2,text=f"KGE_DGM_Triangle_{exports[0]}.tmp")
            
            variable_ = tk.StringVar(frame3_2, f"{exports[0]}")

            def selection():
                label_ex.config(text=f"KGE_DGM_Triangle_{variable_.get()}.tmp")

            for export in exports:
                tk.Radiobutton(frame3_2, text=export, variable=variable_, value=export, command=selection,).pack(anchor="w", padx=10, pady=5) 


            label_ex.pack(anchor="w", padx=10, pady=10) 
             
            button_2 = tk.Button(frame3_2, text='4. markierte Ebene(n) auswerten und Koordinaten speichern      ', command=exportElementsbyType_3_4)
            button_2.pack(padx=10, pady=10) 
            button_2.config(state='disabled')

            #---------FrameExportKoordianten
            
            frame4_1 = ttk.Frame(master=frameExportKoordinaten)
            frame4_2 = ttk.Frame(master=frameExportKoordinaten)
            frame4_1.pack(pady=20)
            frame4_2.pack(pady=20)

            # define columns
            columns_2 = ('level_name_2', 'level_id_2', 'anzahl_elemente_2', 'lines_2', 'linestrings_2', 'shapes_2', 'arcs_2', 'ellipses_2', 'bsplines_2', 'complex_ele_2')

            tree_2 = ttk.Treeview(frame4_1, columns=columns_2, show='headings')            

            # define headings
            tree_2.column('level_name_2',anchor='center', stretch=False, width=200, minwidth=200)
            tree_2.heading('level_name_2', text='Level Name')
            tree_2.column('level_id_2',anchor='center', stretch=False, width=100, minwidth=75)
            tree_2.heading('level_id_2', text='Level ID')
            tree_2.column('anzahl_elemente_2',anchor='center', stretch=False, width=125, minwidth=100)
            tree_2.heading('anzahl_elemente_2', text='Anz. Elemente')
            tree_2.column('lines_2',anchor='center', stretch=False, width=110, minwidth=100)
            tree_2.heading('lines_2', text=' lines')
            tree_2.column('linestrings_2',anchor='center', stretch=False, width=110, minwidth=100)
            tree_2.heading('linestrings_2', text=' linestrings')
            tree_2.column('shapes_2',anchor='center', stretch=False, width=110, minwidth=100)
            tree_2.heading('shapes_2', text=' shapes') 
            tree_2.column('arcs_2',anchor='center', stretch=False, width=110, minwidth=100)
            tree_2.heading('arcs_2', text=' arcs') 
            tree_2.column('ellipses_2',anchor='center', stretch=False, width=110, minwidth=100)
            tree_2.heading('ellipses_2', text=' ellipses')
            tree_2.column('bsplines_2',anchor='center', stretch=False, width=110, minwidth=100)
            tree_2.heading('bsplines_2', text=' bsplines')
            tree_2.column('complex_ele_2',anchor='center', stretch=False, width=110, minwidth=100)
            tree_2.heading('complex_ele_2', text=' complex')                                                
 
            # ebenen eintragen
            for ebene in ebenen:
                tree_2.insert('', tk.END, values=ebene)

            tree_2.grid(row=0, column=0, sticky='nsew')

            # add a scrollbar
            scrollbar_2 = ttk.Scrollbar(frame4_1, orient=tk.VERTICAL, command=tree_2.yview)
            tree_2.configure(yscrollcommand=scrollbar_2.set)
            scrollbar_2.grid(row=0, column=1, sticky='ns') 
            
            scrollbar_hor = ttk.Scrollbar(frame4_1, orient=tk.HORIZONTAL, command=tree_2.xview)
            tree_2.configure(xscrollcommand=scrollbar_hor.set)
            scrollbar_hor.grid(row=1, column=0, sticky='ew')            

            button_5 = tk.Button(frame4_1, text='1. Zeichnung scannen (lines  linestrings shapes)', command=select_2_ElementsbyType_3_4) 
            button_5.grid(row=2, column=0, ipadx=60)

            label_ex5 = tk.Label(frame4_1,text='2. Ebene(n) markieren, danach Elementtypen ausschliessen  ')
            label_ex5.grid(row=3, column=0,pady = 10)
            
            export_lines_var = tk.BooleanVar(value=True)
            export_linestrings_var = tk.BooleanVar(value=True)
            export_arcs_var = tk.BooleanVar(value=True)
            export_bsplines_var = tk.BooleanVar(value=True)
            export_shapes_var = tk.BooleanVar(value=True)
            export_complex_var = tk.BooleanVar(value=True)
            
            checkbox_ex_li = ttk.Checkbutton(frame4_1,text='lines auswerten?',variable=export_lines_var)
            checkbox_ex_ls = ttk.Checkbutton(frame4_1,text='linestrings auswerten?',variable=export_linestrings_var)
            checkbox_ex_ar = ttk.Checkbutton(frame4_1,text='arcs auswerten?',variable=export_arcs_var)
            checkbox_ex_bs = ttk.Checkbutton(frame4_1,text='bsplines  auswerten  ?',variable=export_bsplines_var)
            checkbox_ex_sh = ttk.Checkbutton(frame4_1,text='shapes auswerten ?',variable=export_shapes_var)
            checkbox_ex_co = ttk.Checkbutton(frame4_1,text='complex chains/shapes auswerten ?',variable=export_complex_var)
            
            checkbox_ex_li.grid(sticky='w', row=4, column=0)
            checkbox_ex_ls.grid(sticky='ns', row=4, column=0)
            checkbox_ex_sh.grid(sticky='e', row=4, column=0)
            checkbox_ex_ar.grid(sticky='w', row=5, column=0)
            checkbox_ex_bs.grid(sticky='ns', row=5, column=0)
            checkbox_ex_co.grid(sticky='e', row=5, column=0)
            
            # entlang der elemente weitere zwischenpunkte berechnen im ca. abstand von ....

            weitere_punkte_var = tk.BooleanVar()
            checkbox = ttk.Checkbutton(frame4_1,text='3. weitere Zwischenpunkte berechnen ? ',variable=weitere_punkte_var)
            checkbox.grid(sticky='w', row=7, column=0)
            
            intervall_werte = ('0.05','0.10','0.20', '0.25', '0.50', '1.00', '2.00', '2.50', '5.00', '10.00', '12.50', '20.00', '25.00', '50.00', '100.00')

            label_intervall = tk.Label(frame4_1,text=f"Punktabstand ist ca. : {intervall_werte[6]}")
            label_intervall.grid(sticky='ns', row=7, column=0)

            options_abst = tk.StringVar()
            menu_abst = tk.OptionMenu(frame4_1, options_abst, *intervall_werte, command=callback)
            menu_abst.grid(sticky='e', row=7, column=0)
            options_abst.set(intervall_werte[6])
            
            label_hinweis_ellipsen = tk.Label(frame4_1,text=f"     Hinweis zu 3. : Kreise und Ellipsen vorher brechen mit    trim break bypoint   ")
            label_hinweis_ellipsen.grid(sticky='w', row=8, column=0) 
            
            #button_5_2 = tk.Button(frame4_1, text='>trim break bypoint< senden', command=trim_ellipse_one_point) 
            #button_5_2.grid(sticky='e', row=5, column=0)
            
            #trim_ellipse_one_point

            button_6 = tk.Button(frame4_2, text='4. markierte Ebene(n) auswerten und Koordinaten berechnen      ', command=export_2_ElementsbyType_3_4)
            button_6.pack(padx=10, pady=10) 
            button_6.config(state='disabled')

 
            button_7 = tk.Button(frame4_2, text='  bei Bedarf Koordinatenspeicher leeren      ', command=koordinaten_leeren)
            button_7.pack(padx=10, pady=10, ipadx=85) 

            label_koor_datei = tk.Label(frame4_2,text=f"im Koordinatenspeicher sind : {len(punkte_set)} Punktkoordinaten")
            label_koor_datei.pack(padx=10, pady=5)
                         

            label_ex_datei = tk.Label(frame4_2,text=f"KGE_DGM_Triangle_exportierte_Koordinaten.txt")
            label_ex_datei.pack(padx=10, pady=10)
            
            button_8 = tk.Button(frame4_2, text='5.  Koordinaten speichern                    ', command=exportKoordinaten)
            button_8.pack(padx=10, pady=10, ipadx=105) 
            button_8.config(state='disabled') 

            button_9 = tk.Button(frame4_2, text='0.  externe Koordinatendatei(en) hinzufuegen ', command=importKoordinaten)
            button_9.pack(padx=10, pady=10, ipadx=78) 
            button_9.config(state='active')    
            
            #-------- frame Volumen

            frame5_1 = ttk.Frame(master=frameVolumen)
            #frame5_2 = ttk.Frame(master=frameVolumen)
            frame5_1.pack(pady=20)
            #frame5_2.pack(pady=20) 
            
            #- define columns
            columns_3 = ('level_name_3', 'level_id_3', 'anzahl_maschen_1', 'grundflaeche_1', 'volumen_1',
             'z_aus_1', 'z_min_1', 'z_max_1', 'x_min_1', 'x_max_1', 'y_min_1','y_max_1', 'x_mittel_1', 'y_mittel_1', 'z_mittel_1')

            tree_3 = ttk.Treeview(frame5_1, columns=columns_3, show='headings', selectmode="browse")  # nur ein DGM kann selektiert werden
            
            # define headings
            tree_3.column('level_name_3',anchor='center', stretch=False, width=200, minwidth=200)
            tree_3.heading('level_name_3', text='Level Name')
            tree_3.column('level_id_3',anchor='center', stretch=False, width=100, minwidth=75)
            tree_3.heading('level_id_3', text='Level ID')
            tree_3.column('anzahl_maschen_1',anchor='center', stretch=False, width=120, minwidth=100)
            tree_3.heading('anzahl_maschen_1', text='Anz. Maschen')
            tree_3.column('grundflaeche_1',anchor='center', stretch=False, width=120, minwidth=100)
            tree_3.heading('grundflaeche_1', text='Grundflaeche')
            tree_3.column('volumen_1',anchor='center', stretch=False, width=120, minwidth=100)
            tree_3.heading('volumen_1', text='Volumen gegen 0')
            tree_3.column('z_aus_1',anchor='center', stretch=False, width=50, minwidth=50)
            tree_3.heading('z_aus_1', text='mittlere Hoehe') 
            tree_3.column('z_min_1',anchor='center', stretch=False, width=50, minwidth=50)
            tree_3.heading('z_min_1', text='z min')
            tree_3.column('z_max_1',anchor='center', stretch=False, width=50, minwidth=50)
            tree_3.heading('z_max_1', text='z max')                                                
            tree_3.column('x_min_1',anchor='center', stretch=False, width=50, minwidth=50)
            tree_3.heading('x_min_1', text='x min')
            tree_3.column('x_max_1',anchor='center', stretch=False, width=50, minwidth=50)
            tree_3.heading('x_max_1', text='x max')
            tree_3.column('y_min_1',anchor='center', stretch=False, width=50, minwidth=50)
            tree_3.heading('y_min_1', text='y min')
            tree_3.column('y_max_1',anchor='center', stretch=False, width=50, minwidth=50)
            tree_3.heading('y_max_1', text='y max')
            tree_3.column('x_mittel_1',anchor='center', stretch=False, width=50, minwidth=50)
            tree_3.heading('x_mittel_1', text='x mittel')
            tree_3.column('y_mittel_1',anchor='center', stretch=False, width=50, minwidth=50)
            tree_3.heading('y_mittel_1', text='y mittel')
            tree_3.column('z_mittel_1',anchor='center', stretch=False, width=50, minwidth=50)
            tree_3.heading('z_mittel_1', text='z mittel')                        
                                                
            # ebenen eintragen
            for ebene in ebenen:
                tree_3.insert('', tk.END, values=ebene)

            tree_3.grid(row=0, column=0, sticky='nsew')

            # add scrollbars
            scrollbar_3 = ttk.Scrollbar(frame5_1, orient=tk.VERTICAL, command=tree_3.yview)
            tree_3.configure(yscrollcommand=scrollbar_3.set)
            scrollbar_3.grid(row=0, column=1, sticky='ns') 

            scrollbar_hor_3 = ttk.Scrollbar(frame5_1, orient=tk.HORIZONTAL, command=tree_3.xview)
            tree_3.configure(xscrollcommand=scrollbar_hor_3.set)
            scrollbar_hor_3.grid(row=1, column=0, sticky='ew')             
            
            progressbar_1 = ttk.Progressbar(frame5_1, orient=tk.HORIZONTAL, length=700, mode="determinate")
            progressbar_1.grid(row=2, column=0, sticky='e')  
            
            button_10 = tk.Button(frame5_1, text='1. Zeichnung scannen (Maschen)', command=select_3_ElementsbyType)
            button_10.grid(row=2, column=0, ipadx=60, sticky='w')
            
            frame5_1_1 = ttk.Frame(master=frame5_1)
            frame5_1_1.grid(sticky='w', row=3, column=0) 
            
            label_vol = tk.Label(frame5_1_1,text='       2. Ebene markieren      danach -> ')
            label_vol.grid(row=0, column=0,pady = 10, sticky='ns')
 
            label_abr_hoehe = tk.Label(frame5_1_1, text=' Abrechnungshoehe eingeben :  ')
            label_abr_hoehe.grid(row=0, column=1, sticky='ns')                      
            
            abr_hoehe_eingabe = tk.Entry(frame5_1_1)
            abr_hoehe_eingabe.insert(0,"0.000")
            abr_hoehe_eingabe.grid(row=0, column=2, sticky='e')
            
            utm_var = tk.BooleanVar(value=False) 
            checkbox_utm = ttk.Checkbutton(frame5_1,text='  utm massstab beruecksichtigen?   ',variable=utm_var) 
            checkbox_utm.grid(sticky='w', row=4, column=0)
            
            frame5_1_2 = ttk.Frame(master=frame5_1)
            frame5_1_2.grid(sticky='ns', row=4, column=0)
            
            label_radius = tk.Label(frame5_1_2, text='Schmiegungskugel mittl. Kruemmungs-Radius in km :  ')
            label_radius.grid(row=0, column=0, sticky='w')
            
            radius_eingabe = tk.Entry(frame5_1_2)
            radius_eingabe.insert(0,"6382")
            radius_eingabe.grid(row=0, column=1, sticky='e')            
              
            frame5_1_3 = ttk.Frame(master=frame5_1)
            frame5_1_3.grid(sticky='e', row=4, column=0) 
            
            label_geoid_un = tk.Label(frame5_1_3, text='Geoid Undulation in m :  ')
            label_geoid_un.grid(row=0, column=0, sticky='w') 
            
            geoid_un_eingabe = tk.Entry(frame5_1_3)
            geoid_un_eingabe.insert(0,"43.0")
            geoid_un_eingabe.grid(row=0, column=1, sticky='e') 
            
            utm_hoehen_label = tk.Label(frame5_1, text=" Beachten und einstellen bei   utm : Umkehrung der Hoehenreduktion auf --> ")
            utm_hoehen_label.grid(row=5, column=0, sticky='ns') 
            
            hoehen_werte = ('Hoehen-Mittel je Seite','Hoehen-Mittel je Masche','Hoehen-Mittel Gebiet') 
            
            options_utm_hoehe = tk.StringVar()
            menu_utm_hoehe = tk.OptionMenu(frame5_1, options_utm_hoehe, *hoehen_werte)
            menu_utm_hoehe.grid(sticky='e', row=5, column=0)
            options_utm_hoehe.set(hoehen_werte[1])                                            
            
            button_11 = tk.Button(frame5_1, text='3. Volumen berechnen ', command=dgm_volumen_berechnen_2)
            button_11.grid(row=6, column=0, ipadx=60, sticky='e') 
            button_11.config(state='disabled')                                                                                                                                   

            button_12 = tk.Button(frame5_1, text='Berechnung loeschen ', command=berechnung_loeschen)
            button_12.grid(row=6, column=0, ipadx=60, sticky='w')
            
            button_13 = tk.Button(frame5_1, text='Kopfzeilen eintragen ', command=berechnung_kopf)
            button_13.grid(row=6, column=0, ipadx=60, sticky='ns')  
            
            file_text_vol = scrolledtext.ScrolledText(frame5_1, wrap=tk.WORD, height=15, width=90)
            file_text_vol.grid(row=7, column=0, pady=20, sticky='ns')  
                        
            button_14 = tk.Button(frame5_1, text='Protokoll speichern ', command=berechnung_speichern)
            button_14.grid(row=8, column=0, ipadx=60, sticky='ns') 
            
            #-------- frame Hoehenlinien
            frame6_1 = ttk.Frame(master=frameHoehenlinien)
            frame6_2 = ttk.Frame(master=frameHoehenlinien)
            frame6_1.pack(pady=20)
            frame6_2.pack(pady=20) 
            
            #- define columns
            columns_4 = ('level_name_4', 'level_id_4', 'anzahl_maschen_4',
             'z_min_4', 'z_max_4', 'diff_z_4')

            tree_4 = ttk.Treeview(frame6_1, columns=columns_4, show='headings', selectmode="browse")  # nur ein DGM kann selektiert werden
            
            # define headings 
            tree_4.column('level_name_4',anchor='center', stretch=False, width=200, minwidth=200)
            tree_4.heading('level_name_4', text='Level Name')
            tree_4.column('level_id_4',anchor='center', stretch=False, width=120, minwidth=75)
            tree_4.heading('level_id_4', text='Level ID')
            tree_4.column('anzahl_maschen_4',anchor='center', stretch=False, width=150, minwidth=100)
            tree_4.heading('anzahl_maschen_4', text='Anz. Maschen') 
            tree_4.column('z_min_4',anchor='center', stretch=False, width=150, minwidth=50)
            tree_4.heading('z_min_4', text='z min')
            tree_4.column('z_max_4',anchor='center', stretch=False, width=150, minwidth=50)
            tree_4.heading('z_max_4', text='z max')                                                
            tree_4.column('diff_z_4',anchor='center', stretch=False, width=150, minwidth=50)
            tree_4.heading('diff_z_4', text='z differenz') 

            # ebenen eintragen
            for ebene in ebenen:
                tree_4.insert('', tk.END, values=ebene)

            tree_4.grid(row=0, column=0, sticky='nsew')
                                                                                 
            # add scrollbars
            scrollbar_4 = ttk.Scrollbar(frame6_1, orient=tk.VERTICAL, command=tree_4.yview)
            tree_4.configure(yscrollcommand=scrollbar_4.set)
            scrollbar_4.grid(row=0, column=1, sticky='ns') 

            scrollbar_hor_4 = ttk.Scrollbar(frame6_1, orient=tk.HORIZONTAL, command=tree_4.xview)
            tree_4.configure(xscrollcommand=scrollbar_hor_4.set)
            scrollbar_hor_4.grid(row=1, column=0, sticky='ew')             
            
            progressbar_4 = ttk.Progressbar(frame6_1, orient=tk.HORIZONTAL, length=600, mode="determinate")
            progressbar_4.grid(row=2, column=0, sticky='e') 
            
            button_15 = tk.Button(frame6_1, text='1. Zeichnung scannen (Maschen)', command=select_4_ElementsbyType)
            button_15.grid(row=2, column=0, ipadx=30, sticky='w')
            
            frame6_1_1 = ttk.Frame(master=frame6_1)
            frame6_1_1.grid(sticky='w', row=3, column=0) 
            
            label_hoehl = tk.Label(frame6_1_1,text='       2. Ebene markieren      danach -> ')
            label_hoehl.grid(row=0, column=0,pady = 10, sticky='ns')            

            label_hoehl_intervall = tk.Label(frame6_1_1, text=' Hoehenlinien Grundintervall auswaehlen :  ')
            label_hoehl_intervall.grid(row=0, column=1, sticky='ns')                      
            
            h_intervall_werte = ('0.05','0.10','0.20', '0.25', '0.50', '1.00', '2.00', '2.50', '5.00', '10.00', '12.50', '20.00', '25.00', '50.00', '100.00')

            options_h_intervall = tk.StringVar()
            menu_intervall = tk.OptionMenu(frame6_1_1, options_h_intervall, *h_intervall_werte, command = callback_hoehe_hl1)
            menu_intervall.grid(row=0, column=2, sticky='e')
            options_h_intervall.set(h_intervall_werte[5])

            label_hoehl_ebene = tk.Label(frame6_1_1, text=' Ebene fuer alle Hoehenlinien auswaehlen :  ')
            label_hoehl_ebene.grid(row=0, column=3, sticky='ns')
            
            selected_level_hl = tk.StringVar()
            selected_level_hl.set(levelList[0])
            hl_level_option = tk.OptionMenu(frame6_1_1, selected_level_hl, *levelList) 
            hl_level_option.grid(row=0, column=4, sticky='e')
            
            stricharten_werte = ('0', '1', '2', '3', '4', '5', '6', '7')
            strichbreiten_werte = ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10')
            intervall_vielfaches = ('1', '2', '4', '5', '8', '10', '16', '20', '25', '40', '50', '100') 
            
            frame7_1_1 = ttk.Frame(master=frame6_1)
            frame7_1_1.grid(sticky='w', row=4, column=0) 
            
            label_fabe_hl_0 = tk.Label(frame7_1_1,text='Fuer Grundintervall einstellen: Farbe:')
            label_fabe_hl_0.grid(row=0, column=0,pady = 10, sticky='w') 
            
            hl_0_Color = 9
           
            color_button_hl_0 = tk.Button(frame7_1_1, height=1, width=10, text = hl_0_Color, bg=farben_[hl_0_Color], command=select_tab_h0)
            color_button_hl_0.grid(row=0, column=1,pady = 10, sticky='e')                       

            label_strichart_hl_0 = tk.Label(frame7_1_1, text='     Strichart :  ')
            label_strichart_hl_0.grid(row=0, column=2, sticky='w')                      
            
            options_strichart_hl_0 = tk.StringVar()
            options_strichart_hl_0.set(stricharten_werte[0]) 
            strichart_hl_0_option = tk.OptionMenu(frame7_1_1, options_strichart_hl_0, *stricharten_werte) 
            strichart_hl_0_option.grid(row=0, column=3, sticky='e')
            

            label_strichbreite_hl_0 = tk.Label(frame7_1_1, text='     Strichbreite :  ')
            label_strichbreite_hl_0.grid(row=0, column=4, sticky='w')             

            options_strichbreite_hl_0 = tk.StringVar()
            options_strichbreite_hl_0.set(strichbreiten_werte[0]) 
            strichbreite_hl_0_option = tk.OptionMenu(frame7_1_1, options_strichbreite_hl_0, *strichbreiten_werte) 
            strichbreite_hl_0_option.grid(row=0, column=5, sticky='e')         
            
            # ---------------------------------------------------                                                           
                     
            frame8_1_1 = ttk.Frame(master=frame6_1)
            frame8_1_1.grid(sticky='w', row=5, column=0) 
            
            label_fabe_hl_2 = tk.Label(frame8_1_1,text='Fuer       2. Intervall einstellen: Farbe:')
            label_fabe_hl_2.grid(row=0, column=0,pady = 10, sticky='w') 
            
            hl_2_Color = 10
           
            color_button_hl_2 = tk.Button(frame8_1_1, height=1, width=10, text = hl_2_Color, bg=farben_[hl_2_Color], command=select_tab_h2)
            color_button_hl_2.grid(row=0, column=1,pady = 10, sticky='e')                       

            label_strichart_hl_2 = tk.Label(frame8_1_1, text='     Strichart :  ')
            label_strichart_hl_2.grid(row=0, column=2, sticky='w')                      
            
            options_strichart_hl_2 = tk.StringVar()
            options_strichart_hl_2.set(stricharten_werte[0]) 
            strichart_hl_2_option = tk.OptionMenu(frame8_1_1, options_strichart_hl_2, *stricharten_werte) 
            strichart_hl_2_option.grid(row=0, column=3, sticky='e')
            
            label_strichbreite_hl_2 = tk.Label(frame8_1_1, text='     Strichbreite :  ')
            label_strichbreite_hl_2.grid(row=0, column=4, sticky='w')             

            options_strichbreite_hl_2 = tk.StringVar()
            options_strichbreite_hl_2.set(strichbreiten_werte[1]) 
            strichbreite_hl_2_option = tk.OptionMenu(frame8_1_1, options_strichbreite_hl_2, *strichbreiten_werte) 
            strichbreite_hl_2_option.grid(row=0, column=5, sticky='e')  
            
            label_intervall_multi_hl_2 = tk.Label(frame8_1_1, text='     Grundintervall mal:  ')
            label_intervall_multi_hl_2.grid(row=0, column=6, sticky='w')  
            
            options_intervall_hl_2 = tk.StringVar()
            options_intervall_hl_2.set(intervall_vielfaches[3])
            intervall_multi_hl_2_option = tk.OptionMenu(frame8_1_1, options_intervall_hl_2, *intervall_vielfaches, command=callback_hoehe_hl2)
            intervall_multi_hl_2_option.grid(row=0, column=7, sticky='e')

            label_intervall_h2 = tk.Label(frame8_1_1,text="   = " + str(int(options_intervall_hl_2.get()) * float(options_h_intervall.get())) + "  ")
            label_intervall_h2.grid(sticky='ns', row=0, column=8)              
            
            zeichne_intervall_hl_2_var = tk.BooleanVar(value=False) 
            checkbox_hl_2_cad = ttk.Checkbutton(frame8_1_1,text='  in CAD zeichnen?', variable=zeichne_intervall_hl_2_var)  
            checkbox_hl_2_cad.grid(row=0, column=9, sticky='e')      

            # ---------------------------------------------------
                     
            frame9_1_1 = ttk.Frame(master=frame6_1)
            frame9_1_1.grid(sticky='w', row=6, column=0) 
            
            label_fabe_hl_3 = tk.Label(frame9_1_1,text='Fuer       3. Intervall einstellen: Farbe:')
            label_fabe_hl_3.grid(row=0, column=0,pady = 10, sticky='w') 
            
            hl_3_Color = 11
           
            color_button_hl_3 = tk.Button(frame9_1_1, height=1, width=10, text = hl_3_Color, bg=farben_[hl_3_Color], command=select_tab_h3)
            color_button_hl_3.grid(row=0, column=1,pady = 10, sticky='e')                       

            label_strichart_hl_3 = tk.Label(frame9_1_1, text='     Strichart :  ')
            label_strichart_hl_3.grid(row=0, column=2, sticky='w')                      
            
            options_strichart_hl_3 = tk.StringVar()
            options_strichart_hl_3.set(stricharten_werte[0]) 
            strichart_hl_3_option = tk.OptionMenu(frame9_1_1, options_strichart_hl_3, *stricharten_werte) 
            strichart_hl_3_option.grid(row=0, column=3, sticky='e')
            
            label_strichbreite_hl_3 = tk.Label(frame9_1_1, text='     Strichbreite :  ')
            label_strichbreite_hl_3.grid(row=0, column=4, sticky='w')             

            options_strichbreite_hl_3 = tk.StringVar()
            options_strichbreite_hl_3.set(strichbreiten_werte[2]) 
            strichbreite_hl_3_option = tk.OptionMenu(frame9_1_1, options_strichbreite_hl_3, *strichbreiten_werte) 
            strichbreite_hl_3_option.grid(row=0, column=5, sticky='e')  
            
            label_intervall_multi_hl_3 = tk.Label(frame9_1_1, text='     Grundintervall mal:  ')
            label_intervall_multi_hl_3.grid(row=0, column=6, sticky='w')  
            
            options_intervall_hl_3 = tk.StringVar()
            options_intervall_hl_3.set(intervall_vielfaches[5])
            intervall_multi_hl_3_option = tk.OptionMenu(frame9_1_1, options_intervall_hl_3, *intervall_vielfaches, command=callback_hoehe_hl3)
            intervall_multi_hl_3_option.grid(row=0, column=7, sticky='e')

            label_intervall_h3 = tk.Label(frame9_1_1,text="   = " + str(int(options_intervall_hl_3.get()) * float(options_h_intervall.get())) + "  ")
            label_intervall_h3.grid(sticky='ns', row=0, column=8)              
            
            zeichne_intervall_hl_3_var = tk.BooleanVar(value=False) 
            checkbox_hl_3_cad = ttk.Checkbutton(frame9_1_1,text='  in CAD zeichnen?', variable=zeichne_intervall_hl_3_var)  
            checkbox_hl_3_cad.grid(row=0, column=9, sticky='e')      

            # ---------------------------------------------------            
                     
            frame10_1_1 = ttk.Frame(master=frame6_1)
            frame10_1_1.grid(sticky='w', row=7, column=0) 
            
            label_fabe_hl_4 = tk.Label(frame10_1_1,text='Fuer       4. Intervall einstellen: Farbe:')
            label_fabe_hl_4.grid(row=0, column=0,pady = 10, sticky='w') 
            
            hl_4_Color = 12
           
            color_button_hl_4 = tk.Button(frame10_1_1, height=1, width=10, text = hl_4_Color, bg=farben_[hl_4_Color], command=select_tab_h4)
            color_button_hl_4.grid(row=0, column=1,pady = 10, sticky='e')                       

            label_strichart_hl_4 = tk.Label(frame10_1_1, text='     Strichart :  ')
            label_strichart_hl_4.grid(row=0, column=2, sticky='w')                      
            
            options_strichart_hl_4 = tk.StringVar()
            options_strichart_hl_4.set(stricharten_werte[0]) 
            strichart_hl_4_option = tk.OptionMenu(frame10_1_1, options_strichart_hl_4, *stricharten_werte) 
            strichart_hl_4_option.grid(row=0, column=3, sticky='e')
            
            label_strichbreite_hl_4 = tk.Label(frame10_1_1, text='     Strichbreite :  ')
            label_strichbreite_hl_4.grid(row=0, column=4, sticky='w')             

            options_strichbreite_hl_4 = tk.StringVar()
            options_strichbreite_hl_4.set(strichbreiten_werte[3]) 
            strichbreite_hl_4_option = tk.OptionMenu(frame10_1_1, options_strichbreite_hl_4, *strichbreiten_werte) 
            strichbreite_hl_4_option.grid(row=0, column=5, sticky='e')  
            
            label_intervall_multi_hl_4 = tk.Label(frame10_1_1, text='     Grundintervall mal:  ')
            label_intervall_multi_hl_4.grid(row=0, column=6, sticky='w')  
            
            options_intervall_hl_4 = tk.StringVar()
            options_intervall_hl_4.set(intervall_vielfaches[10])
            intervall_multi_hl_4_option = tk.OptionMenu(frame10_1_1, options_intervall_hl_4, *intervall_vielfaches, command=callback_hoehe_hl4)
            intervall_multi_hl_4_option.grid(row=0, column=7, sticky='e')

            label_intervall_h4 = tk.Label(frame10_1_1,text="   = " + str(int(options_intervall_hl_4.get()) * float(options_h_intervall.get())) + "  ")
            label_intervall_h4.grid(sticky='ns', row=0, column=8)              
            
            zeichne_intervall_hl_4_var = tk.BooleanVar(value=False) 
            checkbox_hl_4_cad = ttk.Checkbutton(frame10_1_1,text='  in CAD zeichnen?', variable=zeichne_intervall_hl_4_var)  
            checkbox_hl_4_cad.grid(row=0, column=9, sticky='e')      

            # --------------------------------------------------- 
            frame11_1_1 = ttk.Frame(master=frame6_1)
            frame11_1_1.grid(sticky='w', row=8, column=0)             
            
            hl_berechnen_button = tk.Button(frame11_1_1, text="3. Hoehenlinien fuer Grundintervall berechnen", command=hoehenlinien_berechnen_2)
            hl_berechnen_button.grid(row=0, column=0, pady = 10, sticky='w')  
            
            label_anz_hoehenlinien = tk.Label(frame11_1_1, text=' Anzahl Hoehenlinien (segmente): 0 ')
            label_anz_hoehenlinien.grid(row=0, column=1, sticky='w') 
            
            # -----------------------------------------------------
            
            frame12_1_1 = ttk.Frame(master=frame6_1)
            frame12_1_1.grid(sticky='w', row=9, column=0)
            
            button_cad_hl = tk.Button(frame12_1_1, text='4. Hoehenlinien in CAD zeichnen', command=hoehenlinien_zeichnen)  
            button_cad_hl.grid(row=0, column=0, pady = 10, sticky='w')  
            button_cad_hl.config(state='disabled')
            
            label_dummy_hl = tk.Label(frame12_1_1, text='vorher einstellen ---->')
            label_dummy_hl.grid(row=0, column=1, sticky='w') 
            
            hl_werte = ('HL als lines zeichnen','HL als linestrings zeichnen','HL als bsplines zeichnen') 
            
            options_hl_art = tk.StringVar()
            menu_hl_werte = tk.OptionMenu(frame12_1_1, options_hl_art, *hl_werte)
            menu_hl_werte.grid(sticky='e', row=0, column=2)
            options_hl_art.set(hl_werte[1]) 
            
            label_dummy2_hl = tk.Label(frame12_1_1, text='             ')
            label_dummy2_hl.grid(row=0, column=3, sticky='w')                                                 
            
            button2_undo_mark = tk.Button(frame12_1_1, text='undo mark senden', command=undo_mark)  
            button2_undo_mark.grid(row=0, column=4, pady = 10, sticky='e')  
            button2_undo_mark.config(state='disabled') 
            
            #-------- frame Umring
            frame7_1 = ttk.Frame(master=frameUmring)
            frame7_2 = ttk.Frame(master=frameUmring)
            frame7_1.pack(pady=20)
            frame7_2.pack(pady=20) 
            
            #- define columns
            columns_5 = ('level_name_5', 'level_id_5', 'anzahl_maschen_5',
             'z_min_5', 'z_max_5', 'diff_z_5')

            tree_5 = ttk.Treeview(frame7_1, columns=columns_5, show='headings', selectmode="browse")  # nur ein DGM kann selektiert werden

            # define headings 
            tree_5.column('level_name_5',anchor='center', stretch=False, width=200, minwidth=200)
            tree_5.heading('level_name_5', text='Level Name')
            tree_5.column('level_id_5',anchor='center', stretch=False, width=120, minwidth=75)
            tree_5.heading('level_id_5', text='Level ID')
            tree_5.column('anzahl_maschen_5',anchor='center', stretch=False, width=150, minwidth=100)
            tree_5.heading('anzahl_maschen_5', text='Anz. Maschen') 
            tree_5.column('z_min_5',anchor='center', stretch=False, width=150, minwidth=50)
            tree_5.heading('z_min_5', text='z min')
            tree_5.column('z_max_5',anchor='center', stretch=False, width=150, minwidth=50)
            tree_5.heading('z_max_5', text='z max')                                                
            tree_5.column('diff_z_5',anchor='center', stretch=False, width=150, minwidth=50)
            tree_5.heading('diff_z_5', text='z differenz') 

            # ebenen eintragen
            for ebene in ebenen:
                tree_5.insert('', tk.END, values=ebene)

            tree_5.grid(row=0, column=0, sticky='nsew')
            
            # add scrollbars
            scrollbar_5 = ttk.Scrollbar(frame7_1, orient=tk.VERTICAL, command=tree_5.yview)
            tree_5.configure(yscrollcommand=scrollbar_5.set)
            scrollbar_5.grid(row=0, column=1, sticky='ns') 

            scrollbar_hor_5 = ttk.Scrollbar(frame7_1, orient=tk.HORIZONTAL, command=tree_5.xview)
            tree_5.configure(xscrollcommand=scrollbar_hor_5.set)
            scrollbar_hor_5.grid(row=1, column=0, sticky='ew')             
            
            progressbar_5 = ttk.Progressbar(frame7_1, orient=tk.HORIZONTAL, length=600, mode="determinate")
            progressbar_5.grid(row=2, column=0, sticky='e') 
            
            button_25 = tk.Button(frame7_1, text='1. Zeichnung scannen (Maschen)', command=select_5_ElementsbyType)
            button_25.grid(row=2, column=0, ipadx=30, sticky='w') 
            
            frame13_1_1 = ttk.Frame(master=frame7_1)
            frame13_1_1.grid(sticky='w', row=3, column=0) 
            
            label_umring = tk.Label(frame13_1_1,text='       2. Ebene markieren, danach einstellen -->   Hoehe_Grund, Farben, Ebenen und Optionen ')
            label_umring.grid(row=0, column=0,pady = 10, sticky='ns')
            
            label_um_un_hoehe = tk.Label(frame13_1_1, text='   Hoehe_Grund eingeben :  ')
            label_um_un_hoehe.grid(row=0, column=1, sticky='ns')                      
            
            um_un_hoehe_eingabe = tk.Entry(frame13_1_1)
            um_un_hoehe_eingabe.insert(0,"0.000")
            um_un_hoehe_eingabe.grid(row=0, column=2, sticky='e')             
            
            # ---------------------------------------------------                                                           
                     
            frame14_1_1 = ttk.Frame(master=frame7_1)
            frame14_1_1.grid(sticky='w', row=5, column=0) 
            
            label_farbe_um_ob = tk.Label(frame14_1_1,text='Fuer  Umring  oben einstellen: Farbe:')
            label_farbe_um_ob.grid(row=0, column=0,pady = 10, sticky='w') 
            
            um_ob_Color = 9
           
            color_button_um_ob = tk.Button(frame14_1_1, height=1, width=10, text=um_ob_Color, bg=farben_[um_ob_Color], command=select_tab_um_ob)
            color_button_um_ob.grid(row=0, column=1,pady = 10, sticky='e')
            
            label_um_ob_ebene = tk.Label(frame14_1_1, text=' Ebene :  ')
            label_um_ob_ebene.grid(row=0, column=3, sticky='ns')
            
            selected_level_um_ob = tk.StringVar()
            selected_level_um_ob.set(levelList[0])
            um_ob_level_option = tk.OptionMenu(frame14_1_1, selected_level_um_ob, *levelList) 
            um_ob_level_option.grid(row=0, column=4, sticky='e') 
            
            label_um_ob_dummy_1 = tk.Label(frame14_1_1, text='   ')
            label_um_ob_dummy_1.grid(row=0, column=8, sticky='ns')            
            
            zeichne_um_ob_var = tk.BooleanVar(value=False) 
            checkbox_um_ob_cad = ttk.Checkbutton(frame14_1_1,text='  in CAD zeichnen?', variable=zeichne_um_ob_var)  
            checkbox_um_ob_cad.grid(row=0, column=9, sticky='e')                        
            
            # --------------------------------------------------- 
            
            frame15_1_1 = ttk.Frame(master=frame7_1)
            frame15_1_1.grid(sticky='w', row=6, column=0) 
            
            label_farbe_um_un = tk.Label(frame15_1_1,text='Fuer  Umring unten einstellen: Farbe:')
            label_farbe_um_un.grid(row=0, column=0,pady = 10, sticky='w') 
            
            um_un_Color = 9
           
            color_button_um_un = tk.Button(frame15_1_1, height=1, width=10, text=um_un_Color, bg=farben_[um_un_Color], command=select_tab_um_un)
            color_button_um_un.grid(row=0, column=1,pady = 10, sticky='e') 
            
            label_um_un_ebene = tk.Label(frame15_1_1, text=' Ebene :  ')
            label_um_un_ebene.grid(row=0, column=3, sticky='ns')
            
            selected_level_um_un = tk.StringVar()
            selected_level_um_un.set(levelList[0])
            um_un_level_option = tk.OptionMenu(frame15_1_1, selected_level_um_un, *levelList) 
            um_un_level_option.grid(row=0, column=4, sticky='e') 
            
            label_um_un_dummy_1 = tk.Label(frame15_1_1, text='   ')
            label_um_un_dummy_1.grid(row=0, column=8, sticky='ns')            
            
            zeichne_um_un_var = tk.BooleanVar(value=False) 
            checkbox_um_un_cad = ttk.Checkbutton(frame15_1_1,text='  in CAD zeichnen?', variable=zeichne_um_un_var)  
            checkbox_um_un_cad.grid(row=0, column=9, sticky='e')                        
            
            # --------------------------------------------------- 
            
            frame16_1_1 = ttk.Frame(master=frame7_1)
            frame16_1_1.grid(sticky='w', row=7, column=0) 
            
            label_farbe_um_sei = tk.Label(frame16_1_1,text='Fuer Umring Seiten einstellen: Farbe:')
            label_farbe_um_sei.grid(row=0, column=0,pady = 10, sticky='w') 
            
            um_sei_Color = 9
           
            color_button_um_sei = tk.Button(frame16_1_1, height=1, width=10, text=um_sei_Color, bg=farben_[um_sei_Color], command=select_tab_um_sei)
            color_button_um_sei.grid(row=0, column=1,pady = 10, sticky='e')
            
            label_um_sei_ebene = tk.Label(frame16_1_1, text=' Ebene :  ')
            label_um_sei_ebene.grid(row=0, column=3, sticky='ns')
            
            selected_level_um_sei = tk.StringVar()
            selected_level_um_sei.set(levelList[0])
            um_sei_level_option = tk.OptionMenu(frame16_1_1, selected_level_um_sei, *levelList) 
            um_sei_level_option.grid(row=0, column=4, sticky='e')

            label_um_sei_dummy_1 = tk.Label(frame16_1_1, text='   ')
            label_um_sei_dummy_1.grid(row=0, column=8, sticky='ns')            
            
            zeichne_um_sei_var = tk.BooleanVar(value=False) 
            checkbox_um_sei_cad = ttk.Checkbutton(frame16_1_1,text='  in CAD zeichnen?', variable=zeichne_um_sei_var)  
            checkbox_um_sei_cad.grid(row=0, column=9, sticky='e') 
            
            label_dummy_um = tk.Label(frame16_1_1, text=' Seiten zeichnen als ->')
            label_dummy_um.grid(row=0, column=10, sticky='w') 
            
            um_sei_werte = ('Dreiecke ','Vierecke') 
            
            options_sei_art = tk.StringVar()
            menu_sei_werte = tk.OptionMenu(frame16_1_1, options_sei_art, *um_sei_werte)
            menu_sei_werte.grid(sticky='e', row=0, column=11)
            options_sei_art.set(um_sei_werte[1])             
                                                
            # --------------------------------------------------- 
            frame17_1_1 = ttk.Frame(master=frame7_1)
            frame17_1_1.grid(sticky='w', row=9, column=0)             
            
            um_berechnen_button = tk.Button(frame17_1_1, text="3. Umring berechnen", command=umring_berechnen)
            um_berechnen_button.grid(row=0, column=0, pady = 10, sticky='w')  
            
            label_anz_umring = tk.Label(frame17_1_1, text=' Anzahl Elemente Umring (polylines): 0 ')
            label_anz_umring.grid(row=0, column=1, sticky='w') 
            
            # ---------------------------------------------------
            frame18_1_1 = ttk.Frame(master=frame7_1)
            frame18_1_1.grid(sticky='w', row=10, column=0)
            
            button_cad_um = tk.Button(frame18_1_1, text='4. Umringe in CAD zeichnen', command=umring_zeichnen)  
            button_cad_um.grid(row=0, column=0, pady = 10, sticky='w')  
            button_cad_um.config(state='disabled') 
            
            label_um_dummy_3 = tk.Label(frame18_1_1, text='   ')
            label_um_dummy_3.grid(row=0, column=1, sticky='ns')            
            
            button3_undo_mark = tk.Button(frame18_1_1, text='undo mark senden', command=undo_mark)  
            button3_undo_mark.grid(row=0, column=4, pady = 10, sticky='e')  
            button3_undo_mark.config(state='disabled')
            
            # ---------------------------------------------------
            frame19_1_1 = ttk.Frame(master=frame7_1)
            frame19_1_1.grid(sticky='w', row=11, column=0)

            button1_export_dgm = tk.Button(frame19_1_1, text='DGM auf gewaehlter Ebene als TIN im Format LandXML exportieren', command=export_dgm_to_lamdxml)  
            button1_export_dgm.grid(row=0, column=0, pady = 10, sticky='w')  
        
            #--- weiter FrameHaupt
            # Level
            # levelList,levelListIDs = GetLevelList()   

            level_label = tk.Label(frameHaupt, text="DGM auf Ebene: ")
            selected_level = tk.StringVar()
            selected_level.set(levelList[0])
            level_option = tk.OptionMenu(frameHaupt, selected_level, *levelList) 
            level_label.pack(side= "left", padx=20, pady=10) 
            level_option.pack(side= "left", padx=20, pady=10)


            button_quit = tk.Button(frameHaupt,text='Programm beenden', command=root.destroy)
            button_quit.pack(side = "right", padx=20, pady=10)

            button_undo_mark = tk.Button(frameHaupt, text='undo mark senden', command=undo_mark)  
            button_undo_mark.pack(side = "right", padx=20, pady=10)
            button_undo_mark.config(state='disabled')          

            open_button = tk.Button(frameHaupt, text="Im Programm fortfahren", command=main)
            open_button.pack(side = "right", padx=20, pady=10)

            root.mainloop()
           


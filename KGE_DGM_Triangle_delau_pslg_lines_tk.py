#!/usr/bin/python3
# -*- coding: utf-8 -*-

from MSPyBentley import *
from MSPyBentleyGeom import *
from MSPyECObjects import *
from MSPyDgnPlatform import *
from MSPyDgnView import *
from MSPyMstnPlatform import *


'''
falls scipy oder triangle noch nicht installiert ist, MicroStation 2024 beenden und scipy/triangle installieren
- Eingabeaufforderung öffnen
- Verzeichnis wechseln:
    cd C:/ProgramData/Bentley/PowerPlatformPython/python
- richtiges Verzeichnis kontrollieren ! , danach eingeben:
    python -m pip install scipy
- falls gewünscht ebenfalls triangle installieren:
    python -m pip install triangle
                   https://pypi.org/project/triangle/
    documentation: https://rufat.be/triangle/
                   https://www.cs.cmu.edu/~quake/triangle.html
- Microstation neu starten

## sonstiges: https://people.sc.fsu.edu/~jburkardt/c_src/triangle_shewchuk/triangle_shewchuk.html


'''


import sys

try:
  import triangle
except:
  print("triangle konnte nicht geladen werden. Bitte installieren")
  sys.exit()

import re
import numpy as np

from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import time
from tkinter import filedialog
from tkinter import scrolledtext
from tkinter import colorchooser
from tkinter import Frame

import tkinter as tk
from tkinter import ttk


def open_file_dialog():
    file_ = filedialog.askopenfilename(title="1 von 4: Waehle Umringkanten - line segments", filetypes=[("Text files", "*.tmp"), ("All files", "*.*")])
    return file_

def open_file_dialog_2():
    file_ = filedialog.askopenfilename(title="2. von 4. Waehle Innenkanten - line segments", filetypes=[("Text files", "*.tmp"), ("All files", "*.*")])
    return file_

def open_file_dialog_3():
    file_ = filedialog.askopenfilename(title="4. von 4. Waehle Lochkoordinaten - line in holes", filetypes=[("Text files", "*.tmp"), ("All files", "*.*")])
    return file_

def daten_einlesen(datei):
    liste_local = []
    with open(datei) as f:
        for zeile in f:
            liste_local.append(zeile.strip())      #  am Ende wird \n entfernt sowie Leerzeichen am Enfang und Ende
        return liste_local

def open_files_dialog():
    files_ = filedialog.askopenfilenames(title="3. von 4. Waehle eine oder mehrere Punktdateien - coordinates", filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
    return files_

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
    file_text.insert(tk.END, 'Grundflaeche des DGM              : {0:.3f}'.format(summe_grund_fl)  + '\n')
    file_text.insert(tk.END, 'Oberflaeche  des DGM              : {0:.3f}'.format(summe_ober_fl)  + '\n')
    file_text.insert(tk.END, 'Volumen ausgeglichen, Hoehe  0.00 : {0:.3f}'.format(volumen_ueber_0)  + '\n')
    file_text.insert(tk.END, 'mittlere Hoehe                    : {0:.6f}'.format(volumen_ueber_0/summe_grund_fl)  + '\n')
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

    shape_eeh = EditElementHandle()

    # Set Shape coordinates
    points = DPoint3dArray()
        
    points.append (DPoint3d (masche[0][0] * mu_, masche[0][1] * mu_, masche[0][2] * mu_))
    points.append (DPoint3d (masche[1][0] * mu_, masche[1][1] * mu_, masche[1][2] * mu_))
    points.append (DPoint3d (masche[2][0] * mu_, masche[2][1] * mu_, masche[2][2] * mu_))
    points.append (DPoint3d (masche[0][0] * mu_, masche[0][1] * mu_, masche[0][2] * mu_))

    
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
    color_label.config(text = farben_[dgm_color], bg=farben_[dgm_color])
    farbe.config(bg=farben_[dgm_color], text = dgm_color)
    
def farb_pick(k):
    dgm_color = k
    farbe.config(bg=farben_[dgm_color], text = dgm_color)
    w_color.set(k)
    
        
def undo_mark():
    PyCadInputQueue.SendKeyin("undo mark")
    button_undo_mark.config(state='disabled')
    lift_window(root)



if __name__ == '__main__':

    global ACTIVEMODEL
    global dgnModel
    global g_1mu
    global dgm_color
    global selected_level
    global farben_
        
        
    # Get the active DGN model reference
    ACTIVEMODEL = ISessionMgr.ActiveDgnModelRef
    if ACTIVEMODEL != None:
            dgnModel = ACTIVEMODEL.GetDgnModel()
            modelInfo = dgnModel.GetModelInfo()
            g_1mu = modelInfo.GetUorPerStorage()

            dgnfile = dgnModel.GetDgnFile()
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

            frameHaupt.pack(fill='both', expand=True)
            frameFarben.pack(fill='both', expand=True)

            notebook.add(frameHaupt, text='Hauptprogramm')
            notebook.add(frameFarben, text='Farbpicker aus Colortable')
            
            frameFarben2 = Frame(master=frameFarben)
            frameFarben2.pack()

            
            text_ = "zuerst Farb-Nr fuer DGM und DGM auf Ebene einstellen \n"
            text_ = text_ + " dann  -Im Programm fortfahren- klicken \n"
            text_ = text_ + " nacheinander auswaehlen bzw. abbrechen \n"
            text_ = text_ + "1. Koordinatendatei Segmente Umring (lines) \n"
            text_ = text_ + "2. Koordinatendatei Segmente Innen (lines) \n"
            text_ = text_ + "3. eine oder mehrere Koordinaten Dateien Punkte innen \n"
            text_ = text_ + "4. Koordinatendatei Linien innerhalb Loecher"


            hinweise_label = tk.Label(frameHaupt, text=text_)
            hinweise_label.pack(padx=20, pady=20)


            file_text = scrolledtext.ScrolledText(frameHaupt, wrap=tk.WORD, height=25, width=100)
            file_text.pack(padx=20, pady=20)

            
            color_label = tk.Label(frameHaupt, text = "Farb-Nr waehlen fuer das DGM :")
            color_label.pack()

            dgmColor = 9
            w_color = tk.Scale(frameHaupt, from_=0, to=254, length=1200, tickinterval=15, orient=tk.HORIZONTAL, command=farbe_aendern)
            w_color.set(dgmColor)
            w_color.pack()
            
            color_label = tk.Label(frameHaupt, text = farben_[dgmColor], bg=farben_[dgmColor])
            color_label.pack(padx=20, pady=20)

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
            
            farbe = tk.Button(frameFarben, text=dgmColor, bg=farben_[dgmColor], height=1, width=6, padx=2, pady=2)
            farbe.pack()



            # Level
            levelList,levelListIDs = GetLevelList()   

            level_label = tk.Label(frameHaupt, text="DGM auf Ebene: ")
            selected_level = tk.StringVar()
            selected_level.set(levelList[0])
            level_option = tk.OptionMenu(frameHaupt, selected_level, *levelList) 
            level_label.pack(side= "left", padx=20, pady=20) 
            level_option.pack(side= "left", padx=20, pady=20)


            button_quit = tk.Button(frameHaupt,text='Programm beenden', command=root.destroy)
            button_quit.pack(side = "right", padx=20, pady=20)

            button_undo_mark = tk.Button(frameHaupt, text='undo mark senden', command=undo_mark)  
            button_undo_mark.pack(side = "right", padx=20, pady=20)
            button_undo_mark.config(state='disabled')          

            open_button = tk.Button(frameHaupt, text="Im Programm fortfahren", command=main)
            open_button.pack(side = "right", padx=20, pady=20)

            root.mainloop()
            #'''

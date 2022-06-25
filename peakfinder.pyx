#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Created on Fri Aug 14 20:10:06 2020
# @author: Fernando Camussi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
PeakFinder, un algortimo para detectar los picos de un blob. Basado en el
método propuesto por Wu et al. [1] para determinar la cantidad de lóbulos en
las células de pavimento.
[1] LobeFinder: A Convex Hull-Based Method for Quantitative Boundary Analyses of
Lobed Plant Cells Tzu-Ching Wu, Samuel A. Belteton, Jessica Pack,
Daniel B. Szymanski, David M. Umulis
Plant Physiology Aug 2016, 171 (4) 2331-2342; DOI: 10.1104/pp.15.00972
"""

import numpy as np
import cv2 as cv
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
import math
from libc.stdlib cimport malloc, free
from libc.math cimport sqrt, ceil


def peakfinder(image, th_dist_peak=0.02):
    """ Busca los picos en un blob

    Argumentos:
        image -- imagen binaria conteniendo el blob
        th_dist_peak -- umbral de distancia entre picos

    Retorna: puntos de la imagen correspondiente a los picos
    """

    if not (len(image.shape) == 2 and len(np.unique(image)) == 2):
        raise Exception('Imagen no binaria')

    contour,_ = cv.findContours(image, cv.RETR_TREE, cv.CHAIN_APPROX_NONE)
    contour = contour[0]

    # Separamos contour en cnt_x y cnt_y
    cnt_x = contour[:,0][:,0]
    cnt_y = contour[:,0][:,1]

    # Suavizado del contorno
    cnt_x = __smooth(cnt_x, 5)
    cnt_y = __smooth(cnt_y, 5)

    # Cierre del contorno
    cnt_x = np.append(cnt_x, cnt_x[0])
    cnt_y = np.append(cnt_y, cnt_y[0])

    # Normalización entre [-1,1]
    min_cnt_x = min(cnt_x)
    max_cnt_x = max(cnt_x)
    min_cnt_y = min(cnt_y)
    max_cnt_y = max(cnt_y)
    cnt_x = ((cnt_x-min_cnt_x) / (max_cnt_x-min_cnt_x))*2 - 1
    cnt_y = ((cnt_y-min_cnt_y) / (max_cnt_y-min_cnt_y))*2 - 1

    # Cápsula convexa del contorno
    hull = cv.convexHull(np.float32(np.vstack((cnt_x,cnt_y)).T))
    hull_x = hull[:,0][:,0]
    hull_y = hull[:,0][:,1]

    # Interpolación spline cúbica del contorno
    x = np.linspace(-1, 1, len(cnt_x))
    xx = np.linspace(-1, 1, len(cnt_x)*10)
    fx = interp1d(x, cnt_x, kind='cubic')
    cnt_x = fx(xx)
    y = np.linspace(-1, 1, len(cnt_y))
    yy = np.linspace(-1, 1, len(cnt_y)*10)
    fy = interp1d(y, cnt_y, kind='cubic')
    cnt_y = fy(yy)

    # Rotación de los puntos de la cápsula convexa para comenzar desde el
    # primer punto correspondiente a un lóbulo
    s = 0
    while (hull_x[s]-hull_x[s-1])**2 + \
     (hull_y[s]-hull_y[s-1])**2 <= th_dist_peak:
        s -= 1
    hull_x = np.roll(hull_x, -s)
    hull_y = np.roll(hull_y, -s)

    # Cierro la capsula convexa
    hull_x = np.append(hull_x, hull_x[0])
    hull_y = np.append(hull_y, hull_y[0])

    # Calcula la mínima distancia para n_PQ_pts puntos por cada segmento PQ
    # en la cápsula convexa a el contorno y se almacenan los valores en fmin
    hull_x = np.double(hull_x)
    hull_y = np.double(hull_y)
    cnt_x = np.double(cnt_x)
    cnt_y = np.double(cnt_y)
    fmin = np.array([], dtype=np.double)
    imin = np.array([], dtype=int)
    cdef long int hull_len, cnt_len
    hull_len = len(hull_x)
    cnt_len = len(cnt_x)
    cdef double *hull_x_ = <double *>malloc(hull_len * sizeof(double))
    cdef double *hull_y_ = <double *>malloc(hull_len * sizeof(double))
    cdef double *cnt_x_ = <double *>malloc(cnt_len * sizeof(double))
    cdef double *cnt_y_ = <double *>malloc(cnt_len * sizeof(double))
    cdef long int i, j, k
    cdef long int n_PQ_pts, ind
    cdef double x0, y0, d, z, m1, m2
    cdef double *hx = NULL
    cdef double *hy = NULL

    for i in range(hull_len):
        hull_x_[i] = hull_x[i]
        hull_y_[i] = hull_y[i]
    for i in range(cnt_len):
        cnt_x_[i] = cnt_x[i]
        cnt_y_[i] = cnt_y[i]
    del hull_x, hull_y, cnt_x, cnt_y

    for i in range(hull_len-1):
        # (x0,y0)=P-Q es un vector paralelo al segmento PQ donde
        # P = (hull_x(i),hull_y(i)) y Q = (hull_x(i+1),hull_y(i+1))
        x0 = hull_x_[i+1] - hull_x_[i]
        y0 = hull_y_[i+1] - hull_y_[i]
        # d es la distancia euclidiana entre los puntos P y Q
        d = sqrt(x0**2 + y0**2)
        # n_PQ_pts es la cantidad de puntos intermedios a considerar en el
        # segmento PQ
        n_PQ_pts = <long int> ceil(d/0.01)
        # hx y hy son puntos intermedios en el segmento PQ (linspace)
        hx = <double *>malloc(n_PQ_pts * sizeof(double))
        hy = <double *>malloc(n_PQ_pts * sizeof(double))
        if n_PQ_pts > 1:
            for j in range(n_PQ_pts):
                hx[j] = hull_x_[i] + j * (hull_x_[i+1] - hull_x_[i]) \
                    / (n_PQ_pts - 1)
                hy[j] = hull_y_[i] + j * (hull_y_[i+1] - hull_y_[i]) \
                    / (n_PQ_pts - 1)
            else:
                hx[j] = hull_x_[i]
                hy[j] = hull_y_[i]
        for j in range(n_PQ_pts):
            # Calculo el valor de z / x*x0+y*y0=z (ecuación de la recta
            # perpendicular a PQ que pasa por el punto intermedio (x,y))
            z = hx[j]*x0+hy[j]*y0
            # Busco el punto del contorno perteneciente a la recta / que la
            # distancia a la cápsula convexa sea mínima
            m1 = abs(z-(cnt_x_[0]*x0+cnt_y_[0]*y0)) + \
                (hx[j]-cnt_x_[0])**2 + (hy[j]-cnt_y_[0])**2
            ind = 0
            for k in range(1, cnt_len):
                m2 = abs(z-(cnt_x_[k]*x0+cnt_y_[k]*y0)) + \
                    (hx[j]-cnt_x_[k])**2 + (hy[j]-cnt_y_[k])**2
                if m2 < m1:
                    m1 = m2
                    ind = k
            # fmin: función de mínima distancia entre la cápsula convexa y el
            # contorno (sin calcular la raíz cuadrada)
            fmin = np.append(fmin, (hx[j]-cnt_x_[ind])**2 + \
                                   (hy[j]-cnt_y_[ind])**2)
            # imin: índice del contorno correspondiente a la función fmin
            imin = np.append(imin, ind)
        free(hx)
        free(hy)

    # Busco mínimos locales en fmin
    local_min,_ = find_peaks(-fmin, prominence=0.01)

    # Se cuenta cada punto en la cápsula convexa dependiendo de la distancia
    # entre puntos vecinos
    # Si la distancia entre un punto de la cápsula convexa y un punto del
    # contorno correspondiente a un mínimo local es <= a th_dist_peak,
    # entonces es removido del conjunto de mínimos locales
    hull_peaks = 0
    hull_peak_label = np.array([], int)
    hull_peak_index = np.array([], int)
    for i in range(hull_len):
        for j in local_min:
            if (hull_x_[i]-cnt_x_[imin[j]])**2 + \
             (hull_y_[i]-cnt_y_[imin[j]])**2 <= th_dist_peak:
                local_min = np.delete(local_min,
                                      np.where(local_min == j)[0][0])
        if (hull_x_[i]-hull_x_[(hull_len+i-1)%hull_len])**2 + \
         (hull_y_[i]-hull_y_[(hull_len+i-1)%hull_len])**2 > th_dist_peak:
            hull_peaks += 1
        hull_peak_label = np.append(hull_peak_label, hull_peaks)

    for i2 in np.unique(hull_peak_label):
        ind2 = np.where(hull_peak_label == i2)[0]
        hull_peak_index = np.append(hull_peak_index, ind2[len(ind2)//2])

    # La suma de la cantidad del conjunto de mínimo locales resultante y de la
    # cantidad de lóbulos contados a partir de la cápsula convexa es la
    # cantidad de lóbulos total
    point_peaks = np.empty([0,2], int)
    for i in hull_peak_index[:-1]:
        x = int(((hull_x_[i]+1)/2)*(max_cnt_x-min_cnt_x)+min_cnt_x)
        y = int(((hull_y_[i]+1)/2)*(max_cnt_y-min_cnt_y)+min_cnt_y)
        point_peaks = np.append(point_peaks, [[x,y]], axis=0)
    for i in imin[local_min]:
        x = int(((cnt_x_[i]+1)/2)*(max_cnt_x-min_cnt_x)+min_cnt_x)
        y = int(((cnt_y_[i]+1)/2)*(max_cnt_y-min_cnt_y)+min_cnt_y)
        point_peaks = np.append(point_peaks, [[x,y]], axis=0)

    free(hull_x_)
    free(hull_y_)
    free(cnt_x_)
    free(cnt_y_)

    return point_peaks


def __smooth(data, w_size):
    window = np.ones(w_size)/w_size
    return np.convolve(data, window, mode='valid')

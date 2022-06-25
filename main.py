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

import cv2 as cv
from peakfinder import peakfinder

image = cv.imread('ejemplo.png', cv.IMREAD_COLOR)
points = peakfinder(cv.cvtColor(image, cv.COLOR_BGR2GRAY))

for p in points:
    cv.circle(image, p, 2, (0,0,255), -1)

cv.namedWindow('PeakFinder', flags=cv.WINDOW_GUI_NORMAL)
cv.imshow('PeakFinder', image)
cv.waitKey(0)
cv.destroyAllWindows()

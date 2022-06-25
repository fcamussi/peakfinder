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

from setuptools import setup
from Cython.Build import cythonize


setup(
    ext_modules=cythonize(["peakfinder.pyx"],
                          language_level = "3")
)

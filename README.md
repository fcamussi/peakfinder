# PeakFinder

Un algortimo para detectar los picos de un blob. Basado en el
método propuesto por Wu et al. [1] para determinar la cantidad de lóbulos en
las células de pavimento.

[1] LobeFinder: A Convex Hull-Based Method for Quantitative Boundary Analyses of
Lobed Plant Cells Tzu-Ching Wu, Samuel A. Belteton, Jessica Pack,
Daniel B. Szymanski, David M. Umulis
Plant Physiology Aug 2016, 171 (4) 2331-2342; DOI: 10.1104/pp.15.00972

![screenshot1](https://user-images.githubusercontent.com/75378876/176070307-2e598b71-c12f-466f-9d01-b205adb69d4e.png)

## Requisitos

* Cython >= 0.29.25
* OpenCV >= 4.5.2
* Scipy  >= 1.7.3

## Documentación

```
peakfinder(image, th_dist_peak=0.02)
```

Busca los picos en un blob

Argumentos:
* image: imagen binaria conteniendo el blob
* th_dist_peak: umbral de distancia entre picos
* valor de retorno: puntos de la imagen correspondiente a los picos

### Compilación

```
python setup.py build_ext --inplace
```

### Ejemplo
```
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
```

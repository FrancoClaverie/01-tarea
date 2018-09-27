#Problema 2

#Primero que nada, se importan los modulos a utilizar

import numpy as np 
import matplotlib.pyplot as pyplot
import matplotlib as mp 

#Parte 1

#Se debe graficar los datos que se encuentran en el archivo dado, cuidando las unidades de medida
#Para esto se lee el archivo y se transforma en un array con numpy.loadtxt('nombre del archivo')

archivo = np.loadtxt('firas_monopole_spec_v1.txt')

#Del archivo, se tiene que cada columna significa:

#Columna 1: Frecuencia [Hz]
#Columna 2: Espectro del monopolo medido por el instrumento FIRAS del satelite COBE [MJy/sr], donde 1 [MJy] = 10^−20 [Wm^−2Hz^−1] 
#Columna 3: Espectro residual del monopolo [kJy/sr]
#Columna 4: Incertidumbre del espectro medido [kJy/sr]
#Columna 5: Modelo del espectro de la galaxia en los polos galácticos [kJy/sr]

#De estas columnas nos interesan las columnas 1 (frecuencia), 2 (esectro medido) y 4 (incertidumbre/error) para la
#realizacion del grafico.


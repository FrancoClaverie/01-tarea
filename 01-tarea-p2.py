#Problema 2

#Primero que nada, se importan los modulos a utilizar

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mp 

#Parte 1

#Se debe graficar los datos que se encuentran en el archivo dado, cuidando las unidades de medida
#Para esto se lee el archivo y se transforma en un array con numpy.loadtxt('nombre del archivo')

archivo = np.loadtxt('firas_monopole_spec_v1.txt')

#Del archivo, se tiene que cada columna significa:

#Columna 1: Frecuencia [cm^-1]
#Columna 2: Espectro del monopolo medido por el instrumento FIRAS del satelite COBE [MJy/sr], donde 1 [MJy] = 10^−20 [Wm^−2Hz^−1] 
#Columna 3: Espectro residual del monopolo [kJy/sr]
#Columna 4: Incertidumbre del espectro medido [kJy/sr]
#Columna 5: Modelo del espectro de la galaxia en los polos galácticos [kJy/sr]

#De estas columnas nos interesan las columnas 1 (frecuencia), 2 (esectro medido) y 4 (incertidumbre/error) para la
#realizacion del grafico.

frecuencia = archivo[:,0]  #primera columna, es de tipo np.ndarray
espectro = archivo[:,1]    #segunda columna
error = archivo[:,3]/1000  #cuarta columna

#Se tiene que el espectro medido esta en [MJy/sr] y el error en [kJy/sr], por lo que se divide el error por mil

#Como el error es muy pequeño por la precision del FIRAS, se multiplica por 400 (factor sugerido como ejemplo en el enunciado)

error_400 = 400*error

#Grafico:
"""
plt.plot(frecuencia, espectro)
plt.title('Espectro de monopolo medido por FIRAS')
plt.xlabel('Frecuencia [$cm^{-1}$]')
plt.ylabel('Espectro del monopolo [MJy/sr]')
plt.errorbar(frecuencia, espectro, xerr=0, yerr=error_400)

plt.savefig('espectro_monopolo_1.png')
plt.show()
"""

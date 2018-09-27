#Problema 2

#Primero que nada, se importan los modulos a utilizar

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mp 
import math as m
from astropy import constants as astro   #modulo con las constantes que se usan en la funcion de Planck

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

frecuencia = archivo[:,0]*100*astro.c.value          #primera columna, es de tipo np.ndarray
espectro = archivo[:,1]                              #segunda columna
error = archivo[:,3]/1000                            #cuarta columna

#Se tiene que el espectro medido esta en [MJy/sr] y el error en [kJy/sr], por lo que se divide el error por mil, tambien la frecuencia
#esta en [cm^-1], por lo que se ultiplica por 100 y por c para quedar en [1/s = Hz]

#Como el error es muy pequeño por la precision del FIRAS, se multiplica por 400 (factor sugerido como ejemplo en el enunciado)

error_400 = 400*error

#Grafico:
"""
plt.plot(frecuencia, espectro)
plt.title('Espectro de monopolo medido por FIRAS')
plt.xlabel('Frecuencia [Hz]')
plt.ylabel('Espectro del monopolo [MJy/sr]')
plt.errorbar(frecuencia, espectro, xerr=0, yerr=error_400)

plt.savefig('espectro_monopolo_1.png')
plt.show()
"""

#Parte 2

#P = *integral_{0}^{infinito}(x^3/(exp(x)-1))dx = K(T)*I

integral_analitica = (np.pi**4)/15 

#Para resolver la integral, se utiliza el cambio de variable sugerido y = arctan(x) y asi poder cambiar lo limites de integracion

#y = arctan(x); x = tan(y)
#dy = (1/(1+x^2))dx; dx = (1+x^2)dy

#x = 0 entonces y = arctan(0) = 0
#x = infinito entonces lim_{x = infinito}arctan(x) = pi/2

#P = K(T)*integral_{0}^{pi/2}(tan^3(y))(1+tan^2(y))/(exp(tan(y)) - 1)dy = K(T)*integral_{0}^{pi/2}f(y)dy

#Se tiene que f(y) tiene numerador y denominador infinitos en y = 0 y solo denominador infinito en en y = pi/2, por lo que
#se usara el metodo del valor medio en ambos extremos. Se resolvera con el metodo de Simpson 1/3 para lo demas

#Definimos la función dentro de la integral

def f_integral(x):
    out = (np.tan(x)**3 + np.tan(x)**5)/(np.exp(np.tan(x)) - 1)
    return out

#Para calcular con dicho metodo, se debe de escoger un h tal que la integral se acerque a su valor analitico, por lo que
#se refinarra el valor de la integral al ir aumentando al doble el numero de terminos que se dividira (b - a) = (pi/2 - 0) = pi/2

#Sabemos que la integral tiene un valor analitico igual a (np.pi**4)/15, por lo que se generara un particion cada vez 
# la mitad de tamaño que la anterior hasta estar cerca de ese valor. Se debe definir entonces una tolerancia

x_0 = 0
x_n = np.pi/2


def valor_integral(a,b, N=10 , toler=1e-10):
    integral = 0
    while np.abs(integral - integral_analitica) > toler:
        integral = 0
        h = (b-a)/N
        for i in range(N+1):
            x = a + h*i    #valor a evaluar la funcion f_integral
            if x == a:
                integral += h*f_integral(x + h/2)   #valor medio en el inicio 
            elif x == b:
                integral += h*f_integral(x - h/2)   #valor medio en el final 
            elif i % 2 == 0:   #si el numero de termino es par
                integral += 2*f_integral(x)
            else: #si el numero de termino es impar y no se esta evaluando en los extremos
                integral += 4*f_integral(x)
        integral += f_integral(a+h) + f_integral(b-h)  #se le suman los terminos inicial y final del metodo de Simpson 1/3
        integral *= h/3        
        N *= 2  #se duplica el numero de divisiones, disminuye a la mitad h
    print('N =', N)
    return integral









#Problema 2

#Primero que nada, se importan los modulos a utilizar

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mp 
import math as m
from astropy import constants as astro   #modulo con las constantes que se usan en la funcion de Planck


#Funcion de Planck: B = ((2*h*v**3)/c**2)/(np.exp((h*v)/(k_B*T) - 1)

#c: Rapidez de la luz en el vacio
c = astro.c.value
#h: Constante de Planck
h = astro.h.value
#v: Frecuencia
#k_B: Constante de Boltzmann
k_B = astro.k_B.value
#T: Temperatura

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

frecuencia = archivo[:,0]*100*c          #primera columna, es de tipo np.ndarray
espectro = archivo[:,1]                              #segunda columna
error = archivo[:,3]/1000                            #cuarta columna

#Se tiene que el espectro medido esta en [MJy/sr] y el error en [kJy/sr], por lo que se divide el error por mil. La frecuencia
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

#P = (2*h/c)*((k_B*T/h)**4)*integral_{0}^{infinito}(x^3/(exp(x)-1))dx = K(T)*I

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


def valor_integral(a,b, N=10 , toler=1e-7):
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
    return integral



#Parte 3

#De la parte 2: P = (2*h/c**2)*((k_B*T/h)**4)*valor_integral, donde integral ya fue calculado

#Se quiere calcular T, para lo cual se calcula la integral de lo graficado en la parte 1

#Para calcular la integral del espectro (del archivo) respecto a la frecuencia se usara el metodo del trapecio. Hay que
#tener en cuenta que para poder igualarlo a P, se utilizaran dimensiones del sistema internacional, asi que hay que
#hacer las conversiones pertinentes. Ya se habia convertido la frecuencia a Hz, pero falta pasar el espectro de MJy/sr a Wm^−2Hz^−1/sr

#Recordar: [MJy] = 10^−20 [Wm^−2Hz^−1] , donde el espectro es la distribucion de energia por unidad de frecuencia

espectro_2 = espectro*10**(-20)

#Teniendo ya bien las dimensiones, se integra con el metodo del trapecio. Para usarlo no es necesario tener un h constante,
#la distancia entre cada termino a evaluar puede ser variable y adaptarse a los puntos si es que ya se dispone de ellos.
#h vendra dado por la distancia entre 2 frecuenciaa consecutivas, siendo el largo de las bases de los trapecios el valor del
#espectro en dichas frecuencias

integral_espectro = 0

for i in range(len(frecuencia)-1):
    integral_espectro += (espectro_2[i]+espectro_2[i+1])*(frecuencia[i+1]-frecuencia[i])/2

#Igualando este resultado con P: integral_espectro = (2*h/c**2)*((k_B*T/h)**4)*valor_integral

#Despejando T: T = (h/k_B)*((integral_espectro*c**2)/(valor_integral*2*h))**(1/4)

T = (h/k_B)*((integral_espectro*c**2)/(valor_integral(x_0,x_n)*2*h))**(1/4)

"""
print('T =', T, '[K]')
"""

#Da T = 2.684 K, siendo que el valor real es 2.725 K, esto ocurre por arrastre de errores de aproximacion


#Parte 4

#Para graficar lo echo en la parte 1 junto a los resultados obtenidos, es conveniente definir la funcion de Planck dependiente
#de la frecuencia y de la temperatura, escrita al principio.

def planck(v,T):
    output = ((2*h*v**(3))/c**(2))/(np.exp((h*v)/(k_B*T)) - 1)
    return output


#Grafico:
"""
plt.plot(frecuencia, espectro_2, label='Medicion FIRAS (S.I.)')              #grafico parte 1 pero en S.I.
plt.plot(frecuencia, planck(frecuencia,T), label='T calculada')       #T es la temperatura calculada en la parte 3
plt.plot(frecuencia, planck(frecuencia,2.725), label='2.725 K')       #2.725 K es la temperatura a la que se queria llegar con T

plt.title('Espectro segun FIRAS y segun Planck')
plt.xlabel('Frecuencia [Hz]')
plt.ylabel('Espectro del monopolo [$W m^{−2} Hz^{−1}/sr$]')

plt.legend()
plt.savefig('espectro_monopolo_2.png')
plt.show()
"""

#Como de grafico de las mediciones de FIRAS son muy cercanos a los de la funcion de Planck en T=2.725 K, es conveniente realizar
#un grafico de la diferencia de estos.

diferencia_espectro = [m.fabs(espectro_2[i] - planck(frecuencia[i],2.725)) for i in range(len(frecuencia))]

"""
plt.title("Diferencia del espectro de FIRAS y el de Planck")

plt.plot(frecuencia, diferencia_espectro)

plt.yscale('log')

plt.xlabel('Frecuencia [Hz]', fontsize=15)
plt.ylabel('Medición FIRAS - Planck en T = 2.725 [K]', fontsize=10)

plt.savefig('diferencia_espectro.png')
plt.show()
"""

#Parte 5

import scipy.integrate as sc

#sc.trapz(y,x) realiza el algortimo del metodo del trapecio sobre la funcion y en los puntos dados por x, 
#por lo que se usara en la integral de la parte 3, es decir, en los datos del archivo de FIRAS.

integral_parte3 = sc.trapz(espectro_2,frecuencia)
print('integral_parte3 =',integral_parte3)

#sc.quad(func,a,b) calcula la integral de la funcion func desde a hasta b. Este se usara para la integral de P (parte 2).
#Para esto, se define la funcion de la integral sin el cambio de variable.

def func(x):
    y = (x**3)/(np.exp(x) - 1)
    return y

integral_parte2 = sc.quad(func, x_0, np.inf)
print('integral_parte2 =', integral_parte2)














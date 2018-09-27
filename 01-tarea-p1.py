##Problema 1

#Mi RUT: 19.687.559-K

#Se pide comparar la precision de 2 metodos de derivacion discreta, uno de error de orden 4 (O(h^4))
#descrito en el enunciado y otro de orden 1 (O(h)) ensenado en clases. Esto para f(x) = -cos(x); f'(x) = sin(x), para
#x = 1.XXX, donde XXX son los 3 ultimos numeros de mi RUT antes del digito verificador (XXX = 559). Primero se 
#realizara la comparacion usando numeros float32 y luego float64 (este computador no puede usar float128), 
#comparando tambien el resultados entre ellos.

#f'(x) =(-f(x+2*h)+8*f(x+h)-8*f(x-h)+f(x-2*h))/(12*h) + O(h^4)

#f'(x) = (f(x+h) - f(x))/h + O(h)

#Antes que nada, se deben importar los modulos a utilizar

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import math as m 

#Definimos ahora las derivadas, para float32 y float64 respectivamente

#Float32:

def derivada_orden_1_float32(x,h):
    """ x en radianes
    """
    y = np.float32(x)
    e = np.float32(h)
    df_dy_1 =  (np.float32(-np.cos(y+e)) - (np.float32(-np.cos(y))))/h
    return np.float32(df_dy_1)

def derivada_orden_4_float32(x,h):
    """ x en radianes
    """
    y = np.float32(x)
    e = np.float32(h)
    df_dy_4 = (-(np.float32(-np.cos(y+2*e)))+8*(np.float32(-np.cos(y+e)))-8*(np.float32(-np.cos(y-e)))+(np.float32(-np.cos(y-2*e))))/(12*e)
    return np.float32(df_dy_4)

#Float64 (default):

def derivada_orden_1(x,h):
    """ x en radianes
    """
    df_dx_1 =  (-np.cos(x+h) - (-np.cos(x)))/h
    return df_dx_1

def derivada_orden_4(x,h):
    """ x en radianes
    """
    df_dx_4 = (-(-np.cos(x+2*h))+8*(-np.cos(x+h))-8*(-np.cos(x-h))+(-np.cos(x-2*h)))/(12*h)
    return df_dx_4


k = 1.559  #x con el que se probara la derivada

seno_float32 = np.float32(m.sin(k)) #valor de la derivada con float32
seno = m.sin(k)  #valor real al que se quiere llegar con la aproximacion con foat 64

#Se debe escoger un rango apropiado de h pequeños, tal
#que la aproximacion sea cercana al valor real. Para esto, se probara con potencias
#negativas de 10.

H = np.logspace(-0, -18, 15, base=10.)   #vector de h

#Primero se haran los graficos de las derivadas para cada h en H, y se comparara con el valor de sin(1.559)

#Grafico orden 1

"""
plt.title("Valor de la derivada (error de orden 1)")

plt.plot(H, derivada_orden_1_float32(k,H), label="float32")
plt.plot(H, derivada_orden_1(k,H), label="float64")

plt.xscale('log')

plt.axhline(seno, color='r', label='valor real')

plt.xlabel('h', fontsize=20)
plt.ylabel('$\sin(1.559)_{aprox}$', fontsize=10)

plt.legend()
plt.savefig('derivada_orden_1.png')
plt.show()
"""

#Grafico orden 4

"""
plt.title("Valor de la derivada (error de orden 4)")

plt.plot(H, derivada_orden_4_float32(k,H), label="float32")
plt.plot(H, derivada_orden_4(k,H), label="float64")

plt.xscale('log')

plt.axhline(seno, color='r', label='valor real')

plt.xlabel('h', fontsize=20)
plt.ylabel('$\sin(1.559)_{aprox}$', fontsize=10)

plt.legend()
plt.savefig('derivada_orden_4.png')
plt.show()
"""

#Otros graficos que pueden ser utiles son los de la diferencia del valor real de la derivada
#con el aproximado (valor absoluto) para ambos floats.

diferencia_orden_1_float32 = [np.float32(m.fabs(derivada_orden_1_float32(k,h) - seno_float32)) for h in H]
diferencia_orden_1 = [m.fabs(derivada_orden_1(k,h) - seno) for h in H]

diferencia_orden_4_float32 = [np.float32(m.fabs(derivada_orden_4_float32(k,h) - seno_float32)) for h in H]
diferencia_orden_4 = [m.fabs(derivada_orden_4(k,h) - seno) for h in H]

#Grafico orden 1

"""
plt.title("Diferencia con la derivada (error de orden 1)")

plt.plot(H, diferencia_orden_1_float32, label="float32")
plt.plot(H, diferencia_orden_1, label="float64")

plt.xscale('log')
plt.yscale('log')

plt.xlabel('h', fontsize=20)
plt.ylabel('$|\sin_(1.559){aprox} - \sin(1.559)_{real}|$', fontsize=10)

plt.legend()
plt.savefig('diferencia_orden_1.png')
plt.show()
"""

#Grafico orden 4


plt.title("Diferencia con la derivada (error de orden 4)")

plt.plot(H, diferencia_orden_4_float32, label="float32")
plt.plot(H, diferencia_orden_4, label="float64")

plt.xscale('log')
plt.yscale('log')

plt.xlabel('h', fontsize=20)
plt.ylabel('$|\sin_(1.559){aprox} - \sin(1.559)_{real}|$', fontsize=10)

plt.legend()
plt.savefig('diferencia_orden_4.png')
plt.show()

##Problema 1

#Mi RUT: 19.687.559-K

#Se pide comparar la precision de 2 metodos de
#derivacion discreta, uno de error de orden 4 (O(h^4))
#descrito en el enunciado y otro de orden 1 (O(h)) ensenado
#en clases. Esto para f(x) = -cos(x); f'(x) = sin(x), para
#x = 1.XXX, donde XXX son los 3 ultimos numeros de mi RUT
#antes del digito verificador (XXX = 559). Primero se 
#realizara la comparacion usando numeros float32 y luego
#float64 (este computador no puede usar float128)

#f'(x) =(-f(x+2*h)+8*f(x+h)-8*f(x-h)+f(x-2*h))/(12*h) + O(h^4)

#f'(x) = (f(x+h) - f(x))/h + O(h)

#Antes que nada, se deben importar los modulos a utilizar
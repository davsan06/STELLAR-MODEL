
# coding: utf-8

# In[137]:


import numpy as np
import scipy as sp
import pandas as pd

# Datos de la estrella (DATOS DEL MODELO DE PRUEBA)
X = 0.85
Y = 0.10
Z = 1 - X - Y
Mtot = 5.1

#Peso molecular medio
mu = 1/(2*X + (3/4)*Y + (1/2)*Z)

#Tolerancia en nuestras comparaciones
Erelmax = 0.0001

#Valores iniciales (DATOS DEL MODELO DE PRUEBA)// Estos son los valores que debemos optimizar con el fin de minimizar el error 
# relativo en el empalme de las integraciones
#Rtot = 10.980000000000002
#Ltot = 22.385
Tcentral = 1.5

Rini = 0.9*Rtot

#Paso de integracion
h = -Rini/100
#Array del radio
radio = np.arange(Rini,0,h)
radio = np.append(radio,0)
#Constantes que aparecen en el calculo
A1 = 1.9022*mu*Mtot
A2 = 10.645*(Mtot/(mu*Z*(1+X)*Ltot))**(1/2)
Cp = 8.084*mu
Ct = 0.01679*Z*(1+X)*mu**2
Cm = 0.01523*mu
Ct_prima = 3.324*mu
A3 = 0.4*8.08435*mu*Mtot

#Arrays donde iremos guardando los resultados de cada capa iesima
temperatura = np.zeros_like(radio)
presion = np.zeros_like(radio)
masa = np.zeros_like(radio)
luminosidad = np.zeros_like(radio)
gradiente_presion = np.zeros_like(radio)
gradiente_temperatura = np.zeros_like(radio)
gradiente_masa = np.zeros_like(radio)
gradiente_luminosidad = np.zeros_like(radio)
epsilon_ciclo_PP = np.zeros_like(radio)
epsilon_ciclo_CN = np.zeros_like(radio)
presion_estimada = np.zeros_like(radio)
delta1_presion = np.zeros_like(radio)
delta2_presion = np.zeros_like(radio)
temperatura_estimada = np.zeros_like(radio)
delta1_temperatura = np.zeros_like(radio)
presion_calculada = np.zeros_like(radio)
temperatura_calculada = np.zeros_like(radio) 
delta1_temperatura = np.zeros_like(radio)
masa_calculada = np.zeros_like(radio)
delta1_masa = np.zeros_like(radio)
delta1_luminosidad = np.zeros_like(radio)
delta2_luminosidad = np.zeros_like(radio)
luminosidad_calculada = np.zeros_like(radio)
parametro_transporte = np.zeros_like(radio)


# In[138]:


#--------------------------------------- GENERACION DE ENERGIA ---------------------------------------------------------


# In[139]:


def ciclo_PP(temperatura):
    
    if ( temperatura[i] > 0.4 and temperatura[i] < 0.6 ):
        
        nuPP = 6
        epsilon1PP = 10**(-6.84)
        
    elif ( temperatura[i] >= 0.6 and temperatura[i] < 0.95 ):
        
        nuPP = 5
        epsilon1PP = 10**(-6.04)
        
    elif ( temperatura[i] >= 0.95 and temperatura[i] < 1.2 ):
        
        nuPP = 4.5
        epsilon1PP = 10**(-5.56)
        
    elif ( temperatura[i] >= 1.2 and temperatura[i] < 1.65 ):
        
        nuPP = 4
        epsilon1PP = 10**(-5.02)
        
    elif ( temperatura[i] >= 1.65 and temperatura[i] < 2.4 ):
        
        nuPP = 3.5
        epsilon1PP = 10**(-4.40)
        
    else :
        
        nuPP = 1
        epsilon1PP = 0
    
    epsilon_ciclo_PP[i] = epsilon1PP*X*X*(10*temperatura[i])**nuPP
    return nuPP, epsilon1PP


# In[140]:


def ciclo_CN(temperatura):
    
    if ( temperatura[i] > 1.2 and temperatura[i] < 1.6 ):
        
        nuCN = 20
        epsilon1CN = 10**(-22.2)
    
    elif ( temperatura[i] >= 1.6 and temperatura[i] < 2.25 ):
        
        nuCN = 18
        epsilon1CN = 10**(-19.8)
        
    elif ( temperatura[i] >= 2.25 and temperatura[i] < 2.75 ):
        
        nuCN = 16
        epsilon1CN = 10**(-17.1)
        
    elif ( temperatura[i] >= 2.75 and temperatura[i] < 3.6 ):
        
        nuCN = 15
        epsilon1CN = 10**(-15.6)
        
    elif ( temperatura[i] >= 3.6 and temperatura[i] < 5 ):
        
        nuCN = 13
        epsilon1CN = 10**(-12.5)
        
    else :
        
        nuCN = 1
        epsilon1CN = 0
        
    epsilon_ciclo_CN[i] = epsilon1CN*X*(Z/3)*(10*temperatura[i])**nuCN
    return nuCN, epsilon1CN


# In[141]:


def generacion_energia(temperatura):
    
    nuPP, epsilon1PP = ciclo_PP(temperatura)
    nuCN, epsilon1CN = ciclo_CN(temperatura)
    
    if (epsilon_ciclo_PP[i] > epsilon_ciclo_CN[i]):
        
        Cl = 0.01845*epsilon1PP*X*X*(10**nuPP)*(mu**2)
        ciclo = "PP"
        nu = nuPP
        
    elif (epsilon_ciclo_PP[i] <= epsilon_ciclo_CN[i]):
        
        Cl = 0.01845*epsilon1CN*X*(Z/3)*(10**nuCN)*(mu**2)
        ciclo = "CN"
        nu = nuCN
        
    return nu, Cl, ciclo   


# In[142]:


#----------------------------------------ENVOLTURA RADIATIVA -----------------------------------------------------------


# In[143]:


# En esta parte del codigo consideramos masa y luminosidad constante. Funcion valida para valores iniciales y A.1.1.
def paso1(radio):
    
    temperatura[i] = A1*(1/radio[i] - 1/Rtot)
    presion[i] = A2*temperatura[i]**4.25
    masa[i] = Mtot
    luminosidad[i] = Ltot
    
    nu, Cl, ciclo = generacion_energia(temperatura)
    
    gradiente_presion[i] = -Cp*(presion[i]/temperatura[i])*(masa[i]/radio[i]**2)
    gradiente_temperatura[i] = -Ct*(presion[i]**2/temperatura[i]**8.5)*(luminosidad[i]/radio[i]**2)
    gradiente_masa[i] = Cm*(presion[i]/temperatura[i])*radio[i]**2
    gradiente_luminosidad[i] = Cl*presion[i]**2*temperatura[i]**(nu-2)*radio[i]**2


# In[144]:


def paso2_a11(radio):
    
    delta1_presion[i-1] = h*gradiente_presion[i-1] - h*gradiente_presion[i-2]
    delta2_presion[i-1] = h*gradiente_presion[i-1] - 2*h*gradiente_presion[i-2] + h*gradiente_presion[i-3]
    presion_estimada[i] = presion[i-1] + h*gradiente_presion[i-1] + (1/2)*delta1_presion[i-1] + (5/12)*delta2_presion[i-1]
    
    delta1_temperatura[i-1] = h*gradiente_temperatura[i-1] - h*gradiente_temperatura[i-2]
    temperatura_estimada[i] = temperatura[i-1] + h*gradiente_temperatura[i-1] + (1/2)*delta1_temperatura[i-1]


# In[145]:


def paso4_a11(radio):
    
    gradiente_presion[i] = -Cp*(presion_estimada[i]/temperatura_estimada[i])*(masa[i]/radio[i]**2)
    delta1_presion[i] = h*gradiente_presion[i] - h*gradiente_presion[i-1]
    presion_calculada[i] = presion[i-1] + h*gradiente_presion[i] - (1/2)*delta1_presion[i]


# In[146]:


def paso5_a11(radio):
    
    while( abs(presion_calculada[i] - presion_estimada[i])/presion_calculada[i] > Erelmax ):
        
        presion_estimada[i] = presion_calculada[i]
        paso4_a11(radio)
        
    presion[i] = presion_calculada[i]


# In[147]:


def paso7_a11(radio):
    
    gradiente_temperatura[i] = -Ct*(presion_calculada[i]**2/temperatura_estimada[i]**8.5)*(Ltot/radio[i]**2)
    delta1_temperatura[i] = h*gradiente_temperatura[i] - h*gradiente_temperatura[i-1]
    temperatura_calculada[i] = temperatura[i-1] + h*gradiente_temperatura[i] - (1/2)*delta1_temperatura[i]


# In[148]:


def paso8_a11(radio):
    
    while( abs(temperatura_calculada[i] - temperatura_estimada[i])/temperatura_calculada[i] > Erelmax ):
        
        temperatura_estimada[i] = temperatura_calculada[i]
        
        paso4_a11(radio)
        paso5_a11(radio)
        paso7_a11(radio)
        
    temperatura[i] = temperatura_calculada[i]


# In[149]:


def paso3_a11(radio):
    
    gradiente_masa[i] = Cm*(presion_calculada[i]/temperatura_calculada[i])*radio[i]**2
    delta1_masa[i] = h*gradiente_masa[i] - h*gradiente_masa[i-1]
    masa_calculada[i] = Mtot + h*gradiente_masa[i] - (1/2)*delta1_masa[i]
    
    masa[i] = masa_calculada[i]


# In[150]:


def paso1_a12(radio):
    
    temperatura[i] = A1*(1/radio[i] - 1/Rtot)
    presion[i] = A2*temperatura[i]**4.25
    masa[i] = masa_calculada[i]
    luminosidad[i] = Ltot
    
    #nu, Cl, ciclo = generacion_energia(temperatura)
    
    gradiente_presion[i] = -Cp*(presion[i]/temperatura[i])*(masa[i]/radio[i]**2)
    gradiente_temperatura[i] = -Ct*(presion[i]**2/temperatura[i]**8.5)*(luminosidad[i]/radio[i]**2)
    gradiente_masa[i] = Cm*(presion[i]/temperatura[i])*radio[i]**2
    #gradiente_luminosidad[i] = Cl*presion[i]**2*temperatura[i]**(nu-2)*radio[i]**2


# In[151]:


def paso3_a12(radio):
    
    gradiente_masa[i] = Cm*(presion_estimada[i]/temperatura_estimada[i])*radio[i]**2
    delta1_masa[i] = h*gradiente_masa[i] - h*gradiente_masa[i-1]
    masa_calculada[i] = masa[i-1] + h*gradiente_masa[i] - (1/2)*delta1_masa[i]
    
    masa[i] = masa_calculada[i] 


# In[152]:


def paso4_a12(radio):
    
    gradiente_presion[i] = -Cp*(presion_estimada[i]/temperatura_estimada[i])*(masa_calculada[i]/radio[i]**2)
    delta1_presion[i] = h*gradiente_presion[i] - h*gradiente_presion[i-1]
    presion_calculada[i] = presion[i-1] + h*gradiente_presion[i] - (1/2)*delta1_presion[i]


# In[153]:


def paso5_a12(radio):
    
    while( abs(presion_calculada[i] - presion_estimada[i])/presion_calculada[i] > Erelmax ):
        
        presion_estimada[i] = presion_calculada[i]
        paso3_a12(radio)
        paso4_a12(radio)
        
    presion[i] = presion_calculada[i]


# In[154]:


def paso8_a12(radio):
    
    while( abs(temperatura_calculada[i] - temperatura_estimada[i])/temperatura_calculada[i] > Erelmax ):
        
        temperatura_estimada[i] = temperatura_calculada[i]
        
        paso3_a12(radio)
        paso4_a12(radio)
        paso5_a12(radio)
        paso7_a11(radio)
            
    temperatura[i] = temperatura_calculada[i] 


# In[155]:


def paso6_a12(radio):

    nu, Cl, ciclo = generacion_energia(temperatura)
    
    gradiente_luminosidad[i] = Cl*(presion_calculada[i]**2)*(temperatura_calculada[i])**(nu-2)*radio[i]**2
    delta1_luminosidad[i] = h*gradiente_luminosidad[i] - h*gradiente_luminosidad[i-1]
    delta2_luminosidad[i] = h*gradiente_luminosidad[i] - 2*h*gradiente_luminosidad[i-1] + h*gradiente_luminosidad[i-2]
    luminosidad_calculada[i] = Ltot + h*gradiente_luminosidad[i] - (1/2)*delta1_luminosidad[i] - (1/12)*delta2_luminosidad[i]
    
    luminosidad[i] = luminosidad_calculada[i]


# Funciones para el algoritmo A.1.3.

# In[156]:


def paso9_a13(radio):
    
    parametro_transporte[i] = (temperatura_calculada[i]/presion_calculada[i])*(gradiente_presion[i]/gradiente_temperatura[i])


# In[157]:


def paso1_a13(radio):
    
    temperatura[i] = A1*(1/radio[i] - 1/Rtot)
    presion[i] = A2*temperatura[i]**4.25
    masa[i] = masa_calculada[i]
    luminosidad[i] = luminosidad_calculada[i]
    
    nu, Cl, ciclo = generacion_energia(temperatura)
    
    gradiente_presion[i] = -Cp*(presion[i]/temperatura[i])*(masa[i]/radio[i]**2)
    gradiente_temperatura[i] = -Ct*(presion[i]**2/temperatura[i]**8.5)*(luminosidad[i]/radio[i]**2)
    gradiente_masa[i] = Cm*(presion[i]/temperatura[i])*radio[i]**2
    gradiente_luminosidad[i] = Cl*presion[i]**2*temperatura[i]**(nu-2)*radio[i]**2


# In[158]:


def paso6_a13(radio):
    
    temperatura[i] = temperatura_estimada[i]
    
    nu, Cl, ciclo = generacion_energia(temperatura)
    
    gradiente_luminosidad[i] = Cl*(presion_calculada[i]**2)*(temperatura_estimada[i])**(nu-2)*radio[i]**2
    delta1_luminosidad[i] = h*gradiente_luminosidad[i] - h*gradiente_luminosidad[i-1]
    delta2_luminosidad[i] = h*gradiente_luminosidad[i] - 2*h*gradiente_luminosidad[i-1] + h*gradiente_luminosidad[i-2]
    luminosidad_calculada[i] = luminosidad[i-1] + h*gradiente_luminosidad[i] - (1/2)*delta1_luminosidad[i] - (1/12)*delta2_luminosidad[i]
    
    luminosidad[i] = luminosidad_calculada[i]


# In[159]:


def paso7_a13(radio):
    
    gradiente_temperatura[i] = -Ct*(presion_calculada[i]**2/temperatura_estimada[i]**8.5)*(luminosidad_calculada[i]/radio[i]**2)
    delta1_temperatura[i] = h*gradiente_temperatura[i] - h*gradiente_temperatura[i-1]
    temperatura_calculada[i] = temperatura[i-1] + h*gradiente_temperatura[i] - (1/2)*delta1_temperatura[i]


# In[160]:


def paso8_a13(radio):
    
    while( abs(temperatura_calculada[i] - temperatura_estimada[i])/temperatura_calculada[i] > Erelmax ):
        
        temperatura_estimada[i] = temperatura_calculada[i]
        
        paso3_a12(radio)
        paso4_a12(radio)
        paso5_a12(radio)
        paso6_a13(radio)
        paso7_a13(radio)
        
    temperatura[i] = temperatura_calculada[i]


# Funciones para el algoritmo A.2.

# In[161]:


def paso1_a2(radio):
    
    temperatura[i] = A3*(1/radio[i] - 1/Rtot)
    presion[i] = k_prima*temperatura[i]**2.5
    masa[i] = masa_calculada[i]
    luminosidad[i] = luminosidad_calculada[i]
    
    nu, Cl, ciclo = generacion_energia(temperatura)
    
    gradiente_presion[i] = -Cp*(presion[i]/temperatura[i])*(masa[i]/radio[i]**2)
    gradiente_temperatura[i] = -Ct*(presion[i]**2/temperatura[i]**8.5)*(luminosidad[i]/radio[i]**2)
    gradiente_masa[i] = Cm*(presion[i]/temperatura[i])*radio[i]**2
    gradiente_luminosidad[i] = Cl*presion[i]**2*temperatura[i]**(nu-2)*radio[i]**2
    
    


# In[162]:


def paso2_a2(radio):
    
    delta1_temperatura[i-1] = h*gradiente_temperatura[i-1] - h*gradiente_temperatura[i-2]
    temperatura_estimada[i] = temperatura[i-1] + h*gradiente_temperatura[i-1] + (1/2)*delta1_temperatura[i-1]
    
    presion_estimada[i] = k_prima*(temperatura_estimada[i]**2.5)


# In[163]:


def paso7_a2(radio):
    
    if(radio[i] == 0):
        
        temperatura_calculada[i] = temperatura_estimada[i]
        presion_calculada[i] = presion_estimada[i]
        
    else:
        
        gradiente_temperatura[i] = -Ct_prima*(masa_calculada[i]/radio[i]**2)
        delta1_temperatura[i] = h*gradiente_temperatura[i] - h*gradiente_temperatura[i-1]
        temperatura_calculada[i] = temperatura[i-1] + h*gradiente_temperatura[i] - (1/2)*delta1_temperatura[i]


# In[164]:


def paso8_a2(radio):
    
    while( abs(temperatura_calculada[i] - temperatura_estimada[i])/temperatura_calculada[i] > Erelmax ):
        
        temperatura_estimada[i] = temperatura_calculada[i]
        
        presion_estimada[i] = k_prima*(temperatura_estimada[i]**2.5)
        paso3_a12(radio)
        paso7_a2(radio)
    
    temperatura[i] = temperatura_calculada[i]
    presion_calculada[i] = k_prima*(temperatura_calculada[i]**2.5)
    presion[i] = presion_calculada[i]


# In[165]:


def paso6_a2(radio):
    
    nu, Cl, ciclo = generacion_energia(temperatura)
    
    gradiente_luminosidad[i] = Cl*(presion_calculada[i]**2)*(temperatura_calculada[i]**(nu-2))*radio[i]**2
    delta1_luminosidad[i] = h*gradiente_luminosidad[i] - h*gradiente_luminosidad[i-1]
    delta2_luminosidad[i] = h*gradiente_luminosidad[i] - 2*h*gradiente_luminosidad[i-1] + h*gradiente_luminosidad[i-2]
    luminosidad_calculada[i] = luminosidad[i-1] + h*gradiente_luminosidad[i] - (1/2)*delta1_luminosidad[i] - (1/12)*delta2_luminosidad[i]
    
    luminosidad[i] = luminosidad_calculada[i]


# # EJECUCION DEL CODIGO

# In[166]:


#INTEGRACION DESDE LA SUPERFICIE

#print("E", "fase", "i", "radio", "presion", "temperatura", "luminosidad", "masa")

i = 0

while (i <=2 ):
    
    paso1(radio)
    #print("--", "INICIO", i, "%.7f" % radio[i], "%.7f" % presion[i],"%.7f" % temperatura[i], "%.6f" % luminosidad[i], "%.6f" % masa[i] )
    
    i = i + 1
    
for r in radio[2:len(radio)]:
    
#Algoritmo A.1.1. Masa y Luminosidad constantes

    paso1(radio)
    paso2_a11(radio)
    paso4_a11(radio)
    paso5_a11(radio)
    paso7_a11(radio)
    paso8_a11(radio)
    paso3_a11(radio)
    
    while( abs(masa[i] - Mtot)/Mtot < Erelmax and abs( luminosidad[i] - Ltot)/Ltot < Erelmax ):
        
        paso1(radio)
        paso2_a11(radio)
        paso4_a11(radio)
        paso5_a11(radio)
        paso7_a11(radio)
        paso8_a11(radio)
        paso3_a11(radio)
    
        #print("--", "A.1.1.", i, "%.7f" % radio[i], "%.7f" % presion[i],"%.7f" % temperatura[i], "%.6f" % luminosidad[i], "%.6f" % Mtot )
    
        i = i + 1 


# In[167]:


#Algoritmo A.1.2. Masa variable y Luminosidad constante

for r in radio[i:len(radio)]:
    
    paso1_a12(radio)
    paso2_a11(radio)
    paso3_a12(radio)
    paso4_a12(radio)
    paso5_a12(radio)
    paso7_a11(radio)    
    paso8_a12(radio)
    paso6_a12(radio)
        
    while( abs(masa[i] - Mtot)/Mtot > Erelmax and abs( luminosidad[i] - Ltot)/Ltot < Erelmax ):
        
        paso1_a12(radio)
        paso2_a11(radio)
        paso3_a12(radio)
        paso4_a12(radio)
        paso5_a12(radio)
        paso7_a11(radio)
        paso8_a12(radio)
        paso6_a12(radio)
    
        #print("--", "A.1.2.", i, "%.7f" % radio[i], "%.7f" % presion[i],"%.7f" % temperatura[i], "%.6f" % Ltot, "%.6f" % masa[i] )
        
        i = i + 1
        
#Algoritmo A.1.3. Masa y luminosidad variables
        
    paso9_a13(radio)
    
    while( parametro_transporte[i] > 2.5 ):
        
        paso1_a13(radio)
        paso2_a11(radio)
        paso3_a12(radio)
        paso4_a12(radio)
        paso5_a12(radio)
        #print(temperatura[i], temperatura_estimada[i], temperatura_calculada[i])
        paso6_a13(radio)
        nu, Cl, ciclo = generacion_energia(temperatura)
        #print(temperatura[i],nu, Cl, ciclo)
        paso7_a13(radio)
        paso8_a13(radio)
        paso9_a13(radio)
        
        #print(ciclo, "A.1.3.", i, "%.7f" % radio[i], "%.7f" % presion[i],"%.7f" % temperatura[i], "%.6f" % luminosidad[i], "%.6f" % masa[i], "%.6f" % parametro_transporte[i] )
        
        i = i + 1


# In[168]:


#Algoritmo A.2. Nucleo Convectivo
    
k_prima = presion_calculada[i]/temperatura_calculada[i]**2.5

if( parametro_transporte[i] <= 2.5 and abs( masa[i] - Mtot)/Mtot > Erelmax and abs( luminosidad[i] - Ltot)/Ltot > Erelmax ):
        
    while( radio[i] >= 0 and i < 96):
#Paramos en la capa 99 para que no de error     
        paso1_a2(radio)
        paso2_a2(radio)
        paso3_a12(radio)
        paso7_a2(radio)
        paso8_a2(radio)
        paso6_a2(radio)
        
        nu, Cl, ciclo = generacion_energia(temperatura)
            
        #print(ciclo, "CONVEC", i, "%.7f" % radio[i], "%.7f" % presion[i],"%.7f" % temperatura[i], "%.6f" % luminosidad[i], "%.6f" % masa[i], "%.6f" % parametro_transporte[i] )
        
        i = i + 1


# Deseamos conocer en que capa se produce la transicion radiación -> convección.
# 
# Indice de transicion en integracion up-down ----> ind_tr[0]

# In[169]:


for x in parametro_transporte:
    
    if x != 0 and x < 2.5:
        
        ind_tr = np.where(parametro_transporte == x)[0]
        #print(ind_tr)


# In[170]:


ind_tr[0] - 1


# Valores en la frontera (desde arriba)

# In[171]:


def valoresfrontera_down(parametro_transporte, radio, presion, temperatura, luminosidad, masa):
    
    radio_transicion_down = np.interp(2.5, [parametro_transporte[ind_tr[0]-1],parametro_transporte[ind_tr[0]]],[radio[ind_tr[0]-1], radio[ind_tr[0]]])
    presion_transicion_down = np.interp(radio_transicion_down, [radio[ind_tr[0]-1],radio[ind_tr[0]]], [presion[ind_tr[0]-1], presion[ind_tr[0]]])
    temperatura_transicion_down = np.interp(radio_transicion_down, [radio[ind_tr[0]-1], radio[ind_tr[0]]], [temperatura[ind_tr[0]-1], temperatura[ind_tr[0]]])
    luminosidad_transicion_down = np.interp(radio_transicion_down, [radio[ind_tr[0]-1], radio[ind_tr[0]]], [luminosidad[ind_tr[0]-1], luminosidad[ind_tr[0]]])
    masa_transicion_down = np.interp(radio_transicion_down, [radio[ind_tr[0]-1], radio[ind_tr[0]]], [masa[ind_tr[0]-1], masa[ind_tr[0]]])
    
    return radio_transicion_down, presion_transicion_down, temperatura_transicion_down, luminosidad_transicion_down, masa_transicion_down
    
    #print(radio_transicion_down)


# In[172]:


valores_down = valoresfrontera_down(parametro_transporte, radio, presion, temperatura, luminosidad, masa)
#print("transicion_down","%.2f" % valores_down[0],"%.2f" %  valores_down[1],"%.2f" %  valores_down[2],"%.2f" %  valores_down[3],"%.2f" %  valores_down[4])
radio_transicion_down = valores_down[0]
#print(radio_transicion_down)


# Guardamos los resultados obtenidos en la integracion desde la superficie

# In[173]:


temperatura_down = np.zeros_like(radio)
presion_down = np.zeros_like(radio)
masa_down = np.zeros_like(radio)
luminosidad_down = np.zeros_like(radio)
radio_down = np.zeros_like(radio)

#------------------------------------------------------------------------------------------------------------

temperatura_down = temperatura
presion_down = presion
masa_down = masa
luminosidad_down = luminosidad
radio_down = radio


# # INTEGRACION DESDE EL CENTRO##

# Invertimos el radio, redefinimos el paso de integracion (positivo) y reiniciamos los arrays con los que trabajamos

# In[174]:


radio = radio[::-1]
radio
#------------------------------------------------------------------------------------------------------------
h = Rini/100
#--------------------------------------------------------------------------------------------------------------
temperatura = np.zeros_like(radio)
presion = np.zeros_like(radio)
masa = np.zeros_like(radio)
luminosidad = np.zeros_like(radio)
gradiente_presion = np.zeros_like(radio)
gradiente_temperatura = np.zeros_like(radio)
gradiente_masa = np.zeros_like(radio)
gradiente_luminosidad = np.zeros_like(radio)
epsilon_ciclo_PP = np.zeros_like(radio)
epsilon_ciclo_CN = np.zeros_like(radio)
presion_estimada = np.zeros_like(radio)
delta1_presion = np.zeros_like(radio)
delta2_presion = np.zeros_like(radio)
temperatura_estimada = np.zeros_like(radio)
delta1_temperatura = np.zeros_like(radio)
presion_calculada = np.zeros_like(radio)
temperatura_calculada = np.zeros_like(radio) 
delta1_temperatura = np.zeros_like(radio)
masa_calculada = np.zeros_like(radio)
delta1_masa = np.zeros_like(radio)
delta1_luminosidad = np.zeros_like(radio)
delta2_luminosidad = np.zeros_like(radio)
luminosidad_calculada = np.zeros_like(radio)
parametro_transporte = np.zeros_like(radio)
#---------------------------------------------------------------------------------------------------------------


# In[175]:


radio


# Condiciones iniciales

# In[176]:


Mcentral = 0
Lcentral = 0
Tcentral = 2 #ponemos la que se nos da, podria ponerse la ultima temperatura de la integracion desde la superfici


# Funciones para el algoritmo B

# In[177]:


def paso1_b(radio):
    
    masa[i] = (Cm/3)*k_prima*Tcentral**1.5*radio[i]**3
    temperatura[i] = Tcentral - 0.008207*(mu**2)*(k_prima)*(Tcentral**1.5)*radio[i]**2
    nu, Cl, ciclo = generacion_energia(temperatura)
    #print(nu, Cl, ciclo)
    luminosidad[i] = (Cl/3)*k_prima**2*(Tcentral**(3+nu))*radio[i]**3
    presion[i] = k_prima*temperatura[i]**2.5
    
    gradiente_masa[i] = Cm*k_prima*temperatura[i]**1.5*radio[i]**2
    #gradiente_presion[i] = -Cp*k_prima*temperatura[i]**1.5*(masa[i]/radio[i]**2)
    gradiente_luminosidad[i] = Cl*(k_prima**2)*(temperatura[i]**(3+nu))*radio[i]**2
    gradiente_temperatura[i] = -Ct_prima*(masa[i]/radio[i]**2)


# # EJECUCION DEL CODIGO

# In[178]:


#INTEGRACION DESDE EL CENTRO

#print("E", "fase", "i", "radio", "presion", "temperatura", "luminosidad", "masa")

i = 0

while (i <=2 ):
    
    paso1_b(radio)
    nu, Cl, ciclo = generacion_energia(temperatura)
    
    #print("--", "CONVEC", i, "%.7f" % radio[i], "%.7f" % presion[i],"%.7f" % temperatura[i], "%.6f" % luminosidad[i], "%.6f" % masa[i] )
    
    i = i + 1


# In[179]:


for r in radio[i:len(radio)]:
    
    while(radio[i] <= radio_transicion_down):
        
        paso1_b(radio)
        paso2_a2(radio)
        paso3_a12(radio)
        paso7_a2(radio)
        paso8_a2(radio)
        paso6_a2(radio)
            
        nu, Cl, ciclo = generacion_energia(temperatura)
        
        #print(ciclo, "CONVEC", i, "%.7f" % radio[i], "%.7f" % presion[i],"%.7f" % temperatura[i], "%.6f" % luminosidad[i], "%.6f" % masa[i] )
    
        i = i + 1
        
        


# Queremos conocer el indice de la capa en la cual se produce la transicion conveccion -> radiacion. Utilizamos la longitud del array de temperatura de forma indiferente, cualquier array hubiera sido utilizado ya que solo queremos conocer su longitud que es igual para todos los parametros.
# 
# Indice de transicion en integracion down-up ---> ind_trans[0]

# In[180]:


ind_trans = len(temperatura) - ind_tr - 1
#print(ind_trans[0])


# Valores en la frontera(desde el centro)

# In[181]:


def valoresfrontera_up(radio, presion, temperatura, luminosidad, masa):
    
    presion_transicion_up = np.interp(radio_transicion_down, [radio[ind_trans[0]-1],radio[ind_trans[0]]], [presion[ind_trans[0]-1], presion[ind_trans[0]]])
    temperatura_transicion_up = np.interp(radio_transicion_down, [radio[ind_trans[0]-1], radio[ind_trans[0]]], [temperatura[ind_trans[0]-1], temperatura[ind_trans[0]]])
    luminosidad_transicion_up = np.interp(radio_transicion_down, [radio[ind_trans[0]-1], radio[ind_trans[0]]], [luminosidad[ind_trans[0]-1], luminosidad[ind_trans[0]]])
    masa_transicion_up = np.interp(radio_transicion_down, [radio[ind_trans[0]-1], radio[ind_trans[0]]], [masa[ind_trans[0]-1], masa[ind_trans[0]]])
    
    return presion_transicion_up, temperatura_transicion_up, luminosidad_transicion_up, masa_transicion_up


# In[182]:


valores_up = valoresfrontera_up(radio, presion, temperatura, luminosidad, masa)
#print("transicion_up","%.2f" % radio_transicion_down ,"%.2f" % valores_up[0],"%.2f" %  valores_up[1],"%.2f" %  valores_up[2],"%.2f" %  valores_up[3])


# Guardamos los resultados obtenidos en la integracion desde el centro

# In[183]:


temperatura_up = np.zeros_like(radio)
presion_up = np.zeros_like(radio)
masa_up = np.zeros_like(radio)
luminosidad_up = np.zeros_like(radio)
radio_up = np.zeros_like(radio)

#------------------------------------------------------------------------------------------------------------

temperatura_up = temperatura
presion_up = presion
masa_up = masa
luminosidad_up = luminosidad
radio_up = radio


# # AJUSTE DE LAS SOLUCIONES A UN RADIO INTERMEDIO

# In[184]:


valores_up = valoresfrontera_up(radio, presion, temperatura, luminosidad, masa)
#print("transicion_down","%.4f" % valores_down[0],"%.4f" %  valores_down[1],"%.4f" %  valores_down[2],"%.4f" %  valores_down[3],"%.4f" %  valores_down[4])
#print("transicion_up","%.4f" % radio_transicion_down ,"%.4f" % valores_up[0],"%.4f" %  valores_up[1],"%.4f" %  valores_up[2],"%.4f" %  valores_up[3])

radio_transicion_down, presion_transicion_down, temperatura_transicion_down, luminosidad_transicion_down, masa_transicion_down = valores_down
presion_transicion_up, temperatura_transicion_up, luminosidad_transicion_up, masa_transicion_up = valores_up

Err_rel_presion = abs(presion_transicion_down - presion_transicion_up)/presion_transicion_down*100
Err_rel_temperatura = abs( temperatura_transicion_down - temperatura_transicion_up)/temperatura_transicion_down*100
Err_rel_luminosidad = abs( luminosidad_transicion_down - luminosidad_transicion_up)/luminosidad_transicion_down*100
Err_rel_masa = abs(masa_transicion_down - masa_transicion_up)/masa_transicion_up*100

Err_rel_total = (Err_rel_presion**2 + Err_rel_temperatura**2 + Err_rel_luminosidad**2 + Err_rel_masa**2)**0.5

#print("Err.relat(%) -->","%.2f" % Err_rel_presion,"%.2f" % Err_rel_temperatura,"%.2f" % Err_rel_luminosidad,"%.2f" % Err_rel_masa,"%.2f" % Err_rel_total)


# ## CÁLCULO DE LA TEMPERATURA CENTRAL QUE MINIMIZA EL ERROR RELATIVO TOTAL

# In[185]:


temps = np.arange(1,2.5,0.001)
#Error_total = np.zeros(len(temps))
Error_total = []

for Tcentral in temps:
    
    l = 0
    i = 0
    
    #print(Tcentral)
    
    while (i <=2 ):
    
        paso1_b(radio)
        nu, Cl, ciclo = generacion_energia(temperatura)
    
        #print("--", "CONVEC", i, "%.7f" % radio[i], "%.7f" % presion[i],"%.7f" % temperatura[i], "%.6f" % luminosidad[i], "%.6f" % masa[i] )
    
        i = i + 1
        
    for r in radio[i:len(radio)]:
    
        while(radio[i] <= radio_transicion_down):
        
            paso1_b(radio)
            paso2_a2(radio)
            paso3_a12(radio)
            paso7_a2(radio)
            paso8_a2(radio)
            paso6_a2(radio)
            
            nu, Cl, ciclo = generacion_energia(temperatura)
        
            #print(ciclo, "CONVEC", i, "%.7f" % radio[i], "%.7f" % presion[i],"%.7f" % temperatura[i], "%.6f" % luminosidad[i], "%.6f" % masa[i] )
    
            i = i + 1
    
    ind_trans = len(temperatura) - ind_tr - 1
    
    temperatura_up = np.zeros_like(radio)
    presion_up = np.zeros_like(radio)
    masa_up = np.zeros_like(radio)
    luminosidad_up = np.zeros_like(radio)
    radio_up = np.zeros_like(radio)

#------------------------------------------------------------------------------------------------------------

    temperatura_up = temperatura
    presion_up = presion
    masa_up = masa
    luminosidad_up = luminosidad
    radio_up = radio
    
    valores_up = valoresfrontera_up(radio, presion, temperatura, luminosidad, masa)

    radio_transicion_down, presion_transicion_down, temperatura_transicion_down, luminosidad_transicion_down, masa_transicion_down = valores_down
    presion_transicion_up, temperatura_transicion_up, luminosidad_transicion_up, masa_transicion_up = valores_up

    Err_rel_presion = abs(presion_transicion_down - presion_transicion_up)/presion_transicion_down*100
    Err_rel_temperatura = abs( temperatura_transicion_down - temperatura_transicion_up)/temperatura_transicion_down*100
    Err_rel_luminosidad = abs( luminosidad_transicion_down - luminosidad_transicion_up)/luminosidad_transicion_down*100
    Err_rel_masa = abs(masa_transicion_down - masa_transicion_up)/masa_transicion_up*100

    Err_rel_total = (Err_rel_presion**2 + Err_rel_temperatura**2 + Err_rel_luminosidad**2 + Err_rel_masa**2)**0.5
    
    Error_total.append(Err_rel_total)
    
    #print(Tcentral, Err_rel_total)
    
    l = l + 1
    
t_opt = temps[Error_total.index(min(Error_total))]
error_opt = min(Error_total)
print("La temperatura central optima es %.3f" % t_opt)
print("y se minimiza el error hasta el %.3f por ciento" % error_opt)


# Ahora calculamos los parametros requeridos de temperatura, presion, luminosidad y masa en cada capa de la estrella cuya TEMPERATURA CENTRAL hemos encontrado que minimiza el error relativo total

# In[186]:


#radio = radio[::-1]
#radio
#------------------------------------------------------------------------------------------------------------
h = Rini/100
#--------------------------------------------------------------------------------------------------------------
temperatura = np.zeros_like(radio)
presion = np.zeros_like(radio)
masa = np.zeros_like(radio)
luminosidad = np.zeros_like(radio)
gradiente_presion = np.zeros_like(radio)
gradiente_temperatura = np.zeros_like(radio)
gradiente_masa = np.zeros_like(radio)
gradiente_luminosidad = np.zeros_like(radio)
epsilon_ciclo_PP = np.zeros_like(radio)
epsilon_ciclo_CN = np.zeros_like(radio)
presion_estimada = np.zeros_like(radio)
delta1_presion = np.zeros_like(radio)
delta2_presion = np.zeros_like(radio)
temperatura_estimada = np.zeros_like(radio)
delta1_temperatura = np.zeros_like(radio)
presion_calculada = np.zeros_like(radio)
temperatura_calculada = np.zeros_like(radio) 
delta1_temperatura = np.zeros_like(radio)
masa_calculada = np.zeros_like(radio)
delta1_masa = np.zeros_like(radio)
delta1_luminosidad = np.zeros_like(radio)
delta2_luminosidad = np.zeros_like(radio)
luminosidad_calculada = np.zeros_like(radio)
parametro_transporte = np.zeros_like(radio)
#---------------------------------------------------------------------------------------------------------------


# In[187]:


radio


# In[188]:


Mcentral = 0
Lcentral = 0
Tcentral = t_opt #ponemos la que temperatura central que minimiza el error


# In[189]:


#INTEGRACION DESDE EL CENTRO

#print("E", "fase", "i", "radio", "presion", "temperatura", "luminosidad", "masa")

i = 0

while (i <=2 ):
    
    paso1_b(radio)
    nu, Cl, ciclo = generacion_energia(temperatura)
    
    #print("--", "CONVEC", i, "%.7f" % radio[i], "%.7f" % presion[i],"%.7f" % temperatura[i], "%.6f" % luminosidad[i], "%.6f" % masa[i] )
    
    i = i + 1


# In[190]:


for r in radio[i:len(radio)]:
    
    while(radio[i] <= radio_transicion_down):
        
        paso1_b(radio)
        paso2_a2(radio)
        paso3_a12(radio)
        paso7_a2(radio)
        paso8_a2(radio)
        paso6_a2(radio)
            
        nu, Cl, ciclo = generacion_energia(temperatura)
        
        #print(ciclo, "CONVEC", i, "%.7f" % radio[i], "%.7f" % presion[i],"%.7f" % temperatura[i], "%.6f" % luminosidad[i], "%.6f" % masa[i] )
    
        i = i + 1


# In[191]:


ind_trans = len(temperatura) - ind_tr - 1


# In[192]:


valores_up = valoresfrontera_up(radio, presion, temperatura, luminosidad, masa)
#print("transicion_up","%.2f" % radio_transicion_down ,"%.2f" % valores_up[0],"%.2f" %  valores_up[1],"%.2f" %  valores_up[2],"%.2f" %  valores_up[3])


# In[193]:


temperatura_up = np.zeros_like(radio)
presion_up = np.zeros_like(radio)
masa_up = np.zeros_like(radio)
luminosidad_up = np.zeros_like(radio)
radio_up = np.zeros_like(radio)

#------------------------------------------------------------------------------------------------------------

temperatura_up = temperatura
presion_up = presion
masa_up = masa
luminosidad_up = luminosidad
radio_up = radio


# In[194]:


valores_up = valoresfrontera_up(radio, presion, temperatura, luminosidad, masa)
#print("transicion_down","%.4f" % valores_down[0],"%.4f" %  valores_down[1],"%.4f" %  valores_down[2],"%.4f" %  valores_down[3],"%.4f" %  valores_down[4])
#print("transicion_up","%.4f" % radio_transicion_down ,"%.4f" % valores_up[0],"%.4f" %  valores_up[1],"%.4f" %  valores_up[2],"%.4f" %  valores_up[3])

radio_transicion_down, presion_transicion_down, temperatura_transicion_down, luminosidad_transicion_down, masa_transicion_down = valores_down
presion_transicion_up, temperatura_transicion_up, luminosidad_transicion_up, masa_transicion_up = valores_up

Err_rel_presion = abs(presion_transicion_down - presion_transicion_up)/presion_transicion_down*100
Err_rel_temperatura = abs( temperatura_transicion_down - temperatura_transicion_up)/temperatura_transicion_down*100
Err_rel_luminosidad = abs( luminosidad_transicion_down - luminosidad_transicion_up)/luminosidad_transicion_down*100
Err_rel_masa = abs(masa_transicion_down - masa_transicion_up)/masa_transicion_up*100

Err_rel_total = (Err_rel_presion**2 + Err_rel_temperatura**2 + Err_rel_luminosidad**2 + Err_rel_masa**2)**0.5

#print("Err.relat(%) -->","%.2f" % Err_rel_presion,"%.2f" % Err_rel_temperatura,"%.2f" % Err_rel_luminosidad,"%.2f" % Err_rel_masa,"%.2f" % Err_rel_total)


# ## CÁLCULO DE LAS CAPAS "EXTRA"

# In[195]:


h = -Rini/100
radio_extra = np.arange(Rtot, Rini-h, h)

temperatura_extra = np.zeros_like(radio_extra)
presion_extra = np.zeros_like(radio_extra)
masa_extra = np.zeros_like(radio_extra)
luminosidad_extra = np.zeros_like(radio_extra)


# In[196]:


def paso1_extra(radio_extra):
    
    temperatura_extra[i] = A1*(1/radio_extra[i] - 1/Rtot)
    presion_extra[i] = A2*temperatura_extra[i]**4.25
    masa_extra[i] = Mtot
    luminosidad_extra[i] = Ltot
    
    #nu, Cl, ciclo = generacion_energia(temperatura)
    
    #gradiente_presion[i] = -Cp*(presion[i]/temperatura[i])*(masa[i]/radio[i]**2)
    #gradiente_temperatura[i] = -Ct*(presion[i]**2/temperatura[i]**8.5)*(luminosidad[i]/radio[i]**2)
    #gradiente_masa[i] = Cm*(presion[i]/temperatura[i])*radio[i]**2
    #gradiente_luminosidad[i] = Cl*presion[i]**2*temperatura[i]**(nu-2)*radio[i]**2


# In[197]:


#INTEGRACION DESDE LA SUPERFICIE PARA LAS CAPAS "EXTRA"

#print("E", "fase", "i", "radio", "presion", "temperatura", "luminosidad", "masa")

i = -len(radio_extra)

while (i < 0 ):
    
    paso1_extra(radio_extra)
    #print("--", "^^^^^^", i, "%.2f" % radio_extra[i], "%.4f" % presion_extra[i],"%.4f" % temperatura_extra[i], "%.3f" % luminosidad_extra[i], "%.3f" % masa_extra[i] )
    
    i = i + 1


# # PRESENTACION FINAL DEL MODELO

# In[198]:


temperatura = np.append(temperatura_extra, temperatura_down[0:81])
temperatura = np.append(temperatura, temperatura_up[19::-1])

presion = np.append(presion_extra, presion_down[0:81])
presion = np.append(presion, presion_up[19::-1])

luminosidad = np.append(luminosidad_extra, luminosidad_down[0:81])
luminosidad = np.append(luminosidad, luminosidad_up[19::-1])

masa = np.append(masa_extra, masa_down[0:81])
masa = np.append(masa, masa_up[19::-1])

radio = np.append(radio_extra, radio_down)


# In[199]:


headers_text = np.array(['Radio', 'Temperatura', 'Presion', 'Luminosidad', 'Masa'])

parametros = {'Radio': radio, 'Temperatura': temperatura, 'Presion': presion, 'Luminosidad': luminosidad, 'Masa': masa} 

parametros_table = pd.DataFrame(parametros, columns = headers_text)
#parametros_table


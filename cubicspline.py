# TFY41xx Fysikk vaaren 2021.
#
# Programmet tar utgangspunkt i hoeyden til de 8 festepunktene.
# Deretter beregnes baneformen y(x) ved hjelp av 7 tredjegradspolynomer, 
# et for hvert intervall mellom to festepunkter, slik at baade banen y, 
# dens stigningstall y' = dy/dx og dens andrederiverte
# y'' = d2y/dx2 er kontinuerlige i de 6 indre festepunktene.
# I tillegg velges null krumning (andrederivert) 
# i banens to ytterste festepunkter (med bc_type='natural' nedenfor).
# Dette gir i alt 28 ligninger som fastlegger de 28 koeffisientene
# i de i alt 7 tredjegradspolynomene.

# De ulike banene er satt opp med tanke paa at kula skal 
# (1) fullfoere hele banen selv om den taper noe mekanisk energi underveis;
# (2) rulle rent, uten aa gli ("slure").

# Vi importerer noedvendige biblioteker:
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

# Horisontal avstand mellom festepunktene er 0.200 m
h = 0.200
xfast=np.asarray([0,h,2*h,3*h,4*h,5*h,6*h,7*h])


# Vi begrenser starthøyden (og samtidig den maksimale høyden) til
# å ligge mellom 250 og 300 mm
ymax = 300
# yfast: tabell med 8 heltall mellom 50 og 300 (mm); representerer
# høyden i de 8 festepunktene
yfast=np.array([0.257,0.2,0.142,0.154,0.127,0.122,0.164,0.143])
#konverter fra m til mm
yfast =yfast/1000
# inttan: tabell med 7 verdier for (yfast[n+1]-yfast[n])/h (n=0..7); dvs
# banens stigningstall beregnet med utgangspunkt i de 8 festepunktene.
inttan = np.diff(yfast)/h
attempts=1
# while-løkken sjekker om en eller flere av de 3 betingelsene ovenfor
# ikke er tilfredsstilt; i så fall velges nye festepunkter inntil
# de 3 betingelsene er oppfylt
while (yfast[0] < yfast[1]*1.04 or
       yfast[0] < yfast[2]*1.08 or
       yfast[0] < yfast[3]*1.12 or
       yfast[0] < yfast[4]*1.16 or
       yfast[0] < yfast[5]*1.20 or
       yfast[0] < yfast[6]*1.24 or
       yfast[0] < yfast[7]*1.28 or
       yfast[0] < 0.250 or
       np.max(np.abs(inttan)) > 0.4 or
       inttan[0] > -0.2):
          yfast=np.asarray(np.random.randint(0, ymax, size=8))
          
          #konverter fra m til mm
          yfast =yfast/1000
          
          inttan = np.diff(yfast)/h
          attempts=attempts+1

# Omregning fra mm til m:
# xfast = xfast/1000
# yfast = yfast/1000

# Når programmet her har avsluttet while-løkka, betyr det at
# tallverdiene i tabellen yfast vil resultere i en tilfredsstillende bane. 

#Programmet beregner deretter de 7 tredjegradspolynomene, et
#for hvert intervall mellom to nabofestepunkter.


#Med scipy.interpolate-funksjonen CubicSpline:
yfast=np.array([0.257,0.2,0.142,0.154,0.127,0.122,0.164,0.143])
cs = CubicSpline(xfast, yfast, bc_type='natural')

xmin = 0.000
xmax = 1.401
dx = 0.001

x = np.arange(xmin, xmax, dx)   

#funksjonen arange returnerer verdier paa det "halvaapne" intervallet
#[xmin,xmax), dvs slik at xmin er med mens xmax ikke er med. Her blir
#dermed x[0]=xmin=0.000, x[1]=xmin+1*dx=0.001, ..., x[1400]=xmax-dx=1.400, 
#dvs x blir en tabell med 1401 elementer
Nx = len(x)
y = cs(x)       #y=tabell med 1401 verdier for y(x)
dy = cs(x,1)    #dy=tabell med 1401 verdier for y'(x)
d2y = cs(x,2)   #d2y=tabell med 1401 verdier for y''(x)

c=2/5
g=9.81



v=np.sqrt(2*g*(y[0]-y)/(1+(c)))




k=d2y/(1+dy**2)**(3/2)

a=(v**2)*k

beta=np.arctan(dy)

M=1

N=M*(g*np.cos(beta)+a)




f=c*M*g*np.sin(beta)/(1+c)

#abs(frisksjon)/Normalkraft:

    #Hvor stor friksjonkoeffisient kreves det for at kula skal rulle rent 
    #dette er friksjonskoeffisienten
absplot=abs(f/N)


#v*npcos(beta)

vx=v*np.cos(beta)

vx_average=np.zeros(Nx-1)

#metode 2
for i in range (Nx-1):
    vx_average[i]=(1/2)*(vx[i]+vx[i+1])


dt=dx/vx_average

t=np.zeros(Nx)

for i in range (Nx-1):
    t[i+1]=t[i]+dt[i]
    


file = open("video1Data.txt", "r")    #Åpnar fila

x_e = []                              
y_e = []                              #Lagar tomme lister
t_e = []
v_e = []

count = 0
for line in file:                   #For-løkke som går gjennom fila, linje for linje
    row = line.split()              #Splittar kvar linje inn i dei tre kolonnene

    if count > 1:
        t_e.append(float(row[0]))
        x_e.append(float(row[1]))  # Gjer om tekst til tal og legg dette
        y_e.append(float(row[2]))  # til i dei forskjellige listene
        try:
            v_e.append(float(row[3]))
        except Exception:
            pass

    count = count+1


file.close()                        #Lukkar fila!

print(t_e)
print(x_e)
print(y_e)

#Numerisk derivasjon

x_e = np.array(x_e)
y_e = np.array(y_e)
t_e = np.array(t_e)
v_e = np.array(v_e)

dy_e = np.zeros(len(x_e)-2)

for i in range(len(dy_e)):
    dy_e[i] = (y_e[i+2]-y_e[i])/(x_e[i+2]-x_e[i])

d2y_e = np.zeros(len(x_e)-4)

for i in range(len(d2y_e)):
    d2y_e[i] = (dy_e[i+2]-dy_e[i])/(x_e[i+3]-x_e[i+1])

#Merk: De deriverte har nå to færre elementer! De mangler den første og siste verdien
#For å bruke disse i beregninger må vi derfor slice arraysene.

    
k_e=d2y_e/(1+dy_e[1:-1]**2)**(3/2)

a_e=(v_e[2:-2]**2)*k_e

beta_e=np.arctan(dy_e[1:-1])

N_e=M*(g*np.cos(beta_e)+a_e)

##f_e=-c*M*a_e

f_e=c*M*g*np.sin(beta_e)/(1+c)

absplot_e=abs(f_e/N_e)

#Eksempel: Plotter banens form y(x)
baneform = plt.figure('y(x)',figsize=(6,3))
plt.plot(x,y,label="teoretisk")

plt.plot(xfast,yfast,'*',label="opphengningspunkter")
#
plt.plot(x_e,y_e,label="eksperimentell")
#
plt.legend()
#plt.title('Banens form')
plt.xlabel('$x$ [m]',fontsize=20)
plt.ylabel('$y(x)$ [m]',fontsize=20)
plt.ylim(0.0,0.40)
plt.grid()
plt.savefig("baneform.pdf", bbox_inches='tight') #denne er riktig
plt.show()
#Figurer kan lagres i det formatet du foretrekker:
#baneform.savefig("baneform.pdf", bbox_inches='tight')
#baneform.savefig("baneform.png", bbox_inches='tight')
#baneform.savefig("baneform.eps", bbox_inches='tight')


#baneform, fart(tid/pos), horisontal pos(t), friksjonskoeffisient, både teoretisk og eks sammen, plt.legend, ingen overskrift, tall på akser og aksetitler skal være lett å lese, figurer skal lagres i pdf, bruker linje som ligger i koden "plt.savefig...pdf

#fart
baneform = plt.figure('y(x)',figsize=(6,3))
plt.plot(x,v,label="teoretisk")
#
plt.plot(x_e,v_e,label="eksperimentell")
plt.legend()
#
#plt.title('Banens fart')
plt.xlabel('$Posisjon(x) [m]$',fontsize=20)
plt.ylabel('$Hastighet(v) [m/s]$',fontsize=20)
plt.grid()
plt.savefig("banefart(pos).pdf", bbox_inches='tight')
plt.show()


#// Krummning


baneform = plt.figure('y(x)',figsize=(6,3))
plt.plot(x,k)
#
plt.plot(x_e[5:-2],k_e[3:])
#
plt.title('Banens krumming')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$k(x)$ (m)',fontsize=20)
plt.grid()
plt.show()


# Sentripetalakselerasjon 
baneform = plt.figure('y(x)',figsize=(6,3))
plt.plot(x,a)
#
plt.plot(x_e[2:-2],a_e)
#
plt.title('Banens akselerasjon')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$a(x)$ (m)',fontsize=20)
plt.grid()
plt.show()



# Normalkraft 
baneform = plt.figure('y(x)',figsize=(6,3))
plt.plot(x,N,label="teoretisk")
#
plt.plot(x_e[2:-2],N_e,label="eksperimentell")
plt.legend()
#
#plt.title('Banens Normalkraft')
plt.xlabel('$Posisjon (x)$ [m]',fontsize=20)
plt.ylabel('$Normalkraft(x)$ [m])',fontsize=20)
plt.grid()
plt.savefig("Normalkraft.pdf", bbox_inches='tight')
plt.show()


#Friksjonskoeffisient
baneform = plt.figure('y(x)',figsize=(6,3))
plt.plot(x,absplot,label="teoretisk")
#
plt.plot(x_e[2:-2],absplot_e,label="eksperimentell")
plt.legend()
#
#plt.title('Banens friksjonskoeffisient abs(f/N)')
plt.xlabel('$Posisjon(x)$ [m]',fontsize=20)
plt.ylabel('$Friksjonskoeffisient$ (μ)',fontsize=20)
plt.grid()
plt.savefig("Friksjonskoeffisient.pdf", bbox_inches='tight')
plt.show()


#plotter horisontal posisjon over tid

baneform = plt.figure('y(x)',figsize=(6,3))
plt.plot(t,x,label="teoretisk")
#
plt.plot(t_e,x_e,label="eksperimentell")
#
plt.legend()
#plt.title('Banens horisontal posisjon over tid')
plt.xlabel('$Tid$ (t) [s]',fontsize=20)
plt.ylabel('$Posisjon(x)$ [m]',fontsize=20)
plt.grid()
plt.savefig("HorisontalPos.pdf", bbox_inches='tight')
plt.show()

#plotter fart over tid

baneform = plt.figure('y(x)',figsize=(6,3))
plt.plot(t,v,label="teoretisk")
#
plt.plot(t_e,v_e,label="eksperimentell")
plt.legend()
#
#plt.title('Banens fart over tid')
plt.xlabel('$Tid (t) [s]$',fontsize=20)
plt.ylabel('$Hastighet(v) [m/s]$',fontsize=20)
plt.grid()
plt.savefig("banefart(t).pdf", bbox_inches='tight')
plt.show()


print('Antall forsøk',attempts)
print('Festepunkthøyder (m)',yfast)
print('Banens høyeste punkt (m)',np.max(y))

print('NB: SKRIV NED festepunkthøydene når du/dere er fornøyd med banen.')
print('Eller kjør programmet på nytt inntil en attraktiv baneform vises.')

















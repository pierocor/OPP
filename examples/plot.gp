set terminal png
set output "morse.png"


eps=0.24
sig=3.4
D=0.24
r0=3.9
a=1.4

l(x)=4*eps*((sig/x)**12 -  (sig/x)**6  )

m(x)=D*(1 - exp(-a*(x-r0)) )**2 -D

set yr [: 0.5]

plot [0:10] l(x) title "Lennard-Jones", m(x) title "Morse"

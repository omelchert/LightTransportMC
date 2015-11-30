set terminal postscript enhanced eps color lw 2 28 font "Times-Roman"
set output "fluenceRate.eps"


LW=1.5
PS=1.25
set style line 1 lw LW lt 1 ps PS pt 6 lc rgb "red"
set style line 2 lw LW lt 2 ps PS pt 8 lc rgb "blue"

set key samplen 1.

set format y "10^{%L}"
set xlabel "{/Times-Italic r}  [cm]" 
set ylabel "{/Times-Italic h(r)}  [W cm^{-3}]"

set label 1 "{/Symbol m}_a = 2   [cm^{-1}]\n\n{/Symbol m}_s = 20 [cm^{-1}]\n\n{/Times-Italic N}_{/Symbol g} = 10^4" at graph 0.5,0.7


ma=2
ms=20
p(x) = (3/(4*pi*x))*ma*exp(-sqrt(3*ma*(ma+ms))*x)*(ma+ms)

set logs y

p\
'mc_ma2ms20_np10000.dat' every 2 u 1:2 w p ls 1 lw 1. t "Monte Carlo simulation",\
p(x) w l ls 2 t "Diffusion theory"


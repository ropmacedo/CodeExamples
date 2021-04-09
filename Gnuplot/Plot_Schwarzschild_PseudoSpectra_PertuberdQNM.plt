reset

set terminal x11


set terminal epslatex color
set output "PseudoSpectra+PertQNM_Schwarzschild.tex"

set multiplot

set size noratio 1.07, 0.26
set origin -0.14, 0.93


set xlabel "\\large{$n$}" offset 0, 1.5
set ylabel "{$\\log\\left(\\dfrac{\\kappa_n}{\\kappa_0}\\right)$}" rotate by 0 offset -0.5
set xrange [0:12]
set xtics 4
#set logscale y
set yrange [0:60]
set ytics 30
#set format y "$10^{%L}$"


pl "Data/Schwarzschild/CondNumbQNM0_spin_-2_l_2_N_200_Prec_500.dat" u 1:(log10($4/$5)) t"" w p pt 7 lc -1


##Countour Plot
set xlabel  "\\large{${\\rm Re}(\\omega_n)$}" offset 0,-0.1
set ylabel "\\large{${\\rm Im}(\\omega_n)$}" rotate by 0 offset -2.
set format y "% h"
unset logscale y

zmin=-50.
zmax=0.

#set palette rgbformulae 23,28,3 

#### COLOURS #####################
##	'#440154' # dark purple		##
##	'#472c7a' # purple			##
##	'#3b518b' # blue			##
##	'#2c718e' # blue			##
##	'#21908d' # blue-green		##
##	'#27ad81' # green			##
##	'#5cc863' # green			##
##	'#aadc32' # lime green		##
##	'#fde725' # yellow			##
##################################
set palette defined ( zmin     "#aadc32", (zmax+zmin)/2   "#21908d", zmax    "#0c0887" )

set pm3d 
set pm3d map
set cbrange[zmin:zmax]

set colorbox
set cntrparam level incremental zmin, 4, zmax
unset key
set cntrlabel onecolor


set xrange [-10.:10.]
set yrange [0.:12.]
set xtics 2
set ytics 1.


set size ratio -1 1.,1.
set origin 0.,0.



set contour base
unset surface


spl "Data/Schwarzschild/PseudoSpectraEnergyNorm0_spin_-2_l_2_N_200_Prec_500Re_s_[-6,0]_Im_s_[0,6].dat" u 2:(-$1):(log10(sqrt($3))) t"" w l lt -1 lw 1.5 lc rgb "white"
spl "Data/Schwarzschild/PseudoSpectraEnergyNorm0_spin_-2_l_2_N_200_Prec_500Re_s_[-6,0]_Im_s_[0,6].dat" u (-$2):(-$1):(log10(sqrt($3))) t"" w l lt -1 lw 1.5 lc rgb "white"

spl "Data/Schwarzschild/PseudoSpectraEnergyNorm0_spin_-2_l_2_N_200_Prec_500Re_s_[-6,0]_Im_s_[6,10].dat" u 2:(-$1):(log10(sqrt($3))) t"" w l lt -1 lw 1.5 lc rgb "white"
spl "Data/Schwarzschild/PseudoSpectraEnergyNorm0_spin_-2_l_2_N_200_Prec_500Re_s_[-6,0]_Im_s_[6,10].dat" u (-$2):(-$1):(log10(sqrt($3))) t"" w l lt -1 lw 1.5 lc rgb "white"

spl "Data/Schwarzschild/PseudoSpectraEnergyNorm0_spin_-2_l_2_N_200_Prec_500Re_s_[-10,-6]_Im_s_[0,6].dat" u 2:(-$1):(log10(sqrt($3))) t"" w l lt -1 lw 1.5 lc rgb "white"
spl "Data/Schwarzschild/PseudoSpectraEnergyNorm0_spin_-2_l_2_N_200_Prec_500Re_s_[-10,-6]_Im_s_[0,6].dat" u (-$2):(-$1):(log10(sqrt($3))) t"" w l lt -1 lw 1.5 lc rgb "white"

spl "Data/Schwarzschild/PseudoSpectraEnergyNorm0_spin_-2_l_2_N_200_Prec_500Re_s_[-10,-6]_Im_s_[6,10].dat" u 2:(-$1):(log10(sqrt($3))) t"" w l lt -1 lw 1.5 lc rgb "white"
spl "Data/Schwarzschild/PseudoSpectraEnergyNorm0_spin_-2_l_2_N_200_Prec_500Re_s_[-10,-6]_Im_s_[6,10].dat" u (-$2):(-$1):(log10(sqrt($3))) t"" w l lt -1 lw 1.5 lc rgb "white"

spl "Data/Schwarzschild/PseudoSpectraEnergyNorm0_spin_-2_l_2_N_200_Prec_500Re_s_[-12,-10]_Im_s_[0,10].dat" u 2:(-$1):(log10(sqrt($3))) t"" w l lt -1 lw 1.5 lc rgb "white"
spl "Data/Schwarzschild/PseudoSpectraEnergyNorm0_spin_-2_l_2_N_200_Prec_500Re_s_[-12,-10]_Im_s_[0,10].dat" u (-$2):(-$1):(log10(sqrt($3))) t"" w l lt -1 lw 1.5 lc rgb "white"


set surface
unset contour
unset pm3d




set key at -0.8, 14. samplen 3.
spl "Data/Schwarzschild/CleanQNM_spin_-2_l_2_NoFuncFreqPert_eps_0.0001_ksig_1._N_250_Prec_500.dat" u 2:(-$1):(1) t"$||\\delta \\tilde V_d||_{_E}\\!=\\!10^{-8},k=1$" w p pt 4 ps 1.6 lw 5 lc rgb "blue"

set key at -1., 13.2 samplen 1.9
spl "Data/Schwarzschild/CleanQNM_spin_-2_l_2_NoFuncFreqPert_eps_0.00000001_ksig_20._N_300_Prec_500.dat" u 2:(-$1):(1) t"$||\\delta \\tilde V_d||_{_E}\\!=\\!10^{-8},k = 20$" w p pt 12 ps 1.2 lw 4 lc rgb "blue"

set key at 10.3, 14. samplen 0.1
spl "Data/Schwarzschild/CleanQNM_spin_-2_l_2_NoFuncFreqPert_eps_0.00000001_ksig_60._N_300_Prec_500.dat" u 2:(-$1):(1) t"$||\\delta \\tilde V_d||_{_E}\\!=\\!10^{-8},k = 60$" w p pt 9 ps 1.4 lw 5 lc rgb "blue"

spl "Data/Schwarzschild/Spectra_spin_-2_l_2_N_200_Prec_500.dat" u 2:(-$1):(1) t"" w p pt 7 ps 1.1 lc rgb "red"

set key at 10.3, 13.2 samplen 0.1
spl "Data/Schwarzschild/CleanQNM_spin_-2_l_2_NoFuncFreqPert_eps_0.0001_ksig_20._N_250_Prec_500.dat" u 2:(-$1):(1) t"$||\\delta \\tilde V_d||_{_E}\\!=\\!10^{-4},k = 20$" w p pt 11 ps 1.4 lw 5 lc rgb "blue"


unset multiplot
#set terminal x11


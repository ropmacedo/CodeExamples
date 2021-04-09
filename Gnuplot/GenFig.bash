#!/bin/bash

########## INSTRUCTION ###################
#  1) Select Gnuplot script one wants to run
#  2) Observer file naming the output within script: e.g. NAME.tex
#  3) Find (or create) in folder Latex_Sources/ specific file with extension NAME_latex.tex
#  4) Certify that line 52 in NAME_latex.tex has the correct NAME.tex as input
#  5) Copy NAME_latex.tex to folder were this script is saved
#  6) Run ./GenFig
##########################################

rm -f *.log
rm -f *.aux
rm -f *.dvi
rm -f *.pdf

#gnuplot Plot_PolschTeller_PertPseudoSpectra_eps_0._N100.plt
#gnuplot Plot_PolschTeller_PertPseudoSpectraZoom_eps_0._N100.plt
#gnuplot Plot_PolschTeller_PertPseudoSpectra_eps_0.00000001_N100.plt
#gnuplot Plot_PolschTeller_PertPseudoSpectra_eps_0.0000000000000001_N100.plt
#gnuplot Plot_PolschTeller_PertPseudoSpectra_eps_1._N100.plt
#gnuplot Plot_PolschTeller-QNM-CondNum.plt
#gnuplot Plot_PolschTeller_PerturbedQNM.plt
#gnuplot Plot_PolschTeller_QNM_Convergence.plt
#gnuplot Plot_PolschTeller_PseudoSpectra_PertuberdQNM.plt
#gnuplot Plot_Schwarzschild-QNM-CondNum.plt
gnuplot Plot_Schwarzschild_PseudoSpectra_PertuberdQNM.plt
#gnuplot Plot_Schwarzschild_QNM-PolarAxial.plt
#gnuplot Plot_Schwarzschild_PseudoSpectra.plt
#gnuplot Plot_Schwarzschild_PseudoSpectraZoom.plt
#gnuplot Plot_Schwarzschild_PerturbedQNM.plt
#gnuplot Plot_WaveSphere-QNM-CondNumb.plt
#gnuplot Plot_WaveSphere_PseudoSpectra.plt
#gnuplot Plot_WaveSphere_PseudoSpectraZoom.plt

for f in *_latex.tex
do
 echo "******** Processing Latex $f********************************"
 latex $f

 f_dvi=${f%.*}.dvi
 

 echo "******** Processing dvipdf $f_dvi ******************************"
 dvipdf $f_dvi
 
 f_pdf=${f_dvi%.*}.pdf

 echo "******** Processing pdfcrop $f_pdf *****************************"
 pdfcrop $f_pdf
 
 echo "********************************************************************************"
 echo ""
done



for f_crop in *-crop.pdf

do
  f_name=${f_crop%_*}.pdf
  echo "******** removing crop and latex names from $f_name *****************************"
  mv $f_crop $f_name
  #okular $f_name &
  cp $f_name ../$f_name
done

rm -f *.log
rm -f *.aux
rm -f *.dvi
rm -f *.out
rm -f *_latex.pdf




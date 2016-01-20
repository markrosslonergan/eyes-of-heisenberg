#!/bin/bash 
cd /home/mark/projects/lepto/lepto++/tests/dispersion_test 
gnuplot << EOF

set terminal png  nocrop enhanced size 900,640 font "arial,16" 
set output 'dispersion_test_wNF.png'
set key inside right top vertical Right noreverse enhanced box title 'M/T'
set title "Dispersion Equations" 
set title  font ",20" norotate
plot 'MoT_2.00_wNF.dat' u 1:2 w l lw 3 title '2', 'MoT_1.50_wNF.dat' w l lw 3 title '1.5', 'MoT_1.20_wNF.dat' w l lw 3 title '1.2', 'MoT_1.00_wNF.dat' w l lw 3 title '1', 'MoT_0.90_wNF.dat' w l lw 3 title '0.9', 'MoT_0.70_wNF.dat' w l lw 3 title '0.7'

EOF

cd -

#!/opt/homebrew/bin/gnuplot

set terminal postscript enhanced color size 8in, 5.59in font 'Helvetica'
set key bottom left

set output 'compare_prism_0.eps'
set ylabel "Excitation Energy (eV)"
set xlabel "Fundamental Coupling Strength (a.u.)"
plot 'compare_prism.txt' u 1:(($3-$2)*27.211) w lp lt 1 pt 6 lc rgb 'black' lw 5 t 'CISS-PF', \
'compare_prism.txt' u 1:(($4-$2)*27.211) w lp lt 3 pt 6 lc rgb 'black' lw 5 notitle, \
'compare_prism.txt' u 1:(($6-$5)*27.211) w lp lt 1 pt 6 lc rgb '#DC143C' lw 4 t 'CISS-JC', \
'compare_prism.txt' u 1:(($7-$5)*27.211) w lp lt 3 pt 6 lc rgb '#DC143C' lw 4 notitle, \
'compare_prism.txt' u 1:(($9-$8)*27.211) w lp lt 1 pt 6 lc rgb '#FF7F50' lw 4 t 'CIS-PF', \
'compare_prism.txt' u 1:(($10-$8)*27.211) w lp lt 3 pt 6 lc rgb '#FF7F50' lw 4 notitle, \
'compare_prism.txt' u 1:(($12-$11)*27.211) w lp lt 1 pt 6 lc rgb '#5F9EA0' lw 4 t 'CIS-JC', \
'compare_prism.txt' u 1:(($13-$11)*27.211) w lp lt 3 pt 6 lc rgb '#5F9EA0' lw 4 notitle

set key top left
set output 'compare_prism_1.eps'
set ylabel "Energy (Hartrees)"
set xlabel "Fundamental Coupling Strength (a.u.)"
plot 'compare_prism.txt' u 1:(($3-$2)*27.211) w lp lt 1 pt 6 lc rgb '#DC143C' lw 5 t 'CISS-PF E_{LP}', \
'compare_prism.txt' u 1:(($4-$2)*27.211) w lp lt 3 pt 6 lc rgb '#8A2BE2' lw 5 t 'CISS-PF E_{UP}', \
#'compare_prism.txt' u 1:6 w lp lt 1 pt 6 lc rgb '#DC143C' lw 4 t 'CISS-JC', \
#'compare_prism.txt' u 1:7 w lp lt 3 pt 6 lc rgb '#DC143C' lw 4 notitle, \
#'compare_prism.txt' u 1:9 w lp lt 1 pt 6 lc rgb '#FF7F50' lw 4 t 'CIS-PF', \
#'compare_prism.txt' u 1:10 w lp lt 3 pt 6 lc rgb '#FF7F50' lw 4 notitle, \
#'compare_prism.txt' u 1:12 w lp lt 1 pt 6 lc rgb '#5F9EA0' lw 4 t 'CIS-JC', 
#'compare_prism.txt' u 1:13 w lp lt 3 pt 6 lc rgb '#5F9EA0' lw 4 notitle



set output 'compare_prism_2.eps'
set ylabel "Rabi Splitting (eV)"
set xlabel "Fundamental Coupling Strength (a.u.)"
plot 'compare_prism.txt' u 1:(($4-$3)*27.211) w lp lt 1 pt 6 lc rgb 'black' lw 4 t 'CISS-PF', \
'compare_prism.txt' u 1:(($7-$6)*27.211) w lp lt 1 pt 6 lc rgb '#DC143C' lw 4 t 'CISS-JC', \
'compare_prism.txt' u 1:(($10-$9)*27.211) w lp lt 1 pt 6 lc rgb '#FF7F50' lw 4 t 'CIS-PF', \
'compare_prism.txt' u 1:(($13-$12)*27.211) w lp lt 1 pt 6 lc rgb '#5F9EA0' lw 4 t 'CIS-JC'

#plot 'hf+_0.0_all.dat' u ($0*0.05 + 0.5):1 w l  lt 1      lc rgb 'black'   lw 4 t 'coherent-state basis', \
#     'hf+_0.0_all.dat' u ($0*0.05 + 0.5):2 w  p lt 1 pt 6 lc rgb '#DC143C' lw 4 t 'n = 2', \
#     'hf+_0.0_all.dat' u ($0*0.05 + 0.5):3 w  p lt 1 pt 6 lc rgb '#FF7F50' lw 4 t 'n = 5', \
#     'hf+_0.0_all.dat' u ($0*0.05 + 0.5):4 w  p lt 1 pt 6 lc rgb '#5F9EA0' lw 4 t 'n = 10', \
#     'hf+_0.0_all.dat' u ($0*0.05 + 0.5):5 w  p lt 1 pt 6 lc rgb '#6495ED' lw 4 t 'n = 15', \
#     'hf+_0.0_all.dat' u ($0*0.05 + 0.5):6 w  p lt 1 pt 6 lc rgb '#8A2BE2' lw 4 t 'n = 20'


#plot 'hf+_0.0_all.dat' u ($0*0.05 + 0.5):1 w l  lt 1      lc rgb 'black'   lw 4 t 'coherent-state basis', \
#     'hf+_0.0_all.dat' u ($0*0.05 + 0.5):2 w  p lt 1 pt 6 lc rgb '#DC143C' lw 4 t 'n = 2', \
#     'hf+_0.0_all.dat' u ($0*0.05 + 0.5):3 w  p lt 1 pt 6 lc rgb '#6495ED' lw 4 t 'n = 5', \
#     'hf+_0.0_all.dat' u ($0*0.05 + 0.5):4 w  p lt 1 pt 6 lc rgb '#8A2BE2' lw 4 t 'n = 10', \
#     'hf+_0.0_all.dat' u ($0*0.05 + 0.5):5 w  p lt 1 pt 6 lc rgb '#6495ED' lw 4 t 'n = 15', \
#     'hf+_0.0_all.dat' u ($0*0.05 + 0.5):6 w  p lt 1 pt 6 lc rgb '#8A2BE2' lw 4 t 'n = 20'

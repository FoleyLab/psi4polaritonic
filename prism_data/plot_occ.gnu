#!/opt/homebrew/bin/gnuplot

set terminal postscript enhanced color size 8in, 5.59in font 'Helvetica'
set key top left
set output 'Photon_occupation.eps'
set ylabel "Photon Occupation"
set xrange [0:0.1]
set xlabel "{/Symbol l} (atomic units)"
plot 'photon_occupation.txt' u 1:5 w lp lt 1 pt 6 lc rgb 'black' lw 4 t 'Total', \
'photon_occupation.txt' u 1:2 w lp lt 1 pt 6 lc rgb '#DC143C' lw 4 t '0^{th} Order', \
'photon_occupation.txt' u 1:3 w lp lt 1 pt 6 lc rgb '#FF7F50' lw 4 t '1^{st} Order', \
'photon_occupation.txt' u 1:4 w lp lt 1 pt 6 lc rgb '#5F9EA0' lw 4 t '2^{nd} Order'

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

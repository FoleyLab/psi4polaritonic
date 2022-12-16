#!/opt/homebrew/bin/gnuplot

set terminal postscript enhanced color size 8in, 5.59in font 'Helvetica'
set key bottom left

set output 'compare_prism_mghp_pes_lz_0.01.eps'
set ylabel "Energy (Hartree)"
set xlabel "Bondlength (Angstroms)"
set yrange [-199.77: -199.57]
plot 'prism_mghp_pes_lambda_0.01.txt' u 1:3 w lp lt 1 pt 6 lc rgb 'black' lw 5 t 'CISS-PF', \
'prism_mghp_pes_lambda_0.01.txt' u 1:4 w lp lt 1 pt 6 lc rgb 'black' lw 5 notitle, \
'prism_mghp_pes_lambda_0.01.txt' u 1:6 w lp lt 2 pt 6 lc rgb '#DC143C' lw 5 t 'CISS-JC', \
'prism_mghp_pes_lambda_0.01.txt' u 1:7 w lp lt 2 pt 6 lc rgb '#DC143C' lw 5 notitle, \
'prism_mghp_pes_lambda_0.01.txt' u 1:9 w lp lt 2 pt 6 lc rgb '#FF7F50' lw 5 t 'CIS-PF', \
'prism_mghp_pes_lambda_0.01.txt' u 1:10 w lp lt 2 pt 6 lc rgb '#FF7F50' lw 5 notitle, \
'prism_mghp_pes_lambda_0.01.txt' u 1:12 w lp lt 2 pt 6 lc rgb '#5F9EA0' lw 5 t 'CIS-JC', \
'prism_mghp_pes_lambda_0.01.txt' u 1:13 w lp lt 2 pt 6 lc rgb '#5F9EA0' lw 5 notitle



set output 'compare_prism_mghp_pes_lz_0.05.eps'
set ylabel "Energy (Hartree)"
set xlabel "Bondlength (Angstroms)"
set yrange [-199.77: -199.57]
plot 'prism_mghp_pes_lambda_0.05.txt' u 1:3 w lp lt 1 pt 6 lc rgb 'black' lw 5 t 'CISS-PF', \
'prism_mghp_pes_lambda_0.05.txt' u 1:4 w lp lt 1 pt 6 lc rgb 'black' lw 5 notitle, \
'prism_mghp_pes_lambda_0.05.txt' u 1:6 w lp lt 1 pt 6 lc rgb '#DC143C' lw 5 t 'CISS-JC', \
'prism_mghp_pes_lambda_0.05.txt' u 1:7 w lp lt 1 pt 6 lc rgb '#DC143C' lw 5 notitle, \
'prism_mghp_pes_lambda_0.05.txt' u 1:9 w lp lt 1 pt 6 lc rgb '#FF7F50' lw 5 t 'CIS-PF', \
'prism_mghp_pes_lambda_0.05.txt' u 1:10 w lp lt 1 pt 6 lc rgb '#FF7F50' lw 5 notitle, \
'prism_mghp_pes_lambda_0.05.txt' u 1:12 w lp lt 1 pt 6 lc rgb '#5F9EA0' lw 5 t 'CIS-JC', \
'prism_mghp_pes_lambda_0.05.txt' u 1:13 w lp lt 1 pt 6 lc rgb '#5F9EA0' lw 5 notitle


set output 'compare_prism_mghp_polariton_energies.eps'
set ylabel "Energy (Hartree)"
set xlabel "Fundamental Coupling Strength (a.u.)"
set yrange [-199.75: -199.6]
plot 'prism_mghp_polariton_energies_vs_lambda.txt' u 1:3 w lp lt 1 pt 6 lc rgb 'black' lw 5 t 'CISS-PF', \
'prism_mghp_polariton_energies_vs_lambda.txt' u 1:4 w lp lt 1 pt 6 lc rgb 'black' lw 5 notitle, \
'prism_mghp_polariton_energies_vs_lambda.txt' u 1:6 w lp lt 1 pt 6 lc rgb '#DC143C' lw 5 t 'CISS-JC', \
'prism_mghp_polariton_energies_vs_lambda.txt' u 1:7 w lp lt 1 pt 6 lc rgb '#DC143C' lw 5 notitle, \
'prism_mghp_polariton_energies_vs_lambda.txt' u 1:9 w lp lt 1 pt 6 lc rgb '#FF7F50' lw 5 t 'CIS-PF', \
'prism_mghp_polariton_energies_vs_lambda.txt' u 1:10 w lp lt 1 pt 6 lc rgb '#FF7F50' lw 5 notitle, \
'prism_mghp_polariton_energies_vs_lambda.txt' u 1:12 w lp lt 1 pt 6 lc rgb '#5F9EA0' lw 5 t 'CIS-JC', \
'prism_mghp_polariton_energies_vs_lambda.txt' u 1:13 w lp lt 1 pt 6 lc rgb '#5F9EA0' lw 5 notitle


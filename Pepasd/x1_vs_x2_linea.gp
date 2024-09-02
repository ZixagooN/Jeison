set term pdf
set title 'x_{1} vs x_{2} de t [0-9]'
set output 'grafica linea x1 vs x2.pdf'
set xlabel 'x_{1}'
set ylabel 'x_{2}'
plot 'resultados.dat' u 1:2 w lp pt 7 linecolor rgb 'purple' lw 3 title 'x1 vs x2' at 5.0, 5.0

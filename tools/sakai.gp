w = 8.5
h = 8.5

borders_width = 1.0
set terminal pdfcairo size w cm, h cm color dl 0.2
set output 'plot.pdf'


# MULTIPLOT SETTING
set multiplot

#layout
# A B C
# D E F

pathA = "mdc_pM_eta0.255_mu-0.850.dat"

reset
N = 501
set grid front; unset grid;
#unset colorbox
set tics scale 0.5
set xrange [0:N-1]
set yrange [0:N-1]
set xtics ('0' 0, 'pi/2' N/2, 'pi' N-1) offset 0.5,0
set ytics ('0' 0, 'pi/2' N/2, 'pi' N-1) offset 0.5,0
set ylabel "k_y" #offset 2.2,0 
set xlabel "k_x" #offset 2.2,0 
set palette defined (0 '#FFFFFF',0.5 '#F0EEEE',0.55 '#C9BBBB',0.6 '#DD6666', 0.7 '#990000',1 '#770000')

set colorbox# user origin graph 0.75,0.5 size 0.02,square/h/2-0.02
set cbrange [0:0.4]
set cbtics 0.2 #offset -0.8,0


set size square
# PlotA
plot pathA matrix notitle with image

unset ylabel




set cbrange [0:1]
set size square
set palette defined (0 'black', 1 'white')
set terminal gif animate
set output "ternary_serial_equal_ie.gif"
unset border
unset xtics
unset ytics
n=10000
i=0
load "animate.gp"
set output

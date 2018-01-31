
set terminal png size 1280,720
set output 'OMP-Time.png'

set xl "No of processes"
set yl "Time (secs)"
 set title "node cn02-01"

set yr [0:95]

set xr [0.9:14.1]
set key outside
set grid
plot 'omp-1-108.dat' u 1:2 w lp ls 1 title "108 atomic updates",\
  'omp-1-2916.dat' u 1:2 w lp ls 2 title "2916 atomic updates",\
  'omp-2-108.dat' u 1:2 w lp ls 3 title "108 no 3d law",\
  'omp-2-2916.dat' u 1:2 w lp ls 4 title "2916 no 3d law"

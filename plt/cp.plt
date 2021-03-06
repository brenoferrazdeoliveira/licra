FT= ` sed -n '/#define FT / p' ../licra.h | cut -f3 -d' ' `
Nx= ` sed -n '/#define Nx / p' ../licra.h | cut -f3 -d' ' `
Ny= ` sed -n '/#define Ny / p' ../licra.h | cut -f3 -d' ' `

if (FT==0){
	set terminal png size 567, 556 crop
	ext="png"
	unset xtics
	unset ytics
	unset colorbox
}
if (FT==1){
	set terminal pdfcairo font "Times, 10" size 8.5cm, 6.8cm
	ext="pdf"
	set xtics  offset  0.0, 0.2 (0, Nx/4, Nx/2, 3*Nx/4, Nx)
	set ytics  offset  0.2, 0.0 (0, Ny/4, Ny/2, 3*Ny/4, Ny)
}

set xrange  [0:Nx]
set yrange  [0:Ny]
set cbrange [0:1.0]
set cbtics offset -0.8, 0.0 autofreq 0.2
set border
set size ratio -1
set nokey
set palette rgbformulae 21, 22, 23

n= ` sed -n '/#define NF / p' ../licra.h | cut -f3 -d' ' `
i=0
while (i <= n ){
	set output sprintf("cp-%d.%s", i, ext)
	plot sprintf("../dat/cp-%d.dat", i) \
	     u ($1):($2):($3) matrix w image t''
	i= i+ 1
}

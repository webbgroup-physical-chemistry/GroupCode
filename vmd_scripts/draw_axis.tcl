proc usual_colors { plrx plrz plrz nrmx nrmy nrmz } { 
# Usage: usual_colors <polar x vector> <polar y vector> <polar z vector>
# <normal x vector> <normal y vector> <normal z vector>
#
# We use the polar axis and normal to the plane to calculate the plane
# because it's easier that way.
# Get the vectors from using g_interface with the -v flag

# Get rid of any reps already present
for {set i 0} {$i < 5000} {incr i} {
	mol delrep 0 top
    graphics top delete $i
}
# Make the pretty figures we usually use!
set gtpase [atomselect top "resid > 100 and not resname GNP and not resname MG"]
mol selection "[$gtpase text]"
mol addrep top
mol modstyle 0 top NewCartoon
mol modcolor 0 top colorid 10
mol modmaterial 0 top Transparent
set ral [atomselect top "resid < 101"]
mol selection "[$ral text]"
mol addrep top
mol modstyle 1 top NewCartoon
mol modcolor 1 top colorid 3 
set mg [atomselect top "resname MG"]
mol selection "[$mg text]"
mol addrep top
mol modstyle 2 top vdw
mol modcolor 2 top colorid 7 
set gnp [atomselect top "resname GNP"]
mol selection "[$gnp text]"
mol addrep top
mol modstyle 3 top licorice
mol modcolor 3 top colorid 6 
set nitrile_CA [atomselect top "name CA and (resid 15 or resid 17 or resid 24 or resid 25 or resid 26 or resid 28 or resid 29 or resid 30 or resid 51)"]
mol selection "[$nitrile_CA text]"
mol addrep top
mol modstyle 4 top vdw
mol modcolor 4 top colorid 1
set mutations [atomselect top "name CA and (resid 135 or resid 136)"]
mol selection "[$mutations text]"
mol addrep top
mol modstyle 5 top vdw
mol modcolor 5 top colorid 0

color Display Background white
color Display Foreground black
axes location Off
material change opacity Transparent 0.20

## I just want to look at the atoms used to make the plane...
#set CAgtpase [atomselect top "index 2015 or index 2147 or index 2167 or index 2173 or index 2187 or index 2206 or index 2221 or index 2244 or index 2265 or index 2289"]
#mol selection "[$CAgtpase text]"
#mol addrep top
#mol modstyle 6 top vdw
#mol modcolor 6 top colorid 10   
#set CAral [atomselect top "index 198 or index 234 or index 344 or index 358 or index 370 or index 384 or index 401 or index 422 or index 444 or index 739 or index 778"]
#mol selection "[$CAral text]"
#mol addrep top
#mol modstyle 7 top vdw
#mol modcolor 7 top colorid 3

#for {set i 0} {$i < 8} {incr i} {
#	mol showrep top $i 0
#}

# Here we start to build the lines...
set w 3
set plr { plrx plry plrz }
set nrm { nrmx nrmy nrmz }

set com1 [measure center $ral weight mass]
set com2 [measure center $gtpase weight mass]
set comt [vecscale .5 [vecadd $com1 $com2]]

set polar_plus [vecadd $comt [vecscale 30 $plr]]
set polar_mins [vecadd $comt [vecscale -30 $plr]]
draw color black
draw line $comt $polar_plus width $w
draw line $comt $polar_mins width $w

set crxprd [veccross $plr $nrm]
set plane_plus [vecadd $comt [vecscale 30 $crxprd]]
set plane_mins [vecadd $comt [vecscale -30 $crxprd]]
draw color black
draw line $comt $plane_plus width $w
draw line $comt $plane_mins width $w

# camera stuff
#display resetview
#display distance -2.0
#display height 3.0
#display scale 1.0

}

proc draw_box { cdx cdy cdz xsize ysize zsize color style } {

set minx [expr $cdx-$xsize/2]
set maxx [expr $cdx+$xsize/2]
set miny [expr $cdy-$ysize/2]
set maxy [expr $cdy+$ysize/2]
set minz [expr $cdz-$zsize/2]
set maxz [expr $cdz+$zsize/2]

draw color $color
draw line "$minx $miny $minz" "$maxx $miny $minz" width 5 style $style
draw line "$minx $miny $minz" "$minx $maxy $minz" width 5 style $style
draw line "$minx $miny $minz" "$minx $miny $maxz" width 5 style $style

draw line "$maxx $miny $minz" "$maxx $maxy $minz" width 5 style $style
draw line "$maxx $miny $minz" "$maxx $miny $maxz" width 5 style $style

draw line "$minx $maxy $minz" "$maxx $maxy $minz" width 5 style $style
draw line "$minx $maxy $minz" "$minx $maxy $maxz" width 5 style $style

draw line "$minx $miny $maxz" "$maxx $miny $maxz" width 5 style $style
draw line "$minx $miny $maxz" "$minx $maxy $maxz" width 5 style $style

draw line "$maxx $maxy $maxz" "$maxx $maxy $minz" width 5 style $style
draw line "$maxx $maxy $maxz" "$minx $maxy $maxz" width 5 style $style
draw line "$maxx $maxy $maxz" "$maxx $miny $maxz" width 5 style $style
}

proc rotate_axis {vec deg {molid top}} {
    # get the current matrix
    lassign [molinfo $molid get rotate_matrix] curr
    # the transformation matrix
    set r [trans axis $vec $deg]
    # get the new matrix
    set m [transmult $r $curr]
    # and apply it to the molecule
    molinfo $molid set rotate_matrix "{ $m }"
}

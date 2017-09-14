set all [atomselect top all frame last]

set zmax 300.

set exc_res [[atomselect top "name OH2 and abs(z) > $zmax"] get residue]
foreach res $exc_res {
	set sel [atomselect top "residue $res"]
	set z [lindex [$sel get z] 0]
	puts "$z $res [$sel num]"
	if {$z > 0} {
		set dz [expr $zmax - $z - 50.]
			} else {
		set dz [expr -$z - $zmax + 50.]
	}
	$sel moveby "0 0 $dz"
}

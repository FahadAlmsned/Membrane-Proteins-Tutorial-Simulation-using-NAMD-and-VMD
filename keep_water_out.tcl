#!/usr/local/bin/namd2
# keep water from entering gaps between protein and lipids
# need membrane be in the XY plane!
# (C) A. Aksimentiev (alek@ks.uiuc.edu) 

# This script will prevent water molecules in your namd simulations 
# from entering into the hydrophobic part of the membrane. This can be very 
# useful for initial equilibration of a protein/lipid/water system.  
# It works with recent NAMD binaries (version 2.4 and higher), which 
# has clearconfig command.

# Brief description:
# At initialization, the arrays of namd atom indices of all water oxygens 
# and of all lipid c21 carbons are created. They will be used later to 
# determine which water molecules to push out of the membrane and to compute 
# membrane dimensions. At every simulation step the pushing forces are applied 
# to a minor number of atoms, which doesn't affect NAMD performance. The 
# membrane dimensions are recalculated every N steps (regulated by the global 
# parameter $lipidCheckFreq) as well as the indices of the 
# water molecules that will be pushed during the next N steps (regulated by 
# the global parameter $waterCheckFreq). The membrane should be in XY plane. 
# To use this script one has to add the following lines to the namd 
# configuration script (without "#" symbol): 
# 
#
#  set waterCheckFreq              100
#  set lipidCheckFreq              100
#  set allatompdb                  a1b2c5_init.pdb
#  tclForces                       on
#  tclForcesScript                 keep_water_out.tcl
#
# # lipidCheckFreq must be an integer multiple of waterCheckFreq .



###################################################################
# user definitions begin
###################################################################

print "Starting Tcl forces"

# define force constant per group (Kcal/(mol*A))
set fconst 0.1
set pressure 1.0
set fconstUp [expr $fconst*$pressure]
set fconstDown [expr -$fconst*$pressure]

# define water and lipid names
set watResName "TIP3"
set watAtomName "OH2"
set lipResName "POPC"
set lipAtomName "C1 "

# define exceptions

set watexcept {0}


###################################################################
# user definitions stop here!
###################################################################



###################################################################
# preprocessing for calcforces
###################################################################

# define all water oxigens and lipid C21 :

set waters_list   {}
set c21plus_list  {}
set c21minus_list {}

set inStream [open $allatompdb r]
foreach line [split [read $inStream] \n] {
    set string1 [string range $line 0 3]
    set string2 [string range $line 6 10]
    set string3 [string range $line 17 20]	
    set string4 [string range $line 13 15]
    set string5 [string range $line 46 53]
    set string6 [string range $line 72 75]
    set string7 [string range $line 22 25]
   
    if { ([string equal $string1 {ATOM}] || \
	      [string equal $string1 {HETA}] ) && \
	     [string equal $watResName $string3] &&\
	     [string equal $watAtomName $string4] } {
	
	lappend waters_list "[string trim $string6]\
			    [string trim $string7] $watAtomName"
    }  
    if { ([string equal $string1 {ATOM}] || \
	      [string equal $string1 {HETA}] ) && \
	     [string equal $lipResName $string3] &&\
	     [string equal $lipAtomName $string4] } {
	if { [string trim $string5] >= 0 } {
	    lappend c21plus_list  "[string trim $string6]\
			    [string trim $string7] $lipAtomName"
	} else {
	    lappend c21minus_list "[string trim $string6]\
			    [string trim $string7] $lipAtomName"
	}
    } 
}
close $inStream

# make list of indices

set waters   {}
set c21plus  {}
set c21minus {}
foreach atomrecord $c21plus_list {
    foreach {segname resid atom} $atomrecord  { break }
    set atomindex [atomid $segname $resid $atom]
    lappend c21plus $atomindex
    addatom  $atomindex
}
foreach atomrecord $c21minus_list {
    foreach {segname resid atom} $atomrecord  { break }
    set atomindex [atomid $segname $resid $atom]
    lappend c21minus $atomindex
    addatom  $atomindex
}
foreach atomrecord $waters_list {
    foreach {segname resid atom} $atomrecord  { break }
    set atomindex [atomid $segname $resid $atom]
    set flag 0
    foreach excep $watexcept {
	if ($atomindex==$excep) {
	    set flag 1
	}
    }
    if ($flag==0) {
	lappend waters $atomindex
	addatom  $atomindex
    }
}

set c21plus [concat $c21plus] 
set c21minus [concat $c21minus]
set waters   [concat $waters]

if {([llength $c21plus] > 0) && ([llength $c21minus] > 0)} {
    set push 1
} else {
    print "WARNING: membrane has not been detected"
    set push 0
}

# initialize printing counter (independent on step counter)

set  pushCount $waterCheckFreq
set checkCount $lipidCheckFreq
set printcount 0
set waterstopushUp   {}
set waterstopushDown {}



###################################################################
# this procedure is executed at each time step
###################################################################

print "Starting calcforces..."
proc calcforces {} {
    global fconstUp fconstDown  stateread fstate fcount
    global stepcount lipidCheckFreq waterCheckFreq
    global pushCount checkCount printcount
    global waters c21plus c21minus zplus zminus
    global waterstopushUp waterstopushDown pressure
    global push

    if {$push == 1} {

##-------------------------  apply forces  ----------------------------###

#    print "  Up: $waterstopushUp"
#    print "Down: $waterstopushDown"
    foreach i $waterstopushUp  {
	set f [list 0.0 0.0 $fconstUp ]
#	print "Push up atom $i at z $z with force $f "
	addforce $i $f
    }

    foreach i $waterstopushDown {
	set f [list 0.0 0.0 $fconstDown ]
#	print "Push down atom $i at z $z with force $f "
	addforce $i $f
    }

###------- get atom indeces from NAMD before recalculation ----###

    if { $pushCount == [expr $waterCheckFreq -1] } {

	loadcoords coord	
#	print "Reconfiguring I ..."

	clearconfig
#	reconfig

	foreach atom $c21plus {
	    addatom  $atom
	}
	foreach atom $c21minus {
	    addatom  $atom
	}
	foreach atom $waters {
	    addatom  $atom
	}
	

    }
###------------ recalculate membrane size and waters to push -----------###

    if { $checkCount == $lipidCheckFreq } {

	loadcoords coord

	set zplus 0.0
	foreach index $c21plus {
	    foreach {x y z} $coord($index) { break }
#	    print "Z($index): $z"
	    set zplus [expr $zplus + $z]
	}
	set zplus [expr $zplus/double([llength $c21plus]) ]

	
	set zminus 0.0
	foreach index $c21minus {
	    foreach {x y z} $coord($index) { break }
#	    print "Z($index): $z"
	    set zminus [expr $zminus + $z]
	}
	set zminus [expr $zminus/double([llength $c21minus ] ) ]


	print "membrane dimensions at step $printcount: {$zplus $zminus}"
	
	set checkCount 0
    }
    

    if { $pushCount == $waterCheckFreq } {

	set waterstopushUp   {}
	set waterstopushDown {}
	set zHalf [expr  ($zplus - $zminus)/2.0 ]

	foreach index $waters {
	    foreach {x y z} $coord($index) { break }
###	    print "$z"
	    if { $z >= $zminus && $z <= $zplus } {
#		print "$z in $index"
		if { [expr $zplus - $z] <= $zHalf } {
		    lappend waterstopushUp $index
		} else {
		    lappend waterstopushDown $index
		}
	    }
	}
	
	print "Waters to push Up:   $waterstopushUp"
	print "Waters to push Down: $waterstopushDown"

	print "Reassign waters to push at step: $printcount"

	set pushCount 0

	clearconfig 
#	reconfig
#	print "Reconfiguring II ..."
	
	set waterstopush [concat $waterstopushUp $waterstopushDown]
	
	foreach atom $waterstopush {
	    addatom $atom
	}

    }


    incr printcount
    incr pushCount
    incr checkCount

    }

    return
}

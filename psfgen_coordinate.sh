#!/bin/bash

/home/shabir/Downloads/NAMD3.1/psfgen << ENDMOL

resetpsf

# Load topology files
topology toppar/top_all36_prot.rtf
topology toppar/top_all36_carb.rtf
topology toppar/top_all36_lipid.rtf
topology toppar/top_all36_na.rtf
topology toppar/top_all36_cgenff.rtf
topology toppar/toppar_all36_carb_imlab.str
topology toppar/toppar_water_ions.str
topology ARA.str

# Aliases for residue/atom name mismatches
alias residue ARA ARA

# ------------------- Protein -------------------
segment PRT {
    pdb protein.pdb
    first none
    last none
}

# Apply patches for special protonation or disulfides
patch NTER PRT:87
patch CTER PRT:450


patch ASPP PRT:142
patch ASPP PRT:171
patch ASPP PRT:191
patch ASPP PRT:247
patch ASPP PRT:265
patch ASPP PRT:280
patch ASPP PRT:282
patch ASPP PRT:336
patch ASPP PRT:338
patch ASPP PRT:380
patch ASPP PRT:387
patch ASPP PRT:545
patch ASPP PRT:587
patch ASPP PRT:630
patch GLUP PRT:99
patch GLUP PRT:165
patch GLUP PRT:322
patch GLUP PRT:500
patch GLUP PRT:520
patch GLUP PRT:558

patch DISU PRT:148 PRT:188
patch DISU PRT:543 PRT:562
patch DISU PRT:585 PRT:604
patch DISU PRT:628 PRT:647

# Coordinates
coordpdb protein.pdb PRT

# ------------------- Ligand -------------------
segment ARA {
    auto none
    pdb ARA.pdb
}
coordpdb ARA.pdb ARA

# ------------------- Waters -------------------
# Define residue and atom aliases for TIP3 waters
pdbalias residue HOH TIP3
pdbalias atom TIP3 O OH2

# TIP3A
segment TIP3A {
    auto none
    pdb TIP3A.pdb
}
coordpdb TIP3A.pdb TIP3A

# TIP3B
segment TIP3B {
    auto none
    pdb TIP3B.pdb
}
coordpdb TIP3B.pdb TIP3B

# TIP3C
segment TIP3C {
    auto none
    pdb TIP3C.pdb
}
coordpdb TIP3C.pdb TIP3C

# ------------------- Ions -------------------
segment SOD {
    auto none
    pdb sodium.pdb
}
coordpdb sodium.pdb SOD

segment CLA {
    auto none
    pdb chloride.pdb
}
coordpdb chloride.pdb CLA

# ------------------- Final -------------------
# Guess missing Hs
guesscoord

# Regenerate angles and dihedrals once — only after all segments
regenerate angles dihedrals

# Write final output
writepdb packmol_system.pdb
writepsf packmol_system.psf

ENDMOL

echo "✅ PSF and PDB built!"

# === 2️⃣ Step: Run VMD to write CHARMM .crd ===
vmd -dispdev text -eofexit << EOF
mol new packmol_system.psf
mol addfile packmol_system.pdb

# Define the procedure inline
proc writecharmmcoor {filename usemolid outtype} {
    set numatoms [molinfo \$usemolid get numatoms]
    set all [atomselect top "all"]

    if {[string match \$outtype "normal"]==1} {
        if { \$numatoms > 99999 } {
            puts "Using expanded format: number of atoms > 99999"
            set outtype "expanded"
        }
        set maxseg 0
        foreach {segname} [lsort -unique [\$all get segname]] {
            set current [string length \$segname]
            if { \$current > \$maxseg } { set maxseg \$current }
        }
        if { \$maxseg > 4 } {
            puts "Using expanded format: segname > 4 chars"
            set outtype "expanded"
        }
        set maxres 0
        foreach {resname} [lsort -unique [\$all get resname]] {
            set current [string length \$resname]
            if { \$current > \$maxres } { set maxres \$current }
        }
        if { \$maxres > 4 } {
            puts "Using expanded format: resname > 4 chars"
            set outtype "expanded"
        }
    }

    \$all delete

    set output [open \$filename "w"]
    puts \$output "* CHARMM coordinates generated from VMD"

    if {[string match \$outtype "normal"]==1} {
        puts \$output "[format "%5i" \$numatoms]"
    } else {
        puts \$output "[format "%10i" \$numatoms]  EXT"
    }

    set weighting "0"
    set countres "1"
    for {set i 0} {\$i < \$numatoms} {incr i} {
        set selection [atomselect \$usemolid "index \$i"]
        set segmentid [\$selection get segname]
        set residueid [\$selection get resid]

        if {\$i > 0} {
            if {\$prevresid != \$residueid || \$prevsegmentid != \$segmentid} {
                set countres [expr "\$countres + 1"]
            }
        }
        set resno "\$countres"
        set prevresid \$residueid
        set prevsegmentid \$segmentid

        if {[string match \$outtype "normal"]==1} {
            puts \$output "[format "%5i" [expr "\$i + 1"]][format "%5i" \$resno] [format "%-4s" [\$selection get resname]] [format "%-4s" [\$selection get name]][format "%10.5f" [\$selection get x]][format "%10.5f" [\$selection get y]][format "%10.5f" [\$selection get z]] [format "%-4s" \$segmentid] [format "%-4s" \$residueid][format "%10.5f" \$weighting]"
        } else {
            puts \$output "[format "%10i" [expr "\$i + 1"]][format "%10i" \$resno]  [format "%-8s" [\$selection get resname]]  [format "%-8s" [\$selection get name]][format "%20.10f" [\$selection get x]][format "%20.10f" [\$selection get y]][format "%20.10f" [\$selection get z]]  [format "%-8s" \$segmentid]  [format "%-8s" \$residueid][format "%20.10f" \$weighting]"
        }

        \$selection delete
    }
    close \$output
    puts "✅ CHARMM .crd written: \$filename"
}

# Run it:
writecharmmcoor packmol_system.crd 0 normal

quit
EOF

echo "✅ All done! PSF, PDB, and CHARMM CRD are ready!"

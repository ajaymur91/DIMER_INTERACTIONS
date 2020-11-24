# Script for dimer interactions
# Requires anaconda
# See ./readme_conda 

#!/bin/bash
# Stop if you encounter error
set -e
export GMX_MAXBACKUP=-1     # Overwrites
export PLUMED_MAXBACKUP=-1  # Unlimited backups
###############################

#source gmx_plumed/bin/activate 
# Edit inputs below

#Ion 
Ion=EC
Ion_q=0

#Solvent
Solv=EC
Solv_q=0


Tot_q=$(($Ion_q+$Solv_q))

#Solute_central_Atom
CA=C4

#Solvent_binding_Atom 
SA=O2

#Number of Solv
n=4

#Inner shell radius (nm)
R0=0.4

#FUNCTIONAL and BASIS SET
FUNCTIONAL=MP2
BASIS_SET=aug-cc-pvdz

#Trajectory sampling time (ps) - do not change
nsteps=250000     #dt is set to 0.5 fs in mdp file

###############################
# Take care of situation where Ion=Solv
[[ "$Ion" == "$Solv" ]] && Solv2="empty" || Solv2="$Solv"

cat << EOF > ion.top
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               3               yes             0.5     0.5

#include "itp/atomtypes.itp"
#include "itp/$Ion.itp"
#include "itp/$Solv2.itp"

[ system ]
; name
$Ion-$Solv-GP 

[ molecules ]
;name number
$Ion      1
EOF

###############################

cat << EOF > md.mdp
define                  = -DFLEXIBLE
integrator              = sd        ; leap-frog integrator
nsteps                  = 1000000     ; 2 * 50000 = 100 ps
dt                      = 0.0005     ; 0.5 fs

; Output control
nstxout                 = 2000     ; save coordinates every 1.0 ps
nstvout                 = 2000     ; save velocities every 1.0 ps
nstenergy               = 2000     ; save energies every 1.0 ps
nstcalcenergy           = 1
nstlog                  = 2000       ; update log file every 1.0 ps


; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
cutoff-scheme            = group    ; Buffered neighbor searching
rlist                    = 0   ;gas ph min (0 means no cutoff)
rcoulomb                 = 0   ;gas ph min (0 means no cutoff)
rvdw                     = 0   ;gas ph min (0 means no cutoff)
pbc                      = no
nstlist                  = 0
ns-type                  = simple
constraints              = none      ; bonds involving H are constrained
continuation             = no  ; does the same thing as unconstrained_start

;       Vacuum simulations are a special case, in which neighbor lists and cutoffs are basically thrown out.
;       All interactions are considered (rlist = rcoulomb = rvdw = 0) and the neighbor list is fixed (nstlist = 0).

; Run parameters

; Temperature coupling is on
;tcoupl                  = nose-hoover           ; modified Berendsen thermostat
tc-grps                 = system                ; two coupling groups - more accurate
tau_t                   = 2                   ; time constant, in ps
ref_t                   = 313                   ; reference temperature, one for each group, in K

EOF

###############################

cat << EOF > min.mdp
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1.0          ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
cutoff-scheme            = group    ; Buffered neighbor searching
rlist                    = 0   ;gas ph min (0 means no cutoff)
rcoulomb                 = 0   ;gas ph min (0 means no cutoff)
rvdw                     = 0   ;gas ph min (0 means no cutoff)
pbc                      = no
nstenergy                = 10
nstxout                  = 10
nstlist                  = 0
ns-type                  = simple
continuation             = no  ; does the same thing as unconstrained_start
;unconstrained_start     = yes ; depricated

;       Vacuum simulations are a special case, in which neighbor lists and cutoffs are basically thrown out.
;       All interactions are considered (rlist = rcoulomb = rvdw = 0) and the neighbor list is fixed (nstlist = 0).
EOF

###############################

# Start from Ion.gro and insert waters and generate force-field
cp ion.top system.top
gmx insert-molecules -f struct/"$Ion".gro -ci struct/$Solv.gro -o IonW.gro -box 1 1 1 -nmol "$n" -try 1000 -scale 0.1 &> /dev/null
cat << EOF >> system.top
$Solv   $n
EOF


###############################

# Make index and plumed.dat ( to enforce QCT criterion )
gmx select -f IonW.gro -s IonW.gro -on CA.ndx -select "atomname $CA and resnr 1" &> /dev/null
gmx select -f IonW.gro -s IonW.gro -on SA.ndx -select "atomname $SA and resnr > 1" &> /dev/null

cat << EOF > plumed.dat
CA: GROUP NDX_FILE=CA.ndx NDX_GROUP=atomname_${CA}_and_resnr_1
SA: GROUP NDX_FILE=SA.ndx NDX_GROUP=atomname_${SA}_and_resnr_>_1
cn: COORDINATION GROUPA=CA GROUPB=SA R_0=$R0 NN=1000000
EOF

k=0;
for i in $(tail -n 1 SA.ndx)
do
k=$(($k+1))
echo "d$k: DISTANCE ATOMS=CA,$i" >> plumed.dat
echo -e "UPPER_WALLS ARG=d$k AT=$R0 KAPPA=100000 LABEL=w$k" >> plumed.dat
cat << EOF >> plumed.dat
EOF
done
echo "PRINT ARG=(d[1-$n]),cn FILE=COLVAR STRIDE=1000" >> plumed.dat

###############################

# Minimize
echo -e "\n Run Minimization $Ion - $n $Solv \n"
gmx grompp -f min.mdp -c IonW.gro -p system.top -o min.tpr  &> /dev/null
gmx mdrun -deffnm min -nsteps 1000000 -plumed plumed.dat  &> /dev/null

# Make box bigger (does not really matter - but do it anyway to be safe)
gmx editconf -f min.gro -box 10 10 10 -o min2.gro &> /dev/null

###############################

# Md (100 ps)
 echo -e "\n Run MD - $Ion - $n $Solv \n"
 gmx grompp -f md.mdp -c min2.gro -p system.top -o md.tpr &> /dev/null
 gmx mdrun -deffnm md -nsteps $nsteps &> /dev/null

###############################

#Extract gro files from Trajectory
rm -rf gro/*
mkdir -p gro
k=1
j=1

# Begin from the middle of the trajectory
b=$(gmx check -f md.trr 2> /dev/stdout | grep "Step   " | awk '{print $2*$3/2}')

for i in $(tail -n 1 SA.ndx)
do
	gmx select -f md.gro -s md.tpr -on "$k".ndx -select "resid 1 || same resnr as atomnr $i" &> /dev/null
	echo "Extract configs -> gro/$p"
	gmx trjconv -f md.trr -b $b -s md.tpr -vel no -o gro/"$k"-.gro -sep -n "$k".ndx -dt 2 &> /dev/null
	k=$(($k+1))
done

###############################

# Make .xyz (coords in Angstrom. For making gaussian input files)
#rm -rf gro2/*
#cp -r gro gro2
#
#for k in $(ls gro/*.gro | cut -d '.' -f1| cut -d '/' -f2)
#do
#sed 's/[^ ]\+/\L\u&/g' gro/$k.gro > gro2/$k.gro
#done

rm -rf xyz/*
mkdir -p xyz
cd gro 
obabel *.gro -oxyz -m
mv *.xyz ../xyz/
cd ..

###############################

# Make gaussian input files (eg run on phoenix)
rm -rf com/*
mkdir -p com
f1=$(sed -n '2p' struct/$Ion.gro | awk '{print $1}')
f2=$(sed -n '2p' struct/$Solv.gro | awk '{print $1}')
rm -rf fragment

# gaussian dimer calc requires fragment definition (eg see com/*.com)
for i in `seq 1 $f1`
do
	echo "(fragment=1)" >> fragment
done
for i in `seq 1 $f2`
do
	echo "(fragment=2)" >> fragment
done


# Iterate over all xyz/*.xyz and make com/*.com
for k in $(ls gro/*.gro | cut -d '.' -f1| cut -d '/' -f2)
do
cat << EOF > com/"$k".com
%chk=$k.chk
%nproc=4
%mem=4GB
# $FUNCTIONAL/$BASIS_SET SP Counterpoise=2

dimer

$Tot_q,1,$Ion_q,1,$Solv_q,1
EOF
tail -n +3 xyz/"$k".xyz | paste fragment /dev/stdin | awk '{printf $2 "" $1 "\t" $3 "\t" $4 "\t" $5 "\n"}'>> com/"$k".com
echo -e "\n \n" >> com/"$k".com
done

###############################
# Extract FF energy (remove intramoleculars and calculate only intermolecular interactions) 

# Make gro/nrexcl.top then remove intramoleculars
cp ion.top 1.top
cat << EOF >> 1.top
$Solv   1
EOF

gmx grompp -f min.mdp -c gro/"$k".gro -p 1.top -o min.tpr -pp gro/nrexcl.top &> /dev/null

# First, edit gro/nrexcl.top to do the following 3 things

#1) Set nrexcl very high (makes lj and couloumb intramoleculars to zero)
awk '/^$/{flag=0} {if(flag==1) $2=50} 1; /  nrexcl/{flag=1}' gro/nrexcl.top > gro/tmp && mv gro/tmp gro/nrexcl.top

#2) delete pairs angles and dihedrals 
sed -i -e '/ pairs /,/^$/{//b' -e '/./d;}' gro/nrexcl.top
sed -i -e '/ angles /,/^$/{//b' -e '/./d;}' gro/nrexcl.top
sed -i -e '/ dihedrals /,/^$/{//b' -e '/./d;}' gro/nrexcl.top

#3) Make bonds strengths zero (we still need bond connectivity for nrexcl stuff)
awk '/^$/{flag=0} {if(flag==1) $4=0; if(flag==1)$5=0} 1; / bonds /{flag=1}' gro/nrexcl.top > gro/tmp && mv gro/tmp gro/nrexcl.top

# Make sp.mdp for FF interactions
cat << EOF > gro/sp.mdp 
integrator               = md
nsteps                   = 0   ;single point calc
cutoff-scheme            = group
rlist                    = 0   ;vaccum s.p. calc
rcoulomb                 = 0   ;vaccum s.p. calc
rvdw                     = 0   ;vaccun s.p. calc
pbc                      = no
nstenergy                = 1
nstxout                  = 1
nstlist                  = 0
ns-type                  = simple
continuation             = yes ; does the same thing as unconstrained_start
;unconstrained_start     = yes ; depricated
EOF

Ngro=$(ls -lth gro/*.gro | wc -l)

# calculate energy
rm -rf results/*FF*.txt
mkdir -p results/
cd gro
for i in $(ls *.gro | cut -d '.' -f1) #`seq 1 $Ngro`
do
	echo "Calc FF Energy $i.gro "
        gmx grompp -f sp.mdp -c $i.gro -p nrexcl.top -o $i.tpr &> /dev/null
        gmx mdrun -deffnm $i -rerun $i.gro &> /dev/null
        echo 4 | gmx energy -f $i.edr -o $i.xvg &> /dev/null
	tail -n 1 $i.xvg | awk '{print $2/4.184}' >> ../results/E_FF_kcal.txt
done
cd ..

######################################

# Clean up
# Remove gromacs backups
rm -rf \#*
# Remove plumed backups if any
rm -rf \#* bck.* step*

# Clean everything for fresh start
#rm -rf com xyz gro results md.* min.* COLVAR 1.ndx 1.top 2.ndx 3.ndx 4.ndx CA.ndx SA.ndx ion.top IonW.gro fragment mdout.mdp md_prev.cpt plumed.dat min2.gro system.top

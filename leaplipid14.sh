# Prepare molecule for simulation
# 1) Run AM1BCC w/ antechamber
# 2) Check for missing parameters
# 3) Prepare leap file from template
# 4) Run tleap to get amber topology

libdir=SimulationWorkflow/1-makeForceField
file=$1
resname=$2
intype=mol2
#outtype=ac
#intype=pdb
outtype=mol2
forcefield=lipid14

module load amber/amber16

# === Run Antechamber ===
antechamber -i $file.$intype -fi $intype -o $file.$outtype -fo $outtype -rn $resname -c bcc -pf y -s 2 -at gaff2

# === Run Parameter Check === 
parmchk2 -i $file.$outtype -f $outtype -o $file.frcmod -s gaff2

# === Create tleap File ===
leapfile=$file.leap.in
sed "s/__filetype__/${outtype}/g" ${libdir}/template.leaplipid14.in > $leapfile
#sed "s/__filetype__/${outtype}/g" ./template.leap.in > $leapfile
sed -i "s/__name__/${file}/g" $leapfile
sed -i "s/__forcefield__/${forcefield}/g" $leapfile

# === Create Amber Topology with tleap ===
# need to make lipid14 atom type substituions manually
python ${libdir}/lipid14subs.py $file.$outtype
tleap -f $leapfile > ${file}lipid14.leap.out
#xleap -f $leapfile > $file.leap.out

# === Convert to Gromacs using parmed ===
: '
python -c "
import parmed
amber = parmed.load_file('${file}.parm7', xyz='${file}.rst7')
amber.save( '${file}.top', parameters='${file}.itp' )
amber.save( '${file}.gro' )
"
'
python ${libdir}/amber2gro.py $file -itp
#python ../hydroxynator.py ${file}.top -s 1 -e 1 -o ${file}_scaled.top



#!/usr/bin/bash

# This script prepares a double-layered 'computational electrophysiology' 
# molecular dynamics system from an existing single-layered GROMACS system.
# Works with GROMACS 5.0.
# 
# Adapt the bottom section to your needs!


#==============================================================================
function func.testquit
{
    if [ "$1" = "0" ] ; then
        local RESULT=OK  # dummy line; do nothing :)
    else
        echo "ERROR: exit code of the last command was $1. Exiting."
        
        exit
    fi
}
#==============================================================================



#==============================================================================
# Extract index group names from an .mdp file.
# Expects the .mdp file identifier string (as e.g. "tc_grps") as second argument 
function func.getGroups
{
    if [ $# -ne 2 ]; then
        echo "ERROR: func.getGroups needs two arguments: 1) an .mdp file name 2) .mdp file identifier string!" >&2
        echo "       It got: '$@'" >&2
        exit 1234
    fi

    local MDPFILE=$1
    local STRING=$2
    
    local GROUPNAMES=$( grep "^$STRING" "$MDPFILE" | awk '{ for(i=3;i<=NF;i++) { print $i } }' )
    func.testquit $?
    echo "$GROUPNAMES"
}
#==============================================================================



#==============================================================================
# Build double-membrane setup from single-layered system!
function func.makeSandwich
{
    local SINGLE=".single"  # suffix for single-layered files
    local DOUBLE=".double"  # suffix for double-layered output files

    echo "=== Checking whether all necessary files can be found ==="

    # Check whether necessary input files can be found
    for FILE in "$GRO.gro" "$MDP.mdp" "$NDX.ndx" "$TOP.top" ; do
        if [ ! -f "$FILE" ] ; then
            echo "ERROR: Can not find input file $FILE"
            return 1
        fi
    done

    # Check whether ion and solvent group names appear in the provided index file:
    for GROUP in "$IONSGROUP" "$SOLVGROUP" "$SPLITGROUP" ; do
        CHECK=$( grep "\[ $GROUP \]" "$NDX.ndx" | awk '{ print $2 }')
        if [ "$CHECK" == "$GROUP" ] ; then
            echo "Found group \"$GROUP\" in index file $NDX.ndx"
        else
            echo "ERROR: Could not find group \"$GROUP\" in index file $NDX.ndx!"
            exit 111
        fi
    done
    
    # Preprocess the single-layered system, to
    # a) see wether that works, 
    # b) get a recent .mdp output file
    echo
    echo "=== Preprocessing single-layered system ==="
    COMMAND="$GMX grompp -quiet -f $MDP.mdp -c $GRO.gro -r $GRO.gro -p $TOP.top -n $NDX.ndx -o $TOP$SINGLE.tpr -po $MDP$SINGLE.mdp -maxwarn 1" 
    echo "RUNNING: $COMMAND"
    $COMMAND > grompp.single.out 2>&1
    func.testquit $?
    echo "Single-layered system passed preprocessing step :)"
    echo
   
    # Make molecules whole in single-layered starting conformation
    echo "Since we are changing periodicity, we need to make all molecules whole *before* duplicating the system:"
    COMMAND="$GMX trjconv -quiet -f $GRO.gro -o $GRO.whole.gro -pbc mol -s $TOP$SINGLE.tpr"
    echo "RUNNING: $COMMAND"
    echo 0 | $COMMAND > trjconv.out 2>&1
    func.testquit $?

    # Duplicate & translate:
    echo
    echo "=== Duplicating and translating ==="
    COMMAND="$GMX genconf -quiet -f $GRO.whole.gro -nbox 1 1 2 -o $GRO$DOUBLE.gro"
    echo "RUNNING: $COMMAND"
    $COMMAND
    func.testquit $?

    # For convenience also provide a .pdb file:
    rm -f $GRO$DOUBLE.pdb
    $GMX editconf -quiet -f $GRO$DOUBLE.gro -o $GRO$DOUBLE.pdb
    func.testquit $?

    cat $TOP.top > $TOP$DOUBLE.top
    func.testquit $?
    # Attach the [ molecules ] section a second time at the end of the topology file:
    cat $TOP.top | sed -e '1,/molecules/d' >> $TOP$DOUBLE.top
    func.testquit $?

    # Now prepare the index file by duplication of all index groups found in the single-layered file:
    rm -rf $NDX$DOUBLE.ndx
    COMMAND="$GMX make_ndx -quiet -twin -n $NDX.ndx -o $NDX$DOUBLE.ndx"
    echo "RUNNING: $COMMAND"
    echo "q\n" | $COMMAND 1> make_ndx.twin.out
    func.testquit $?
    

    # Extract temperature coupling group names from .mdp file.
    # These also need to be duplicated for the sandwich
    TCGROUPS=$( func.getGroups $MDP$SINGLE.mdp "tc_grps" )
    echo
    echo "Found these temperature coupling groups in the single-layered system:"
    echo $TCGROUPS

    COMGROUPS=$( func.getGroups $MDP$SINGLE.mdp "comm_grps" )
    echo
    echo "Found these center of mass motion removal groups in the single-layered system:"
    echo $COMGROUPS

    # Add a group with the complete solvent and ions,
    # also add the temperature coupling and COM-removal groups to the index file by
    # combining the groups of each of the layers
    export GMX_MAXBACKUP=-1
    unique_groups=$(echo "$SOLVGROUP $IONSGROUP $SPLITGROUP $TCGROUPS $COMGROUPS" | tr ' ' '\n' | sort -u)
    for GROUP in $unique_groups; do
        printf "\"%s\" | \"%s_copy\"\nq\n" "$GROUP" "$GROUP" | gmx make_ndx -quiet -n "$NDX$DOUBLE.ndx" -o "$NDX$DOUBLE.ndx" 1>make_ndx.$GROUP.out
        func.testquit $?
    done
    

    export GMX_MAXBACKUP=99        
    for GROUP in $TCGROUPS ; do
        local TCGROUPS2="$TCGROUPS2 ${GROUP}_${GROUP}_copy "
        func.testquit $?
    done
    echo
    echo "Temperature coupling group(s) of the double-layered system will be"
    echo $TCGROUPS2
    echo
    
    for GROUP in $COMGROUPS ; do
        local COMGROUPS2="$COMGROUPS2 ${GROUP}_${GROUP}_copy "
        func.testquit $?
    done
    echo "COM removal group(s) of the double-layered system will be"
    echo $COMGROUPS2
    echo

    # Prepare .mdp file for double-layered system:
    cp $MDP$SINGLE.mdp $MDP$DOUBLE.mdp

    # Expand temperature coupling and COM removal groups:
    sed -i "s/tc_grps.*/tc_grps     = $TCGROUPS2/g"  $MDP$DOUBLE.mdp
    func.testquit $?
    sed -i "s/comm_grps.*/comm_grps = $COMGROUPS2/g" $MDP$DOUBLE.mdp
    func.testquit $?
    
    # Switch on Computational Electrophysiology:
    sed -i 's/swapcoords.*/swapcoords = Z/g' $MDP$DOUBLE.mdp
    func.testquit $?
    echo "swap_frequency = $SWAPFREQ"                      >> $MDP$DOUBLE.mdp
    echo "split_group0   = $SPLITGROUP"                    >> $MDP$DOUBLE.mdp
    echo "split_group1   = ${SPLITGROUP}_copy"             >> $MDP$DOUBLE.mdp
    sed -i 's/swapcoords.*/swapcoords = Z/g' $MDP$DOUBLE.mdp
    func.testquit $?
    echo "iontypes = 2"                                    >> $MDP$DOUBLE.mdp
    echo "iontype0-name  = $CATION"                        >> $MDP$DOUBLE.mdp
    echo "iontype0-in-A  = -1"                             >> $MDP$DOUBLE.mdp # change this number
    echo "iontype0-in-B  = -1"                             >> $MDP$DOUBLE.mdp # change this number
    echo "iontype1-name  = $ANION"                         >> $MDP$DOUBLE.mdp
    echo "iontype1-in-A  = -1"                             >> $MDP$DOUBLE.mdp # change this number
    echo "iontype1-in-B  = -1"                             >> $MDP$DOUBLE.mdp # change this number
    echo "solvent-group  = ${SOLVGROUP}_${SOLVGROUP}_copy" >> $MDP$DOUBLE.mdp
    
    # Preprocess the double-layered system:
    echo "=== Building double-layered system ==="


    COMMAND="$GMX grompp -quiet -f $MDP$DOUBLE.mdp -c $GRO$DOUBLE.gro -r $GRO$DOUBLE.gro -p $TOP$DOUBLE.top -n $NDX$DOUBLE.ndx -o $GRO$DOUBLE.tpr -maxwarn 1"
    echo "RUNNING: $COMMAND"
    $COMMAND > grompp.double.out 2>&1
    func.testquit $?
    
    echo "Wrote .tpr file $GRO$DOUBLE.tpr"
    echo
}
#==============================================================================


#==============================================================================
# Please adapt the following section to your needs:
#==============================================================================
source /usr/local/bin
GMX=gmx
func.testquit $?
export GMX_NO_QUOTES=1
1
# Input files (leave away extension!)
# don't use standard output file names like 'mdout.mdp' as input files,
# they could be overwritten!
MDP=step7_production    # .mdp  (i.e. membrane.mdp, be sure to set swapcoords=no in this parameter file!)
GRO=membrane_md_0_1    # .gro  (membrane.gro)
TOP=topol    # .top  ...
NDX=newest_index    # .ndx
CATION=CAM # residue name
ANION=CLA # residue name

# Index file group names of ions that can be exchanged with solvent molecules:
IONSGROUP="Ions"
SOLVGROUP="TIP3"
# Group name of the index group that contains the compartment-partitioning atoms (i.e. the channel)
SPLITGROUP="MEMB"
SWAPFREQ=500

# Topology include files, force field files:
#export GMXLIB=./myTop/$GMXLIB

# Call the sandwich maker!
#==============================================================================
func.makeSandwich
#==============================================================================




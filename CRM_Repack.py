#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pyrosetta
pyrosetta.init()
from rosetta.protocols.analysis import InterfaceAnalyzerMover


def mutate_residue( pose , mutant_position , mutant_aa ,
        pack_radius , pack_scorefxn = '' ):
    """
    Replaces the residue at  <mutant_position>  in  <pose>  with  <mutant_aa>
        and repack any residues within  <pack_radius>  Angstroms of the mutating
        residue's center (nbr_atom) using  <pack_scorefxn>
    note: <mutant_aa>  is the single letter name for the desired ResidueType
    example:
        mutate_residue(pose, 505, H, 30)
    """

    if pose.is_fullatom() == False:
        IOError( 'mutate_residue only works with fullatom poses' )
    # create a standard scorefxn by default
    if not pack_scorefxn:
        pack_scorefxn = get_fa_scorefxn()
    task = pyrosetta.standard_packer_task( pose )
    aa_bool = pyrosetta.rosetta.utility.vector1_bool()
    mutant_aa = pyrosetta.rosetta.core.chemical.aa_from_oneletter_code( mutant_aa )
    # mutation is performed by using a PackerTask with only the mutant
    #    amino acid available during design

    for i in range( 1 , 21 ):
        # in Python, logical expression are evaluated with priority, thus the
        #    line below appends to aa_bool the truth (True or False) of the
        #    statement i == mutant_aa
        aa_bool.append( i == mutant_aa )

    task.nonconst_residue_task( mutant_position
        ).restrict_absent_canonical_aas( aa_bool )
    # prevent residues from packing by setting the per-residue "options" of
    #    the PackerTask
    center = pose.residue( mutant_position ).nbr_atom_xyz()
    for i in range( 1 , pose.total_residue() + 1 ):
        # only pack the mutating residue and any within the pack_radius
        if not (i == mutant_position or center.distance_squared(
                pose.residue( i ).nbr_atom_xyz() ) <= pack_radius**2):
            task.nonconst_residue_task( i ).prevent_repacking()
        elif i!=mutant_position:
            aa_bool = pyrosetta.rosetta.utility.vector1_bool()
            for j in range( 1 , 21 ):
                aa_bool.append( j == pose.residues[i].aa() )
            task.nonconst_residue_task( i
                        ).restrict_absent_canonical_aas( aa_bool )
    # apply the mutation and pack nearby residues
    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover( pack_scorefxn , task )
    packer.apply( pose )

def Iam_Score(_pose,chains):
    # get interface complementarity between specified chains for a pose 
    iam = InterfaceAnalyzerMover(chains)
    iam.set_pack_separated(True)
    iam.apply(_pose)
    toret = {}
    toret['sc_total'] = iam.get_all_data().sc_value
    
    return toret

"""
Repack WT and Omicron Fixed Model
"""
# set path to CRM PDBs path
PDB_path = './CRM_PDBs/'
PDBs = ['6zgi_RBD_dimer','model_01_RBD_dimer']
pack_range = 30
scores = {}

# WT 
pdb = PDBs[0]
# initialize pyrosetta for the PDB
pose_WT = pyrosetta.rosetta.core.import_pose.pose_from_file('{}{}.pdb'.format(PDB_path,pdb))
scorefxn = pyrosetta.get_fa_scorefxn()
# repack 30A from interface site 505 
chain_mut_num = pose_WT.pdb_info().pdb2pose('A', 505) 
mutate_residue( pose_WT ,  chain_mut_num, 'Y', pack_radius = pack_range , pack_scorefxn = scorefxn )
pose_WT.dump_pdb(pdb + "_repack_" + str(pack_range)+ ".pdb")
#score the interface
score_IAM = Iam_Score( pose_WT , 'A'+'_'+'B' )
scores[pdb] = score_IAM


# for  Omicron
pdb = PDBs[1]
# initialize pyrosetta for the PDB
pose_VOC = pyrosetta.rosetta.core.import_pose.pose_from_file('{}{}.pdb'.format(PDB_path,pdb))
scorefxn = pyrosetta.get_fa_scorefxn()
# repack 30A from interface site 505 
chain_mut_num = pose_VOC.pdb_info().pdb2pose('A', 505) 
mutate_residue( pose_VOC ,  chain_mut_num, 'H', pack_radius = pack_range  , pack_scorefxn = scorefxn )
pose_VOC.dump_pdb(pdb + "_repack_" + str(pack_range)+ ".pdb")
#score the interface
score_IAM = Iam_Score( pose_VOC , 'A'+'_'+'B' )
scores[pdb] = score_IAM





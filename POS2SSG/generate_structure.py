from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from irssg import pos2abr,load_one_ssg

from .deal_cell import std_cell_pos2ssg, prim_cell, extend_cell
from .generate_mag import ssg2magmom

def find_possible_ssg_numbers(structure, ssg_list, Ik_filter=None, tolerance=1e-2):
    """Find possible SSG numbers for a given structure."""
    sga = SpacegroupAnalyzer(structure, symprec=tolerance)
    sgnum = sga.get_space_group_number()
    
    print(f'Atomic space group: {sga.get_space_group_symbol()} ({sgnum})')
    print('Possible SSG numbers:')
    
    found_ssg = []
    
    if Ik_filter is not None:
        # Filter by specific Ik value
        for ssg in ssg_list:
            ssgnum = ssg['ssgNum']
            Gnum = int(ssgnum.split('.')[0])
            Ik = int(ssgnum.split('.')[1])
            if Gnum == sgnum and Ik == Ik_filter:
                print(f'  {ssgnum}')
                found_ssg.append(ssgnum)
    else:
        # Show all possible Ik values
        for Ik in range(1, 13):
            ik_found = False
            for ssg in ssg_list:
                ssgnum = ssg['ssgNum']
                Gnum = int(ssgnum.split('.')[0])
                Ik_val = int(ssgnum.split('.')[1])
                if Gnum == sgnum and Ik_val == Ik:
                    if not ik_found:
                        print(f'  Ik = {Ik}:')
                        ik_found = True
                    print(f'    {ssgnum}')
                    found_ssg.append(ssgnum)
    
    if not found_ssg:
        print('  No matching SSG found for this structure.')
    
    return found_ssg

def generate_ssg_structure(ssgnum_in, cell, magmom_dict, ssg_list, tolerance=1e-2, verbose=False):
    """Generate SSG structure with magnetic moments."""
    print(f'Target spin space group: {ssgnum_in}')
    
    # Find the target SSG in the database
    ssg_found = False
    target_ssg = None
    
    for ssg in ssg_list:
        if ssg['ssgNum'] == ssgnum_in:
            ssg_found = True
            target_ssg = ssg
            break
    
    if not ssg_found:
        print(f'Error: SSG number {ssgnum_in} not found in the database.')
        return False
    
    # Load SSG data
    ssgdic = load_one_ssg(ssgnum_in)
    superCell = ssgdic['superCell']
    
    if verbose:
        print(f"SSG dimension: {ssgdic['dim']}")
        print(f"Supercell matrix: {superCell}")
    
    # Process the cell
    cell_std = std_cell_pos2ssg(cell, tol=tolerance)
    cell_prim = prim_cell(cell_std, tol=tolerance)
    cell_pos2abr, shift = pos2abr(cell_prim)
    
    # Extend cell considering supercell matrix
    cell_sc = extend_cell(cell_pos2abr, superCell)
    
    # Generate magnetic moments
    result = ssg2magmom(ssgdic, cell_sc, magmom_dict, tol=tolerance)
    
    return result
    
    
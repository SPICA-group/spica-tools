#!/usr/bin/env python3
"""
gen_gns_input: Generate GENESIS input file for SPICA/pSPICA simulations

This module reads box size from a PDB file and generates a GENESIS control file
for NPT ensemble MD simulations.
"""

import argparse
import sys
import os


def read_pdb_box(pdb_file):
    """
    Read box size information from CRYST1 line in PDB file.
    
    Parameters
    ----------
    pdb_file : str
        Path to PDB file
        
    Returns
    -------
    tuple
        (box_x, box_y, box_z) in Angstroms, or (None, None, None) if not found
    """
    box_x = box_y = box_z = None
    
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('CRYST1'):
                    items = line.split()
                    if len(items) >= 4:
                        box_x = float(items[1])
                        box_y = float(items[2])
                        box_z = float(items[3])
                        break
    except FileNotFoundError:
        print(f"ERROR: PDB file not found: {pdb_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to read PDB file: {e}", file=sys.stderr)
        sys.exit(1)
    
    return box_x, box_y, box_z


def check_enm_files(toppar_dir='toppar'):
    """
    Check which ENM files exist in the toppar directory.
    
    Parameters
    ----------
    toppar_dir : str
        Directory containing topology and parameter files
        
    Returns
    -------
    list
        List of existing ENM file paths
    """
    enm_files = [
        'enm_bond_index.ndx',
        'enm_bond_parm.prm',
        'enm_angle_index.ndx',
        'enm_angle_parm.prm'
    ]
    
    existing_files = []
    for fname in enm_files:
        fpath = os.path.join(toppar_dir, fname)
        if os.path.exists(fpath):
            existing_files.append(os.path.join(toppar_dir, fname))
    
    return existing_files


def calculate_domain_candidates(box_size, min_cell_size=100.0):
    """
    Calculate reasonable domain decomposition numbers.
    
    Domain number must be:
    - At least 2
    - Multiples of 2 or 3 (divisible by 2 or 3)
    - box_size / domain >= min_cell_size
    
    Parameters
    ----------
    box_size : float
        Box size in one dimension (Angstroms)
    min_cell_size : float
        Minimum cell size (default: 100.0 Angstroms)
        
    Returns
    -------
    list
        List of reasonable domain numbers
    """
    max_domain = int(box_size / min_cell_size)
    
    if max_domain < 2:
        return [2]
    
    candidates = []
    
    # Include all numbers from 2 to max_domain that are divisible by 2 or 3
    for n in range(2, max_domain + 1):
        if n % 2 == 0 or n % 3 == 0:
            candidates.append(n)
    
    return candidates if candidates else [2]


def generate_genesis_input(pdb_file, output_file='npt.inp', pspica=False, toppar_dir='toppar', scale_factors=None):
    """
    Generate GENESIS control file for NPT MD simulation.
    
    Parameters
    ----------
    pdb_file : str
        Input PDB file containing coordinates and box information
    output_file : str
        Output control file name (default: npt.inp)
    pspica : bool
        Use pSPICA parameters (default: False)
    toppar_dir : str
        Directory containing topology and parameter files (default: toppar)
    scale_factors : list of float or None
        Box size scaling factors. None (default), [scale] for uniform scaling,
        or [scale_x, scale_y, scale_z] for per-dimension scaling
    """
    # Read box dimensions from PDB
    box_x0, box_y0, box_z0 = read_pdb_box(pdb_file)
    
    if box_x0 is None or box_y0 is None or box_z0 is None:
        print("ERROR: No CRYST LINE exist", file=sys.stderr)
        print("       CRYST1 line is required in the PDB file.", file=sys.stderr)
        sys.exit(1)
    
    # Process scaling factors
    if scale_factors is None:
        # Default: no scaling
        scale_x = scale_y = scale_z = 1.0
    elif len(scale_factors) == 1:
        # Single value: apply to all dimensions
        scale_x = scale_y = scale_z = scale_factors[0]
    elif len(scale_factors) == 3:
        # Three values: apply to x, y, z respectively
        scale_x, scale_y, scale_z = scale_factors
    else:
        print(f"ERROR: Invalid number of scaling factors: {len(scale_factors)}", file=sys.stderr)
        print("       Provide either 1 value (uniform) or 3 values (x, y, z)", file=sys.stderr)
        sys.exit(1)
    
    # Calculate scaled box sizes
    box_x = box_x0 * scale_x
    box_y = box_y0 * scale_y
    box_z = box_z0 * scale_z
    
    print(f"Box size: {box_x:.2f} x {box_y:.2f} x {box_z:.2f} Angstrom")

    # Calculate domain decomposition candidates
    domain_x_candidates = calculate_domain_candidates(box_x)
    domain_y_candidates = calculate_domain_candidates(box_y)
    domain_z_candidates = calculate_domain_candidates(box_z)
    
    # Default to minimum (2)
    domain_x = 2
    domain_y = 2
    domain_z = 2
    
   
    # Check for ENM files
    enm_files = check_enm_files(toppar_dir)
    
    # Build parfile line
    parfile_list = [f'{toppar_dir}/par_spica_charmm.prm']
    parfile_list.extend(enm_files)
    parfile_line = ', '.join(parfile_list)
    
    if enm_files:
        print(f"\nFound ENM files: {len(enm_files)}")
        for f in enm_files:
            print(f"  - {f}")
    
    # Set pSPICA-specific parameters
    if pspica:
        pme_max_spacing = 2.0
        fast_water = 'YES'
        print("\npSPICA mode: pme_max_spacing=2.0, fast_water=YES")
    else:
        pme_max_spacing = 5.0
        fast_water = 'NO'
        print("\nSPICA mode: pme_max_spacing=5.0, fast_water=NO")
    
    # Format domain candidate strings
    domain_x_str = ', '.join(map(str, domain_x_candidates))
    domain_y_str = ', '.join(map(str, domain_y_candidates))
    domain_z_str = ', '.join(map(str, domain_z_candidates))
    
    # Generate control file content
    content = f"""[INPUT]
topfile = {toppar_dir}/top_spica_charmm.rtf
parfile = {parfile_line}
psffile = gns_spica.psf
pdbfile = gns_spica.pdb
reffile = gns_spica.pdb
# rstfile = equilibrated.rst
 
[OUTPUT]
rstfile = mdout.rst
dcdfile = mdout.dcd

[ENERGY]
forcefield      = SPICA        # [CHARMM]
electrostatic   = PME          # [CUTOFF,PME]
switchdist      = 12.0         # switch distance
cutoffdist      = 15.0         # cutoff distance
pairlistdist    = 18.0         # pair-list distance
pme_nspline     = 4
pme_max_spacing = {pme_max_spacing}
#contact_check   = YES          # avoid atomic crash
                               # minimization or early equilibration
#nonb_limiter    = YES          # minimization or early equilibration
                               # for large ENM system
#minimum_contact = 0.5          # in the case of nonb_limiter = YES

# contact_check, nonb_limiter, minimum_contact can be used for initial unstable structure.

[DYNAMICS]
integrator      = VVER         # [VVER, VRES]
timestep        = 0.01         # timestep (ps)
nsteps          = 10000        # short test
crdout_period   = 10000
eneout_period   = 1000         # energy output period
rstout_period   = 10000      
nbupdate_period = 10 

# In the case of MTS, you can write 
# integrator        = VRES   
# elec_long_period  = 5
# thermostat_period = 10
# barostat_period   = 10
 
[CONSTRAINTS]
rigid_bond      = NO           # constraints all bonds involving hydrogen
fast_water      = {fast_water}

[ENSEMBLE]
ensemble        = NPT          # [NVE,NVT,NPT]
tpcontrol       = BUSSI        # thermostat and barostat
temperature     = 310          # K
pressure        = 1.0          # atm
isotropy        = ISO          # SEMI-ISO for the membrane case

[BOUNDARY]
type            = PBC           # [PBC]
box_size_x      = {box_x:.3f}
box_size_y      = {box_y:.3f}
box_size_z      = {box_z:.3f}
#domain_x = {domain_x}
#domain_y = {domain_y}
#domain_z = {domain_z}
"""
    
    # Write to file
    try:
        with open(output_file, 'w') as f:
            f.write(content)
        print(f"\nGENESIS control file generated: {output_file}")
    except Exception as e:
        print(f"ERROR: Failed to write output file: {e}", file=sys.stderr)
        sys.exit(1)


def get_option_script(args):
    """Parse command-line arguments for cg_spica integration"""
    parser = argparse.ArgumentParser(
        description='Generate GENESIS input file for SPICA/pSPICA CG simulations',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('-pdb', '--pdb', required=True,
                        help='Input PDB file with box information (CRYST1 line)')
    parser.add_argument('-o', '--output', default='npt.inp',
                        help='Output control file name (default: npt.inp)')
    parser.add_argument('-pspica', '--pspica', action='store_true',
                        help='Use pSPICA parameters (pme_max_spacing=2, fast_water=YES)')
    parser.add_argument('--toppar_dir', default='toppar',
                        help='Directory containing topology/parameter files (default: toppar)')
    parser.add_argument('-scale', '--scale', nargs='+', type=float, metavar='SCALE',
                        help='Box size scaling factors. One value applies to all dimensions (e.g., -scale 1.02), '
                             'three values apply to x, y, z respectively (e.g., -scale 1.02 1.02 1.02). '
                             'Default: no scaling (1.0 1.0 1.0)')
    
    return parser.parse_args(args)


def run(args):
    """Main execution function for cg_spica integration"""
    if not os.path.exists(args.pdb):
        print(f"ERROR: PDB file not found: {args.pdb}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.toppar_dir):
        print(f"WARNING: toppar directory not found: {args.toppar_dir}", file=sys.stderr)
        print(f"         Please make sure topology files exist before running GENESIS.", file=sys.stderr)
    
    generate_genesis_input(args.pdb, args.output, args.pspica, args.toppar_dir, args.scale)
    

if __name__ == '__main__':
    args = get_option_script(sys.argv[1:])
    run(args)

import os
import shutil
import argparse
import numpy as np
import MDAnalysis as mda
from scipy import optimize
from MDAnalysis.analysis import align
from MDAnalysis.transformations import rotate
from string import ascii_uppercase as auc
from string import ascii_lowercase as alc

#!/usr/bin/env python3
'''
TODO:
- cite pdbs in readme
- figure out how to install conda from this file
 '''

# define a function to convert string to boolean
def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# define a function to return radius given a set of x, y and z coord
# and center coordinates i.e. xc, yc and zc
def calc_R(xc, yc):
    """Calculate radius given x, y coordinates and center coordinates."""
    return np.sqrt((x-xc)**2 + (y-yc)**2)

# Function to minimize or optimize
# calculate the algebraic distance between the data points and the 
# mean circle centered at point c=(xc, yc, zc) 
def func(c):
    """Calculate algebraic distance between data points and mean circle."""
    Ri = calc_R(*c)
    return Ri - Ri.mean()

def build_circle(topdir, no_subunits, in_structure, z_rotation, xy_rotation, tilt, in_radius, delete_intermediates):
    """
    Main function to build symmetric rings of circular protein assemblies.

    Parameters:
        topdir (str): Top directory containing input files.
        no_subunits (int): Number of subunits to assemble in the ring.
        in_structure (str): Input PDB file with high resolution monomers.
        z_rotation (float): Angle (deg) to rotate each subunit of the ring around its z-axis.
        xy_rotation (float): Angle (deg) to rotate each subunit of the ring around its xy-axis.
        in_radius (float or None): Input radius for the ring. If None, the radius is fitted from the input structure.
        delete_intermediates (bool): Whether to delete rotated individual subunit structures after building the ring.
    """
    ring = mda.Universe(os.path.join(topdir, in_structure))
    # Create intermediates directory if it doesn't exist
    intermediates_dir = os.path.join(topdir, "intermediate_structures")
    os.makedirs(intermediates_dir, exist_ok=True)
    
    # mdanalysis is sometimes a bit difficult with atom seletions
    # we therefore have to write out one subunit first, to then load it in again as the monomer
    # this will take the first subunit of the input structure
    # note that this assumes the first subunit is a high resolution monomer itself and the ring will be built from
    subunit = ring.segments[0].atoms
    subunit.write(os.path.join(intermediates_dir,'monomer.pdb'))
    mono = mda.Universe(os.path.join(intermediates_dir, 'monomer.pdb'))
    mono_a = mono.select_atoms('protein')

    if in_radius is None:
        print("Fitting the radius of the ring from the input structure.")
        # Calculate radius from low-resolution fitted ring structure
        global x, y
        x = ring.atoms.positions[:,0]
        y = ring.atoms.positions[:,1]
        center_estimate = np.mean(x), np.mean(y)
        final_center, _ = optimize.leastsq(func, center_estimate)
        Ri_f = calc_R(*final_center)
        radius = Ri_f.mean()
        print(f"The circular assembly has a radius of about {round(radius)} Angstrom.")
    else:
        # Use a fixed radius for the ring
        radius = in_radius
        print(f"Using a fixed radius of {radius} Angstrom for the circular assembly.")   

    segid_list = auc + alc  # Create a list of segment IDs using uppercase and lowercase letters
                            # Note: this assumes no more than 52 subunits
    segid_list = iter(segid_list)  # Convert to an iterator for sequential access 

    if no_subunits > len(auc + alc):
        raise ValueError(f"Number of subunits ({no_subunits}) exceeds the number of available segment IDs ({len(auc + alc)}). Please reduce the number of subunits or extend the segment ID scheme.")

    for idx in range(no_subunits):
        ring_part = ring.segments[0].atoms
        mono_a = mono.select_atoms('all')
        align.alignto(mono_a, ring_part, select=('name CA'))

        # Rotate the subunit around the z-axis
        angle = 360/no_subunits * idx
        mono_a = rotate.rotateby(angle=(angle + z_rotation), direction=[0,0,1], ag=mono_a)(mono_a)

        # Translate the subunit to the correct position in the ring
        # shift the subunit, to center it at 0,0
        to_zero = -np.concatenate((mono_a.positions[:, :2].mean(axis=0), np.zeros(1)))
        mono_a.atoms.translate(to_zero)
        # now offset it to its location in the ring
        mono_a = mono_a.select_atoms("all")
        angle *= np.pi/180
        mono_a.atoms.translate([radius*np.cos(angle), radius*np.sin(angle), 0])
        
        # Rotate the subunit around axis that connects it to the center
        second_axis = -mono_a.positions[:, :2].mean(axis=0)
        second_axis /= np.linalg.norm(second_axis)
        second_axis = np.concatenate((second_axis, np.zeros(1)))
        mono_a = rotate.rotateby(angle=xy_rotation, direction=second_axis, ag=mono_a)(mono_a)

        # Tilt the subunit in and out of the plane of the ring (this is perpendicular to the second axis)
        tilt_axis = np.array([-second_axis[1], second_axis[0], 0])
        tilt_axis /= np.linalg.norm(tilt_axis)
        mono_a = rotate.rotateby(angle=tilt, direction=tilt_axis, ag=mono_a)(mono_a)

        # write out temporary copy of subunit
        chain_id = next(segid_list)
        mono_a.segments.segids = chain_id
        mono_a.write(os.path.join(intermediates_dir, f"{idx}subunit.pdb"))

    # Assemble the circle
    for monomer in range(0, no_subunits):
        mono = mda.Universe(os.path.join(intermediates_dir, f'{monomer}subunit.pdb'))
        if monomer == 0:
            system = mono.select_atoms("all")
        else:
            system = mda.Merge(system.select_atoms("all"), mono.select_atoms("all"))
    
    system = system.select_atoms("all")
    system.write(os.path.join(topdir, f"{no_subunits}mer_rad{round(radius,2)}A_z{z_rotation}deg_xy{xy_rotation}deg_tilt{tilt}deg.pdb"))
    if delete_intermediates:
        if os.path.exists(intermediates_dir):
            print("Deleting individual subunit structures.")
        shutil.rmtree(intermediates_dir)


def main():
    parser = argparse.ArgumentParser(description='Build circular assemblies of proteins')
    parser.add_argument('--in_structure', required=True, type=str,
                        help='PDB file with high resolution monomers.' \
                        'Expects a crude ring structure from which the radius can be calculated.' \
                        'Alternatatively, a single subunit works with a user decided input radius.') 
    parser.add_argument('--topdir', required=False, default='.',
                        help='Top directory containing input files. Defaults to current directory.')
    parser.add_argument('--no_subunits', type=int, required=True,
                        help='Desired number of intermediates in the ring. Must be an integer greater than 0.')
    parser.add_argument('--z_rotation', type=float, default=0,
                        help='Angle (deg) to rotate each subunit of the ring around its z-axis.')
    parser.add_argument('--xy_rotation', type=float, default=0,
                        help='Angle (deg) to rotate each subunit of the ring around its xy-axis.')
    parser.add_argument('--tilt', type=float, default=0,
                        help='Angle (deg) to tilt each subunit into the plane of the ring.')
    parser.add_argument('--in_radius', type=float, required=False, default=None,
                        help='Input radius for the ring. Disables fitting the radius to the original structure.')
    parser.add_argument('--delete_intermediates', type=str2bool, default=True,
                        help='Delete rotated individual subunit structures after building the ring. Defaults to True.')
    args = parser.parse_args()
    
    build_circle(args.topdir, args.no_subunits, args.in_structure,
                 args.z_rotation, args.xy_rotation, args.tilt,
                 args.in_radius, args.delete_intermediates
                )
    
if __name__ == "__main__":
    main()
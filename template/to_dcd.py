#!/usr/bin/env python

"""This is a simple script to convert an ensemble of output models
   (multiple single-frame RMF files) into DCD (CHARMM/NAMD) trajectories,
   with each bead represented as a single 'atom', in residue number order
   (the same order as in the mmCIF file).

   We use the updated version of MDTools that's bundled with Chimera, so you
   may need to change sys.path accordingly, below.
"""

from __future__ import print_function, division
import glob
import os
import re
import sys
sys.path.append('/opt/chimera-1.11-1.fc24/share/Trajectory/DCD/MDToolsMarch97/')
import md
import RMF
import math
import IMP.atom
import IMP.rmf
import IMP.pmi1.output

class DCDOutput(object):
    """Dump a series of scaffold coordinates from RMFs to a DCD file."""

    def __init__(self, fname, comps, bead_resrange):
        self.fname = fname
        self.comps = comps
        self.bead_resrange = {}
        for comp in comps:
            self.bead_resrange[comp] = bead_resrange[comp]
        n_coords = sum(len(b) for b in bead_resrange.values())
        self._init_dcd(n_coords)
    
    def _count_coords(self, coords):
        return sum(len(b) for b in coords.values())

    def _get_resrange(self, coords):
        for c in coords:
            if c[5] is not None: # return resrange if available
                yield c[5][0], c[5][-1]
            else:
                yield c[4], c[4]

    def dump(self, coords):
        """Dump a single set of coordinates to the DCD file"""
        len_coords = self._count_coords(coords)
        for comp in self.comps:
            resrange = list(self._get_resrange(coords[comp]))
            if resrange != self.bead_resrange[comp]:
                raise ValueError("Residue ranges mismatch for component %s; "
                                 "mmCIF: %s; RMF: %s"
                                 % (comp, self.bead_resrange[comp], resrange))
        assert(len_coords == len(self._ag.atoms))
        for atom, coord in zip(self._ag.atoms,
                               self._get_coords(coords)):
            atom.x, atom.y, atom.z = coord
        self._d.append()

    def _init_dcd(self, n_coords):
        self._ag = md.AtomGroup()
        for i in range(n_coords):
            self._ag.atoms.append(md.Atom())
        self._d = md.DCDWrite(self.fname, self._ag)

    def _get_coords(self, coords):
        """Get all bead coordinates in the same order as in the mmCIF file"""
        for comp in self.comps:
            for coord in coords[comp]:
                # Everything should be non-atomic for now
                assert(coord[1] is None)
                yield coord[0]


class RMFReader(object):
    """Read scaffold coordinates from a set of RMF files"""

    def __init__(self, rmf_dir):
        self.rmf_dir = rmf_dir

    def get_all_rmfs(self):
        # Assume files are named <int>-<int>.rmf3
        namere = re.compile('(\d+)\-(\d+)\.rmf3')
        rmf_by_num = {}
        for f in os.listdir(self.rmf_dir):
            m = namere.match(f)
            if m:
                rmf_by_num[int(m.group(1))] = f
        for num in sorted(rmf_by_num.keys()):
            yield os.path.join(self.rmf_dir, rmf_by_num[num])

    def read(self, rmf, comps):
        """Read coordinates from the RMF"""
        rmf.set_current_frame(RMF.FrameID(0))
        m = IMP.Model()
        mhs = IMP.rmf.create_hierarchies(rmf, m)
        m.update()
        assert(len(mhs) == 1)
        return self._get_coords(m, mhs, comps)

    def _get_coords(self, m, mhs, comps):
        name = 'dcd-output'
        o = IMP.pmi1.output.Output(atomistic=True)
        o.dictionary_pdbs[name] = mhs[0]
        o._init_dictchain(name, mhs[0], multichar_chain=True)

        coords = {}
        for mol in mhs[0].get_children():
            comp = mol.get_name()
            if comp in comps:
                # Force get_particle_infos to only look at this one chain
                o.dictionary_pdbs[name] = mol
                c, _ = o.get_particle_infos_for_pdb_writing(name)
                coords[comp] = self._exclude_coordinates(c, comp)
        return coords

    def _exclude_coordinates(self, coords, component):
        def not_excluded(coord, seqrange):
            all_indexes = coord[5] if coord[5] else [coord[4], coord[4]]
            return all_indexes[0] < seqrange[0] or all_indexes[-1] > seqrange[1]

        excludes = {'Mlp1': (717, 1875),
                    'Mlp2': (691,1679),
                    'Nup42': (364,430),
                    'Gle1': (121,538)}
        exclude = excludes.get(component.split('@')[0], None)
        if exclude is None:
            return coords
        if exclude:
            return [c for c in coords if not_excluded(c, exclude)]


class CifParser(object):
    """Read an mmCIF file.
       Read reference coordinates plus the order of components from the mmCIF
       file. This is a pretty simple parser, and will need to be updated
       if the formatting of the mmCIF file changes."""

    def _parse_struct_asym(self, fh):
        """Get a mapping from asym_ids to IMP component names"""
        asym_from_comp = {}
        for line in fh:
            if line.startswith('_'): continue
            if line.startswith('#'):
                return asym_from_comp
            asym_id, entity_id, imp_name = line.rstrip('\r\n').split()
            asym_from_comp[imp_name] = asym_id

    def _parse_obj_site(self, fh, asym_from_comp, model_num):
        """Get correct component ordering and bead residue ranges"""
        line_end = ' %s\n' % model_num
        comps = []
        beads = {}
        current_asym = None
        comp_from_asym = {}
        for comp, asym in asym_from_comp.items():
            comp_from_asym[asym] = comp
        for line in fh:
            if line.startswith('_'): continue
            if line.startswith('#'):
                return comps, beads
            if line.endswith(line_end):
                sphere_obj = line.split()
                asym_id = sphere_obj[4]
                comp = comp_from_asym[asym_id]
                if asym_id != current_asym:
                    current_asym = asym_id
                    comps.append(comp)
                    beads[comp] = []
                beads[comp].append((int(sphere_obj[2]), int(sphere_obj[3])))
        return comps, beads

    def parse(self, fh, model_num='1'):
        asym_from_comp = {}
        for line in fh:
            if line == '_struct_asym.id\n':
                asym_from_comp = self._parse_struct_asym(fh)
            elif line == '_ihm_sphere_obj_site.ordinal_id\n':
                return self._parse_obj_site(fh, asym_from_comp, model_num)

def parse_args():
    if len(sys.argv) != 4:
        print("Usage: %s path-to-nup82.cif rmf_dir dcd-fname\n\n"
              % sys.argv[0], file=sys.stderr)
        # e.g. ./to_dcd.py template/nup82.cif /netapp/sali/kimsj/NPC/nup82-master/chef13_final_refinement/analysis_allEM2D_chef9/kmeans_10000_2/cluster.0_370 cluster0.dcd
        sys.exit(1)
    return sys.argv[1:]

def main():
    cif_fname, rmf_dir, dcd_fname = parse_args()

    c = CifParser()
    with open(cif_fname) as fh:
        comps, beads = c.parse(fh)

    d = DCDOutput(dcd_fname, comps, beads)

    reader = RMFReader(rmf_dir)
    for num, rmf in enumerate(reader.get_all_rmfs()):
        print("Handling RMF %d, %s" % (num + 1, rmf))
        r = RMF.open_rmf_file_read_only(rmf)
        if r.get_number_of_frames() == 0:
            print("    Skipped empty RMF file\n")
            continue
        coords = reader.read(r, comps)
        d.dump(coords)

if __name__ == '__main__':
    main()

import nanome
import tempfile
import os
import sys
import math

from nanome.api.ui import Menu
from nanome.util import async_callback, Logs, Process, Vector3
from nanome.api.structure import Complex
from nanome._internal._structure import _Bond

BASE_DIR = os.path.join(os.path.dirname(__file__))

class HAADPlugin(nanome.AsyncPluginInstance):

    def start(self):
        self.menu = Menu()
        self.menu.title = 'HAAD Plugin'
        self.menu.width = 1
        self.menu.height = 1

        self.temp_dir = tempfile.TemporaryDirectory()
        self.input_file = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb', dir=self.temp_dir.name)

        msg = 'Hello Nanome!'
        node = self.menu.root.create_child_node()
        self.label = node.add_new_label(msg)
        Logs.message(msg)

    @async_callback
    async def on_run(self):
        complex_indices = [comp.index for comp in await self.request_complex_list()]
        deep = await self.request_complexes(complex_indices)
        result = await self.add_hydrogens(complexes=deep)
        await self.update_structures_deep(result)

    @async_callback
    async def add_hydrogens(self, request=None, complexes=None):
        if request:
            complexes = request.get_args()

        for complex in complexes:
            complex.io.to_pdb(self.input_file.name)

            # remember atoms by position
            atom_by_position = dict()
            for atom in complex.atoms:
                atom_by_position[get_position_key(atom)] = atom

            # compute all hydrogens
            result_complex, hydrogen_ids = await call_HAAD(self.input_file.name)
            if not result_complex:
                continue

            # add hydrogens to original complex
            self.match_and_update(atom_by_position, result_complex, hydrogen_ids)

        if request:
            request.send_response(complexes)

        return complexes


    def match_and_update(self, atom_by_position, result_complex, hydrogen_ids):
        """
        Match and add hydrogens to original complex using positions to match.

        :param atom_by_position: dict mapping position key to atom in source complex
        :type atom_by_position: dict
        :param result_complex: Output complex from hydrogens calculation
        :type result_complex: :class:`nanome.structure.Complex`
        :param hydrogen_ids: List of hydrogen atom ids to use
        :type hydrogen_ids: list of int
        """
        atoms = list(result_complex.atoms)
        for id in hydrogen_ids:
            h_atom = atoms[id]
            h_atom.symbol = "H"
            
            bonded_atom, dist = self.get_closest_heavy_atom_in_residue(result_complex, h_atom)

            if not bonded_atom:
                continue
            
            if dist > 2.0:
                continue


            bond = self.add_bond(h_atom, bonded_atom)

            bonded_atom_key = get_position_key(bonded_atom)

            if bonded_atom_key is None or bonded_atom_key not in atom_by_position:
                Logs.warning(f'H {h_atom.serial} bonded with unknown atom {bonded_atom.symbol} at {bonded_atom.position}')
                continue

            source_atom = atom_by_position[bonded_atom_key]

            new_atom = h_atom._shallow_copy()
            new_atom._display_mode = source_atom._display_mode
            new_atom.is_het = source_atom.is_het
            new_atom.selected = source_atom.selected

            # copy bfactor and occupancy for surface coloring
            new_atom.bfactor = source_atom.bfactor
            new_atom.occupancy = source_atom.occupancy

            new_bond = bond._shallow_copy()
            new_bond.atom1 = source_atom
            new_bond.atom2 = new_atom

            residue = source_atom.residue
            residue.add_atom(new_atom)
            residue.add_bond(new_bond)

    def add_bond(self, atom1, atom2):
        new_bond = _Bond._create()
        new_bond._kind = nanome.util.enums.Kind.CovalentSingle
        new_bond._atom1 = atom1
        new_bond._atom2 = atom2
        atom2.residue._add_bond(new_bond)
        return new_bond

    def get_closest_heavy_atom_in_residue(self, complex, atom):
        result = None
        minD = 10000.0
        for a in atom.residue.atoms:
            if a != atom and a.symbol != "H":
                d = Vector3.distance(a.position, atom.position)
                if d < minD:
                    result = a
                    minD = d
        return result, minD

def get_position_key(atom):
    """
    Get atom position as tuple of integers to be used as lookup key.
    Rounds the coordinates to 4 decimal places before multiplying by 50
    to get unique integer-space coordinates, and avoid floating point errors.

    :param atom: Atom to get key for
    :type atom: :class:`nanome.structure.Atom`
    :return: Position key tuple
    :rtype: (int, int, int)
    """
    return tuple(map(lambda x: int(50 * round(x, 4)), atom.position))


HAAD_PATH = ""
if sys.platform == 'linux':
    HAAD_PATH = os.path.join(BASE_DIR, 'bin/linux/haad')
elif sys.platform == 'darwin':
    HAAD_PATH = os.path.join(BASE_DIR, 'bin/darwin/haad')
elif sys.platform == 'win32':
    HAAD_PATH = os.path.join(BASE_DIR, 'bin/win32/haad.exe')

async def call_HAAD(pdb_path):

    output_pdb = pdb_path + ".h"

    p = Process()
    p.executable_path = HAAD_PATH
    p.args = [pdb_path]
    p.on_error = Logs.error
    p.on_output = Logs.debug

    # if error, notify user and return
    exit_code = await p.start()
    if exit_code:
        self.send_notification(enums.NotificationTypes.error, 'Error computing hydrogens, haad failed')
        return
    
    fix_haad_chains(output_pdb)

    hydrogenated = Complex.io.from_pdb(path=output_pdb)
    new_hydrogens = []
    id = 0
    with open(output_pdb) as f:
        lines = f.readlines()
        for line in lines:
            if len(line) > 56 and line[56] != " " and int(line[56]) < 2:
                new_hydrogens.append(id)
            id+=1

    return (hydrogenated, new_hydrogens)

def fix_haad_chains(path):
    id = 0
    with open(path) as f:
        lines = f.readlines()
        output_lines = lines[:]
        for line in lines:
            if len(line) > 56 and line[56] != " " and int(line[56]) < 2:
                if id > 0:
                    #search in previous lines
                    for i in range(id - 1, -1, -1):
                        if len(lines[i]) > 22:
                            chain = lines[i][21]
                            if chain != " ":
                                l = lines[i]
                                break
                    #Same residue
                    if line[17:17+4] == l[17:17+4] and line[23:23+4] == l[23:23+4]:
                        output_lines[id] = output_lines[id][:21] + l[21] + output_lines[id][22:-2] + "1.00 1.00            H \n"
            id+=1

    with open(path, "w") as f:
        for l in output_lines:
            f.write(l)
    
    print(path)


def main():
    plugin = nanome.Plugin('HAAD Plugin', 'Hydrogenation plugin using HAAD', 'other', False)
    plugin.set_plugin_class(HAADPlugin)
    plugin.run()


if __name__ == '__main__':
    main()

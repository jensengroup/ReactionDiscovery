import copy
import numpy as np

from rdkit import Chem

from .xyz2mol_local import xyz2mol, int_atom, AC2mol, get_proto_mol

class Compound():
    """Compound object """

    def __init__(self, name=None, symbols=None, coordinates=None, 
                 charge=None, multiplicity=None, charged_fragments=None):

        self.name = name

        # information about compund
        if symbols is not None:
            self.natoms = len(symbols)

        self.atomic_symbols = symbols
        self.multiplicity = multiplicity
        self.charge = charge
        self.charged_fragments = charged_fragments

        self.coordinates = coordinates
        
        # placeholders for ase and rdkit objects
        self.rdkit_mol = None
        self.ase_mol = None

        # placeholder for calculated results
        self.properties = {}

    @classmethod
    def from_xyz(cls, filename, name="mol", charge=0, multiplicity=1, 
            charged_fragments=True):
        """Initialize ´Compound´from xyz file  """

        with open(filename, 'r') as f:
            natoms = int(f.readline())

            atom_coords = np.zeros((natoms,3))
            atom_symbols = []

            # skip header n.b. set name form this line.
            f.readline() 

            # Read each line, and use natoms to avoid reading beyond last line.
            for line_num, line in enumerate(f):
                 
                if line_num == natoms:
                    break

                line_values = line.split()

                atom_symbol = line_values[0]
                atom_pos = [float(x) for x in line_values[1:]]
                
                if len(atom_pos) == 3:
                    atom_coords[line_num] = atom_pos
                    atom_symbols.append(atom_symbol)

        obj = Compound(name=name, symbols=atom_symbols, coordinates=atom_coords,
                       charge=charge, multiplicity=multiplicity,
                       charged_fragments=charged_fragments)

        return obj

    @classmethod
    def from_gaussian(cls, filename, name='mol', charge=0, multiplicity=1,
            charged_fragments=True):
        """ """
        pt = Chem.GetPeriodicTable() # Periodic Table info
        
        with open(filename, 'r') as f:
            content = f.read()
        
        lines = content.split('Standard orientation')[-1].split('\n')

        # First 5 lines are header
        del lines[:5]

        atom_symbols = []
        atom_coords = []
        for line in lines:
            line = line.strip()
            
            #if only - in line it is the end
            if set(line).issubset(set(['-', ' '])):
                break

            tmp_line = line.split()
            if not len(tmp_line) == 6:
                raise RuntimeError('Length of line does not match structure!')

            # read atoms and positions:
            try:
                atom_symbol = pt.GetElementSymbol(int(tmp_line[1]))
                atom_position = list(map(float, tmp_line[3:]))
            except:
                raise ValueError('Expected a line with three integers and three floats.')
            
            atom_symbols.append(atom_symbol)
            atom_coords.append(atom_position)

        obj = Compound(name=name, symbols=atom_symbols, coordinates=atom_coords,
                charge=charge, multiplicity=multiplicity,
                charged_fragments=charged_fragments)
        
        return obj

    @classmethod
    def from_xtbout(cls, filename, name='mol', charge=0, multiplicity=1,
            charged_fragments=True):
        """ Read output from xTB 6.3.2 """

        with open(filename, 'r') as f:
            content = f.read()

        temp_items = content.split('final structure:')[1:]

        atom_coords = []
        atom_symbols = []
        for item_i in temp_items:
            lines = [ line for line in item_i.split('\n') if len(line) > 0]
            del lines[:3] # first 2 lines are header lines

            for line in lines:
                line = line.strip()
                # is line is empty - mol block ends.
                if set(line).issubset(set([' '])):
                    break

                tmp_line = line.split()
                if not len(tmp_line) == 4:
                    raise RuntimeError('Length of line does not match structure!')
            
                # read atoms and positions
                try:
                    atom_symbols.append(str(tmp_line[0]))
                    atom_position = list(map(float, tmp_line[1:]))
                    #atom_position = [x*0.529177249 for x in atom_position] # bohr to AA.
                    atom_coords.append(atom_position)
                except: 
                    raise ValueError('Expected a line with one string and three floats.')
        
        obj = Compound(name=name, symbols=atom_symbols, coordinates=atom_coords,
                charge=charge, multiplicity=multiplicity,
                charged_fragments=charged_fragments)
            
        return obj


    @classmethod
    def from_rdkit(cls, rdkit_mol, name='mol', charge=0, multiplicity=1,
                    charged_fragments=True, rdkit_conf=True):
        """Intialize ´Compound´ from RDKit object
        
        rdkit_mol : RDKit Mol object
        rdkit_conf : bool
            if True uses coordinates from RDKit conformer (-1).
        """
        if rdkit_conf:
            try:
                conf = rdkit_mol.GetConformer()
                coords = conf.GetPositions()
            except:
                RuntimeError("RDKit doesn't have a conformer. Switch to rdkit_conf=False.")
        else: 
            coords = None

        atom_symbols = [atom.GetSymbol() for atom in rdkit_mol.GetAtoms()]
        obj = Compound(name=name, symbols=atom_symbols, coordinates=coords,
                       charge=charge, multiplicity=multiplicity,
                       charged_fragments=charged_fragments)

        # copy rdkit_mol to self.rdkit_mol 
        # TODO: is it possible to reomove old instance.
        obj.rdkit_mol = copy.deepcopy(rdkit_mol)
        
        return obj

        
    # @classmethod
    # def from_ase(cls, ase_atoms, name='mol', charged_fragments=True):
    #     """Initialize ´Compound´ from ASE Atoms 
    #     """
        
    #     charge = int(np.sum(ase_atoms.get_initial_charges()))
    #     multiplicity = int(np.sum(ase_atoms.get_initial_magnetic_moments()))

    #     obj = Compound(name=name, symbols=ase_atoms.symbols, 
    #                    coordinates=ase_atoms.positions,
    #                    charge=charge, multiplicity=multiplicity, 
    #                    charged_fragments = charged_fragments)

    #     # copy ase_atoms to self.ase_mol 
    #     # TODO: is it possible to reomove old instance.
    #     obj.ase_mol = copy.deepcopy(ase_atoms)
        
    #     return obj

    @classmethod
    def from_ac_matrix(cls, I, symbols, charge=0, multiplicity=1, name="mol",
                       charged_fragments=True, addHs=False, use_atom_maps=True):
        """Initialize ´Compound´ from AC matrix."""
        
        atomic_num_list = [int_atom(symbol) for symbol in symbols]
        proto_mol = get_proto_mol(atomic_num_list)

        rdkit_mol = AC2mol(proto_mol, I, atomic_num_list, charge, 
                           allow_charged_fragments=charged_fragments,
                           use_graph=True, use_atom_maps=use_atom_maps)

        Chem.Kekulize(rdkit_mol, clearAromaticFlags=True)
        if addHs:
            Chem.AddHs(rdkit_mol)
        
        obj = Compound(symbols=symbols, charge=charge, multiplicity=multiplicity,
                       name=name, charged_fragments=charged_fragments)
                       
        obj.rdkit_mol = rdkit_mol

        return obj

    def get_rdkit_mol(self, charged_fragments=True, huckel=False, use_atom_maps=False, embed_chiral=True):
        """Return rdkit mol object. 

        Problem: if charged fragment is set. It should use that value.
        """
        if getattr(self, 'rdkit_mol') is not None:
            return self.rdkit_mol
        else:
            atomic_numbers = [int_atom(atom) for atom in self.atomic_symbols]
            rdkit_mol = xyz2mol(atomic_numbers, self.coordinates, charge=self.charge,
                                allow_charged_fragments=charged_fragments, 
                                use_graph=True, use_huckel=huckel, embed_chiral=embed_chiral, 
                                use_atom_maps=use_atom_maps)
            self.rdkit_mol = rdkit_mol
            
            return rdkit_mol
    
    # def get_ase_atoms(self, multiplicity_change=0):
    #     """Return ASE Atoms object.

    #     multiplicity_change: int
    #         Allows one to change the multiplicity. Fx. change multiplicity from
    #         2S+1 to 2S by setting multiplicity_change=-1.
    #     """
    #     ase_atoms = Atoms(symbols=self.atomic_symbols, positions=self.coordinates, 
    #                       charges=[self.charge]+[0]*(self.natoms-1), 
    #                       magmoms=[self.multiplicity + multiplicity_change]+[0]*(self.natoms-1))
    #     self.ase_atoms = ase_atoms
        
    #     return ase_atoms

    def write_xyz(self, to_file=True):
        xyz_string = f"{self.natoms}\n{self.name}\n"
        for symbol, pos in zip(self.atomic_symbols, self.coordinates):
            xyz_string += symbol + "  " + " ".join([str(x) for x in pos]) + "\n"
        
        if to_file:
            with open(f'{self.name}.xyz', 'w') as f:
                f.write(xyz_string)
        else:
            return xyz_string

    @staticmethod
    def symbols2numbers(symbols):
        return [int_atom(symbol) for symbol in symbols]
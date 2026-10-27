import re
import struct
from pathlib import Path
import numpy as np
from ase import Atoms


BOHR_A=0.52917720859
MASS_AU=15.9994*1822.88848426455
TIME_AU_S=2.41888432650478e-17
C_CM_S=2.99792458e10


def restart_geometry(text):
    coord = re.search(r'&COORD\s*\n(.*?)&END COORD', text, re.S)[1]
    rows = [line.split() for line in coord.splitlines() if line.strip()]
    cell_text = re.search(r'&CELL\s*\n(.*?)&END CELL', text, re.S)[1]
    cell = [np.fromstring(re.search(r'(?m)^\s*' + label + r'\s+(.+)$', cell_text)[1], sep=' ')
            for label in 'ABC']
    assert all(len(v) == 3 for v in cell)
    assert all(len(row) == 4 for row in rows)
    return Atoms([r[0] for r in rows], positions=np.array([r[1:] for r in rows], float), cell=cell, pbc=True)

def replica_coordinates(text):
    blocks = re.findall(r'&REPLICA\s*\n(.*?)&END REPLICA', text, re.S)
    return np.array([np.array([line.split() for line in
        re.search(r'&COORD\s*\n(.*?)&END COORD', block, re.S)[1].strip().splitlines()], float)
        for block in blocks])

def binary_matrix(path):
    data = Path(path).read_bytes()
    assert len(data) == 96, 'Three 24-byte double-precision Fortran records expected'
    rows = []
    for offset in range(0, 96, 32):
        assert struct.unpack('<I', data[offset:offset+4])[0] == 24
        assert struct.unpack('<I', data[offset+28:offset+32])[0] == 24
        rows.append(np.frombuffer(data[offset+4:offset+28], '<f8').copy())
    return np.array(rows)

def signed_frequencies(matrix):
    values = np.linalg.eigvalsh(matrix) * 1e-6
    return np.sign(values) * np.sqrt(np.abs(values)) / (2*np.pi*TIME_AU_S*C_CM_S)

def displaced_hessian(text):
    pattern = (r'VIB\| REPLICA Nr\.\s+\d+- Energy and Forces for particle:\s+(\d+)\s+'
               r'coordinate: ([XYZ]) ([+-]) D[XYZ]\s*\nV?\s*VIB\|\s+Total energy:\s+([-\d.]+)'
               r'\s*\n\s*VIB\|\s+ATOM[^\n]*\n\s*VIB\|\s+O\s+([^\n]+)')
    entries = re.findall(pattern, text)
    assert len(entries) == 6
    forces, energies = {}, []
    for atom, axis, sign, energy, vector in entries:
        assert atom == '65'
        key = (axis, sign)
        assert key not in forces
        forces[key] = np.fromstring(vector, sep=' ')
        assert forces[key].shape == (3,)
        energies.append(float(energy))
    hessian = np.column_stack([-(forces[(axis,'+')]-forces[(axis,'-')]) / 0.02 for axis in 'XYZ'])
    return hessian, forces, energies

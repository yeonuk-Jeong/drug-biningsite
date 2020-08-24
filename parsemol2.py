import pandas as pd
import numpy as np
from biopandas.mol2.mol2_io import split_multimol2

COLUMN_NAMES = (
 'subst_id',
 'subst_name',
 'root_atom',
 'subst_type',
 'dic_type',
 'chain',
 'sub_type',
 'inter_bonds',
 'status',
 'sub_type2',
 'chain2',
 'residue_num'
)

COLUMN_TYPES = (int, str, int, str, int, str, str, int, str, str, str, int)

class parsemol2(object):
    
    def __init__(self):
        self._df = None
        self.mol2_text = ''
        self.header = ''
        self.code = ''        
        self.mol2_path = ''
    
    @property    
    def df(self):
        """Acccesses the pandas DataFrame"""
        return self._df
    
    def _load_mol2(self, mol2_lines, mol2_code, columns):
        """Load mol2 contents into assert_raise_message instance"""
        if columns is None:
            col_names = COLUMN_NAMES
            col_types = COLUMN_TYPES
        else:
            col_names, col_types = [], []
            for i in range(len(columns)):
                col_names.append(columns[i][0])
                col_types.append(columns[i][1])

        try:
            self.mol2_text = ''.join(mol2_lines)
            self.code = mol2_code
        except TypeError:
            mol2_lines = [m.decode() for m in mol2_lines]
            self.mol2_text = ''.join(mol2_lines)
            self.code = mol2_code.decode()
            
       
        self._df = self._construct_df(mol2_lines, col_names, col_types)
        
    def read_mol2(self, path, columns=None):
 
        mol2_code, mol2_lines = next(split_multimol2(path))
        self._load_mol2(mol2_lines, mol2_code, columns)
        self.mol2_path = path
        return self
    
    def _construct_df(self, mol2_lines, col_names, col_types):
        """Construct DataFrames from list of PDB lines."""
        return self._atomsection_to_pandas(self._get_atomsection(mol2_lines),
                                           col_names=col_names,
                                           col_types=col_types)
    @staticmethod
    def _get_atomsection(mol2_lst):
        """Returns atom section from mol2 provided as list of strings"""
        started = False
        last_idx_plus1 = -1
        for idx, s in enumerate(mol2_lst):
            if s.startswith('@<TRIPOS>SUBSTRUCTURE'):
                first_idx = idx + 1
                started = True
            elif started and s.startswith('@<TRIPOS>'):
                last_idx_plus1 = idx
                break
                
        if last_idx_plus1 == -1 : return mol2_lst[first_idx:]
        return mol2_lst[first_idx:last_idx_plus1]

    @staticmethod
    def _atomsection_to_pandas(mol2_seq_lst, col_names, col_types):
        
        bining_list = [lst.split() for lst in mol2_seq_lst]
        for i in range(len(bining_list)):
            if len(bining_list[i]) == 11:
                str1 = bining_list[i][-1]
                bining_list[i][-1] = filter(str.isalpha, str1)
                bining_list[i].append(filter(str.isdigit, str1))

        df = pd.DataFrame(bining_list,columns=col_names)

    
        for i in range(df.shape[1]):
            df[col_names[i]] = df[col_names[i]].astype(col_types[i])

        return df
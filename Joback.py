"""
                    Joback Group Contribution method

author: Edgar Ivan Sanchez Medina
email: sanchez@mpi-magdeburg.mpg.de
Date: October 2020

SMARTS visualization: https://smartsview.zbh.uni-hamburg.de/
SMARTS info: https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
SMARTS from: https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg07446.html

Updates:
    23.10.2020 - Tm and Tb functions
"""
from rdkit import Chem
import pandas as pd
import numpy as np

class Joback:
    '''
    Joback group contribution method.
    
    Multiple linear regression techniques were employed to determine the group
    contributions for each structurally-dependent parameter.
    
    "High accuracy is not claimed, but the proposed methods are often as 
    accurate as or more accurate than techniques in common use today." [2]
    
    "Note we have assumed no
    interaction between groups, and structurally-dependent parameters 0:) are
    thereby determined simply by summing the number frequency of each group
    times its group contribution. Clearly this is only a first-order 
    approximation, but it is our experience that there rarely exist sufficient 
    data to obtain reliable contributions between so-called "next-nearest 
    neighbors." Also, introducing group interactions often leads to a 
    significantly more complex estimation method."[2]
    
    "Estimations of the normal boiling point and, in particular, the normal freezing
    point are not accurate and should be considered as only very approximate. There
    is a strong dependence of the actual conformation of the molecule on various 
    properties" (T_m and H_f).[2]
    
    References
    ----------
    .. [1] Joback, K.G., 1984. A unified approach to physical property 
            estimation using multivariate statistical techniques (Doctoral 
            dissertation, Massachusetts Institute of Technology).
    .. [2] Joback, K.G. and Reid, R.C., 1987. Estimation of pure-component 
            properties from group-contributions. Chemical Engineering 
            Communications, 57(1-6), pp.233-243.
    '''
    ##########################
    # --- Initialization --- #
    ##########################    
    def __init__(self, SMILES):
        mol = Chem.MolFromSmiles(SMILES); self.flag=''
        # Check smiles
        if type(mol) != Chem.rdchem.Mol:
            raise Exception('SMILES: '+SMILES+' was NOT recognized by RDKit')
        num_atoms = len(mol.GetAtoms())
        
        # Count groups in molecule
        self.data = pd.read_csv('Joback_data.csv'); self.count=np.zeros(41) 
        atoms_matched=0; self.groups_selected={}
        for i, smarts in enumerate(self.data['SMARTS']):
            group         = Chem.MolFromSmarts(smarts)
            matches       = mol.GetSubstructMatches(group)
            self.count[i] = len(matches)
            for match in matches:
                atoms_matched += len(match)
            if self.count[i] > 0:
                self.groups_selected[self.data['Group'][i]] = self.count[i]
        
        # Count atoms for checking
        if num_atoms != atoms_matched:
            self.flag = f'{atoms_matched} atoms matched out of {num_atoms} atoms'
        
    ##################################
    # --- Flag of atoms matching --- #
    ##################################
    def flag(self):
        return self.flag
    
    ###########################
    # --- Count of groups --- #
    ###########################
    def groups(self):
        return self.groups_selected
    
    ###############################
    # --- Temperature boiling --- #
    ###############################
    def Tb(self):
        '''
        Estimates normal boiling temperature.
        
         .. math::
            T_b = 198.2 + \sum_i {T_{b,i}}
            
        For obtaining the parameters:
            * 438 compounds were used.
            * AAE: 12.9 K
            * std: 17.9 K
            * APE: 3.6 %

        Returns
        -------
        Tb : float
            Estimated normal boiling point [K]
        '''
        count = self.count; data= self.data
        
        Tb_values = data['Tb'].to_numpy()   
        where_are_NaNs = np.isnan(Tb_values)
        Tb_values[where_are_NaNs] = 0       # Coeff values
        count = self.count                  # Count of groups
        sum_Tb = np.dot(Tb_values, count)   # Dot product
         
        Tb = 198.2 + sum_Tb
        return Tb
    
    ###############################
    # --- Temperature melting --- #
    ###############################
    def Tm(self):
        '''
        Estimates normal melting temperature.
        
         .. math::
            T_m = 122.5 + \sum_i {T_{m,i}}
            
        For obtaining the parameters:
            * 388 compounds were used.
            * AAE: 22.6 K
            * std: 24.7 K
            * APE: 11.2 %

        Returns
        -------
        Tb : float
            Estimated normal boiling point [K]
        '''
        count = self.count; data= self.data
        
        Tm_values = data['Tm'].to_numpy()   
        where_are_NaNs = np.isnan(Tm_values)
        Tm_values[where_are_NaNs] = 0       # Coeff values
        count = self.count                  # Count of groups
        sum_Tm = np.dot(Tm_values, count)   # Dot product
         
        Tm = 122.5 + sum_Tm
        return Tm
        





















from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

proton = 1.00784

class ChemReader:
    
    """
    >>> reader = ChemReader('pentamer.esi.json', maxlength=5)
    
    >>> reader.sequence('*AB', '[NH]C([C+]=O)CCSCCC([O])=O')
    684.226076164
    >>> reader.sequence('*ABC', '[NH]C([C+]=O)CCSCCC([O])=O')
    1044.361583652
    >>> reader.sequence('*ABCD', '[NH]C(CCSCCC([O])=O)C(NCCO[C+]=O)=O')
    1463.53420492
    
    >>> reader.match_weight(394.13)
    [(394.132433488, '*A')]
    >>> reader.match_weight(437.17)
    [(437.174632648, '*A')]
    >>> reader.match_weight(727.26773)
    [(727.268275324, '*AB')]
    >>> reader.match_weight(1087.40323)
    [(1087.4037828120001, '*ABC')]
    >>> reader.match_weight(1376.50163)
    [(1376.50217652, '*ABCD')]
    >>> reader.match_weight(1419.54)
    [(1419.54437568, '*ABCD')]
    >>> reader.match_weight(530.24)
    [(530.235853, 'E*')]
    >>> reader.match_weight(617.26733)
    [(617.2678814, 'E*')]
    >>> reader.match_weight(862.3759)
    [(862.3764458680001, 'DE*')]
    >>> reader.match_weight(1222.5114)
    [(1222.511953356, 'CDE*')]
    >>> reader.match_weight(1463.53366)
    [(1463.53420492, '*ABCD')]
    >>> reader.match_weight(1512.61)
    [(1512.605596032, 'BCDE*'), (1512.629405476, 'DDE*')]
    >>> reader.match_weight(1599.64)
    [(1599.637624432, 'BCDE*')]
    >>> reader.match_weight(1992.76168)
    [(1992.7622478560004, '*ABCDE*')]
    
    >>> reader.match_weight(996.88448)
    [(1991.7544078560004, '*ABCDE* (m/2 + 0 protons)')]
    
    >>> reader.match_weight(1856.64)  # '*ABCD?*'
    [(1856.624456312, 'ABCC*')]

    >>> reader = ChemReader('octamer.esi.json', maxlength=8)
    >>> reader.sequence('*ABCDEACB*')
    3012.085657204
    
    >>> reader.match_weight(441.15124)
    [(441.151789028, '*A')]
    >>> reader.match_weight(398.10904)
    [(398.109589868, '*A')]
    >>> reader.match_weight(731.24488)
    [(731.245431704, '*AB')]
    >>> reader.match_weight(1091.38039) # *ABC
    [(1091.3809391920001, '*ABC'), (1091.368363128, 'AD*'), (1091.368363128, 'AD*')]
    >>> reader.match_weight(1423.52098) # *ABCD
    [(1423.52153206, '*ABCD'), (1423.5089559960002, 'ADD*'), (1423.5089559960002, 'ADD*')]
    >>> reader.match_weight(1835.72417) # *ABCDE
    [(1835.7094690440003, '*CCDDD'), (1835.7247251840001, '*ABCDE'), (1835.71214912, 'ADDE*'), (1835.71214912, 'ADDE*')]
    >>> reader.match_weight(2201.84912) # *ABCDEA -> *AABCDE
    [(2201.8344118480004, '*ACCDDD'), (2201.849667988, '*AABCDE'), (2201.837091924, 'AADDE*'), (2201.8370919240006, 'AADDE*'), (2201.858221292, 'CCCCE*')]
    >>> reader.match_weight(2561.98462) # *ABCDEAC -> *AABCCDE
    [(2561.969919336, '*ACCCDDD'), (2561.9851754759998, '*AABCCDE'), (2561.9725994120004, 'AACDDE*'), (2561.9725994120004, 'AACDDE*'), (2561.9878555520004, 'AAABEE*'), (2561.99372878, 'CCCCCE*')]
    >>> reader.match_weight(1507.55178) # *ABCDEACB* (m/2)
    [(3012.085657204, '*AABBCCDE* (m/2 + 1 protons)')]

    >>> reader.match_weight(408.12575)
    [(408.126302552, 'B*')]
    >>> reader.match_weight(768.26126)
    [(768.26181004, 'BC*')]
    >>> reader.match_weight(1134.3862)
    [(1134.386752844, 'ABC*')]
    >>> reader.match_weight(1546.5894)
    [(1546.589945968, 'ABCE*')]
    >>> reader.match_weight(1878.72999) # DEACB* -> ABCDE*
    [(1878.734561596, '*AABEE'), (1878.730538836, 'ABCDE*')]
    >>> reader.match_weight(2238.8655) # CDEACB* -> ABCCDE*
    [(2238.8700690840005, '*AABCEE'), (2238.8660463240003, 'ABCCDE*')]
    >>> reader.match_weight(2528.95914) # BCDEACB* -> ABBCCDE*
    [(2528.9637117600005, '*AABBCEE'), (2528.9596890000003, 'ABBCCDE*'), (2528.9682423040003, 'CCDDDD*')]
    >>> reader.match_weight(2615.99117) # BCDEACB* -> ABBCCDE*
    [(2615.9957401600004, '*AABBCEE'), (2616.0109963, '*AAAADDE'), (2615.9917174, 'ABBCCDE*')]

    >>> reader = ChemReader('pentamer.maldi.json', maxlength=8)
    >>> reader.match_weight(459.214057, eps=0.06)
    [(459.15656192800003, '*A')]
    >>> reader.match_weight(749.252737, eps=0.06)
    [(749.250204604, '*AB')]
    >>> reader.match_weight(1109.415772, eps=0.06)
    [(1109.3857120920002, '*ABC')]
    >>> reader.match_weight(1534.628544, eps=0.06)
    [(1534.587525312, 'BCDE*'), (1534.6113347560001, 'DDE*')]
    >>> reader.match_weight(1398.466924, eps=0.06)
    [(1398.4841058, '*ABCD'), (1398.4497487360002, 'BBC*')]
    >>> reader.match_weight(1331.597447, eps=0.08)
    [(1331.525911036, 'CDE*')]
    """
    
    def __init__(self, components, maxlength):
        
        components = eval(open(components, 'r').read())

        self.start = Chem.MolFromSmiles(components['start'])
        self.stop = Chem.MolFromSmiles(components['stop'])
        self.backbone = Chem.MolFromSmiles(components['backbone'])
        self.multiplyCharged = components['multiplyCharged']
        self.charge = components['charge']

        self.start_fragments = [
            Chem.MolFromSmiles(m) 
            for m in components['start_fragments']
        ]

        self.stop_fragments = [
            Chem.MolFromSmiles(m) 
            for m in components['stop_fragments']
        ]

        self.flags = {
            letter:Chem.MolFromSmiles(m) 
            for letter, m in components['flags'].items()
        }
        
        self.maxlength = maxlength
        self.generate_fragments()
        
    def __str__(self):
        
        return 'start: {} ({})\nbackbone: {} ({})'.format(
            Chem.MolToSmiles(self.start),
            Chem.rdMolDescriptors.CalcExactMolWt(self.start),
            Chem.MolToSmiles(self.backbone),
            Chem.rdMolDescriptors.CalcExactMolWt(self.backbone),
        )
        
    def sequence(self, flags, extra=None):
        
        start, stop = False, False
        if flags.startswith('*'):
            start = True
            flags = flags[1:]
        if flags.endswith('*'):
            stop = True
            flags = flags[:-1]
        
        weight = 0
        
        if start:
            weight += Chem.rdMolDescriptors.CalcExactMolWt(self.start)
            
        if stop:
            weight += Chem.rdMolDescriptors.CalcExactMolWt(self.stop)

        for flag in flags:
            weight += Chem.rdMolDescriptors.CalcExactMolWt(self.flags[flag])
            
        n = len(flags)
        if not start:
            n -= 1
        if not stop:
            n -= 1
            
        weight += n * Chem.rdMolDescriptors.CalcExactMolWt(self.backbone)
        
        if extra is not None:
            weight += Chem.rdMolDescriptors.CalcExactMolWt(Chem.MolFromSmiles(extra))
            
        return weight
    
    def generate_fragments(self):
        
        self.missing_start = []
        self.missing_stop = []
        self.missing_none = []
        self.missing_both = []
        
        flags = sorted(self.flags)
        
        fragments = flags[:]
        
        for _ in range(self.maxlength - 1):
            
            self.missing_start.extend(
                (self.sequence(fragment + '*'), fragment + '*') 
                for fragment in fragments
            )
            
            self.missing_stop.extend(
                (self.sequence('*' + fragment), '*' + fragment) 
                for fragment in fragments
            )
            
            self.missing_both.extend(
                (self.sequence(fragment), fragment) 
                for fragment in fragments
            )
            
            fragments = [
                fragment + flag
                for fragment in fragments
                for flag in flags[flags.index(fragment[-1]):]
            ]

        self.missing_start.extend(
            (self.sequence(fragment + '*'), fragment + '*') 
            for fragment in fragments
        )
        
        self.missing_stop.extend(
            (self.sequence('*' + fragment), '*' + fragment) 
            for fragment in fragments
        )
            
        self.missing_both.extend(
            (self.sequence(fragment), fragment) 
            for fragment in fragments
        )
            
        self.missing_none.extend(
            (self.sequence('*' + fragment + '*'), '*' + fragment + '*') 
            for fragment in fragments
        )
        
        self.missing_start = sorted(self.missing_start)
        self.missing_stop = sorted(self.missing_stop)
        self.missing_none = sorted(self.missing_none)
        
    def match_weight(self, weight, eps=0.02):
        
        matches = []

        # find matches with start sequences
        stop_fragments = [
            (f, Chem.rdMolDescriptors.CalcExactMolWt(f)) 
            for f in self.stop_fragments
        ]
        for w, s in self.missing_stop:
            for f, fw in stop_fragments:
                if w + fw - eps <= weight + proton - self.charge <= w + fw + eps:
                    matches.append((w + fw + self.charge - proton, s))
                    
        # find matches with stop sequences
        start_fragments = [
            (f, Chem.rdMolDescriptors.CalcExactMolWt(f)) 
            for f in self.start_fragments
        ]
        for w, s in self.missing_start:
            for f, fw in start_fragments:
                if w + fw - eps <= weight + proton - self.charge <= w + fw + eps:
                    matches.append((w + fw + self.charge - proton, s))
                    
        # find matches with internal sequences
        for w, s in self.missing_start:
            for f1, f1w in stop_fragments:
                for f2, f2w in start_fragments:
                    if w + f1w + f2w - eps <= weight + proton - self.charge <= w + f1w + f2w + eps:
                        matches.append((w + f1w + f2w + self.charge - proton, s))
                    
        # find matches with complete sequences
        for w, s in self.missing_none:
            if w - eps <= weight - self.charge <= w + eps:
                matches.append((w + self.charge, s))
                    
        # find matches with complete sequences
        if self.multiplyCharged:
            for w, s in self.missing_none:
                for protons in range(2):
                    if w - eps <= 2 * (weight - proton) - protons * proton <= w + eps:
                        matches.append((w, s + ' (m/2 + {} protons)'.format(protons)))
                    
        return matches
                            
if __name__ == '__main__':
    import doctest
    doctest.testmod()
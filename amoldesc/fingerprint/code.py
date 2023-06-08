from rdkit import Chem,DataStructs
from rdkit.Chem import MACCSkeys, AllChem, rdMolDescriptors
from rdkit.Chem.AtomPairs import Pairs,Torsions

def fingerprint_generation(smile):
    value = {}
    
    molval=Chem.MolFromSmiles(smile)

    maccs = MACCSkeys.GenMACCSKeys(molval)
    pval_maccs=maccs.ToBitString()

    rdk = Chem.RDKFingerprint(molval)
    pval_rdk=rdk.ToBitString()

    pairs = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(molval)
    pval_pairs=pairs.ToBitString()

    tts = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(molval)
    pval_tts=tts.ToBitString()

    morgan = AllChem.GetMorganFingerprintAsBitVect(molval,2)
    pval_morgan=morgan.ToBitString()

    value['smile']=smile
    value['maccs']=pval_maccs
    value['rdk']=pval_rdk
    value['pairs']=pval_pairs
    value['tts']=pval_tts
    value['morgan']=pval_morgan

    return value

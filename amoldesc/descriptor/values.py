from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors,rdMolDescriptors
from rdkit.Chem import Crippen
from rdkit.Chem import AllChem
from rdkit.Chem.EState import EState

def get_properties_and_analysis(smile):

    values={}
    mol = Chem.MolFromSmiles(smile)

    Draw.MolToFile(mol, 'static/images/mol.png')

    k=Descriptors.MolWt(mol)
    values['Molecular Weight']=k
    g=Descriptors.HeavyAtomMolWt(mol)
    values['Number of Heavy Atoms']=g
    b=Descriptors.FractionCSP3(mol)
    values['FractionCSP3']=b
    i=Descriptors.NumRotatableBonds(mol)
    values['Number of Rotatable Bonds']=i
    a=Descriptors.NumHAcceptors(mol)
    values['NumHAcceptors']=a
    l=Descriptors.NumHDonors(mol)
    values['NumHDonors']=l
    o=Descriptors.MolMR(mol)
    values['Molar Refractivity']=o
    C=Descriptors.TPSA(mol)
    values['TPSA']=C

    # Liphocity
    lipho={}
    x=Descriptors.MolLogP(mol)
    lipho['Log P']=x
    MLogP = Crippen.MolLogP(mol)/Chem.Descriptors.MolWt(mol)
    lipho['Mog P']=MLogP
    aLogP = Descriptors.MolLogP(mol)
    lipho['aLog P']=aLogP
    XLogp3 = Crippen.MolLogP(mol, True)
    lipho['XLog P3']=XLogp3
    iLogP = (Descriptors.MolLogP(mol) - 0.74) * (Descriptors.NumRotatableBonds(mol) - 0.007) * (
                                Descriptors.MolWt(mol) < 5000) + 0.22
    lipho['iLog P']=iLogP
    WLOGP = sum(EState.EStateIndices(mol))
    lipho['WLog P']=WLOGP
    LOGD = Descriptors.MolLogP(mol) - 0.74 * (mol.GetNumHeavyAtoms() ** 0.5) - 0.47
    lipho['Log D']=LOGD
    LOGS = (- 0.048 * (rdMolDescriptors.CalcTPSA(mol)) - 0.104 * (Descriptors.MolLogP(mol))-0.295)
    lipho['Log S']=LOGS

 # Drug Likeness
    drug={}
    if(l<=5 and a<=10)and(k<=500 and x<=5):
        drug['Lipinski']='Yes; 0 violation'
    else:
        drug['Lipinski']='No;  violation'
    
    if(-0.4>=x<=5.6 and 40>=o<=130)and(160>=k<=480 and 70>=x<=20):
        drug['Ghose']='No'
    else:
        drug['Ghose']='Yes'
    
    if(i<=10)and(C<=140):
        drug['Veber']='Yes'
    else:
        drug['Veber']='No'

    if(x<=5.88)and(C<=131.6):
        drug['Egan']='Yes'
    else:
        drug['Egan']='No'
    
    m=Descriptors.NumAliphaticRings(mol)
    values['NumAliphaticRings']=m
    t=Descriptors.NumHeteroatoms(mol)
    values['NumHeteroatoms']=t

    if(200>=k<=600 and -2>=x<=5)and(C<=150 and m<=7) and(t>1 and i<=15) and(a<=10 and l<=10):
        drug['Muegge']='Yes'
    else:
        drug['Muegge']='No , Violation Found'

    if(1.3>=aLogP<=4.1 and 70>=o<=110)and (230>=k<=390 and 30>=g<=55):
        drug['CMC 50 Like rule']='Yes'
    else:
        drug['CMC 50 Like rule']='No , Violation Found'



    
    # Medicinal Chemistry
    medicinal_chem={}
    if(250>=k<=350 and x<=3.5 and i<=7):
        medicinal_chem['Leadlikeness'] = 'Yes'
    else:
        medicinal_chem['Leadlikeness'] = 'No,Violation Found'

    # Excretion
    ex={}
    Cl =  0.025 * (Descriptors.MolWt(mol)) ** 0.75 * 10 ** (0.107 * (Descriptors.MolLogP(mol)))
    ex['Cl']=Cl
    CLint =  0.025 * Descriptors.MolWt(mol) ** 0.75
    ex['CLint']=CLint
    THalf = 0.693 * (Descriptors.MolWt(mol)) ** 0.5 / (10 ** (0.006*(Descriptors.MolLogP(mol)))+1)
    ex['THalf']=THalf
    Tox =Descriptors.MolLogP(mol) > 5.0
    ex['Tox']=Tox
    BCF = 0.176 * (Descriptors.MolLogP(mol)) - 0.00358 * (Descriptors.MolWt(mol))+1.351
    ex['BCF']=BCF


    # Absorption
    ab={}
    p_gp = Descriptors.MolLogP(mol) > 0.0
    HIA = 10 ** (0.022 * Descriptors.TPSA(mol) - 0.675 * Descriptors.MolLogP(mol) - 0.005 * Descriptors.MolWt(mol) + 0.861)
    PPB = Descriptors.MolWt(mol)<500
    ab['pgp']=p_gp
    ab['HIA']=HIA
    ab['PPB']=PPB

    # Distribution
    dis={}
    BBB = Descriptors.MolLogP(mol) < -0.3
    Fraction_Unbound = 0.74 * Descriptors.MolLogP(mol) - 0.007 * Descriptors.MolWt(mol) - 0.27 * Descriptors.NumRotatableBonds(mol) - 0.42 * Descriptors.NumHAcceptors(mol)-1.12
    dis['BBB']=BBB
    dis['Fraction_Unbound']=Fraction_Unbound
    result={'search_term':smile,'prop':values, 'lipho':lipho, 'drug':drug, 'med':medicinal_chem,
            'ex':ex,'ab':ab,'dis':dis}
    
    return result
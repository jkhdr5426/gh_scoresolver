import xml.etree.ElementTree as ET

def parse_pubchem_xml(pubchem_xml_path):
    tree = ET.parse(pubchem_xml_path)
    root = tree.getroot()
    
    # Extract atom information
    atoms = []
    atom_elements = root.findall('.//PC-Atoms_element/PC-Element')
    atom_ids = root.findall('.//PC-Atoms_aid/PC-Atoms_aid_E')
    
    for atom_id, element in zip(atom_ids, atom_elements):
        atoms.append({'id': atom_id.text, 'element': element.attrib['value']})

    # Extract bond information
    bonds = []
    bond_from = root.findall('.//PC-Bonds_aid1/PC-Bonds_aid1_E')
    bond_to = root.findall('.//PC-Bonds_aid2/PC-Bonds_aid2_E')
    bond_order = root.findall('.//PC-Bonds_order/PC-BondType')
    
    for b_from, b_to, b_order in zip(bond_from, bond_to, bond_order):
        bonds.append({'from': b_from.text, 'to': b_to.text, 'order': b_order.attrib['value']})

    return atoms, bonds

def create_openmm_xml(atoms, bonds):
    root = ET.Element('ForceField')

    atom_types = ET.SubElement(root, 'AtomTypes')
    for i, atom in enumerate(atoms):
        element = atom['element']
        mass = '12.01' if element == 'C' else ('16.00' if element == 'O' else '18.998')
        ET.SubElement(atom_types, 'Type', attrib={
            'name': f'{element}{i+1}', 'class': element, 'element': element, 'mass': mass
        })

    residues = ET.SubElement(root, 'Residues')
    residue = ET.SubElement(residues, 'Residue', attrib={'name': 'PFOA'})

    for i, atom in enumerate(atoms):
        ET.SubElement(residue, 'Atom', attrib={
            'name': f'{atom["element"]}{atom["id"]}', 'type': f'{atom["element"]}{i+1}', 'charge': '0.0'
        })

    for bond in bonds:
        ET.SubElement(residue, 'Bond', attrib={'from': f'{bond["from"]}', 'to': f'{bond["to"]}'})

    xml_str = ET.tostring(root, encoding='unicode')
    with open('deprotonated_pfoa.xml', 'w') as f:
        f.write(xml_str)

    print(xml_str)

# Path to your PubChem XML file
pubchem_xml_path = '/mnt/data/Conformer3D_COMPOUND_CID_9554.xml'

# Parse PubChem XML and create OpenMM XML
atoms, bonds = parse_pubchem_xml(pubchem_xml_path)
create_openmm_xml(atoms, bonds)

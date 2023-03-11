'''
xml_reader.py reads a metabolite xml file from the website http://www.hmdb.ca and turns it into a 
tabular dataset
'''

import xml.etree.ElementTree as ET 

'''
params: 
    tree [ElementTree]: the ElementTree for the basic subtree that will be parsed into a dict
    (basic tree means that the tree has several repeating tags with identical subtags, including
    a subtag 'name')

returns: the subdictionary for the basic tree
'''
def basic_tree_to_dict(tree, tags_to_keys):
    dict = {}

    for child in tree:
        child_dict = {}
        name = ''

        for child_tag in child:
            if child_tag.tag == '{http://www.hmdb.ca}name':
                name = child_tag.text
            elif child_tag.tag in tags_to_keys:
                child_dict[tags_to_keys[child_tag.tag]] = child_tag.text
        
        dict[name] = child_dict
        
    return dict


'''
params: 
    biological_properties [ElementTree]: the ElementTree for the biological_properties section of
    the metabolite XML file being parsed

returns: the biological properties subdictionary based on the biological_properties section of the 
xml file
'''
def bio_properties_dict(biological_properties):
    dict = {}

    for bio_tag in biological_properties:
        if bio_tag.tag == '{http://www.hmdb.ca}cellular_locations':
            locations = []
            for location in bio_tag:
                locations.append(location.text)
        
            dict['Cellular_Locations'] = locations
        elif bio_tag.tag == '{http://www.hmdb.ca}pathways':
            tags_to_keys = {
                '{http://www.hmdb.ca}smpdb_id': 'SMPDB_ID',
                '{http://www.hmdb.ca}kegg_map_id': 'Kegg_Map_ID',
            }

            dict['Pathways'] = basic_tree_to_dict(bio_tag, tags_to_keys)
        
    return dict

'''
params: 
    diseases [ElementTree]: the ElementTree for the diseases section of the metabolite XML file
    being parsed

returns: the associated diseases and disorders subdictionary based on the biological_properties
section of the xml file
'''
def diseases_dict(diseases):
    dict = {}

    for disease in diseases:
        disease_dict = {}
        name = ''

        for disease_tag in disease:
            if disease_tag.tag == '{http://www.hmdb.ca}name':
                name = disease_tag.text
            elif disease_tag.tag == '{http://www.hmdb.ca}omim_id':
                disease_dict['OMIM_ID'] = disease_tag.text
            elif disease_tag.tag == '{http://www.hmdb.ca}references':
                references = {}
                i = 1
                for reference in disease_tag: 
                    references[i] = {}
                    for reference_tag in reference: 
                        if reference_tag.tag == '{http://www.hmdb.ca}pubmed_id':
                            references[i]['PubMed_ID'] = reference_tag.text
                        elif reference_tag.tag == '{http://www.hmdb.ca}reference_text':
                            references[i]['Reference_Text'] = reference_tag.text
                    i += 1
                
                disease_dict['References'] = references

        dict[name] = disease_dict
    return dict


'''
params: 
    ontology [ElementTree]: the ElementTree for the ontology section of
    the metabolite XML file being parsed

returns: the ontology subdictionary based on the ontology section of the 
xml file
'''
def ontology_dict(ontology):
    dict = {}
    desired_roots = ['Physiological effect', 'Process', 'Role']

    for root in ontology:
       root_name = ''
       for root_tag in root:
           if root_tag.tag == '{http://www.hmdb.ca}term':
                root_name = root_tag.text
       
       if root_name not in desired_roots:
           continue
       
       for root_tag in root:
           if root_tag.tag == '{http://www.hmdb.ca}descendants':
                descendants_dict = {}
                for descendant in root_tag:
                    descendants_dict = parse_root_descendants(descendants_dict, descendant, root_name)

                dict[root_name] = descendants_dict
        
    return dict


'''
params: 
    descendants_dict [dict]: the current dictionary of descendants
    descedant [ElementTree]: the descedant node of the root
    parent_node [str]: the string name of the parent node (found in the term tag)

returns: the completed dictionary of the ontology root descendants
'''
def parse_root_descendants(descendants_dict, descendant, parent_node):
    if descendant == None:
        return descendants_dict
    
    descendant_name = ''
    descendant_dict = {'parent_node': parent_node}
    for tag in descendant: 
        if tag.tag == '{http://www.hmdb.ca}term':
            descendant_name = tag.text
        elif tag.tag == '{http://www.hmdb.ca}definition':
            descendant_dict['annotation'] = tag.text
        elif tag.tag == '{http://www.hmdb.ca}descendants':
            for child in tag: 
                descendants_dict = parse_root_descendants(descendants_dict, child, descendant_name)
        
        descendants_dict[descendant_name] = descendant_dict
    
    return descendants_dict


'''
params: 
    xml_file [str]: file name/pathway of the xml file being parsed

returns: the dictionary of the csv file 
'''
def xml_to_dict(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    metabolites = {}

    for metabolite in root:      # list of metabolites
        name = ''
        metabolite_dict = {}

        for meta_tag in metabolite:
            # print(metaTag.tag)
            if meta_tag.tag == '{http://www.hmdb.ca}name':
                name = meta_tag.text
            elif meta_tag.tag == '{http://www.hmdb.ca}accession':
                metabolite_dict['HMDBP_ID'] = meta_tag.text
            elif meta_tag.tag == '{http://www.hmdb.ca}drugbank_id':
                metabolite_dict['DrugBank_ID'] = meta_tag.text
            elif meta_tag.tag == '{http://www.hmdb.ca}ontology':
                metabolite_dict['Ontology'] = ontology_dict(meta_tag)
            elif meta_tag.tag == '{http://www.hmdb.ca}biological_properties':
                metabolite_dict['Biological_Properties'] = bio_properties_dict(meta_tag)
            elif meta_tag.tag == '{http://www.hmdb.ca}diseases':
                metabolite_dict['Diseases'] = diseases_dict(meta_tag)
                pass
            elif meta_tag.tag == '{http://www.hmdb.ca}protein_associations':
                tags_to_keys = {
                    '{http://www.hmdb.ca}protein_accession': 'HMDBP_ID',
                    '{http://www.hmdb.ca}uniprot_id': 'Uniprot_ID',
                    '{http://www.hmdb.ca}gene_name': 'Gene'
                }

                metabolite_dict['Enzymes'] = basic_tree_to_dict(meta_tag, tags_to_keys)

        metabolites[name] = metabolite_dict
    
    return metabolites
        

if __name__ == '__main__':
    urine_metabolites = "urine_metabolites.xml"
    serum_metabolites = "serum_metabolites.xml"
    all_metabolites = "hmdb_metabolites.xml"

    xml_file = urine_metabolites
    
    dict = xml_to_dict(xml_file)
    #print(dict['1-Methylhistidine'])       # test for urine metabolites
    print(dict[list(dict.keys())[1]])

    
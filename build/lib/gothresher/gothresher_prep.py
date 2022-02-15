import sys
import argparse
import configparser
import networkx as nx
import pickle as cp
import os


# Some filenames
# FILE_ALTERNATE_ID_TO_ID_MAPPING = "data/alt_to_id.graph"
# FILE_CAFA_ID_TO_UNIPROT_ID_MAP = "data/CAFAIDTOUniprotIDMap.txt"
# FILE_MFO_ONTOLOGY_GRAPH = "data/mf.graph"
# FILE_BPO_ONTOLOGY_GRAPH = "data/bp.graph"
# FILE_CCO_ONTOLOGY_GRAPH = "data/cc.graph"
# FILE_MFO_ONTOLOGY_ANCESTORS_GRAPH = "data/mf_ancestors.map"
# FILE_BPO_ONTOLOGY_ANCESTORS_GRAPH = "data/bp_ancestors.map"
# FILE_CCO_ONTOLOGY_ANCESTORS_GRAPH = "data/cc_ancestors.map"

# Iddo: fixed paths to incorporate ONTO_DIR (passed as an argument from argparse)
# FILE_ALTERNATE_ID_TO_ID_MAPPING = "alt_to_id.graph"
# FILE_CAFA_ID_TO_UNIPROT_ID_MAP = "CAFAIDTOUniprotIDMap.txt"
# FILE_MFO_ONTOLOGY_GRAPH = "mf.graph"
# FILE_BPO_ONTOLOGY_GRAPH = "bp.graph"
# FILE_CCO_ONTOLOGY_GRAPH = "cc.graph"
# FILE_MFO_ONTOLOGY_ANCESTORS_GRAPH = "mf_ancestors.map"
# FILE_BPO_ONTOLOGY_ANCESTORS_GRAPH = "bp_ancestors.map"
# FILE_CCO_ONTOLOGY_ANCESTORS_GRAPH = "cc_ancestors.map"
# ROOT_BPO = 'GO:0008150'
# ROOT_CCO = 'GO:0005575'
# ROOT_MFO = 'GO:0003674'


def findAllAncestors(g, node):
    ancestors = [node]
    for immediate_ancestor in g.successors(node):
        # print(immediate_ancestor)
        ancestors.append(immediate_ancestor)
        ancestor_temp = findAllAncestors(g, immediate_ancestor)
        ancestors.extend(ancestor_temp)
    return list(set(ancestors))


def parseGOTerms(inputfile):
    fhr = open(inputfile, "r")
    mf_g = nx.DiGraph()
    cc_g = nx.DiGraph()
    bp_g = nx.DiGraph()
    allGOterms = fhr.read().split("[Term]")
    GO_term = ""
    namespace = ""
    relation = ""
    alt_id_to_id_mapping = dict()
    mf = []
    bp = []
    cc = []
    for term in allGOterms[1:]:
        split_term = term.split("\n")
        # alt_ids=[]
        for line in split_term:
            if "id:" in line and "GO:" in line and "alt_id" not in line:
                GO_term = "GO:" + line.split("GO:")[-1].strip()
            if "namespace: biological_process" in line:
                namespace = "bp"
                bp.append(GO_term)
            elif "namespace: cellular_component" in line:
                namespace = "cc"
                cc.append(GO_term)
            elif "namespace: molecular_function" in line:
                namespace = "mf"
                mf.append(GO_term)

    mf_g.add_nodes_from(mf)
    bp_g.add_nodes_from(bp)
    cc_g.add_nodes_from(cc)
    for term in allGOterms[1:]:
        split_term = term.split("\n")
        # alt_ids=[]
        for line in split_term:
            if "id:" in line and "GO:" in line and "alt_id" not in line:
                GO_term = "GO:" + line.split("GO:")[-1].strip()
            if "GO:" in line and "alt_id" in line:
                alt_id = "GO:" + line.split("GO:")[-1].strip()
                alt_id_to_id_mapping[alt_id] = GO_term
            if "namespace: biological_process" in line:
                namespace = "bp"
            elif "namespace: cellular_component" in line:
                namespace = "cc"
            elif "namespace: molecular_function" in line:
                namespace = "mf"

            # Iddo: this looks dangerous!
            if ("is_a:" in line or "relationship: part_of" in line) and "GO" in line:
                if "is_a" in line:
                    relation = "is_a"
                else:
                    relation = "part_of"
                try:
                    parent_GO_term = "GO:" + line.split("GO:")[1][:7]
                except IndexError:
                    print(f"Error in obo file {line}")
                    exit()

                # print(GO_term,line)
                if namespace == "bp" and parent_GO_term in bp:
                    bp_g.add_edge(GO_term, parent_GO_term, weight=relation)
                elif namespace == "cc" and parent_GO_term in cc:
                    cc_g.add_edge(GO_term, parent_GO_term, weight=relation)
                elif namespace == "mf" and parent_GO_term in mf:
                    mf_g.add_edge(GO_term, parent_GO_term, weight=relation)

    return mf_g, bp_g, cc_g, alt_id_to_id_mapping


def parseCommandLineArguments():
    """
    Iddo: Fix this docstring
    Please do not use this function yet.
    """
    parser = argparse.ArgumentParser(prog="gothresher_prep.py")

    parser.add_argument("--config", "-c", help="INI file. gothresher.ini", default=os.getcwd()+"/gothresher.ini")
    parser.add_argument("--input", "-i", help="path of GO file to be processed", required=True)
#    parser.add_argument("--ontodir","-d",help="GO files directory (default ./onto)", 
#                        default=os.getcwd()+"/onto")
    # parser.add_argument("--output","-o",help="Path of the output file",default="data/parseGO.parsed")

    # parser.add_argument("--keep_obsolete", "-k", help="Enter 1 if you want to keep obsolete terms", default=0)
    args = parser.parse_args()
    return args


def findAllCommonAncestorsAndDisjointCommonAncestors(mf_g, bp_g, cc_g, go1, go2):
    mf, bp, cc = 0, 0, 0
    if mf_g.has_node(go1):
        all_paths_go1 = nx.all_simple_paths(mf_g, source=go1, target=ROOT_MFO)
        mf += 1
    elif bp_g.has_node(go1):
        all_paths_go1 = nx.all_simple_paths(bp_g, source=go1, target=ROOT_BPO)
        bp += 1
    elif cc_g.has_node(go1):
        all_paths_go1 = nx.all_simple_paths(cc_g, source=go1, target=ROOT_CCO)
        cc += 1
    # all_paths_go1=list(all_paths_go1)
    # all_paths_go1.append(go1)

    if mf_g.has_node(go2):
        all_paths_go2 = nx.all_simple_paths(mf_g, source=go2, target=ROOT_MFO)
        mf += 1
    elif bp_g.has_node(go2):
        all_paths_go2 = nx.all_simple_paths(bp_g, source=go2, target=ROOT_BPO)
        bp += 1
    elif cc_g.has_node(go2):
        all_paths_go2 = nx.all_simple_paths(cc_g, source=go2, target=ROOT_CCO)
        cc += 1

    all_paths_go1 = list(all_paths_go1)
    all_paths_go2 = list(all_paths_go2)

    # print((mf, bp, cc))
    if mf == 2:
        current_graph = mf_g
        root_node = ROOT_MFO
    elif bp == 2:
        current_graph = bp_g
        root_node = ROOT_BPO
    elif cc == 2:
        current_graph = cc_g
        root_node = ROOT_CCO

    common_ancestors = []
    for path1 in all_paths_go1:
        for path2 in all_paths_go2:
            print(path1)
            print(path2)
            print((list(set(path1) & set(path2))))
            print(("*" * 50))
            if list(set(path1) & set(path2)) not in common_ancestors:
                common_ancestors.append(list(set(path1) & set(path2)))
            # common_ancestors=list(set(common_ancestors))
    dca = []
    for each_path in common_ancestors:
        node_to_root_distance = dict()
        for node in each_path:
            # if node is not root_node:
            node_to_root_distance[node] = nx.shortest_path_length(current_graph, source=node, target=root_node)
        print(node_to_root_distance)
        if len(node_to_root_distance) != 0:
            key, _ = max(iter(list(node_to_root_distance.items())), key=lambda x: x[1])
            if key not in dca:
                dca.append(key)

    return common_ancestors, dca


def findAllAncestorsForAllNodesForOntology(ontology,ONTO_DIR):
    # graph = cp.load(open("data/" + ontology + ".graph", "rb"))
    print(ontology)
    graph = cp.load(open(f"{ONTO_DIR}/{ontology}.graph", "rb"))
    graph_ancestors = dict()
    for nodes in graph.nodes():
        graph_ancestors[nodes] = findAllAncestors(graph, nodes)
    return graph_ancestors

def init_globals(gothresh_ini_file="gothresher.ini"):
    config = configparser.ConfigParser()
    config.read(gothresh_ini_file)
    print(config.sections())
    global FILE_ALTERNATE_ID_TO_ID_MAPPING 
    global FILE_CAFA_ID_TO_UNIPROT_ID_MAP 
    global FILE_MFO_ONTOLOGY_GRAPH
    global FILE_BPO_ONTOLOGY_GRAPH
    global FILE_CCO_ONTOLOGY_GRAPH
    global FILE_MFO_ONTOLOGY_ANCESTORS_GRAPH
    global FILE_BPO_ONTOLOGY_ANCESTORS_GRAPH
    global FILE_CCO_ONTOLOGY_ANCESTORS_GRAPH

    global ONTO_DIR
    ONTO_DIR = os.path.abspath(config['DIRS']['ONTO_DIR'])
    FILE_ALTERNATE_ID_TO_ID_MAPPING = "alt_to_id.graph"
    FILE_CAFA_ID_TO_UNIPROT_ID_MAP = "CAFAIDTOUniprotIDMap.txt"
    FILE_MFO_ONTOLOGY_GRAPH = "mf.graph"
    FILE_BPO_ONTOLOGY_GRAPH = "bp.graph"
    FILE_CCO_ONTOLOGY_GRAPH = "cc.graph"
    FILE_MFO_ONTOLOGY_ANCESTORS_GRAPH = f"{ONTO_DIR}/mf_ancestors.map"
    FILE_BPO_ONTOLOGY_ANCESTORS_GRAPH = f"{ONTO_DIR}/bp_ancestors.map"
    FILE_CCO_ONTOLOGY_ANCESTORS_GRAPH = f"{ONTO_DIR}/cc_ancestors.map"

    ROOT_BPO = config['GO']['root_bpo']
    ROOT_CCO = config['GO']['root_cco']
    ROOT_MFO = config['GO']['root_mfo']



    
def main():
    options = parseCommandLineArguments()
    init_globals(options.config)
    mf_g, bp_g, cc_g, alt_id_to_id_mapping = parseGOTerms(options.input)
    #print(mf_g.has_node("GO:0019786"))
    #print(bp_g.has_node("GO:0019786"))
    #print(cc_g.has_node("GO:0019786"))
    # if os.path.isdir("data") == False:
    #    os.makedirs("data")
    # Iddofix
    if os.path.exists(ONTO_DIR):
        if not os.path.isdir(ONTO_DIR):
            raise ValueError(f"{ONTO_DIR} is not a directory")
    else:
        sys.stderr.write(f"{ONTO_DIR} does not exist. Creating...\n")
        os.makedirs(ONTO_DIR)
    # Iddo: fixed all paths to ONTO_DIR root path
    cp.dump(mf_g, open(f"{ONTO_DIR}/{FILE_MFO_ONTOLOGY_GRAPH}", "wb"))
    cp.dump(bp_g, open(f"{ONTO_DIR}/{FILE_BPO_ONTOLOGY_GRAPH}", "wb"))
    cp.dump(cc_g, open(f"{ONTO_DIR}/{FILE_CCO_ONTOLOGY_GRAPH}", "wb"))
    cp.dump(alt_id_to_id_mapping, open(f"{ONTO_DIR}/{FILE_ALTERNATE_ID_TO_ID_MAPPING}", "wb"))

    mf_g = cp.load(open(f"{ONTO_DIR}/{FILE_MFO_ONTOLOGY_GRAPH}", "rb"))
    bp_g = cp.load(open(f"{ONTO_DIR}/{FILE_BPO_ONTOLOGY_GRAPH}", "rb"))
    cc_g = cp.load(open(f"{ONTO_DIR}/{FILE_CCO_ONTOLOGY_GRAPH}", "rb"))

    mf_ancestors = findAllAncestorsForAllNodesForOntology("mf",ONTO_DIR)
    bp_ancestors = findAllAncestorsForAllNodesForOntology("bp",ONTO_DIR)
    cc_ancestors = findAllAncestorsForAllNodesForOntology("cc",ONTO_DIR)

    cp.dump(mf_ancestors, open(FILE_MFO_ONTOLOGY_ANCESTORS_GRAPH, "wb"))
    cp.dump(bp_ancestors, open(FILE_BPO_ONTOLOGY_ANCESTORS_GRAPH, "wb"))
    cp.dump(cc_ancestors, open(FILE_CCO_ONTOLOGY_ANCESTORS_GRAPH, "wb"))

    # print(findAllAncestors(mf_g, "GO:0019786"))
    """ca,dca=findAllCommonAncestorsAndDisjointCommonAncestors(mf_g, bp_g, cc_g, "GO:0035556", "GO:0009966")
    for eachca in ca:
        print(eachca)
    print("DCA")
    print(dca)"""


if __name__ == "__main__":
    main()

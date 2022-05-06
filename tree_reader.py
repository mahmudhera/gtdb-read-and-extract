from Bio import Phylo
from matplotlib import pyplot as plt
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceMatrix

#plt.figure(figsize=(8, 6000), dpi=80)

if __name__ == "__main__":
    tree = Phylo.read("bac120_r207.tree", "newick")
    leaf_nodes = ['RS_GCF_016405145.1', 'RS_GCF_001746855.1', 'GB_GCA_017934115.1', 'RS_GCF_015240065.1', 'GB_GCA_009784895.1', 'RS_GCF_003003235.1', 'RS_GCF_014204945.1', 'GB_GCA_015711725.1', 'GB_GCA_018435605.1']
    for l in leaf_nodes:
        print(l)

    matrix = []
    for i in range(len(leaf_nodes)):
        distance_list = []
        for j in range(i+1):
            distance_list.append( tree.distance(leaf_nodes[i], leaf_nodes[j]) )
        matrix.append(distance_list)

    dm = DistanceMatrix( names=leaf_nodes, matrix=matrix )
    print(dm)

    constructor = DistanceTreeConstructor()
    tree_new = constructor.nj(dm)
    for clade in tree_new.find_clades():
        if "Inner" in clade.name:
            clade.name = ''

    fig, ax = plt.subplots()
    ax.axis("off")
    Phylo.draw(tree_new, do_show=False)
    plt.savefig('phylogenetic_tree_selected.pdf')

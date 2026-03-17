import sys
import re
import os
import os.path
import gzip
import glob
import shutil
import argparse
import subprocess
import requests
import json
import dendropy

num_ancestral_nodes = 3 # number of nodes to step down to find outgroups
num_outgroups_per_ancestral_node = 2 
inclusive_from_mrca = False

gtdb_base_url = "https://data.gtdb.aau.ecogenomic.org/"
gtdb_base_url = "https://data.ace.uq.edu.au/public/gtdb/data/"
bvbrc_base_url="https://www.bv-brc.org/api-internal/";

def download_file(url):
    #https://stackoverflow.com/questions/16694907/download-a-large-file-in-python-with-requests
    local_filename = url.split('/')[-1]
    # NOTE the stream=True parameter below
    print(f"download_file {url} to {local_filename}")
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192): 
                f.write(chunk)
    return local_filename

def download_latest_gtdb_data(url_base):
    url_base += f"releases/latest/"
    retval = requests.get(url_base+"/VERSION.txt")
    version = None
    if retval.status_code == 200:
        m = re.search("v([\d.]+)", retval.text)
        if m:
            version = float(m.group(1))
    else:
        raise Exception(f"download_latest_gtdb_data got status: {retval.status_code}")
    print(f"download_latest_gtdb_data returns {version}")
    if not version:
        raise Exception("cannot read latest version ID from {url_base}")
    local_gtdb_dir = f"gtdb_release_{version}"
    if os.path.exists(local_gtdb_dir): # get rid of existing data if any
        raise Exception("latest directory already exists: gtdb_version_{version}")
    os.mkdir(local_gtdb_dir)
    cur_dir = os.getcwd()
    os.chdir(local_gtdb_dir)
    print(f"local gtdb dir = {local_gtdb_dir}")
    download_file(url_base+"bac120.tree")
    download_file(url_base+"bac120_metadata.tsv.gz")
    download_file(url_base+"ar53.tree")
    download_file(url_base+"ar53_metadata.tsv.gz")
    download_file(url_base+"VERSION.txt")
    os.chdir(cur_dir)
    return local_gtdb_dir;

def integrate_bvbrc_gtdb_data(local_gtdb_dir):
    cur_dir = os.getcwd()
    os.chdir(local_gtdb_dir)
    gtdb_assembly = {} # dict of gtdb_id to genbank assembly id
    assembly_gtdb = {} # dict of assembly to gtdb_id
    gtdb_species = {}
    total_lines = 0
    num_matches = 0
    for metadata_file in ["bac120_metadata.tsv.gz", "ar53_metadata.tsv.gz"]:
        with gzip.open(metadata_file, "rt") as F:
            print("reading "+metadata_file)
            header = F.readline().rstrip().split("\t")
            if not (header[0] == 'accession' and header[17] == 'gtdb_genome_representative' and header[57] == 'ncbi_genbank_assembly_accession'):
                raise Exception("headers not as expected")
            for line in F:
                total_lines += 1
                fields = line.rstrip().split("\t")
                gtdb_id = fields[0]
                is_gtdb_ref = fields[18]
                assembly = fields[57]
                if (is_gtdb_ref == 't'):
                    gtdb_assembly[gtdb_id] = assembly
                    assembly_gtdb[assembly] = gtdb_id
                    species = fields[19].split('__')[-1]
                    gtdb_species[gtdb_id] = species
                    num_matches += 1
                    if num_matches < 5:
                        print(f"match gt={gtdb_id} as={assembly} sp={species}")
        print(f"number of gtdb reference entries = {len(gtdb_assembly)}")
        print(f"total lines was {total_lines}")

    bvbrc_gtdb = {}
    gtdb_bvbrc = {}
    division = {}
    num_matches = 0
    total_lines = 0
    cmd_base = "p3-all-genomes -r assembly_accession -a assembly_accession -e superkingdom,"
    for div in ['Bacteria', 'Archaea']:
        print(f"run {cmd_base}{div}")
        i = 0
        result = subprocess.run(cmd_base+div, shell=True, stdout=subprocess.PIPE, text=True)
        for line in result.stdout.split("\n"):
            (bvbrc_id, assembly, species) = ('na', 'na', 'na')
            fields = line.rstrip().split("\t")
            total_lines += 1
            bvbrc_id = fields[0]
            if len(fields) > 1:
                assembly = fields[1]
                if assembly in assembly_gtdb:
                    gtdb_id = assembly_gtdb[assembly]
                    bvbrc_gtdb[bvbrc_id] = gtdb_id
                    gtdb_bvbrc[gtdb_id] = bvbrc_id
                    division[bvbrc_id] = div
                    num_matches += 1
                    if num_matches < 4:
                        print(f"got match, {bvbrc_id}\t{gtdb_id}\t{assembly}")
            else:
                print("line = "+line)
                print(f"bi={bvbrc_id}\tas={assembly}")

    with open("bvbrc_gtdb_ids.tsv", "w") as F:
        F.write("bvbrc_id\tgtdb_id\tdivision\n")
        for bvbrc_id in sorted(bvbrc_gtdb):
            gtdb_id = bvbrc_gtdb[bvbrc_id]
            F.write(f"{bvbrc_id}\t{gtdb_id}\t{division[bvbrc_id]}\n")
        F.close()

    cmd = "cut -f 1 bvbrc_gtdb_ids.tsv | p3-get-genome-data -a taxon_lineage_ids > bvbrc_taxonomy.tsv"
    print("run "+cmd)
    subprocess.run(cmd, shell=True)
    print(f"number of matches = {num_matches}")

    os.chdir(cur_dir)

def test_data_integrity(data_dir):
    result = True
    needed_files = ['bac120.tree', 'ar53.tree', 'bac120_metadata.tsv.gz', 'ar53_metadata.tsv.gz', 'VERSION.txt', 'bvbrc_gtdb_ids.tsv', 'bvbrc_taxonomy.tsv']
    for f in needed_files:
        if not os.path.isfile(f"{data_dir}/{f}"):
            result = False
            print(f"Needed file {f} not found in {data_dir}")
    if not result:
        raise Exception("Not all needed file found. Exiting.")

def read_gtdb_bvbrc(data_dir):
    gtdb_bvbrc = {}
    bvbrc_division = {}
    with open(data_dir+"/bvbrc_gtdb_ids.tsv") as F:
        for line in F:
            fields = line.rstrip().split("\t")
            bvbrc_id = fields[0]
            gtdb_id = fields[1]
            division = fields[2]
            gtdb_bvbrc[gtdb_id] = bvbrc_id
            bvbrc_division[bvbrc_id] = division
    return gtdb_bvbrc, bvbrc_division

def read_bvbrc_taxonomy(data_dir):
    taxon_genomes = {}
    taxon_level = {}
    taxon_division = {}
    taxon_subtaxa = {}
    taxon_parent = {}
    with open(data_dir+"/bvbrc_taxonomy.tsv") as F:
        header = F.readline()
        for line in F:
            genome_id, taxonomy = line.rstrip().split("\t")
            taxa = taxonomy.split("::")
            for level, taxon in enumerate(taxa):
                #if level < 4:
                #    continue
                if not taxon in taxon_genomes:
                    taxon_genomes[taxon] = set()
                    taxon_level[taxon] = level
                    #taxon_division[taxon] = taxa[1] # 2=Bacteria, 2157=Archea
                    taxon_subtaxa[taxon] = set()
                taxon_genomes[taxon].add(genome_id)
                if (level+1) < len(taxa):
                    taxon_parent[taxa[level+1]] = taxon
                    for level2 in range(level+1, len(taxa)):
                        taxon_subtaxa[taxon].add(taxa[level2])
    return(taxon_genomes, taxon_subtaxa, taxon_parent)

def get_taxon_name_rank(taxon_set):
    get_taxonomy_data_cmd = "p3-get-taxonomy-data -a taxon_name,taxon_rank"
    taxon_name = {}
    taxon_rank = {}
    taxon_list = list(taxon_set)
    stride = 300
    start = 0
    while start < len(taxon_list):
        sublist = taxon_list[start : start+stride]
        print(f"get range {start}: {start+stride}")
        start += stride
        proc = subprocess.Popen(get_taxonomy_data_cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        #for taxon in taxon_set:
        proc.stdin.write("\n".join(sublist)+"\n")
        proc.stdin.close()
        proc.wait()
        for line in proc.stdout:
            (taxon_id, name, rank) = line.rstrip().split("\t")
            taxon_name[taxon_id] = name
            taxon_rank[taxon_id] = rank
            #print(f"tid:{taxon_id}, n:{name}, r:{rank}")
    return taxon_name, taxon_rank

def get_id_from_taxon_name(taxon_name):
    Session = requests.Session()
    Session.headers.update({ 'accept': "application/json" })
    Session.headers.update({ "Content-Type": "application/rqlquery+x-www-form-urlencoded" })

    query = f"eq(taxon_name,{taxon_name})"
    query += "&select(taxon_id,taxon_name)"
    query += "&limit(100000)"
    response = Session.get(bvbrc_base_url+"taxonomy", params=query)
    print(response.url+"\n")
    print(f"response={response}\n")
    retval = ''
    numrows = 0
    for record in response.json():
        numrows += 1
        if record['taxon_name'] == taxon_name:
            print(f"{record['taxon_id']}\t{record['taxon_name']}")
            retval = record['taxon_id']
            break
    print(f"numrows = {numrows}")
    print("returning: "+retval)
    return(retval)

    

def filter_taxon_name(taxon_name, valid_taxon_set=None):
    if not valid_taxon_set:
        valid_taxon_set = set(taxon_name.keys())
    for taxon_id in taxon_name:
        name = taxon_name[taxon_id]
        if re.match("environmental", name):
            valid_taxon_set.discard(taxon_id)
        elif re.match("unclassified", name):
            valid_taxon_set.discard(taxon_id)
        elif re.match("uncultured", name):
            valid_taxon_set.discard(taxon_id)
    return valid_taxon_set

def filter_taxon_rank(taxon_rank, allowed_ranks, valid_taxon_set=None):
    if not valid_taxon_set:
        valid_taxon_set = set(taxon_rank.keys())
    for taxon_id in taxon_rank:
        if not taxon_rank[taxon_id] in allowed_ranks:
            valid_taxon_set.discard(taxon_id)
    return valid_taxon_set

def filter_taxon_size(taxon_genomes, valid_taxon_set=None, min_size=0, max_size=1000000):
    if max_size < min_size:
        raise Exception(f"max_size {max_size} < min_size {min_size}")
    if not valid_taxon_set:
        valid_taxon_set = set(taxon_genomes.keys())
    for taxon_id in taxon_genomes:
        num_genomes = len(taxon_genomes[taxon_id])
        valid = num_genomes >= min_size
        if max_size:
            valid &= (num_genomes <= max_size)
        if not valid:
            valid_taxon_set.discard(taxon_id)
    return valid_taxon_set

def read_gtdb_tree(tree_file, gtdb_bvbrc_map):
    print(f"read tree file: {tree_file}")
    tree = None
    tree = dendropy.Tree.get(path=tree_file, schema='newick')
    tree.is_rooted = True
    num_tips = 0
    num_replacements = 0
    bvbrc_tip_labels = set()
    for leaf in tree.leaf_nodes():
        num_tips += 1
        gtdb_id = leaf.taxon.label.replace(' ', '_') # undo modification done by dendropy
        # Check if the current label is in the map
        if gtdb_id in gtdb_bvbrc_map:
            leaf.taxon.label = gtdb_bvbrc_map[gtdb_id]
            leaf.label = gtdb_bvbrc_map[gtdb_id]
            #print("set label to "+leaf.label)
            bvbrc_tip_labels.add(leaf.label)
            num_replacements += 1
            # It may also be useful to update the node label itself if it's used
            # leaf.label = label_map[leaf.label]
    print(f"    gtdb tree from {tree_file}\n\tnum tips = {num_tips}\n\tnum_replacements = {num_replacements}")
    return tree, bvbrc_tip_labels

def extract_subtree(tree, tip_ids, genome_priority=None, num_outgroups_per_ancestral_node=2, num_ancestral_nodes=2, target_tree_size=0):
    # use dendropy
    #for tip in tip_ids:
    #    node = tree.find_node_with_label(tip)
    #    print(f"find node with label {tip} yields {node}")
    mrca = tree.mrca(taxon_labels = tip_ids)
    print(f"found mrca: {mrca}")
    rep_outgroups = set()
    other_outgroups = set()
    cur_anc = mrca
    for anc_node in mrca.ancestor_iter():
        #print("got anc {}".format(anc_node))
        #child_nodes = anc_node.child_nodes()
        potential_rep_outgroups = set()
        potential_other_outgroups = set()
        less_good_outgroups = set()
        for child in anc_node.child_nodes():
            #print("try anc child {}".format(child))
            if child != cur_anc:
                desc_nodes = child.leaf_nodes()
                #print("child not anc, num leafs={}".format(len(desc_nodes)))
                for node in desc_nodes:
                    bvbrc_id = node.taxon.label
                    if genome_priority and bvbrc_id in genome_priority:
                        if genome_priority[bvbrc_id] > 0:
                            potential_rep_outgroups.add(node.taxon.label)
                        elif genome_priority[bvbrc_id] < 0:
                            less_good_outgroups.add(node.taxon.label)
                        else:
                            potential_other_outgroups.add(node.taxon.label)
                    else:
                        potential_other_outgroups.add(node.taxon.label)
                    #print("adding pot out {}".format(node.taxon.label))
        for i, label in enumerate(potential_rep_outgroups):
            rep_outgroups.add(label)
            if i >= num_outgroups_per_ancestral_node:
                break
        if len(rep_outgroups) > num_outgroups_per_ancestral_node * num_ancestral_nodes:
            break
        for i, label in enumerate(potential_other_outgroups):
            other_outgroups.add(label)
            if i >= num_outgroups_per_ancestral_node:
                break
        cur_anc = anc_node

    tip_ids.update(rep_outgroups)
    num_other_outgroups_to_use = (num_outgroups_per_ancestral_node * num_ancestral_nodes) - len(rep_outgroups)
    for i, label in enumerate(other_outgroups):
        if i >= num_other_outgroups_to_use:
            break
        tip_ids.add(label)

    sys.stderr.write("after adding outgroups, num tips is {}\n".format(len(tip_ids)))
    if target_tree_size and (len(tip_ids) < target_tree_size): # try to fill out tree with additional members to get to target size
        sys.stderr.write(f"add additional members from mrca of representative members\n")
        for node in mrca.leaf_nodes():
            label = node.taxon.label
            if not label in tip_ids:
                tip_ids.add(label)
                if len(tip_ids) >= target_tree_size:
                    break

    subtree = tree.extract_tree_with_taxa_labels(tip_ids)
    return subtree



def main():
    print("parse arguments")
    parser = argparse.ArgumentParser(description="Exctract from the large GTDB phylogeny subtrees corresponding to NCBI taxa.", formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    parser.add_argument("--download_latest_gtdb_release", action='store_true', help="Look online for GTDB latest release and download trees and metadata to 'gtdb_release_XXX'")
    parser.add_argument("--gtdb_data_url", metavar="URL", type=str, default=gtdb_base_url,  help="URL for top of data, must contain 'releases' folder.")
    parser.add_argument("--link_bvbrc_gtdb_ids", action="store_true", help="given GTDB metadata, use Genbank accession ID to link to BVBRC genome IDs")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--gtdb_data_dir", metavar="directory", type=str, help="Optional, searches for latest version: gtdb_release_XXX")
    #parser.add_argument("--division", metavar="bacteria|archaea", choices = ['bacteria', 'archaea'], type=str, default="bacteria", help="Analyze Bacteria or Archaea.")
    parser.add_argument("--target_taxon", metavar="name or ID", type=str, help="Extract this one taxon.")
    parser.add_argument("--target_taxon_file", metavar="file", type=str, help="File with list of taxon IDs for extraction.")
    parser.add_argument("--encompassing_taxon", metavar="name or ID", type=str, help="Consider extracting only subtaxa of this.")
    parser.add_argument("--ranks", metavar="rank1 [rank2 ...]", nargs='*', type=str, default=['genus', 'family', 'order', 'class', 'phylum'], help=" List of ranks to extract.")
    #parser.add_argument("--all_taxa", action='store_true', help="Generate trees for all taxa
    #parser.add_argument("--target_tree_size", metavar="NUM", type=int, default=20, help="Target number of tips on tree.")
    parser.add_argument("--min_taxon_size", metavar="NUM", type=int, default=0, help="Minimum size of taxon.")
    parser.add_argument("--max_taxon_size", metavar="NUM", type=int, default=0, help="Maximum size of taxon (0=no limit).")
    parser.add_argument("--num_ancestral_nodes", metavar="NUM", type=int, default=2, help="How many nodes to step down to find outgroups.")
    parser.add_argument("--outgroups_per_anc_node", metavar="NUM", type=int, default=2,help="How many outgroup tips to select from each ancestral node.")

    args = parser.parse_args()
    print(f"parsed args: args=\n{args}")
    #starttime = time()

    if args.download_latest_gtdb_release:
        args.gtdb_data_dir = download_latest_gtdb_data(args.gtdb_data_url)
        args.link_bvbrc_gtdb_ids = True

    if not args.gtdb_data_dir:
        latest_version = 0
        dirlist = glob.glob("gtdb_release_*")
        for d in dirlist:
            m = re.search("gtdb_release_(.*)", d)
            version = float(m.group(1))
            if version > latest_version:
                latest_version = version
                args.gtdb_data_dir = d
        if not args.gtdb_data_dir:
            raise Exception("No gtdb release directory found (perhaps run with --download_latest_gtdb_release)")

    if args.link_bvbrc_gtdb_ids:
        integrate_bvbrc_gtdb_data(args.gtdb_data_dir)

    test_data_integrity(args.gtdb_data_dir)

    gtdb_bvbrc, genome_division = read_gtdb_bvbrc(args.gtdb_data_dir)
    taxon_genomes, taxon_subtaxa, taxon_parent = read_bvbrc_taxonomy(args.gtdb_data_dir)

    division_tree = {}
    division_genomes = {}
    division_tree_file = {"Bacteria": 'bac120.tree', "Archaea": 'ar53.tree', '2': 'bac120.tree', '2157': 'ar53.tree'}
    taxon_set = set()
    taxon_name = None;
    taxon_rank = None;
    if (args.target_taxon or args.target_taxon_file):
        if args.target_taxon:
            if not args.target_taxon in taxon_parent:
                raise Exception(f"specified target taxon {args.target_taxon} not found.")
            taxon_set.add(args.target_taxon)
        if args.target_taxon_file:
            print("reading taxa IDs from "+args.target_taxon_file)
            with open(args.target_taxon_file) as F:
                for line in F:
                    taxon_id = line.split()[0]
                    if not args.target_taxon in taxon_parent:
                        raise Exception(f"specified target taxon {taxon_id} from file {args.target_taxon_file} not found.")
                    taxon_set.add(taxon_id)
        taxon_name, taxon_rank = get_taxon_name_rank(taxon_set)
    
    else:
        if args.encompassing_taxon:
            if not args.encompassing_taxon in taxon_subtaxa:
                raise Exception(f"specified encompassing taxon {args.encompassing_taxon} not found.")
            taxon_set.update(taxon_subtaxa[args.encompassing_taxon])
            print(f"after restricting to descendant taxa of {args.encompassing_taxon}, num in taxon_set = {len(taxon_set)}")

        else:
            taxon_set = set(taxon_genomes.keys())
            print(f"after adding all taxa, num in taxon_set = {len(taxon_set)}")

        taxon_set = filter_taxon_size(taxon_genomes, taxon_set, min_size=args.min_taxon_size, max_size=args.max_taxon_size)
        print(f"after filter taxon size: num in taxon_set = {len(taxon_set)}")
        taxon_name, taxon_rank = get_taxon_name_rank(taxon_set)
        taxon_set = filter_taxon_name(taxon_name, taxon_set)
        print(f"after filter taxon name: num in taxon_set = {len(taxon_set)}")
        taxon_set = filter_taxon_rank(taxon_rank, args.ranks, taxon_set)
        print(f"after filter taxon rank: num in taxon_set = {len(taxon_set)}")

    newick_directory = args.gtdb_data_dir + "/newick_trees"
    if not os.path.exists(newick_directory):
        os.mkdir(newick_directory)
    phyloxml_directory = args.gtdb_data_dir + "/phyloxml_trees"
    if not os.path.exists(phyloxml_directory):
        os.mkdir(phyloxml_directory)
    taxon_tree_file = {}
    for taxon in taxon_set:
        division = None
        if taxon in taxon_subtaxa['2']:
            division = '2'
        elif taxon in taxon_subtaxa['2157']:
            division = '2157'
        if not division in division_tree:
            tree_file = args.gtdb_data_dir + "/" + division_tree_file[division]
            division_tree[division], division_genomes[division] = read_gtdb_tree(tree_file, gtdb_bvbrc)
        tree = division_tree[division]
        tree_genomes = division_genomes[division]
        available_tree_tips = (taxon_genomes[taxon] & tree_genomes)
        #if len(available_tree_tips) < len(taxon_genomes[taxon]):
        print(f"for taxon {taxon}, {taxon_name[taxon]} {len(available_tree_tips)} of {len(taxon_genomes[taxon])} tips available")
        subtree = extract_subtree(tree, available_tree_tips,num_outgroups_per_ancestral_node=args.outgroups_per_anc_node, num_ancestral_nodes=args.num_ancestral_nodes)
        rank = taxon_rank[taxon]
        name = taxon_name[taxon].replace(' ', '_')
        name = name.replace('(', '-')
        name = name.replace(')', '-')
        newick_file = f"taxon_{taxon}_{rank}_{name}_gtdb_subtree.nwk"
        print(f"save tree as {newick_file}")
        with open(f"{newick_directory}/{newick_file}", 'w') as F:
            subtree.write(file=F, schema="newick", suppress_rooting=True)
        # now run p3x-newick-to-phyloxml
        phyloxml_file = newick_file[:-4]+".phyloxml"
        cmd = f"p3x-newick-to-phyloxml -l genome_id "
        cmd += " --genomefields genome_name"
        if rank == 'genus':
            cmd += ",species"
        elif rank == 'family':
            cmd += ",genus,species"
        elif rank == 'order':
            cmd += ",family,genus,species"
        elif rank == 'class':
            cmd += ",order,family,genus,species"
        elif rank == 'phylum':
            cmd += ",class,order,family,genus,species"
        cmd += f" {newick_directory}/{newick_file}"
        print("run "+cmd)
        subprocess.run(cmd, shell=True)
        if os.path.exists(f"{phyloxml_directory}/{phyloxml_file}"):
            os.remove(f"{phyloxml_directory}/{phyloxml_file}")
        print(f"move {newick_directory}/{phyloxml_file} {phyloxml_directory}")
        shutil.move(f"{newick_directory}/{phyloxml_file}", f"{phyloxml_directory}")
        taxon_tree_file[taxon] = phyloxml_file

    # now write taxon-to-tree_file dict as json file
    with open(args.gtdb_data_dir + "/taxon_tree_dict.json", 'w') as F:
        F.write(json.dumps(taxon_tree_file, indent=2))
        F.close()


if __name__ == "__main__":
    main()

cp $1 $1.0

a=$(basename "$1")
b="${a%.*}"

search_terms=("Did not match ontology id or number for term")
search_terms+=("Invalid Ontology ID")
search_terms+=("Skipping invalid vertex name")
search_terms+=("Skipping invalid subject type in triple")
search_terms+=("Skipping invalid predicate type in triple")
search_terms+=("Skipping invalid object type in triple")

l=0

wc -l ${b}_triples.txt
wc -l ${b}_fnode_triples.txt
wc -l ${b}_bnode_triples.txt

for i in ${!search_terms[@]}; do

    j=$(($l + 1))
    k=$(($j + 1))

    cat $1.$l | grep "${search_terms[$i]}" > $1.$j
    cat $1.$l | grep -v "${search_terms[$i]}" > $1.$k

    l=$(($l + 2))

    echo ${search_terms[$i]}
    wc -l $1.$j

done

cat $1.5 | cut -d " " -f 5 | sort | uniq > $1.5_uniq



# cat load_bioportal_cell_ontology.log | grep 'Skipping non-literal object in triple' > load_bioportal_cell_ontology.log.1
# cat load_bioportal_cell_ontology.log | grep -v 'Skipping non-literal object in triple' > load_bioportal_cell_ontology.log.2

# cat load_bioportal_cell_ontology.log.2 | grep 'Skipping invalid predicate type in triple' > load_bioportal_cell_ontology.log.3
# cat load_bioportal_cell_ontology.log.2 | grep -v 'Skipping invalid predicate type in triple' > load_bioportal_cell_ontology.log.4

# cat load_bioportal_cell_ontology.log.4 | grep 'Skipping invalid vertex name' > load_bioportal_cell_ontology.log.5
# cat load_bioportal_cell_ontology.log.4 | grep -v 'Skipping invalid vertex name' > load_bioportal_cell_ontology.log.6

# cat load_bioportal_cell_ontology.log.6 | grep 'Skipping invalid subject type in triple' > load_bioportal_cell_ontology.log.7
# cat load_bioportal_cell_ontology.log.6 | grep -v 'Skipping invalid subject type in triple' > load_bioportal_cell_ontology.log.8

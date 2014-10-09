md=ms/maintext.md
pdf=ms/probabilistic_measures.pdf

ALL: $(pdf) data

$(pdf): $(md)
	pandoc $< -o $@ --template=ms/template.tex

data:
	wget http://datadryad.org/bitstream/handle/10255/dryad.56778/Interactions_matrices.txt?sequence=1 -O data/raw.txt
	cd data; csplit -s raw.txt '/P/' '{*}' -f hp -n 3
	rm data/hp000
	sed -i -e "$$ d" data/hp*
	sed -i -r 's/\S+(\s+)?//1' data/hp*
	sed -i -e "1d" data/hp*

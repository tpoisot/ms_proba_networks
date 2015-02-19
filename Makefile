md=ms/maintext.md
pdf=ms/probabilistic_measures.pdf
refs=/home/tpoisot/.pandoc/default.bib
csl=/home/tpoisot/.pandoc/styles/methods-in-ecology-and-evolution.csl

ALL: $(pdf) data

$(pdf): $(md)
	pandoc $< -o $@ --template=ms/template.tex --bibliography=$(refs) --csl=$(csl)

$(refs): ms/bib.keys
	python2 ms/extractbib.py ms/bib.keys /home/tpoisot/.pandoc/default.bib $(refs)

ms/bib.keys: $(md)
	grep @[-:_a-zA-Z0-9]* $(md) -oh --color=never | sort | uniq | sed 's/@//g' > ms/bib.keys

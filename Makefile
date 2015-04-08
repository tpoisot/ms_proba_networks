md=ms/maintext.md
pdf=ms/probabilistic_measures.pdf
refs=/home/tpoisot/.pandoc/default.bib
csl=/home/tpoisot/.pandoc/styles/methods-in-ecology-and-evolution.csl
pflags= --template=ms/template.tex --bibliography=$(refs) --csl=$(csl)

ALL: $(pdf) diff.pdf

$(pdf): $(md)
	pandoc $< -o $@ $(pflags)

$(refs): ms/bib.keys
	python ms/extractbib.py ms/bib.keys /home/tpoisot/.pandoc/default.bib $(refs)

ms/bib.keys: $(md)
	grep @[-:_a-zA-Z0-9]* $(md) -oh --color=never | sort | uniq | sed 's/@//g' > ms/bib.keys

diff.pdf: $(md)
	wget -O sub1.md https://raw.githubusercontent.com/tpoisot/ms_proba_networks/master/ms/maintext.md
	pandoc $< -o rev.tex $(pflags)
	pandoc sub1.md -o sub.tex $(pflags)
	latexdiff sub.tex rev.tex > diff.tex
	pdflatex diff
	rm *.tex
	rm sub1.md
	rm diff.aux
	rm diff.log
	rm diff.out

md=ms/maintext.md
pdf=ms/probabilistic_measures.pdf

ALL: $(pdf)

$(pdf): $(md)
	pandoc $< -o $@ --template=ms/template.tex

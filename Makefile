# Clear the suffix list; no suffix rules in this makefile.  See section 7.2.1 of the GNU
# Coding Standards.
.SUFFICES:

# The default target.
protocol.html: protocol.md pandoc.css
	pandoc -o protocol.html protocol.md --css pandoc.css --toc --toc-depth=4 --mathjax

protocol.pdf: protocol.md pandoc.css
	pandoc -o protocol.pdf protocol.md latex.yaml --latex-engine xelatex --toc --toc-depth=4 --highlight-style haddock

.PHONY: clean

clean:
	$(RM) protocol.html protocol.pdf

# vim: tw=90 ts=8 sts=-1 sw=3 noet

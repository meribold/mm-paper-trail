protocol.html: protocol.md pandoc.css
	pandoc -o protocol.html protocol.md --css pandoc.css --toc --toc-depth=4 --mathjax

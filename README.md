# Statistical Power Analysis and Related Issues in Human Genetic Linkage and Association

NB

Nowadays, it is typical to use \epsfig{file=foo,a=b,x=y} rather than \includegraphics[a=b,x=y]{foo}.

The original 2001.dvi is sufficient to regenerate 2001.pdf with better quality text (OCR). Under Fedora 28 it is possible with 

```bash
latex 2001
bibtex 2001
latex 2001
latex 2001
dvips 2001
dvipdf 2001
```
while the much simplified `pdflatex` does not work well with figures in PostScript (.ps) format. File `epsf.tex` is copied here so the same could be done under Ubuntu 18.04. 


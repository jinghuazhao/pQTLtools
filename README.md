# Statistical Power Analysis and Related Issues in Human Genetic Linkage and Association

NB

It is now typical to use \epsfig{file=foo,a=b,x=y} instead of \includegraphics[a=b,x=y]{foo}. Nevertheless the original 2001.dvi (.ps) is sufficient to regenerate 2001.pdf with better quality (OCR) text. Under Fedora 28 it is to reproduce with 

```bash
# .dvi --> .pdf
dvipdf 2001
# .dvi --> .ps and then .pdf with smaller size from ps2pdf
dvips 2001
ps2pdf 2001.ps
# .tex, .bib --> dvi
latex 2001
bibtex 2001
latex 2001
latex 2001
```
while the much simplified `pdflatex` does not recognise figures in PostScript (.ps) format. File `epsf.tex` copied here enables the same to be done with Ubuntu 18.04. 

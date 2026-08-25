@echo off
xelatex -interaction=nonstopmode %1
biber "%~n1"
xelatex -interaction=nonstopmode %1
xelatex -interaction=nonstopmode %1
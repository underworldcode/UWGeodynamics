#!/bin/bash
jupyter nbconvert --to markdown --template=markdown.tpl ../manual/*.ipynb
mv ../manual/*.md ../docs/wiki/.

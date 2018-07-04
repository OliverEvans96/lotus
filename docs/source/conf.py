#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import sphinx_rtd_theme

subprocess.call('cd ..; doxygen', shell=True)

html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

def setup(app):
    app.add_stylesheet("main_stylesheet.css")

extensions = ['breathe', 'sphinx.ext.mathjax', 'sphinx.ext.todo']
breathe_projects = { 'lotus': '../xml' }
breathe_default_project = "lotus"
templates_path = ['_templates']
html_static_path = ['_static']
source_suffix = '.rst'
master_doc = 'index'
project = 'lotus'
copyright = '2018, Oliver Evans'
author = 'Oliver Evans'

html_logo = ''

exclude_patterns = []
highlight_language = 'c++'
pygments_style = 'sphinx'
todo_include_todos = True
todo_link_only = True
htmlhelp_basename = 'lotus-doc'

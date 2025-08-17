# -- Path setup --------------------------------------------------------------

import os
import sys
sys.path.insert(0, os.path.abspath('..'))  # 让 Sphinx 可以找到你的项目模块

# -- Project information -----------------------------------------------------

project = 'MyProject'       # 项目名称
author = 'Your Name'        # 作者
release = '0.1.0'           # 版本号

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',    # 自动生成 API 文档
    'sphinx.ext.viewcode',   # 在文档中显示源代码
]

templates_path = ['_templates']
exclude_patterns = []

language = 'en'  # 文档语言

# -- Options for HTML output -------------------------------------------------

html_theme = 'alabaster'          # 简单默认主题
html_static_path = ['_static']    # 静态文件目录

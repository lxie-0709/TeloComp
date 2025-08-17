# -*- coding: utf-8 -*-
import os
import sys
from unittest.mock import Mock
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

# ----------------------------
# Mock 外部依赖包，避免 RTD 构建失败
# ----------------------------
MOCK_MODULES = [
    'numpy', 'pandas', 'scipy', 'matplotlib', 'pyBigWig', 'h5py', 'sklearn'
]

for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = Mock()

autodoc_mock_imports = MOCK_MODULES

# ----------------------------
# Python 模块路径
# ----------------------------
sys.path.insert(0, os.path.abspath('..'))  # 项目根目录

# ----------------------------
# 项目信息
# ----------------------------
project = 'TeloComp'
author = 'Shoubian Huang && Liang Xie'
release = '0.1.0'
version = release

# ----------------------------
# Sphinx 配置
# ----------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
]

source_suffix = '.rst'
master_doc = 'index'
language = 'en'
exclude_patterns = []

# ----------------------------
# HTML 输出配置
# ----------------------------
# 使用 Sphinx 默认主题
html_theme = 'alabaster'

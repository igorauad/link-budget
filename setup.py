import re
import sys
from pathlib import Path
from setuptools import setup, find_packages

if sys.version_info[0] < 3:
    raise SystemExit("Error: the link-budget tool requires Python 3")

version = re.search(r'^__version__\s*=\s*"(.*)"',
                    open('linkbudget/cli.py').read(), re.M).group(1)

long_description = (Path(__file__).parent / "README.md").read_text()

setup(name="link-budget",
      packages=find_packages(),
      entry_points={"console_scripts": ['link-budget = linkbudget.cli:main']},
      version=version,
      description="Link budget calculator for satellite and radar systems",
      long_description=long_description,
      long_description_content_type='text/markdown',
      author="Igor Freire",
      author_email="igor@blockstream.com",
      url="https://github.com/igorauad/link-budget",
      classifiers=[
          'Programming Language :: Python :: 3',
          "Programming Language :: Python :: 3 :: Only",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
      ],
      install_requires=[
          "numpy>=1.19.4",
          "itur>=0.3.2"
      ],
      python_requires='>=3')

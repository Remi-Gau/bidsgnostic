#!/usr/bin/env python
from __future__ import annotations

import versioneer


import setuptools

setuptools.setup(
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    include_package_data=True,
)

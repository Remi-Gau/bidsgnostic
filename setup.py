#!/usr/bin/env python
from __future__ import annotations

import setuptools
import versioneer

setuptools.setup(
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    include_package_data=True,
)

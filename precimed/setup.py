from __future__ import print_function

from setuptools import setup

setup_kwargs = dict(
    name='precimed',
    packages=['mixer'],
    install_requires=['numpy', 'numdifftools'],  # dependencies
)

if __name__ == '__main__':
    setup(**setup_kwargs)

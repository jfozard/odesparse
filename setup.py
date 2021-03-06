#!/usr/bin/env python

from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('odesparse', parent_package, top_path)

    config.add_library('lsodes',
                       sources=[join('lsodes','*.f')])

    # odepack
    libs = ['lsodes']
    
    config.add_extension('_odesparse',
                         sources=['_odesparsemodule.c'],
                         libraries=libs,
                         depends=['__odesparse.h','multipack.h'])

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())

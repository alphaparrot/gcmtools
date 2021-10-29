from setuptools import setup
import os

setup(
    name='gcmtools',
    version='0.5.7',
    packages=['gcmtools',],
    zip_safe=False,
    install_requires=["numpy","netCDF4","matplotlib","scipy"],
    include_package_data=True,
    author='Adiv Paradise',
    author_email='paradise.astro@gmail.com',
    license='GNU General Public License',
    license_files="LICENSE.txt",
    url='https://github.com/alphaparrot/gcmtools',
    description='GCM Output Analysis Tools',
    long_description_content_type="text/markdown",
    long_description=open('README.md').read(),
    )

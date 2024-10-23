from setuptools import setup, find_packages

setup(
    name='JOnTADS',
    version='0.8',
    author='Qiuhai Zeng',
    author_email='qiuhai.stat@gmail.com',
    description='JOnTADS: a unified caller for TADs and stripes in Hi-C data',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    packages=find_packages(),
    python_requires='>=3.10.5',
    install_requires=[
        'numba==0.56.4', 
        'numpy==1.23.5',
        'scipy==1.12.0',
        'qpsolvers==2.7.3'
    ],
    entry_points={
        'console_scripts': [
            'JOnTADS=JOnTADS.JOnTADS:main',  # "script_name=package.module:function"
        ],
    },
    url='https://github.com/qunhualilab/JOnTADS',
)

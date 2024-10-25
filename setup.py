from setuptools import setup, find_packages


setup(
    name='pathsingle',
    version='0.1.4',
    description='PathSingle is a Python package for pathway analysis of single-cell RNA-seq data.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Dani Livne',
    author_email='zurkin@gmail.com',
    url='https://github.com/zurkin1/pathsingle',
    license='MIT',
    packages=find_packages(include=['pathsingle', 'pathsingle.*', 'pathway_db', 'pathway_db.*']),
    include_package_data=True,  # Include non-Python files specified in MANIFEST.in.
    install_requires=[
        'numpy', 
        'scipy',
        'pandas',
        'scikit-learn',
        'ipywidgets',
        'networkx',
        'plotly',
        'bs4',
        'scanpy',
        'scikit-misc',
        'tqdm',
        'omnipath',
        'umap-learn',
        'magic-impute',
        'torch',
        'torchvision',
        'torchaudio',
        'matplotlib',
        'scprep'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',  # Change as appropriate
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
    python_requires='>=3.6',
)
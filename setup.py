from setuptools import setup, find_packages

setup(
    name='microorganism-snps',
    version='0.1.0',
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=[
        'typer',
        'rich',
        'biopython',
    ],
    entry_points={
        'console_scripts': [
            'microorganism_snps=microorganism_snps.cli:run_app',
        ],
    },
    include_package_data=True,
    zip_safe=False,
)

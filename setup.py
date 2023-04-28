from setuptools import setup, find_packages

with open("requirements.txt", "r") as f:
    requirements = f.read().strip().split("\n")

setup(
    name = 'SnapFISH',
    version = '0.1.0',
    author = 'Lindsay Lee, Hongyu Yu',
    author_email = 'hongyuyu@unc.edu',
    license = 'GNU',
    description = 'A computational pipeline to identify chromatin loops from multiplexed DNA FISH data',
    long_description = open('README.md').read(),
    long_description_content_type = "text/markdown",
    url = 'https://github.com/lindsayhrlee/SnapFISH',
    py_modules = ['run_SnapFISH', 'src'],
    packages = find_packages(),
    install_requires = [requirements],
    python_requires = '>=3.6.8',
    entry_points = '''
        [console_scripts]
        SnapFISH=run_SnapFISH:main
    '''
)
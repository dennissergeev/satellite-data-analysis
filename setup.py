import os
from setuptools import setup

# Get version and release info, which is all stored in umtools/version.py
ver_file = os.path.join('sattools', 'version.py')
with open(ver_file) as f:
    exec(f.read())

# Get version and release info, which is all stored in umtools/version.py
req_file = 'requirements.txt'
REQUIRES = [line.rstrip('\n') for line in open(req_file)]

scripts = ['bin/satview']

opts = dict(name=NAME,
            maintainer=MAINTAINER,
            maintainer_email=MAINTAINER_EMAIL,
            description=DESCRIPTION,
            long_description=LONG_DESCRIPTION,
            url=URL,
            download_url=DOWNLOAD_URL,
            license=LICENSE,
            classifiers=CLASSIFIERS,
            author=AUTHOR,
            author_email=AUTHOR_EMAIL,
            version=VERSION,
            packages=PACKAGES,
            requires=REQUIRES,
            scripts=scripts)


if __name__ == '__main__':
    setup(**opts)

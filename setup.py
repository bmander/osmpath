import setuptools

setuptools.setup(
    name="osmpath",
    version="0.1.0",

    author="Brandon Martin-Anderson",
    author_email="badhill@gmail.com",

    description="Shortest paths using OpenStreetMap.",
    long_description=open('README.md').read(),

    packages=setuptools.find_packages(),

    install_requires=["osmium","numpy","scipy","pyproj"],

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
)

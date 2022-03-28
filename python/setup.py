import setuptools

setuptools.setup(
    name='sphtools',
    version='0.2.0',
    url='',
    classifiers=[
        "Private :: Do not Upload"
    ],
    python_requires='>=3.6',
    packages=setuptools.find_packages(),
    install_requires=['numpy', 'scipy', 'matplotlib', 'pyvistaqt', 'imageio-ffmpeg', 'PyQt5'],
    entry_points={
        "console_scripts": [
            "sphviz = sphtools.visualize:main",
            "sphprocess = sphtools.process_data:main",
        ]
    }
)

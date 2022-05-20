import setuptools

setuptools.setup(
    name="conet-py",
    version="1.0.0",
    author="Tomasz Cakala",
    author_email="tc360950@gmail",
    description="Utility wrapper for CONET",
    url="https://github.com/tc360950/CONET",
    classifiers=[
        "Programming Language :: Python :: 3"
    ],
    install_requires=[
      'numpy>=1.22.3',
      'pandas>=1.4.2',
      'networkx>=2.8'
    ],
    packages=["conet", "conet.generative_model", "conet.data_converter"],
    python_requires=">=3.7",
)

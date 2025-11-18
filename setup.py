from setuptools import find_packages, setup

about: dict = {}
exec(open("polyat/__init__.py", encoding="utf-8").read(), about)

setup(
    name="polyat",
    version=about["__version__"],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[],
    entry_points={
        "console_scripts": [
            "polyat=polyat.polyat:main",
        ],
    },
    description="Tool for summarizing poly-A/T homopolymers in sequencing reads",
    author="Daan Jansen",
    author_email="jansendaan94@gmail.com",
    url="https://github.com/DaanJansen94/polyat",
    python_requires=">=3.8",
)


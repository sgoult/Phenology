from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


ext_modules = [
    Extension(
        name="fast_polarity.polarity",
        sources=[
            "fast_polarity/polarity.pyx"
        ],
        libraries=[
            'm'
        ],
        extra_compile_args=[
            '-ffast-math'
        ]
    )
]

setup(
    name="fast_polarity/polarity",
    cmdclass = {"build_ext": build_ext},
    ext_modules=ext_modules
)

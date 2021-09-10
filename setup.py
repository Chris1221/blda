from distutils.core import setup

setup(name='bulk_lda',
      version='0.2',
      description = "LDA for Bulk ATAC",
      author='Chris Cole',
      author_email='chris.c.1221@gmail.com',
      packages=['bulk_lda'],
      install_requires = [
          'numpy'
      ]
      )
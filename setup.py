from distutils.core import setup

setup(
  name = 'apobec',
  packages = ['apobec'],   
  version = '0.3',      
  license='MIT',        
  description = 'script for SNP counting in deep sequencing data',   
  author = 'Yuriy Babin',                  
  author_email = 'babin.yurii@gmail.com',      
  url = 'https://github.com/babinyurii/apobec', 
  download_url = 'https://github.com/babinyurii/apobec/archive/refs/tags/v_0.3.tar.gz',
  keywords = ['bioinformatics', 'SNP'],   
  install_requires=[            
          'pandas',
          'biopython',
          'matplotlib',
          'numpy',
          'seaborn'
      ],
  classifiers=[
    'Development Status :: 4 - Beta',      
    'Intended Audience :: Developers',      
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   
    'Programming Language :: Python :: 3',      
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
  ],
)

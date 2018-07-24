# my_bigdft

# Ma branche de BigDFT pour mes changements au code de base.

Notes d'installation:

1. mkdir ../build
2. cd ../build
3. S'assurer d'être dans le bon environnement
    - Python 2.7 (alias lpy2 sur briarée et sur le laptop)
    - compilateur openmpi version 1.x , (marche pas sur 3.x ?) (alias bigenv sur briarée, par défaut sur le laptop)
4. cp ../bigdft-master/bigdft/rcfiles/briaree.rc buildrc (marche sur les deux)
5. ../my_bigdft/bigdft-master/bigdft/Installer.py -f buildrc -v build -y

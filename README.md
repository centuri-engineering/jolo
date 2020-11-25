# jolo
Jupyter / voilà interface to OLOGRAM


## Install

Install [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)

```sh
sudo apt-get install bedtools
```



## Dependencies




## OLOGRAM documentation

https://dputhier.github.io/pygtftk/ologram.html



## Tentative spec

- Interface web
  - Trois onglets:
    - OLOGRAM
    - OLOGRAM-MODL
    - MERGE_STATS
      - Input (+)
      - Ajouterai un nom et un tsv
- Proposer des jeux exemple.
  - Jeux artificiel
  - Données réels

- Aide contextuelle sur chaque argument (cf CLI)
  - Zone chargement de fichiers
  - Upload de plusieurs fichiers
  - paramètres

- Comment fait on avec les GTF de référence

- Output
    - Un fichier pdf + un tsv

- On devrait pouvoir construire un format html en sortie
    - Graph interactif
        - bqplot
        - altair

- On voit la commande envoyée par l’utilisateur

- Serveur avec un scheduler (load balancing ?):
- Slurm
    - Définir la mémoire
- Qsub

- Comment l’utilisateur suit son travail
    - s’identifier / on simplement donner son adresse.
- Une URL qu’il pourra consulter avec son résultat.
- Les résultats ont une certaine pérennité
- Migrable facilement


## Ressource management

* https://docs.hpc.cofc.edu/using-the-hpc/scheduling-jobs/jupyterhub
* https://docs.python.org/3/library/resource.html

# Alignement d'embedding par programmation dynamique

Conception d’un programme d’alignement d’embedding par programmation dynamique

Réaliser un programme permettant d’aligner des séquences
protéiques représentées par un embedding. L’article [1] servira de support pour
comprendre l’intérêt des embeddings (notamment les figure 3 et 4C) pour
l’alignement. Les embeddings sont disponibles à l’adresse suivante :
www.dsimb.inserm.fr/~gelly/data/Emb.zip. 
Ils ont été obtenus par la méthode T5 ProtTrans (github.com/agemagician/ProtTrans) et ont une longueur de 1024
variables. L’archive contient des fichiers de protéines dont les séquences ont
été encodé en embedding de longueur 1024. Dans ces fichiers, chaque ligne de
1024 valeurs représente une position de la séquence protéique. Les séquences
de ces protéines peuvent être obtenues à l’adresse suivante :
www.dsimb.inserm.fr/~gelly/data/Sequences.zip . Pour avoir une idée de la
ressemblance entre toutes ces protéines, vous pouvez consulter le fichier
www.dsimb.inserm.fr/~gelly/data/TMSCORES_HOMSTRAD.txt. La première
colonne est l’identifiant de la protéine 1 la deuxième colonne celle de la
protéine 2 et les colonnes 7 et 8 sont les TMscores mesuré sur les structures de
la protéine 1 contre la protéine 2 et l’inverse. Un TMscore supérieur à 0.5 est
considéré comme une ressemblance significative entre protéines. Ces
TMscores pourront dont être utilisé comme référence pour mesurer la capacité
qu’à l’algorithme pour retrouver une ressemblance structurale à partir des
embeddings des séquences.

Le programme devra reprendre les algorithmes Needleman-Wunsch, Smith-
Waterman et semi-global (Glocal ou Gloloc) (c’est-à-dire global sur une des
séquences et global sur l’autre) et aligner les séquences représentées par des
embeddings. Le score de match sera calculé par le dot product (ou la corrélation)
entre les vecteurs d’embedding d’une position de la protéine 1 et une autre
position de la protéine 2. Le programme gérera les pénalités fixe et affine des
brêches (gap). En sortie il donnera l’alignement représenté sous forme de
séquence. Il permettra également de gérer des embedding de taille arbitraire
(donc pas uniquement 1024).

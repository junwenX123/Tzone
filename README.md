# Enrichissement dans la zone en T au sein d'un modèle stochastique de mort cellulaire

Ce dépôt rassemble les expériences de simulation et les analyses statistiques dédiées à l'étude de l'enrichissement spatial des morts cellulaires observées au sein d'une **zone spatiale en forme de T** prédéfinie.

Le modèle intègre :

- des événements d'activation stochastiques ;
- des régions actives spatiales générées par ces activations ;
- des événements de mort cellulaire observés ;
- une protection locale médiée par ERK consécutive à chaque mort cellulaire ;
- une région fixe en forme de T présentant une intensité d'activation basale accrue.

Les simulations reposent sur un **algorithme de Gillespie événementiel** couplé à une méthode d'amincissement (thinning).

L'objectif principal consiste à identifier les régimes paramétriques sous lesquels les morts observées se concentrent de manière significative dans la zone en T, et plus particulièrement à dissocier les effets :

1. du contraste d'activation entre la zone en T et le reste du domaine ;
2. de la localisation spatiale des régions actives ;
3. de l'intensité de mort ;
4. du rayon d'action et de la durée de protection ERK.

---

## 1. Domaine spatial et zone en T

La fenêtre d'observation spatiale est définie par

\[W=[0,L_x]\times[0,L_y],\]

avec pour valeurs par défaut

\[L_x=L_y=10.\]

Le domaine est découpé en tiers selon les deux axes spatiaux.

La zone fixe en T est constituée de :

- la colonne verticale centrale ;
- la section gauche de la rangée horizontale centrale.

Son aire est donc

\[\vert{}T\vert{} = \frac{1}{3}\vert{}W\vert{}+\frac{1}{9}\vert{}W\vert{} = \frac{4}{9}\vert{}W\vert{}.\]

Par conséquent, sous l'hypothèse nulle d'un processus spatialement uniforme, la probabilité qu'une mort survienne dans la zone en T vaut

\[p_0=\frac{\vert{}T\vert{}}{\vert{}W\vert{}}=\frac49.\]

---

## 2. Modèle stochastique

### Activation

Les candidats à l'activation sont générés selon un processus de Poisson dominant d'intensité

\[\lambda_{a,1}\vert{}W\vert{}.\]

Un candidat à l'activation en position \(x\) est accepté selon l'intensité spatiale d'activation

\[\lambda_a(x)= \begin{cases} \lambda_{a,1}, & x \text{ se trouve déjà dans une région active},\\ \lambda_{a,T}, & x\in T \text{ mais hors des régions actives},\\ \lambda_{a,c}, & x\notin T \text{ et hors des régions actives}. \end{cases}\]

La simulation applique ainsi un rejet par amincissement (thinning) basé sur l'intensité dominante \(\lambda_{a,1}\).

Chaque activation acceptée engendre un disque actif dont le rayon aléatoire suit une loi exponentielle

\[R^a\sim \mathrm{Exp}(\beta_{a,R}),\]

et les centres actifs disparaissent à un taux

\[\beta_{a,T}.\]

Dès lors,

\[\mathbb E[R^a]=\frac{1}{\beta_{a,R}}, \qquad \mathbb E[T^a]=\frac{1}{\beta_{a,T}}.\]

---

### Mort cellulaire

Les candidats à la mort cellulaire sont générés uniformément sur \(W\) à un taux total

\[\lambda_d \vert{}W\vert{}.\]

Un événement de mort candidat n'est validé que si sa position vérifie :

- l'appartenance à l'union des régions actives, et
- la non-appartenance aux régions protégées par ERK à cet instant.

Ainsi, l'intensité spatiale effective de mort s'écrit

\[r(t,x) = \lambda_d \mathbf 1_{\{ x\in A_t^{a}\setminus A_t^{\mathrm{ERK}} \}}.\]

---

### Rétroaction négative par ERK

Chaque mort cellulaire observée acceptée génère un disque de protection locale ERK.

Son rayon suit

\[R^d\sim \mathrm{Exp}(\beta_{d,R}),\]

tandis que les zones de protection ERK se dissipent à un taux

\[\beta_{d,T}.\]

Par conséquent,

\[\mathbb E[R^d]=\frac{1}{\beta_{d,R}}, \qquad \mathbb E[T^d]=\frac{1}{\beta_{d,T}}.\]

Un disque ERK plus large ou une durée de vie accrue d'ERK inhibe davantage les morts ultérieures à proximité.

---

## 3. Simulation par algorithme de Gillespie

À l'instant \(t\), quatre classes d'événements entrent en compétition :

\[\begin{aligned} &\text{candidat activation :} &&\lambda_{a,1}\vert{}W\vert{},\\ &\text{candidat mort :} &&\lambda_d\vert{}W\vert{},\\ &\text{expiration d'un centre actif :} &&N_t^a\beta_{a,T},\\ &\text{expiration d'une zone ERK :} &&N_t^p\beta_{d,T}. \end{aligned}\]

Le taux total de Gillespie s'élève ainsi à

\[a_0(t) = \lambda_{a,1}\vert{}W\vert{} + \lambda_d\vert{}W\vert{} + N_t^a\beta_{a,T} + N_t^p\beta_{d,T}.\]

Le temps d'attente jusqu'au prochain événement est tiré selon

\[\Delta t\sim\mathrm{Exp}(a_0(t)),\]

et la nature de l'événement est sélectionnée proportionnellement aux quatre taux ci-dessus.

---

## 4. Phase de chauffe (burn-in) et protocole expérimental

Les simulations principales adoptent le protocole suivant :

- **500 morts observées pour la phase de chauffe (burn-in)** ;
- les **1000 morts observées suivantes pour l'analyse statistique**.

Chaque réplicat exploitable totalise donc

\[500+1000=1500\]

morts observées acceptées, les 500 premières étant systématiquement écartées des statistiques spatiales finales.

Cette phase transitoire atténue la dépendance vis-à-vis de l'état initial vide arbitraire.

Des graines aléatoires indépendantes sont attribuées aux réplicats.

---

## 5. Test statistique de l'enrichissement dans la zone en T

Posons

\[ K = \#\{\text{morts analysées situées dans la zone en T}\}, \]

et désignons par \(N\) le nombre total de morts analysées.

Sous l'hypothèse nulle d'homogénéité spatiale,

\[K\sim\mathrm{Binomiale}(N,p_0), \qquad p_0=\frac{\vert{}T\vert{}}{\vert{}W\vert{}}.\]

Le dépôt met en œuvre le test binomial exact unilatéral

\[H_0:\pi_T=p_0\]

contre

\[H_1:\pi_T>p_0.\]

La p-valeur s'exprime par

\[P_{H_0}(X\ge K), \qquad X\sim\mathrm{Binomiale}(N,p_0).\]

Le seuil de significativité nominal retenu pour les analyses actuelles est

\[\alpha=0.05.\]

Aucune correction de Bonferroni n'est appliquée sur les analyses binomiales agrégées finales.

---

## 6. Ratio de densité dans la zone en T

En complément du test binomial, les simulations évaluent le ratio de densité spatiale

\[ R_T = \frac{ K/\vert{}T\vert{} }{ (N-K)/(\vert{}W\vert{}-\vert{}T\vert{}) }. \]

Interprétation :

- \(R_T>1\) : concentration plus élevée des événements de mort dans la zone en T ;
- \(R_T=1\) : densité spatiale uniforme entre l'intérieur et l'extérieur de la zone en T ;
- \(R_T<1\) : sous-représentation des morts dans la zone en T.

---

# 7. Déroulement des analyses

Le dépôt s'articule autour de trois étapes successives de simulation et de diagnostic.

---

## Étape 1 — Contraste d'activation et localisation spatiale

`step1code.py`

Cette étape correspond à l'exploration paramétrique globale préliminaire.

Elle examine la sensibilité du signal dans la zone en T vis-à-vis des paramètres d'activation :

\[\lambda_{a,T}, \qquad \lambda_{a,c}, \qquad \beta_{a,R}, \qquad \beta_{a,T}.\]

Ces simulations étudient plus particulièrement :

- différents contrastes d'activation \(\lambda_{a,T}/\lambda_{a,c}\) entre la zone en T et le fond ;
- différentes intensités d'activation absolues ;
- des degrés variables de localisation des régions actives ;
- l'impact du rayon et de la durée de vie des zones actives.

Chaque jeu de paramètres est simulé sur plusieurs réplicats indépendants.

### Fichiers de l'Étape 1

| Fichier | Description |
|---|---|
| `step1code.py` | Simulation de Gillespie et balayage des paramètres d'activation et de zone en T |
| `step1result.csv` | Résultats de simulation bruts par réplicat pour l'Étape 1 |
| `step1resultsummary.csv` | Statistiques descriptives agrégées par configuration de paramètres |
| `step1pooled_binomial_analysis.py` | Agrégation des réplicats indépendants et calcul du test binomial exact unilatéral |
| `step1pooled_binomial_summary.csv` | Résultats de l'enrichissement dans la zone en T après agrégation |

---

## Étape 2 — Balayage factoriel des paramètres de mort et d'ERK

`step2code.py`

L'Étape 2 fixe un régime d'activation où l'enrichissement spatial dans la zone en T est déjà clairement caractérisé :

\[\lambda_{a,T}=0.5, \qquad \lambda_{a,c}=0.005, \qquad \beta_{a,R}=2.5, \qquad \beta_{a,T}=3.0.\]

Les paramètres régissant la mort et la boucle ERK varient ensuite de façon systématique.

Le plan factoriel complet explore

\[\lambda_d\in\{0.5,1,2\},\]

\[\beta_{d,R}\in\{1,2,4\},\]

et

\[\beta_{d,T}\in\{0.4,0.8,1.6\}.\]

Soit

\[3\times3\times3=27\]

combinaisons paramétriques de mort et d'inhibition ERK.

Chaque configuration comprend cinq réplicats indépendants.

L'enjeu consiste à tester si l'enrichissement dans la zone en T identifié à l'Étape 1 demeure robuste face aux modulations de la dynamique de mort et de la rétroaction ERK.

### Fichiers de l'Étape 2

| Fichier | Description |
|---|---|
| `step2code.py` | Balayage factoriel complet des paramètres de mort et de rétroaction ERK |
| `step2result.csv` | Résultats bruts des simulations par réplicat |
| `step2resultsummary.csv` | Statistiques synthétiques agrégées par configuration de mort/ERK |
| `step2binomialtestcode.py` | Test binomial unilatéral exact sur données agrégées |
| `step2pooled_binomial_summary.csv` | Résultats finaux du test binomial agrégé |

---

## Étape 3 — Analyse diagnostique de la localisation d'activation

`step3code.py`

L'Étape 3 isole et teste directement le rôle de la **localisation de l'activation**.

Trois régimes d'activation y sont confrontés :

| Régime | \(\beta_{a,R}\) | \(\beta_{a,T}\) | Interprétation |
|---|---:|---:|---|
| `non_local_activation` | 1.0 | 1.2 | disques actifs étendus et de longue durée |
| `medium_activation` | 2.5 | 3.0 | localisation intermédiaire |
| `local_activation` | 5.0 | 3.0 | disques actifs restreints et éphémères |

Les paramètres de bruit de fond d'activation sont fixés à

\[\lambda_{a,T}=0.5, \qquad \lambda_{a,c}=0.005.\]

Pour chacun des trois régimes d'activation, l'expérience fait également varier

\[\lambda_d \in \{0.25,0.5,1,2,4\},\]

croisé avec trois niveaux d'intensité pour la protection ERK :

- protection ERK forte ;
- protection ERK de référence ;
- protection ERK faible.

Chaque régime d'activation fait ainsi l'objet d'un examen sous

\[5\times3=15\]

configurations de mort et d'ERK.

Chaque configuration s'appuie sur **10 réplicats de simulation indépendants**, correspondant à 10 000 morts observées analysées par condition après agrégation.

### Fichiers de l'Étape 3

| Fichier | Description |
|---|---|
| `step3code.py` | Expérience diagnostique croisant localisation d'activation et régimes mort/ERK |
| `step3result.csv` | Résultats des tests binomiaux agrégés pour chaque jeu de paramètres |
| `step3resultsummary.csv` | Synthèse de robustesse pour les trois régimes d'activation |

---

# 8. Principal résultat diagnostique

Les observations de l'Étape 3 démontrent sans ambiguïté que la **localisation spatiale de l'activation constitue le mécanisme prépondérant dans l'émergence d'un enrichissement net au sein de la zone en T**.

| Régime d'activation | Configurations significatives | Fraction moyenne zone en T | Ratio de densité moyen |
|---|---:|---:|---:|
| activation locale | 15 / 15 | 0.907 | 12.21 |
| activation intermédiaire | 15 / 15 | 0.622 | 2.07 |
| activation non locale | 0 / 15 | 0.424 | 0.92 |

Pour rappel, la proportion géométrique de la zone en T sous l'hypothèse nulle d'homogénéité spatiale s'élève à

\[p_0=\frac49\approx0.444.\]

Le régime non local ne présente aucun enrichissement significatif, alors que les régimes intermédiaire et fortement localisé maintiennent une surreprésentation robuste à travers l'ensemble des conditions de mort et d'ERK testées.

Ces mesures confirment que la morphologie spatiale du patron en T est gouvernée au premier ordre par le **degré de confinement spatial du processus d'activation**, l'intensité de mort et la rétroaction ERK modulant principalement l'échelle temporelle et l'amplitude quantitative du signal observé.

---

# 9. Paramètres du modèle

| Paramètre | Signification |
|---|---|
| `lambda_a_1` | intensité d'activation dominante / au sein des régions actives |
| `lambda_a_T` | intensité d'activation basale dans la zone en T |
| `lambda_a_c` | intensité d'activation basale hors de la zone en T |
| `lambda_d` | intensité des événements candidats de mort |
| `beta_a_R` | paramètre d'échelle (taux) pour le rayon des disques d'activation |
| `beta_a_T` | taux d'extinction des régions actives |
| `beta_d_R` | paramètre d'échelle (taux) pour le rayon de protection ERK |
| `beta_d_T` | taux d'extinction des zones de protection ERK |

Pour des variables aléatoires suivant une loi exponentielle,

\[\mathbb E[R^a]=\frac1{\beta_{a,R}}, \qquad \mathbb E[T^a]=\frac1{\beta_{a,T}},\]

et

\[\mathbb E[R^d]=\frac1{\beta_{d,R}}, \qquad \mathbb E[T^d]=\frac1{\beta_{d,T}}.\]

Un paramètre \(\beta_R\) plus **élevé** correspond ainsi à un **rayon moyen plus petit**, tandis qu'un \(\beta_T\) plus grand équivaut à un **temps de persistance moyen plus court**.

---

# 10. Prérequis logiciels

Le code nécessite Python 3 ainsi que les bibliothèques suivantes :

```bash
pip install numpy pandas scipy matplotlib

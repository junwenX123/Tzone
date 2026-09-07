### Conclusion de l'Étape 1

- **Le régime de référence ne produit aucun enrichissement.** Sur données agrégées ($N=5000$) : $K/N = 2187/5000$, soit $\hat\pi_T = 0.4374$, valeur *inférieure* à $p_0 = 0.4444$, avec $p_{\text{agrégée}} = 0.845$. Les contrôles simples (mort plus rapide, rayon ERK plus grand, durée ERK plus courte, activation de $T$ plus forte) restent tous non significatifs au seuil de 5 %.
- **Le contraste d'intensité seul ne suffit pas.** Porter $\lambda_{a,T}/\lambda_{a,c}$ à $10^2$ ($\hat\pi_T = 0.4420$, $p = 0.641$), voire à $10^3$ ou $10^4$, laisse la fraction au voisinage de $p_0$ ou en dessous : si les disques actifs sont trop grands ou trop persistants, le signal local de $T$ est spatialement dilué.
- **La localisation change tout.** À contraste identique, les jeux à activation locale atteignent $\hat\pi_T \approx 0.905$–$0.932$ et $R_T \approx 11.9$–$17.2$. Le balayage conjoint $(\beta_{a,R},\beta_{a,T})$ fait passer $R_T$ de $0.99$ à $23.0$ : l'effet dominant vient de l'augmentation de $\beta_{a,R}$ (rayon moyen plus petit), $\beta_{a,T}$ affinant le signal surtout lorsque le rayon n'est pas encore très petit.

> **À retenir.** L'enrichissement de la zone en T n'est pas causé par une intensité d'activation plus élevée dans $T$ ; il apparaît lorsque le champ d'activation devient suffisamment **local en espace et en temps**.

---

### Conclusion de l'Étape 2

- **Robustesse complète.** Les $27$ configurations mort/ERK sont significatives. La fraction agrégée s'étend de $0.5742$ à $0.6642$, toujours nettement au-dessus de $p_0 = 0.4444$, avec $R_T \in [1.686,\ 2.472]$ : même dans le cas le plus défavorable, la densité de morts dans $T$ vaut $\approx 1.69$ fois la densité extérieure.
- **Extrêmes du plan factoriel.** Le plus fort pour $(\lambda_d,\beta_{d,R},\beta_{d,T}) = (0.5,\ 4.0,\ 1.6)$ : $k_T = 3321/5000$, $\hat\pi_T = 0.6642$, $R_T = 2.472$, $p \approx 9.1\times10^{-215}$. Le plus faible pour $(2.0,\ 1.0,\ 0.4)$ : $k_T = 2871/5000$, $\hat\pi_T = 0.5742$, $R_T = 1.686$, $p \approx 1.0\times10^{-75}$.
- **L'effet dominant est temporel, non spatial.** Augmenter $\lambda_d$ réduit le temps moyen nécessaire pour accumuler 1000 morts analysées ($114.5 \to 75.2 \to 53.3$ pour $\lambda_d = 0.5,\,1,\,2$) sans modifier substantiellement la fraction dans la zone en T ($0.637 \to 0.623 \to 0.608$) ni le ratio de densité ($2.22 \to 2.08 \to 1.96$). Le rayon et la durée de la protection ERK jouent de même sur la force et la vitesse du processus, pas sur l'existence de l'enrichissement.

> **À retenir.** Une fois le régime d'activation fixé dans une configuration modérément localisée, les paramètres de mort et d'ERK **ne suppriment pas** l'enrichissement de la zone en T : ils contrôlent principalement l'échelle de temps du processus de mort observé et n'en modulent qu'accessoirement l'amplitude.

---

### Conclusion de l'Étape 3

À contraste d'activation fixé ($\lambda_{a,1}=5.0$, $\lambda_{a,T}=0.5$, $\lambda_{a,c}=0.005$, soit $\lambda_{a,T}/\lambda_{a,c}=100$) et sur les **mêmes** 15 configurations mort/ERK, 10 réplicats chacune ($\approx 10\,000$ morts analysées par condition) :

| Régime d'activation | $(\beta_{a,R},\beta_{a,T})$ | Configurations significatives | Fraction zone T | $R_T$ moyen |
|---|---:|---:|---:|---:|
| non local | $(1.0\ ;\ 1.2)$ | 0 / 15 | $0.3865$–$0.4474$ (moy. $0.4236$) | $0.921$ |
| intermédiaire | $(2.5\ ;\ 3.0)$ | 15 / 15 | $0.5686$–$0.6574$ (moy. $0.6220$) | $2.071$ |
| local | $(5.0\ ;\ 3.0)$ | 15 / 15 | $0.8950$–$0.9136$ (moy. $0.9068$) | $12.211$ |

- **Le régime non local est indiscernable du modèle nul homogène** ($R_T \le 1.012$, fractions au niveau de $p_0$ ou en dessous), et ce pour *toutes* les conditions de mort et d'ERK testées : l'enrichissement disparaît purement et simplement.
- Le passage non local $\to$ intermédiaire $\to$ local fait croître $R_T$ de $0.92$ à $2.07$ puis à $12.2$ — soit plus d'un ordre de grandeur — alors que **seule la géométrie des zones actives a changé**, à contraste et à paramètres mort/ERK identiques.
- L'interprétation est directe : puisque $\mathbb E[R^a]=1/\beta_{a,R}$ et $\mathbb E[T^a]=1/\beta_{a,T}$, des valeurs plus grandes correspondent à des zones d'activation **plus petites et plus brèves**, qui préservent la structure spatiale du paysage d'activation fixe en forme de T au lieu de la diffuser.

> **À retenir.** **L'enrichissement de la zone en T est principalement contrôlé par la localisation de l'activation.** Le contraste d'intensité est nécessaire mais non suffisant (Étape 1) ; la dynamique de mort et la rétroaction ERK ne règlent que la vitesse et l'amplitude du signal (Étape 2).

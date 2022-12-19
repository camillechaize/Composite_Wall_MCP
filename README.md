> - [Méthode Physique de résolution](#intro)
> - [Fonctionnement du programme Python](#prog)

<h4 id="intro"></h4>

# Introduction au problème

On étudie l’influence des matériaux à changement de phase (**MCP**) sur l’isolation des batiments.

> ![original image](https://cdn.mathpix.com/snip/images/86JViDn_5w5EE3N9ugJBcKIU_BMChFuT9xcGdamfG_c.original.fullsize.png)
> Illustration du problème

## Hypothèses
- Les transferts thermiques sont supposés __unidimensionnels__
- La __température extérieure__ est variable dans le temps, elle est __connue__
- La __température de la pièce__ ne varie pas au cours du temps, elle est __connue__
- Les __échanges aux interfaces air/solide__ sont gérés par la __loi de Newton__
- Le mur est __composite__ et peut comporter __plusieurs MCP__
- Les __propriétés thermodynamiques du MCP__ sont supposées __constantes__ dans sa phase solide ou dans sa phase liquide (la __conductivité__, la __capacité calorifique volumique__, la __masse volumique__ ne __dépendent que__ de l’__état solide/liquide__ du MCP et non de sa température, pression, etc…)

 ## Grandeurs utilisées

| Symbole | Signification | Unité |
|-----|-----|-----|
| $H$ | Enthalpie volumique | $J.m^{-3}$ |
| $h$ | Enthalpie sensible volumique | $J.m^{-3}$ |
| $T$ | Température | $K$ |
| $\lambda$ | Conductivité thermique| $W.m^{-1}.K^{-1}$ |
| $\rho$ | Masse volumique | $kg.m^{-3}$ |
| $c$ | Capacité calorifique massique | $J.kg^{-1}$ |
| $L$ | Chaleur latente massique | $J.kg^{-1}$ |
| $f$ | fraction fondue, liquide | SU |
| $h_{int}$, $h_{ext}$ | Coefficients de convection | $W.m^{-2}.K^{-1}$ |


___

# Modélisation du problème

> On scinde le transfert de chaleur en trois parties:
> - à  __l’intérieur__ du MCP
> - à __l’interface Solide/MCP__
> - à __l’interface Air/MCP__

___

## Intérieur du MCP

Un bilan d’énergie sur un volume de contrôle du mur permet d’écrire:

<h4 id="eqn-1"></h4>

$$
\frac{\partial H}{\partial t}=\vec{\nabla}\left(\lambda_{k} \vec{\nabla} T\right)
\tag{1}
$$

En décomposant l’enthalpie volumique totale $H$ en enthalpie sensible volumique et chaleur latente volumique:

<h4 id="eqn-2"></h4>

$$
H(T)=h(T)+\rho_{L} L \cdot f
\tag{2}
$$

> - $f$ est la fraction liquide (fondue)

L’enthalpie sensible volumique a pour expression:

<h4 id="eqn-3"></h4>

$$
h(T)=\int_{T_{f}}^{T} \rho_{k} \cdot c_{k} \cdot d T
\tag{3}
$$
> - $k$ correspond à la phase du MCP
> - $T_f$ correspond à la température de fusion (et solidification) du MCP
> - _l’enthalpie étant définie à une constante près, ce choix de borne inférieur est arbitraire mais permet de simplifer les équations_

En combinant [(1)](#eqn-1), [(2)](#eqn-2) et [(3)](#eqn-3) on obtient:

$$
\frac{\partial h}{\partial t}=\frac{\partial}{\partial x}\left(\alpha \cdot \frac{\partial h}{\partial x}\right)-\rho_{L} \cdot L \cdot \frac{\partial f}{\partial t}
\tag{4}
$$

> - $\alpha = \frac{\lambda}{\rho\cdot c}$

Que l’on discrétise selon:

<h4 id="eqn-5"></h4>

$$
h_{i}^{t+1}=h_{i}^{t}+\alpha R \cdot\left(h_{i+1}^{t+1}-2 h_{i}^{t+1}+h_{i-1}^{t+1}\right)+p_{L} L \cdot\left(f^{t}-f^{t+1}\right)
\tag{5}
$$

Que l’on écrit (_pour résoudre matriciellement_) sous la forme:

<h4 id="eqn-6"></h4>

$$
a_{i-1} h_{i-1}^{t+1}+a_{i} h_{i}^{t+1}+a_{i+1} h_{i+1}^{t+1}=Q
\tag{6}
$$

> - $a_{i-1}=a_{i+1}=-\alpha R$
> - $Q=h_{i}^{t}+\rho_{L} L \cdot\left(f_{i}^{t}-f_{i}^{t+1}\right)$

___


## Interface Solide/MCP
> - ![original image](https://cdn.mathpix.com/snip/images/KrEI-tVByxkjCWQOzmEwIHanUpClseNnoRI_RKBKG3A.original.fullsize.png)
>  _Illustration de deux MCP en contact (les noeuds fantômes ne sont pas représentés)_

> - Les solides étant aussi des matériaux à changement de phase (simplements utilisés en _phase constante_) il est alors logique de s’intéresser à une interface généralle __MCP/MCP__

### Conditions d’interface:
En traduisant la __continuité du flux de chaleur__ et la __continuité de la température__ à __l’interface__, on obtient:

$$
\left.\left.\lambda_{1} \cdot \frac{\partial T_{1}}{\partial x}\right)_{x_{int}}=\lambda_{2} \cdot \frac{\partial T_{2}}{\partial x}\right)_{x_{int}}
$$

$$
T_{\text {int }}=T_{1}\left(x_{\text {int }}\right)=T_{2}\left(x_{\text {int }}\right)
$$

Que l’on discrétise en:

$$
\lambda_{1} \cdot \frac{\overline{T_{1}}(N+1)-T_{1}(N)}{d x}=\lambda_{2} \cdot \frac{T_{2}(1)-\overline{T_{2}}(0)}{d x}
$$

$$
T_{i n t}=\frac{T_{1}(N)+\overline{T}_{1}(N+1)}{2}=\frac{T_{2}(1)+\overline{T}_{2}(0)}{2}
$$

> - ![original image](https://cdn.mathpix.com/snip/images/9QTZGQGU5VY1lZlDuFz4JWIhsPQIyKTYAVG7Vf-XltI.original.fullsize.png)
> _Illustration du noeud fantôme du matériau de gauche situé  « à la place » du noeud 1 du matériau de droite_

> - $\overline{T_{1}}(N+1)$ et $\overline{T_{2}}(0)$ sont les noeuds fantômes (_seul_ $\overline{T_{1}}(N+1)$ _est représenté sur l’illustration_)
> - Les noeuds fantômes permettent de considérer en $N$ (resp. en $1$) que le $\text{MCP}_1$ (resp. $\text{MCP}_2$) est entouré de noeuds constitués uniquement de  $\text{MCP}_1$ (resp. $\text{MCP}_2$) et ainsi pouvoir discrétiser les conditions sur les flux de chaleur et les conditions sur les températures (ces conditions nécessitent des propriétés thermodynamiques uniformes).
> - Ces noeuds n’existant pas, on les élimine des équations par la suite.


On résout pour les noeuds fantômes:

$$
\overline{T_{2}}(0)=T_{1}(N) \cdot\left(\frac{2 \lambda_{1}}{\lambda_{1}+\lambda_{2}}\right)+T_{2}(1) \cdot\left(\frac{\lambda_{2}-\lambda_{1}}{\lambda_{1}+\lambda_{2}}\right)
$$

$$
\overline{T_{1}}(N+1)=T_{1}(N) \cdot\left(\frac{\lambda_{1}-\lambda_{2}}{\lambda_{1}+\lambda_{2}}\right)+T_{2}(1) \cdot\left(\frac{2 \lambda_{2}}{\lambda_{1}+\lambda_{2}}\right)
$$

Que l’on notera (_pour simplifier_):

$$
\overline{T_{2}}(0)=\left(1-\beta_{1 / 2}\right) \cdot T_{1}(N)+\beta_{1 / 2} \cdot T_{2}(1)
$$

$$
\overline{T_{1}}(N+1)=\beta_{2 / 1} \cdot T_{1}(N)+\left(1-\beta_{2 / 1}\right) \cdot T_{2}(1)
$$

> - $\beta_{j / \mu}=\frac{\lambda_{\mu}-\lambda_{j}}{\lambda_{\mu}+\lambda_{j}}$

On fait apparaître l’enthalpie sensible volumique des noeuds fantômes:


$$
\overline{h_{2,0}^{t+1}}=\left(1 - \beta_{1 / 2}\right) \cdot \frac{\rho_{2} \cdot c_{2}}{\rho_{1} \cdot c_{1}} \cdot h_{1, N}^{t+1}+\beta_{1 / 2} \cdot h_{2,1}^{t+1} + \left(1 - \beta_{1 / 2}\right) \cdot \left(T_{f_1} - T_{f_2}\right) \cdot \rho_2 c_2
$$

$$
\overline{h_{1, N+1}^{t+1}}=\beta_{2 / 1} \cdot h_{1, N}^{t+1}+\frac{\rho_{1} \cdot c_{1}}{\rho_{2} \cdot c_{2}} \cdot\left(1 - \beta_{2 / 1}\right) \cdot h_{2,1}^{t+1} + \left(1 - \beta_{2 / 1}\right) \cdot \left(T_{f_2} - T_{f_1}\right) \cdot \rho_1 c_1
$$

> - On s’est servi du fait que $h(T)=\rho c \cdot (T-T_{f_m})$ où $T_{f_m}$ correspond à la température de fusion du matériau (noté $m$) du noeud

Que l’on notera (_pour simplifier_):

$$
\overline{h_{2,0}^{t+1}}=\left(1 - \beta_{1 / 2}\right) \cdot p_{2/1} \cdot h_{1, N}^{t+1}+\beta_{1 / 2} \cdot h_{2,1}^{t+1} + \left(1 - \beta_{1 / 2}\right) \cdot \left(T_{f_1} - T_{f_2}\right) \cdot \rho_2 c_2
$$

$$
\overline{h_{1, N+1}^{t+1}}=\beta_{2 / 1} \cdot h_{1, N}^{t+1}+\left(1 - \beta_{2 / 1}\right) \cdot p_{1/2} \cdot h_{2,1}^{t+1} + \left(1 - \beta_{2 / 1}\right) \cdot \left(T_{f_2} - T_{f_1}\right) \cdot \rho_1 c_1
$$

> - $p_{j / \mu}=\frac{\rho_{j} c_{j}}{\rho_{\mu} c_{\mu}}$


#### Cas où un seul MCP est présent dans le mur

$$
\overline{h_{2,0}^{t+1}}=\left(1 - \beta_{1 / 2}\right) \cdot p_{2/1} \cdot h_{1, N}^{t+1}+\beta_{1 / 2} \cdot h_{2,1}^{t+1}
$$

$$
\overline{h_{1, N+1}^{t+1}}=\beta_{2 / 1} \cdot h_{1, N}^{t+1}+\left(1 - \beta_{2 / 1}\right) \cdot p_{1/2} \cdot h_{2,1}^{t+1}
$$

> - Ces formules ne sont valables que si le mur contient un seul MCP « réel » (_hors les solides constants comme le béton_) et où pour tout matériau $m$ présent dans le mur on a défini $h_m(T)$ comme $h_m(T)=\rho c \cdot (T-T_{f_{ref}})$
> - $T_{f_{ref}}$ correspondant à la température de fusion du seul MCP « réel » du mur.


_On notera qu’en créeant une séparation artificielle sur un même matériau (i.e. en considérant les paramètres thermodynamiques 1 et 2 égaux dans les équations précédentes), on a cohérence des enthalpies fantômes:_

$$
\left\{\begin{array}{l}\overline{h_{2,0}^{t+1}}=h_{1, N}^{t+1} \\ \overline{h_{1, N+1}^{t+1}}=h_{2,1}^{t+1}\end{array}\right.
$$

On reprend l’équation [(5)](#eqn-5) traduisant les transferts, et l’on réinjecte les noeuds fantômes:

#### Pour le noeud de gauche (matériau 1)

$$
-\gamma_{1} \cdot h_{i-1}^{t+1}+\left[1+\gamma_{1}\left(2-\beta_{2 / 1}\right)\right] \cdot h_{i}^{t+1}-\left[\gamma_{1} \cdot p_{1 / 2} \cdot\left(1-\beta_{2 / 1}\right)\right] \cdot h_{i+1}^{t+1}=h_{i}^{t}+\gamma_{1}\cdot\left(1 - \beta_{2 / 1}\right) \cdot \left(\Delta T_{f}\right)_{1}^{2} \cdot p_{1} c_{1}+\eta \cdot \left(f_{i}^{t} - f_{i}^{t+1}\right)
$$

> - $\gamma_{j}=\alpha_{j} R$

Soit en utilisant la forme de [(6)](#eqn-6):
> -  $a_{i-1}=-\gamma_{1}$
> - $a_{i}=1+\gamma_{1}\left(2-\beta_{2 / 1}\right)$
> - $a_{i+1}=-\gamma_{1} \cdot p_{1 / 2} \cdot\left(1-\beta_{2 / 1}\right)$
> - $Q = h_{i}^{t}+\gamma_{1}\cdot\left(1 - \beta_{2 / 1}\right)\cdot \left(\Delta T_{f}\right)_ {1}^{2} \cdot p_{1} c_{1}+\eta \cdot\left(f_i^{t}-f_i^{t+1}\right)$

#### Pour le noeud de droite (matériau 2)

$$
-\left[\gamma_{2} \cdot p_{2 / 1} \cdot\left(1-\beta_{1 / 2}\right)\right]\cdot h_{i-1}^{t+1}+\left[1+\gamma_{2}\left(2-\beta_{1 / 2}\right)\right] \cdot h_{i}^{t+1} -\gamma_{2} \cdot h_{i+1}^{t+1}=h_{i}^{t}+\gamma_{2}\cdot\left(1 - \beta_{1 / 2}\right)\cdot \left(\Delta T_{f}\right)_{2}^{1} \cdot p_{2} c_{2}+\eta \cdot \left(f_{i}^{t} - f_{i}^{t+1}\right)
$$

Soit en utilisant la forme de [(6)](#eqn-6):
> -  $a_{i-1}=-\gamma_{2} \cdot p_{2 / 1} \cdot\left(1-\beta_{1 / 2}\right)$
> - $a_{i}=1+\gamma_{2}\left(2-\beta_{1 / 2}\right)$
> - $a_{i+1}=-\gamma_{2}$
> - $Q = h_{i}^{t}+\gamma_{2}\cdot\left(1 - \beta_{1 / 2}\right)\cdot \left(\Delta T_{f}\right)_ {2}^{1} \cdot p_{2} c_{2}+\eta \cdot \left(f_{i}^{t} - f_{i}^{t+1}\right)$


## Interface Air/MCP
> - Il y a ici deux interfaces différentes Air/MCP:
>   1. Extérieur/Mur: _La température exterieure est fixée_
>   2. Mur/Intérieur: _La température de la pièce est inconnue_
> - Ces deux interfaces suivront chacune une loi de Newton de paramètres respectifs $h_{ext}$ et $h_{int}$

### Interface Extérieur/Mur
> - ![original image](https://cdn.mathpix.com/snip/images/kfsEz-kEN3VTgpCJYcKmEZ_MWvoDZJDhSt-uvCDAlWo.original.fullsize.png)
> _Illustration du noeud fantôme (hachuré)_

En exprimant la continuité du flux de chaleur à l’interface ( $x = x_{int}$ ) :

$$
\overrightarrow{j_{th}(x_{int})}=-\lambda \cdot \left.\frac{\partial T(x)}{\partial x} \right|_{x =x_{int}} \cdot \vec{x} =h_{\text {ext}} \cdot\left(T_{\text {ext}} - T_{\text {interface}}\right)\cdot \vec{x}
$$

> - $T(x)$ correspond exclusivement à la température dans le MCP

et en se servant de la méthode des noeuds fantômes (_en considérant une cellule imaginaire de MCP à gauche de l’interface_), tout en discrétisant, on obtient:

<h4 id="eqn-int-1"></h4>

$$
-\lambda \cdot \frac{T(1)-\overline{T}(0)}{d x}=h_{e x t} \cdot\left(T_{ext} -\frac{T(1)+\overline{T}(0)}{2}\right)
\tag{Int.1}
$$

> - où $T_{interface}$ est la température à la surface du mur, soit la température entre le $0^{eme}$ (_fantôme_) et $1^{er}$ noeud: on considère ici la moyenne des deux noeuds comme étant la température de l’interface.

En isolant $\overline{T}(0)$ :

$$
\overline{T}(0)=\frac{2 \lambda-h_{\text {ext }} d x}{2 \lambda+h_{e x t} d x} \cdot T(1)+\frac{2 h_{e x t} dx}{2 \lambda+h_{e x t} d x} \cdot T_{\text {ext }}
$$

Que l’on notera (_pour simplifier_):

$$
T(0)=\nu \cdot T(1)+(1-\nu) \cdot T_{e x t}
$$

> - $\nu = \frac{2 \lambda-h_{\text {ext }} d x}{2 \lambda+h_{e x t} d x}$

On fait ainsi apparaître l’enthalpie du noeud fantôme:

$$
\overline{h_{0}}=\nu \cdot h_{1}+(1-\nu)\cdot \rho c \cdot (T_{ext} - T_{f})
$$

> - $T_{f}$ est la température de fusion du MCP
> - On rappelle que selon [(3)](#eqn-3): $h(T) = \rho c \cdot(T-T_{f})$

En réinjectant dans l’équation [(5)](#eqn-5) traduisant les transferts:

$$
h_{1}^{t+1}=h_{1}^{t}+\alpha R\left[(\nu-2) \cdot h_{1}^{t+1}+h_{2}^{t+1}+(1-\nu) \cdot \rho c \cdot (T_{ext} - T_{f})\right]+\rho_{L} L\left(f^{t}-f^{t+1}\right)
$$

Que l’on met sous la forme (_pour résoudre matriciellement_):

$$
[1+\gamma \cdot(2-\nu)] \cdot h_{1}^{t+1} - \gamma \cdot h_{2}^{t+1} = h_{1}^{t} + [\gamma \cdot (1 - \nu) \cdot \rho c \cdot (T_{ext} - T_{f})] + \rho_{L} L\left(f^{t}-f^{t+1}\right)
$$

Soit en utilisant la forme de [(6)](#eqn-6):
> - $a_{i-1}=0$
> - $a_{i}=1+\gamma \cdot(2-\nu)$
> - $a_{i+1}=-\gamma$
> - $Q = h_{1}^{t} + [\gamma \cdot (1 - \nu) \cdot \rho c \cdot (T_{ext} - T_{f})] + \eta \cdot \left(f^{t}-f^{t+1}\right)$
> - $\nu = \frac{2 \lambda-h_{\text {ext }} d x}{2 \lambda+h_{e x t} d x}$
> _En rappelant:_ $\gamma = \alpha R$ et $\eta = \rho_{L}L$

### Interface Mur/Intérieur
> - On fixe la température intérieure à $T_{room}$
> - Cette température est valable _loin_ du mur
#### Point de vue du dernier noeud du mur

> - ![original image](https://cdn.mathpix.com/snip/images/wSNqO17KYJ6UyLSnta-9rAC142kbG24P4q-YRh2VfkA.original.fullsize.png)
> _Illustration du noeud fantôme utilisé (hachuré)_

En utilisant la même méthode que celle utilisée pour parvenir à [(Int.1)](#eqn-int-1), on obtient en isolant $\overline{T_{w}}(N+1)$:

$$
\overline{T_{w}}(N+1)=\varphi \cdot T_{w}(N)+(1-\varphi) \cdot T_{\text {room }}
$$

> - $\varphi=\frac{2\lambda-h_{\text {int}} \cdot dx}{2\lambda+h_{\text {int}}\cdot dx}$

En faisant apparaître l’enthalpie volumique du mur:

$$
\overline{h}(N+1)=\varphi \cdot h(N)+(1-\varphi) \cdot \rho c \cdot \left(T_{room} - T_f\right)
$$

> - $\rho$ et $c$ sont les paramètres thermodynamiques du mur.
> - $T_f$ est la température de fusion du mur

Finalement, en réinjectant dans l’équation [(5)](#eqn-5), traduisant les transferts puis en y mettant sous la forme de [(6)](#eqn-6) :

> - $a_{i-1}=-\gamma$
> - $a_{i}=1+\gamma \cdot(2-\varphi)$
> - $a_{i+1}=0$
> - $Q= h_i^{t} +\gamma \cdot(1-\varphi) \cdot \rho c \cdot \left(T_{room} - T_f\right)+\rho_{L} L \cdot\left(f^{t}-f^{t+1}\right)$


___
___


# Résolution Numérique

## Matrice globale
Le schéma numérique utilisé est ici implicite. Sa mise sous forme matricelle se fait en considérant un vecteur contenant les enthalpies volumiques de chaque noeud (_de la gauche du mur jusqu’à sa droite_) au temps $t$:

$$
H^{t}=\left(\begin{array}{c}
h_{1}^{t} \\
\vdots \\
h_{N}^{t}
\end{array}\right)
$$

Pour passer à $t+dt$, il convient alors de résoudre:

$$
\left(\begin{array}{ccccc}
a_{1}^{1} & a_{2}^{1} & 0 & \dots & \dots & 0\\
a_{1}^{2} & a_{2}^{2} & a_{3}^{2} & 0 & \dots & 0 \\
0 & \ddots & \ddots & \ddots & \ddots & \vdots \\
\vdots & \ddots & \ddots & \ddots & \ddots & 0\\
0 & \dots & 0 & a_{N-2}^{N-1} & a_{N-1}^{N-1} & a_{N}^{N-1} \\
0 & \dots & \dots & 0 & a_{N-1}^{N} & a_{N}^{N}
\end{array}\right)\cdot 
\left(\begin{array}{c}
h_{1}^{t+d t} \\
\vdots \\
\vdots \\
\vdots \\
\vdots \\
h_{N}^{t+d t}
\end{array}\right)=
\left(\begin{array}{c}
h_{1}^{t} \\
\vdots \\
\vdots \\
\vdots \\
\vdots \\
h_{N}^{t}
\end{array}\right)+\eta\cdot
\left(\begin{array}{c}
f_{1}^{t}-f_{1}^{t+d t} \\
\vdots \\
\vdots \\
\vdots \\
\vdots \\
f_{N}^{t}-f_{N}^{t+d t}
\end{array}\right)+E
$$

> - $a_{i-1}^{j}$, $a_{i}^{j}$ et $a_{i+1}^{j}$ correspondent aux coefficients trouvés dans les parties précédentes pour un noeud $j$
> - Les lignes correspondant à un noeud « solide constant » dans le vecteur des différences de fractions liquides seront nulles: étant donné que le matériau ne change pas de phase.
> - $E$ correspond aux termes supplémentaires relatifs aux bords/interfaces


Il convient donc de résoudre le système:

$$
A \cdot H^{t+d t}=H^{t} + \eta\Delta F + E
$$

ou:

$$
A \cdot H^{t+d t}=B
$$

> - $B = H^{t} + \eta\Delta F + E$
> - $A$ est une matrice tridiagonale que l’on peut inverser avec l’algorithme de Thomas

## Problème d’inconnues

> - Le nouveau terme de fraction liquide étant aussi une inconnue, on ne peut résoudre ce système  « d’un seul coup »: à chaque pas de temps, il conviendra donc de trouver le __nouveau__ vecteur __fraction liquide__

En reprenant l’équation « de base », pour un pas de temps quelconque __fixé__ (on n’affiche pas $t+dt$, on met simplement $old$ en exposant pour repérer les anciennes valeurs):

$$
a_Wh_W + a_P h_P + a_Eh_E
= h_P^{old} +\rho Lf^{old} −\rho Lf_k\ \ \left[2\right]
$$

> - $a_W \equiv a_{i-1}$
> - $a_P \equiv a_{i}$
> - $a_E \equiv a_{i+1}$

A chaque itération (dans un seul pas de temps), on applique un solveur TDMA à cette équation (où $f_k$ correspond à la $k^{eme}$ évaluation de la fraction liquide au nœud P). On met ensuite à jour le champ de fraction liquide (i.e. le tableau, vecteur entier des fractions liquides).
Cette __mise à jour__ est __cruciale__ dans cette méthode. Après la 1ère itération, on a :

$$
a_P h_P = −a_Eh_E − a_Wh_W + h_P^{old} +\rho Lf_l^{old} −\rho Lf_0 \ \
\left[3\right]
$$

Donc, après la $(k+1)^{eme}$ itération, on a :

$$
a_P h_P = −a_Eh_E − a_Wh_W + h_P^{old} +\rho Lf^{old} −\rho Lf_k \ \
\left[4\right]
$$

Si un changement de phase se passe au $P^{eme}$ nœud (i.e. $0 < f_k < 1$ ), alors la $k^{eme}$ estimation de la fraction liquide — i.e. $f_k$— doit être mise à jour de sorte que pour la $(k+2)^{eme}$ itération, la partie gauche de [4] soit égale à zéro — grâce au nouveau $f_{k + 1}$ — c’est-à-dire :

$$
0 = −a_Eh_E − a_Wh_W + h_P^{old} +\rho Lf_l^{old} −\rho Lf_k + update\ \left[5\right]
$$

On a donc :

$$
update =\ −a_P h_P 
$$

L’équation étant utilisée pour la $(k+2)^{eme}$ itération étant :

$$
a_P h_P = −a_Eh_E − a_Wh_W + h_\rho^{old} +\rho Lf^{old} −\rho Lf_{k + 1}
$$

Il vient alors :

$$
−\rho Lf_k + update =\ −\rho Lf_{k + 1}
$$

Puis :

$$
f_{k + 1} = f_k +\frac{a_Ph_P}{\rho L} 
$$

Que l’on peut modifier en :

$$
f_{k + 1}
=f_k +\lambda\frac{a_P h_P}{\rho L}\ \ \left[6\right]
$$

Où $\lambda$ est un coefficient de sous-relaxation qui permet d’atteindre la convergence.


En pratique, après la $k^{eme}$ solution du TDMA, on applique cette mise à jour à chaque nœud (pas seulement aux nœuds en changement de phase). Comme cette mise à jour n’est pas adaptée à tous les nœuds, on la suit d’une correction de dépassement/sous-dépassement:

$$
f=\left\{\begin{array}{ll}
0 & \text { si } f_{k+1}<0 \\
1 & \text { si } f_{k+1}>1
\end{array}\right. \ \ \left[7\right]
$$

> - Cette correction a deux avantages :
> - 1. Il n’y a pas besoin de vérifier où et quand appliquer la correction
> - 2. Il y a une totale et correcte prise en charge des situations quand le front passe d’un nœud à
un autre.
>- Le secret de cette méthode réside dans la mise à jour de la fraction liquide depuis le vecteur d’enthalpie sensible : [6] était l’équation utilisée. 

### Optimisation

Ici, on introduit une modification qui permettra des calculs plus rapides.
En effet, on n’utilisait pas une information potentiellement importante : si la fraction liquide au nœud $P$ est strictement dans l’intervalle $]0, 1[$ alors, $h_P = 0$, (selon la définition de l’enthalpie [(3)](#eqn-3)). 
Cette information peut être _forcée_ sur le solveur TDMA en mettant simplement le coefficient $a_P = BIG\left({10}^{15}\right)$ à chaque fois que $0 < f_k < 1$ pour le nœud en question. De cette manière, le solveur TDMA retournera une valeur pour $\left(h_P\right)_{k + 1}$ proche de zéro.
En mettant par la suite $h_P = 0$ on a la nouvelle mise à jour suivante:

$$
\left(f_l\right)_{k + 1}
=\frac{−a_Wh_W − a_Eh_E + h_P^{old}}{\rho L} + f^{old}\ \ \left[8\right]
$$

Qui ne nécessite pas de coefficient de relaxation. Le nouveau critère de convergence sera :

$$
\frac{\sum ABS\left({RES}_P\right)}{c_P}< TOL\ \ \left[9\right]
$$


> - ${RES}_P = a_Wh_w + a_P h_P + a_Eh_E − h_P^{old} −\rho Lf^{old} +\rho\ Lf_k$

---

<h4 id="prog"></h4>

# Fonctionnement du programme Python

## Librairies/Modules utilisés

- numpy
- pathlib
- matplotlib
- math

## Configuration de l'expérience

Cette partie se passe exclusivement dans le dossier `Config`

1. Ajouter des matériaux dans le fichier `materials_config.csv` (ces matériaux n'ont pas besoin d'être utilisés par la suite)
2. Configurer le mur dans le fichier `wall_config.csv`:
  > - chaque ligne correspond à une épaisseur de matériau (le haut correspondant à la gauche du mur)
  > - le nom du matériau d'une couche doit correspondre à un matériau dans `materials_config.csv`
  > - l'épaisseur (en m) doit être indiquée sur la colonne de droite pour chaque épaisseur de matériau

3. Configurer certaines données de l'expérience dans le fichier `experiment_config.ini` (se référer aux commentaires présents)
4. Modifier la loi de température extérieure dans `temperature_config.py`, e.g. si la température extérieure doit être constante:
```python
def update_outside_temperature(t: float) -> float:
    pass
```

## Simulation

Lancer `main.py`, une barre de chargement indiquant l'avancée du calcul devrait apparaître dans la console. Une fois terminé, une fenêtre matplotlib affiche les résultats de l'expérience. 
Pour supprimer le darkmode, accéder à `Plot\standard_plot.py` puis supprimer la première ligne de la première fonction:

```python
def plot_simulation(simulation: Simulation):
    plt.style.use('dark_background')
```

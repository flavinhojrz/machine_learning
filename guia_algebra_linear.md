# Guia Completo de Álgebra Linear
## Matrizes, Sistemas Lineares e Suas Propriedades

> **Para o estudante:** Este guia foi concebido para ser lido sequencialmente. Cada tópico constrói sobre o anterior. Não pule etapas — a intuição é tão importante quanto o cálculo.

---

## Tópico 1 — Matrizes: Definição e Tipos Especiais

### 1.1 O que é uma Matriz? (A Intuição)

Imagine que você precisa organizar as notas de 3 alunos em 4 disciplinas. Você naturalmente faria uma **tabela**. Uma matriz é exatamente isso: uma forma matemática estruturada de organizar informação em **linhas e colunas**.

Na computação, uma matriz é um **array bidimensional**. Em Python com NumPy: `np.array([[1, 2], [3, 4]])` é uma matriz $2 \times 2$.

### 1.2 Definição Formal

Uma matriz $A$ de ordem $m \times n$ é um arranjo retangular de $m \cdot n$ números (reais ou complexos) dispostos em $m$ **linhas** e $n$ **colunas**:

$$A = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix}$$

**Notação:** O elemento na linha $i$ e coluna $j$ é denotado $a_{ij}$ (primeiro índice = linha, segundo = coluna). Escrevemos $A = [a_{ij}]_{m \times n}$.

> **Regra mnemônica:** "$ij$" = "**I**ndo **J**untos" → "linha **i**, coluna **j**".

### 1.3 Tipos Especiais de Matrizes

#### Matriz Quadrada
$m = n$. Tem o mesmo número de linhas e colunas. A **diagonal principal** é formada pelos elementos $a_{11}, a_{22}, \ldots, a_{nn}$.

#### Matriz Identidade $I_n$
A "unidade" das matrizes. É quadrada com $a_{ii} = 1$ (diagonal) e $a_{ij} = 0$ para $i \neq j$. Satisfaz $AI = IA = A$.

$$I_3 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

#### Matriz Diagonal
Quadrada com todos os elementos fora da diagonal iguais a zero: $a_{ij} = 0$ para $i \neq j$.

$$D = \begin{pmatrix} 3 & 0 & 0 \\ 0 & -1 & 0 \\ 0 & 0 & 5 \end{pmatrix}$$

> **Intuição geométrica:** Multiplicar por uma matriz diagonal **escala** cada eixo coordenado de forma independente. $D$ acima estica o eixo $x$ por 3, inverte $y$ e estica $z$ por 5.

#### Matriz Triangular Superior
$a_{ij} = 0$ para todo $i > j$ (zeros abaixo da diagonal):

$$U = \begin{pmatrix} 2 & 5 & 1 \\ 0 & 3 & -4 \\ 0 & 0 & 7 \end{pmatrix}$$

#### Matriz Triangular Inferior
$a_{ij} = 0$ para todo $i < j$ (zeros acima da diagonal):

$$L = \begin{pmatrix} 2 & 0 & 0 \\ 5 & 3 & 0 \\ 1 & -4 & 7 \end{pmatrix}$$

> **Por que importam?** Sistemas triangulares são **trivialmente resolúveis** por substituição direta. A eliminação de Gauss transforma qualquer sistema nesta forma.

#### Matriz Simétrica
Quadrada com $a_{ij} = a_{ji}$ para todos $i, j$, ou seja, $A = A^T$.

$$S = \begin{pmatrix} 1 & 4 & 2 \\ 4 & 0 & -3 \\ 2 & -3 & 5 \end{pmatrix}$$

Matrizes simétricas aparecem naturalmente em **covariância estatística** e em sistemas de equações derivadas de energia potencial.

#### Matriz Nula
Todos os elementos são zero: $a_{ij} = 0\;\forall\,i,j$. Denotada $\mathbf{0}$.

---

## Tópico 2 — Operações Básicas

### 2.1 Adição e Subtração

**Regra de existência:** Só podemos somar (ou subtrair) matrizes de **mesma ordem**.

**Definição:** $(A \pm B)_{ij} = a_{ij} \pm b_{ij}$. Operamos elemento a elemento.

**Exemplo:**
$$\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix} + \begin{pmatrix} 5 & 0 \\ -1 & 2 \end{pmatrix} = \begin{pmatrix} 6 & 2 \\ 2 & 6 \end{pmatrix}$$

**Propriedades da adição:**
- **Comutatividade:** $A + B = B + A$
- **Associatividade:** $(A + B) + C = A + (B + C)$
- **Elemento neutro:** $A + \mathbf{0} = A$
- **Elemento oposto:** $A + (-A) = \mathbf{0}$

### 2.2 Multiplicação por Escalar

Multiplicar uma matriz por um número real (escalar) $k$ é simples: **multiplique cada elemento por $k$**.

$$(\alpha A)_{ij} = \alpha \cdot a_{ij}$$

**Exemplo:**

$$3 \cdot \begin{pmatrix} 1 & -2 \\ 0 & 4 \end{pmatrix} = \begin{pmatrix} 3 & -6 \\ 0 & 12 \end{pmatrix}$$

**Propriedades:**
- $\alpha(A + B) = \alpha A + \alpha B$
- $(\alpha + \beta)A = \alpha A + \beta A$
- $(\alpha \beta)A = \alpha(\beta A)$
- $1 \cdot A = A$

### 2.3 Matriz Transposta

**Intuição:** Transpor uma matriz é como "virar" ela ao longo da diagonal principal — as linhas tornam-se colunas e vice-versa.

**Definição:** Se $A$ é $m \times n$, então $A^T$ é $n \times m$, onde $(A^T)_{ij} = a_{ji}$.

**Exemplo:**

$$A = \begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \end{pmatrix}_{2 \times 3} \implies A^T = \begin{pmatrix} 1 & 4 \\ 2 & 5 \\ 3 & 6 \end{pmatrix}_{3 \times 2}$$

**Propriedades da transposta:**

| Propriedade | Fórmula |
|---|---|
| Involução | $(A^T)^T = A$ |
| Linearidade | $(\alpha A + \beta B)^T = \alpha A^T + \beta B^T$ |
| Produto (ordem invertida) | $(AB)^T = B^T A^T$ |
| Determinante preservado | $\det(A^T) = \det(A)$ |

> A propriedade $(AB)^T = B^T A^T$ é **contra-intuitiva**: a ordem se inverte! Isso ocorre porque a transposição "reflete" as operações matriciais.

---

## Tópico 3 — Produto de Matrizes

### 3.1 Por que o Produto é Definido Assim?

Antes do "como", entenda o "porquê". O produto de matrizes modela a **composição de transformações lineares**. Se $B$ representa uma rotação e $A$ representa uma escala, então $AB$ representa "primeiro rotacione, depois escale". Essa composição **não é comutativa** — a ordem importa.

### 3.2 Definição Formal

Para $A$ ser $m \times n$ e $B$ ser $n \times p$, o produto $C = AB$ é $m \times p$, com:

$$c_{ij} = \sum_{k=1}^{n} a_{ik} \cdot b_{kj}$$

Ou seja: o elemento $(i,j)$ de $C$ é o **produto escalar** da linha $i$ de $A$ com a coluna $j$ de $B$.

> **Regra de existência:** O número de **colunas** de $A$ deve ser igual ao número de **linhas** de $B$. Memorize: $A_{m \times \mathbf{n}} \cdot B_{\mathbf{n} \times p} = C_{m \times p}$.

### 3.3 Exemplo de Produto 2×2

$$A = \begin{pmatrix} 2 & 1 \\ 3 & 0 \end{pmatrix}, \quad B = \begin{pmatrix} 1 & 4 \\ 2 & -1 \end{pmatrix}$$

$$c_{11} = 2 \cdot 1 + 1 \cdot 2 = 4, \quad c_{12} = 2 \cdot 4 + 1 \cdot (-1) = 7$$
$$c_{21} = 3 \cdot 1 + 0 \cdot 2 = 3, \quad c_{22} = 3 \cdot 4 + 0 \cdot (-1) = 12$$

$$AB = \begin{pmatrix} 4 & 7 \\ 3 & 12 \end{pmatrix}$$

### 3.4 Demonstração da Não-Comutatividade

> **O produto de matrizes NÃO é comutativo em geral:** $AB \neq BA$.

Usando as mesmas matrizes acima:

$$BA = \begin{pmatrix} 1 & 4 \\ 2 & -1 \end{pmatrix} \begin{pmatrix} 2 & 1 \\ 3 & 0 \end{pmatrix}$$

$$b_{11}a_{11} \text{ col}: \quad 1 \cdot 2 + 4 \cdot 3 = 14, \quad 1 \cdot 1 + 4 \cdot 0 = 1$$
$$\phantom{b_{11}a_{11} \text{ col}: \quad} 2 \cdot 2 + (-1) \cdot 3 = 1, \quad 2 \cdot 1 + (-1) \cdot 0 = 2$$

$$BA = \begin{pmatrix} 14 & 1 \\ 1 & 2 \end{pmatrix}$$

Comparando:
$$AB = \begin{pmatrix} 4 & 7 \\ 3 & 12 \end{pmatrix} \neq \begin{pmatrix} 14 & 1 \\ 1 & 2 \end{pmatrix} = BA \quad \checkmark$$

### 3.5 Propriedades do Produto

| Propriedade | Fórmula | Observação |
|---|---|---|
| Associatividade | $(AB)C = A(BC)$ | Sempre vale |
| Distributividade à esquerda | $A(B+C) = AB + AC$ | Sempre vale |
| Distributividade à direita | $(A+B)C = AC + BC$ | Sempre vale |
| Escalar | $\alpha(AB) = (\alpha A)B = A(\alpha B)$ | Sempre vale |
| Identidade | $AI = IA = A$ | Sempre vale |
| **NÃO comutatividade** | $AB \neq BA$ em geral | **Atenção!** |
| **NÃO cancelamento** | $AB = AC \not\Rightarrow B = C$ | Se $A \neq \mathbf{0}$ |

---

## Tópico 4 — Matrizes Inversas

### 4.1 Intuição

A inversa de uma matriz $A$ é o análogo matricial de $1/a$ para números. Assim como $a \cdot (1/a) = 1$, temos $A \cdot A^{-1} = I$. Geometricamente, $A^{-1}$ **desfaz** a transformação representada por $A$.

### 4.2 Definição e Critério de Existência

A matriz $A$ (quadrada, $n \times n$) é **inversível** (ou não-singular) se existe $B$ tal que $AB = BA = I_n$. Nesse caso, $B = A^{-1}$.

> **Critério fundamental:** $A$ é inversível **se e somente se** $\det(A) \neq 0$.

Se $\det(A) = 0$, $A$ é chamada **singular** (ou degenerada) e não possui inversa.

### 4.3 Fórmula Direta para Matrizes 2×2

Para $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$ com $\det(A) = ad - bc \neq 0$:

$$A^{-1} = \frac{1}{ad - bc} \begin{pmatrix} d & -b \\ -c & a \end{pmatrix}$$

**A receita:** Troque os elementos da diagonal principal, troque os sinais da diagonal secundária, divida pelo determinante.

**Exemplo:**
$$A = \begin{pmatrix} 3 & 1 \\ 5 & 2 \end{pmatrix}, \quad \det(A) = 6 - 5 = 1$$

$$A^{-1} = \frac{1}{1}\begin{pmatrix} 2 & -1 \\ -5 & 3 \end{pmatrix} = \begin{pmatrix} 2 & -1 \\ -5 & 3 \end{pmatrix}$$

**Verificação:**
$$AA^{-1} = \begin{pmatrix} 3 & 1 \\ 5 & 2 \end{pmatrix}\begin{pmatrix} 2 & -1 \\ -5 & 3 \end{pmatrix} = \begin{pmatrix} 6-5 & -3+3 \\ 10-10 & -5+6 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I \checkmark$$

### 4.4 Introdução ao Método Gauss-Jordan (visão geral)

Para matrizes maiores, a fórmula direta torna-se impraticável. O método **Gauss-Jordan** é mais geral: construímos a matriz aumentada $[A \mid I]$ e aplicamos operações elementares de linha até transformar $A$ na identidade. O lado direito automaticamente torna-se $A^{-1}$:

$$[A \mid I] \xrightarrow{\text{operações}} [I \mid A^{-1}]$$

O detalhamento completo está no Tópico 8.

### 4.5 Propriedades da Inversa

- $(A^{-1})^{-1} = A$
- $(AB)^{-1} = B^{-1}A^{-1}$ (ordem invertida!)
- $(A^T)^{-1} = (A^{-1})^T$
- $(\alpha A)^{-1} = \frac{1}{\alpha} A^{-1}$, para $\alpha \neq 0$
- $\det(A^{-1}) = \frac{1}{\det(A)}$

---

## Tópico 5 — Matrizes Elementares

### 5.1 O que São e Por que Importam

Uma **matriz elementar** $E$ é obtida aplicando **uma única operação elementar de linha** à matriz identidade $I$. O motivo de estudá-las: **toda operação de linha em $A$ pode ser representada como uma multiplicação à esquerda por $E$**.

Isso conecta as operações "mecânicas" da eliminação de Gauss com a teoria algébrica das inversas.

### 5.2 Os Três Tipos

#### Tipo I — Troca de Linhas ($L_i \leftrightarrow L_j$)

Troca as linhas $i$ e $j$ de $I$. Trocar a linha 1 e a linha 2 em $I_3$:

$$E_1 = \begin{pmatrix} 0 & 1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

Multiplicar $A$ à esquerda por $E_1$ troca as linhas 1 e 2 de $A$.

**Inversa:** $(E_1)^{-1} = E_1$ (a própria $E_1$, pois trocar duas vezes restaura a ordem).

#### Tipo II — Escalonamento ($L_i \leftarrow \alpha L_i$, $\alpha \neq 0$)

Multiplica a linha $i$ de $I$ por um escalar $\alpha$. Multiplicar a linha 2 por $-3$:

$$E_2 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & -3 & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

**Inversa:** $(E_2)^{-1}$ é obtida colocando $1/\alpha$ no lugar de $\alpha$:

$$(E_2)^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & -1/3 & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

#### Tipo III — Combinação Linear ($L_i \leftarrow L_i + \alpha L_j$)

Soma $\alpha$ vezes a linha $j$ à linha $i$. Somar $3 \times$ linha 1 à linha 2:

$$E_3 = \begin{pmatrix} 1 & 0 & 0 \\ 3 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

**Inversa:** $(E_3)^{-1}$ usa $-\alpha$:

$$(E_3)^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ -3 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

> **Intuição:** Para desfazer "some 3 vezes a linha 1 à linha 2", basta "subtraia 3 vezes a linha 1 da linha 2".

### 5.3 Conexão com a Eliminação de Gauss

A eliminação de Gauss aplicada a uma matriz $A$ nada mais é do que a multiplicação sequencial por matrizes elementares do Tipo III (principalmente):

$$E_k \cdots E_2 E_1 A = U \quad \text{(forma escalonada)}$$

Isso implica que toda matriz inversível pode ser **fatorada** como:

$$A = L \cdot U$$

onde $L$ é uma matriz triangular inferior (produto das inversas das elementares) e $U$ é triangular superior. Esta é a famosa **fatoração LU**, base de muitos algoritmos numéricos.

---

## Tópico 6 — Sistemas Lineares

### 6.1 Representação Matricial $Ax = b$

Um sistema de $m$ equações e $n$ incógnitas:

$$\begin{cases} a_{11}x_1 + a_{12}x_2 + \cdots + a_{1n}x_n = b_1 \\ a_{21}x_1 + a_{22}x_2 + \cdots + a_{2n}x_n = b_2 \\ \quad\vdots \\ a_{m1}x_1 + a_{m2}x_2 + \cdots + a_{mn}x_n = b_m \end{cases}$$

pode ser escrito compactamente como $Ax = b$, onde:
- $A$ é a **matriz dos coeficientes** ($m \times n$)
- $x$ é o **vetor de incógnitas** ($n \times 1$)
- $b$ é o **vetor dos termos independentes** ($m \times 1$)

### 6.2 Matriz Aumentada

Para resolver o sistema, construímos a **matriz aumentada** $[A \mid b]$:

$$[A \mid b] = \left(\begin{array}{cccc|c} a_{11} & a_{12} & \cdots & a_{1n} & b_1 \\ a_{21} & a_{22} & \cdots & a_{2n} & b_2 \\ \vdots & & \ddots & \vdots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} & b_m \end{array}\right)$$

Todas as operações de linha são aplicadas nesta matriz aumentada.

### 6.3 Classificação dos Sistemas e Teorema de Rouché-Capelli

Definimos duas quantidades-chave:
- $\text{posto}(A) = r$: número de linhas linearmente independentes de $A$
- $\text{posto}([A \mid b]) = r'$: posto da matriz aumentada

> **Teorema de Rouché-Capelli:** O sistema $Ax = b$ é **compatível** (tem solução) se e somente se $\text{posto}(A) = \text{posto}([A \mid b])$.

#### Sistema Possível e Determinado (SPD)
$$\text{posto}(A) = \text{posto}([A \mid b]) = n$$
Solução **única**. Geometricamente: as hiperplanas se intersectam num único ponto.

#### Sistema Possível e Indeterminado (SPI)
$$\text{posto}(A) = \text{posto}([A \mid b]) = r < n$$
**Infinitas soluções** (família paramétrica com $n - r$ parâmetros livres). As hiperplanas se intersectam numa reta, plano ou subespaço de dimensão maior.

#### Sistema Impossível (SI)
$$\text{posto}(A) < \text{posto}([A \mid b])$$
**Sem solução**. As hiperplanas não possuem ponto comum.

**Exemplo de classificação:**
$$[A \mid b] = \left(\begin{array}{cc|c} 1 & 2 & 5 \\ 2 & 4 & 10 \end{array}\right) \xrightarrow{L_2 \leftarrow L_2 - 2L_1} \left(\begin{array}{cc|c} 1 & 2 & 5 \\ 0 & 0 & 0 \end{array}\right)$$

$\text{posto}(A) = \text{posto}([A \mid b]) = 1 < 2 = n$ → **SPI** (infinitas soluções: $x_1 = 5 - 2t$, $x_2 = t$, $t \in \mathbb{R}$).

---

## Tópico 7 — Eliminação de Gauss

### 7.1 O Algoritmo

A eliminação de Gauss transforma a matriz aumentada $[A \mid b]$ em uma forma **triangular superior** (escalonada), que é então resolvida por **substituição retroativa**.

**Passos formais:**
1. Identifique o elemento **pivô** (primeiro não-nulo da coluna atual).
2. Se necessário, troque linhas (**pivoteamento parcial**: escolha o maior pivô em módulo).
3. Elimine todos os elementos **abaixo** do pivô usando operações de Tipo III.
4. Avance para a próxima coluna à direita e repita.
5. Ao atingir a forma triangular, aplique **substituição retroativa**.

### 7.2 Exemplo Completo 3×3

**Sistema:**
$$\begin{cases} 2x + y - z = 8 \\ -3x - y + 2z = -11 \\ -2x + y + 2z = -3 \end{cases}$$

**Matriz aumentada inicial:**
$$\left(\begin{array}{ccc|c} 2 & 1 & -1 & 8 \\ -3 & -1 & 2 & -11 \\ -2 & 1 & 2 & -3 \end{array}\right)$$

---

**Etapa 1:** Pivot na posição (1,1) = 2. Eliminar coluna 1 abaixo do pivô.

$L_2 \leftarrow L_2 + \dfrac{3}{2}L_1$:
$$L_2: \; -3+3 \;=\; 0, \quad -1+\frac{3}{2} \;=\; \frac{1}{2}, \quad 2-\frac{3}{2} \;=\; \frac{1}{2}, \quad -11+12 \;=\; 1$$

$L_3 \leftarrow L_3 + L_1$:
$$L_3: \; -2+2 \;=\; 0, \quad 1+1 \;=\; 2, \quad 2-1 \;=\; 1, \quad -3+8 \;=\; 5$$

$$\left(\begin{array}{ccc|c} 2 & 1 & -1 & 8 \\ 0 & \tfrac{1}{2} & \tfrac{1}{2} & 1 \\ 0 & 2 & 1 & 5 \end{array}\right)$$

---

**Etapa 2:** Pivô na posição (2,2) = 1/2. Eliminar coluna 2 abaixo do pivô.

$L_3 \leftarrow L_3 - 4 \cdot L_2$:
$$L_3: \; 0-0=0, \quad 2-2=0, \quad 1-2=-1, \quad 5-4=1$$

$$\left(\begin{array}{ccc|c} 2 & 1 & -1 & 8 \\ 0 & \tfrac{1}{2} & \tfrac{1}{2} & 1 \\ 0 & 0 & -1 & 1 \end{array}\right)$$

A matriz está na **forma escalonada triangular superior**. ✓

---

**Etapa 3: Substituição retroativa**

**Da linha 3:**
$$-z = 1 \implies z = -1$$

**Da linha 2:**
$$\frac{1}{2}y + \frac{1}{2}(-1) = 1 \implies \frac{y}{2} = \frac{3}{2} \implies y = 3$$

**Da linha 1:**
$$2x + 3 - (-1) = 8 \implies 2x + 4 = 8 \implies x = 2$$

**Solução:** $x = 2,\; y = 3,\; z = -1$.

**Verificação:**
- $2(2) + 3 - (-1) = 4 + 3 + 1 = 8$ ✓
- $-3(2) - 3 + 2(-1) = -6 - 3 - 2 = -11$ ✓
- $-2(2) + 3 + 2(-1) = -4 + 3 - 2 = -3$ ✓

---

## Tópico 8 — Método de Gauss-Jordan e FERL

### 8.1 Forma Escalonada Reduzida por Linhas (FERL)

A eliminação de Gauss produz uma forma **triangular**. O Gauss-Jordan vai além: produz a **Forma Escalonada Reduzida por Linhas (FERL)**, onde os elementos **acima e abaixo** de cada pivô são zerados, e cada pivô vale 1.

**Critérios da FERL:**
1. Todas as linhas nulas, se existirem, estão abaixo das não-nulas.
2. O primeiro elemento não-nulo de cada linha (o **pivô**) vale 1.
3. Cada pivô está à **direita** dos pivôs das linhas superiores.
4. Cada coluna pivô tem **apenas zeros** exceto no pivô.

### 8.2 Exemplo Numérico Completo 3×3 (FERL)

Usando o mesmo sistema do Tópico 7, partindo de onde paramos:

$$\left(\begin{array}{ccc|c} 2 & 1 & -1 & 8 \\ 0 & \tfrac{1}{2} & \tfrac{1}{2} & 1 \\ 0 & 0 & -1 & 1 \end{array}\right)$$

**Etapa 4:** Tornar os pivôs iguais a 1.

$L_1 \leftarrow \dfrac{1}{2}L_1$, $L_2 \leftarrow 2L_2$, $L_3 \leftarrow (-1) L_3$:

$$\left(\begin{array}{ccc|c} 1 & \tfrac{1}{2} & -\tfrac{1}{2} & 4 \\ 0 & 1 & 1 & 2 \\ 0 & 0 & 1 & -1 \end{array}\right)$$

**Etapa 5:** Eliminar acima do pivô da coluna 3.

$L_2 \leftarrow L_2 - L_3$: $\quad [0,\; 1,\; 0 \mid 3]$

$L_1 \leftarrow L_1 + \tfrac{1}{2}L_3$: $\quad [1,\; \tfrac{1}{2},\; 0 \mid \tfrac{7}{2}]$

$$\left(\begin{array}{ccc|c} 1 & \tfrac{1}{2} & 0 & \tfrac{7}{2} \\ 0 & 1 & 0 & 3 \\ 0 & 0 & 1 & -1 \end{array}\right)$$

**Etapa 6:** Eliminar acima do pivô da coluna 2.

$L_1 \leftarrow L_1 - \tfrac{1}{2}L_2$: $\quad [1,\; 0,\; 0 \mid 2]$

$$\left(\begin{array}{ccc|c} 1 & 0 & 0 & 2 \\ 0 & 1 & 0 & 3 \\ 0 & 0 & 1 & -1 \end{array}\right)$$

A solução aparece diretamente: $x=2,\; y=3,\; z=-1$. ✓

### 8.3 Cálculo de $A^{-1}$ via Gauss-Jordan

**Ideia:** Partimos de $[A \mid I]$ e aplicamos operações de linha até obter $[I \mid A^{-1}]$.

**Por quê funciona?** Se $E_k \cdots E_1 A = I$, então $E_k \cdots E_1 = A^{-1}$. As mesmas operações aplicadas a $I$ produzem $A^{-1}$.

**Exemplo com** $A = \begin{pmatrix} 1 & 0 & 2 \\ 0 & 1 & 3 \\ 1 & 2 & 0 \end{pmatrix}$, $\det(A) = -8 \neq 0$.

**Matriz aumentada $[A \mid I_3]$:**

$$\left(\begin{array}{ccc|ccc} 1 & 0 & 2 & 1 & 0 & 0 \\ 0 & 1 & 3 & 0 & 1 & 0 \\ 1 & 2 & 0 & 0 & 0 & 1 \end{array}\right)$$

**Passo 1:** $L_3 \leftarrow L_3 - L_1$:

$$\left(\begin{array}{ccc|ccc} 1 & 0 & 2 & 1 & 0 & 0 \\ 0 & 1 & 3 & 0 & 1 & 0 \\ 0 & 2 & -2 & -1 & 0 & 1 \end{array}\right)$$

**Passo 2:** $L_3 \leftarrow L_3 - 2L_2$:

$$\left(\begin{array}{ccc|ccc} 1 & 0 & 2 & 1 & 0 & 0 \\ 0 & 1 & 3 & 0 & 1 & 0 \\ 0 & 0 & -8 & -1 & -2 & 1 \end{array}\right)$$

**Passo 3:** $L_3 \leftarrow -\dfrac{1}{8}L_3$:

$$\left(\begin{array}{ccc|ccc} 1 & 0 & 2 & 1 & 0 & 0 \\ 0 & 1 & 3 & 0 & 1 & 0 \\ 0 & 0 & 1 & \tfrac{1}{8} & \tfrac{1}{4} & -\tfrac{1}{8} \end{array}\right)$$

**Passo 4:** Eliminar acima. $L_1 \leftarrow L_1 - 2L_3$, $L_2 \leftarrow L_2 - 3L_3$:

$$L_1: \left[1,\; 0,\; 0 \;\Big|\; 1-\tfrac{1}{4},\; -\tfrac{1}{2},\; \tfrac{1}{4}\right] = \left[1,0,0\;\Big|\;\tfrac{3}{4},\;-\tfrac{1}{2},\;\tfrac{1}{4}\right]$$

$$L_2: \left[0,\; 1,\; 0 \;\Big|\; -\tfrac{3}{8},\; 1-\tfrac{3}{4},\; \tfrac{3}{8}\right] = \left[0,1,0\;\Big|\;-\tfrac{3}{8},\;\tfrac{1}{4},\;\tfrac{3}{8}\right]$$

**Resultado:**

$$A^{-1} = \begin{pmatrix} \tfrac{3}{4} & -\tfrac{1}{2} & \tfrac{1}{4} \\[4pt] -\tfrac{3}{8} & \tfrac{1}{4} & \tfrac{3}{8} \\[4pt] \tfrac{1}{8} & \tfrac{1}{4} & -\tfrac{1}{8} \end{pmatrix}$$

---

## Tópico 9 — Determinantes

### 9.1 Intuição Geométrica Primeiro

O determinante de uma matriz $A$ mede o **fator de escala de volume** da transformação linear que $A$ representa:
- Para $2 \times 2$: $|\det(A)|$ é a **área** do paralelogramo formado pelas colunas (ou linhas) de $A$.
- Para $3 \times 3$: $|\det(A)|$ é o **volume** do paralelepípedo formado pelas colunas de $A$.
- O **sinal** indica orientação: positivo = orientação preservada; negativo = orientação invertida.
- $\det(A) = 0$: a transformação **colapsa** o espaço em dimensão menor (as colunas são linearmente dependentes).

### 9.2 Determinante 2×2

$$\det\begin{pmatrix} a & b \\ c & d \end{pmatrix} = ad - bc$$

**Exemplo:**
$$\det\begin{pmatrix} 3 & 2 \\ 1 & 4 \end{pmatrix} = 3 \cdot 4 - 2 \cdot 1 = 12 - 2 = 10$$

### 9.3 Determinante 3×3 — Regra de Sarrus

Para $A = \begin{pmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{pmatrix}$:

**Sarrus:** Repita as duas primeiras colunas à direita da matriz e calcule as diagonais.

$$\det(A) = (a_{11}a_{22}a_{33} + a_{12}a_{23}a_{31} + a_{13}a_{21}a_{32}) - (a_{13}a_{22}a_{31} + a_{11}a_{23}a_{32} + a_{12}a_{21}a_{33})$$

> ⚠️ **Atenção:** A Regra de Sarrus funciona **apenas** para matrizes $3 \times 3$.

**Exemplo:**
$$A = \begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{pmatrix}$$

Diagonais positivas: $(1 \cdot 5 \cdot 9) + (2 \cdot 6 \cdot 7) + (3 \cdot 4 \cdot 8) = 45 + 84 + 96 = 225$

Diagonais negativas: $(3 \cdot 5 \cdot 7) + (1 \cdot 6 \cdot 8) + (2 \cdot 4 \cdot 9) = 105 + 48 + 72 = 225$

$\det(A) = 225 - 225 = 0$

O determinante nulo indica que as linhas (ou colunas) são **linearmente dependentes** (note que $L_3 = 2L_2 - L_1$). Esta matriz é **singular** e não tem inversa.

### 9.4 Expansão em Cofatores (Teorema de Laplace)

Para matrizes $n \times n$, o determinante pode ser calculado expandindo por qualquer linha ou coluna.

**Definições:**
- **Menor** $M_{ij}$: determinante da submatriz obtida removendo a linha $i$ e coluna $j$.
- **Cofator** $C_{ij} = (-1)^{i+j} M_{ij}$. O sinal segue o padrão: $\begin{pmatrix} + & - & + \\ - & + & - \\ + & - & + \end{pmatrix}$

**Expansão pela linha $i$:**

$$\det(A) = \sum_{j=1}^{n} a_{ij} \cdot C_{ij}$$

**Exemplo** (expandindo pela 1ª linha da matriz $A$ anterior com det = 10):

$$A = \begin{pmatrix} 1 & 2 & 1 \\ 0 & 3 & -1 \\ 2 & 1 & 4 \end{pmatrix}$$

$$\det(A) = 1 \cdot \det\begin{pmatrix}3&-1\\1&4\end{pmatrix} - 2 \cdot \det\begin{pmatrix}0&-1\\2&4\end{pmatrix} + 1 \cdot \det\begin{pmatrix}0&3\\2&1\end{pmatrix}$$

$$= 1(12+1) - 2(0+2) + 1(0-6) = 13 - 4 - 6 = 3$$

### 9.5 Propriedades dos Determinantes

| Propriedade | Descrição |
|---|---|
| $\det(I) = 1$ | Identidade tem det 1 |
| $\det(AB) = \det(A)\det(B)$ | Produto de dets |
| $\det(A^T) = \det(A)$ | Transposta não altera |
| $\det(A^{-1}) = 1/\det(A)$ | Inversa inverte o det |
| Troca de linhas | Multiplica det por $-1$ |
| Linha proporcional | $\det = 0$ |
| $L_i \leftarrow L_i + \alpha L_j$ | Não altera o det |
| $\det(\alpha A) = \alpha^n \det(A)$ | Escalonamento |

---

## Tópico 10 — Autovalores e Autovetores

### 10.1 Definição Geométrica e Algébrica

**Intuição:** A maioria dos vetores muda de direção quando multiplicados por uma matriz $A$. Um **autovetor** (ou vetor próprio) é tão especial que, ao ser transformado por $A$, apenas muda de **escala** — não de direção. O fator de escala é o **autovalor** (valor próprio) correspondente.

**Definição algébrica:** $\lambda$ é autovalor e $v \neq \mathbf{0}$ é autovetor associado se:

$$Av = \lambda v$$

Equivalentemente: $Av - \lambda v = \mathbf{0}$, ou $(A - \lambda I)v = \mathbf{0}$.

Para que este sistema homogêneo tenha solução não trivial, é necessário:

$$\det(A - \lambda I) = 0$$

Esta é a **equação característica** de $A$.

### 10.2 Cálculo dos Autovalores

O polinômio $p(\lambda) = \det(A - \lambda I)$ é chamado **polinômio característico**. Suas raízes são os autovalores.

**Exemplo 2×2:**
$$A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$$

$$A - \lambda I = \begin{pmatrix} 2-\lambda & 1 \\ 1 & 2-\lambda \end{pmatrix}$$

$$\det(A - \lambda I) = (2-\lambda)^2 - 1 = \lambda^2 - 4\lambda + 3 = (\lambda-1)(\lambda-3) = 0$$

**Autovalores:** $\lambda_1 = 1$ e $\lambda_2 = 3$.

### 10.3 Cálculo dos Autovetores

Para cada $\lambda_i$, resolvemos o sistema homogêneo $(A - \lambda_i I)v = \mathbf{0}$.

**Para $\lambda_1 = 1$:**
$$A - I = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} \xrightarrow{L_2 \leftarrow L_2 - L_1} \begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}$$

$v_1 + v_2 = 0 \implies v_2 = -v_1$. Escolhendo $v_1 = 1$: $\mathbf{v}_1 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$.

**Para $\lambda_2 = 3$:**
$$A - 3I = \begin{pmatrix} -1 & 1 \\ 1 & -1 \end{pmatrix} \xrightarrow{L_2 \leftarrow L_2 + L_1} \begin{pmatrix} -1 & 1 \\ 0 & 0 \end{pmatrix}$$

$-v_1 + v_2 = 0 \implies v_2 = v_1$. Escolhendo $v_1 = 1$: $\mathbf{v}_2 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$.

**Verificação:**
$$A\mathbf{v}_1 = \begin{pmatrix}2&1\\1&2\end{pmatrix}\begin{pmatrix}1\\-1\end{pmatrix} = \begin{pmatrix}1\\-1\end{pmatrix} = 1 \cdot \mathbf{v}_1 \checkmark$$
$$A\mathbf{v}_2 = \begin{pmatrix}2&1\\1&2\end{pmatrix}\begin{pmatrix}1\\1\end{pmatrix} = \begin{pmatrix}3\\3\end{pmatrix} = 3 \cdot \mathbf{v}_2 \checkmark$$

### 10.4 Diagonalização de Matrizes

Uma matriz $n \times n$ é **diagonalizável** se possui $n$ autovetores linearmente independentes. Nesse caso:

$$A = P D P^{-1}$$

onde:
- $P$: matriz cujas **colunas** são os autovetores $[\mathbf{v}_1 \mid \mathbf{v}_2 \mid \cdots \mid \mathbf{v}_n]$
- $D$: matriz **diagonal** com os autovalores correspondentes na diagonal

**Para o exemplo acima:**
$$P = \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}, \quad D = \begin{pmatrix} 1 & 0 \\ 0 & 3 \end{pmatrix}$$

$$\det(P) = 1+1 = 2 \neq 0 \implies P \text{ é inversível} \checkmark$$

**Por que diagonalizar importa?** Calcular $A^k$ (potências de $A$) torna-se trivial:

$$A^k = P D^k P^{-1}, \quad D^k = \begin{pmatrix} \lambda_1^k & 0 \\ 0 & \lambda_2^k \end{pmatrix}$$

Isso é **exponencialmente mais eficiente** computacionalmente.

### 10.5 Autovalores de Matrizes 3×3

O processo é análogo: o polinômio característico será cúbico. Por exemplo, para uma matriz diagonal $D = \text{diag}(d_1, d_2, d_3)$, os autovalores são imediatamente $d_1, d_2, d_3$ (os elementos da diagonal). Para matrizes gerais 3×3, resolve-se a equação cúbica $\det(A - \lambda I) = 0$ numericamente ou analiticamente.

---

## Tópico 11 — Aplicações Reais

### 11.1 Circuitos Elétricos — Leis de Kirchhoff

Em circuitos com múltiplas malhas, as leis de correntes e tensões de Kirchhoff geram um sistema linear. Para um circuito com 3 malhas independentes, obtemos:

$$\begin{pmatrix} R_{11} & -R_{12} & -R_{13} \\ -R_{21} & R_{22} & -R_{23} \\ -R_{31} & -R_{32} & R_{33} \end{pmatrix} \begin{pmatrix} I_1 \\ I_2 \\ I_3 \end{pmatrix} = \begin{pmatrix} V_1 \\ V_2 \\ V_3 \end{pmatrix}$$

onde $R_{ii}$ é a resistência total da malha $i$ e $R_{ij}$ é a resistência compartilhada entre as malhas $i$ e $j$. A solução $I = A^{-1}V$ fornece todas as correntes. Resolvedores modernos de circuitos (SPICE, por exemplo) fazem exatamente isso para sistemas com milhares de nós.

### 11.2 Regressão Linear — Mínimos Quadrados

Dado um conjunto de $m$ observações $(x_i, y_i)$ e querendo ajustar um modelo linear $y = \beta_0 + \beta_1 x$, montamos o sistema **super-determinado** $X\beta = y$, onde:

$$X = \begin{pmatrix} 1 & x_1 \\ 1 & x_2 \\ \vdots & \vdots \\ 1 & x_m \end{pmatrix}$$

Como $m > 2$ (mais equações que incógnitas), o sistema geralmente não tem solução exata. A solução de mínimos quadrados é:

$$\hat{\beta} = (X^T X)^{-1} X^T y$$

Esta é a **Equação Normal**. A expressão $(X^T X)^{-1} X^T$ é chamada pseudo-inversa de Moore-Penrose. Em Python: `np.linalg.lstsq(X, y)` e `sklearn.linear_model.LinearRegression` implementam exatamente esta fórmula por baixo dos panos.

### 11.3 PageRank do Google

O algoritmo original do Google representa a web como uma matriz de adjacência $A$, onde $a_{ij} = 1/d_j$ se existe link da página $j$ para a página $i$ (e $d_j$ é o número de links saindo de $j$).

O vetor de relevância $r$ satisfaz o sistema de autovalores:

$$r = Ar \quad \Longleftrightarrow \quad Ar = 1 \cdot r$$

O PageRank é o **autovetor** de $A$ associado ao autovalor dominante $\lambda = 1$. O Teorema de Perron-Frobenius garante que este autovetor existe, é único e tem entradas positivas (sob certas hipóteses). Na prática, é calculado pelo **Método da Potência**: $r^{(k+1)} = A r^{(k)}$, que converge para o autovetor dominante.

### 11.4 Transformações em Computação Gráfica

Em gráficos 2D e 3D, **toda transformação linear** é uma matriz. Os transformações usam **coordenadas homogêneas** (dimensão extra para incluir translação), permitindo combinar todas as transformações numa única multiplicação matricial:

| Transformação | Matriz 2D (homog.) |
|---|---|
| Rotação por $\theta$ | $\begin{pmatrix}\cos\theta & -\sin\theta & 0\\\sin\theta & \cos\theta & 0\\0&0&1\end{pmatrix}$ |
| Escala $(s_x, s_y)$ | $\begin{pmatrix}s_x & 0 & 0\\0 & s_y & 0\\0&0&1\end{pmatrix}$ |
| Translação $(t_x, t_y)$ | $\begin{pmatrix}1 & 0 & t_x\\0 & 1 & t_y\\0&0&1\end{pmatrix}$ |

Uma pipeline de renderização 3D (como OpenGL ou DirectX) aplica sequencialmente: **Modelo → Visão → Projeção**, cada uma sendo uma multiplicação matricial $4 \times 4$. A composição de todas é uma única matriz $M = M_{proj} \cdot M_{view} \cdot M_{model}$, calculada uma vez e aplicada a milhões de vértices.

### 11.5 Machine Learning — Análise de Componentes Principais (PCA)

O PCA é um dos algoritmos mais importantes de redução de dimensionalidade e usa diretamente a teoria de autovalores/autovetores.

**O problema:** Dado um dataset com $n$ observações e $p$ variáveis (matriz $X$ de tamanho $n \times p$, com $p$ enorme), queremos encontrar uma representação em $k \ll p$ dimensões que preserve o máximo de variância.

**O algoritmo, passo a passo:**

**Passo 1 — Centralização:** Subtrair a média de cada variável: $\tilde{X} = X - \bar{X}$.

**Passo 2 — Matriz de covariância:** Calcular

$$C = \frac{1}{n-1}\tilde{X}^T \tilde{X}$$

$C$ é uma matriz $p \times p$, **simétrica** e positiva semi-definida.

**Passo 3 — Decomposição espectral:** Calcular os autovalores $\lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_p \geq 0$ e autovetores $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_p$ de $C$. Os autovetores são as **Componentes Principais**.

**Passo 4 — Projeção:** Selecionar os $k$ autovetores com maiores autovalores e construir $W = [\mathbf{v}_1 \mid \cdots \mid \mathbf{v}_k]$. A representação reduzida é:

$$Z = \tilde{X} W \quad \text{(matriz } n \times k\text{)}$$

**Por que funciona?** O autovalor $\lambda_i$ representa a **variância** dos dados na direção $\mathbf{v}_i$. Manter as direções de maior variância é equivalente a manter a informação mais relevante do dataset. A "variância explicada" pela $i$-ésima componente é $\lambda_i / \sum_j \lambda_j$.

**Em Python:**
```python
from sklearn.decomposition import PCA
import numpy as np

# X: matriz de dados (n_samples, n_features)
pca = PCA(n_components=2)
Z = pca.fit_transform(X)
# Z: representação em 2D
# pca.explained_variance_ratio_: % de variância por componente
```

Internamente, `sklearn` usa a **Decomposição em Valores Singulares (SVD)**: $\tilde{X} = U \Sigma V^T$, que é numericamente mais estável do que calcular $C$ explicitamente. Os autovetores de $C$ correspondem às colunas de $V$.

---

## Resumo Visual das Conexões

```
Matrizes & Operações
        │
        ▼
   Produto de Matrizes ──────────────────────────────► Transformações Gráficas
        │
        ▼
   Matrizes Inversas ◄──── Matrizes Elementares ◄──── Eliminação de Gauss
        │                                                      │
        ▼                                                      ▼
   Determinantes ────────────────────────────────► Sistemas Lineares (Ax=b)
        │                                                      │
        ▼                                                      ▼
Autovalores & Autovetores ──────────────────────► PageRank, PCA, Circuitos
        │
        ▼
   Diagonalização ───────────────────────────────► Potências de Matrizes
                                                   Equações Diferenciais
```

---

*Guia elaborado com rigor matemático e foco em aplicações práticas. Use-o como referência durante seus estudos e não hesite em revisitar os exemplos numéricos até que os cálculos se tornem naturais.*

# redes_complexas

Esse repositório foi criado para organizar meus programas da iniciação científica, onde estudo redes complexas.
Pretendo dividir melhor meu código ao longo do tempo.


- v0.1
  - Medidas básicas e medida de posição relativa funcionando para redes com arestas sem peso e sem direção.

- v0.2
  - Adicionado programa para trabalhar com redes espaciais de Waxman.
  - Medida de posição relativa alterada (agora pode-se selecionar o número de níveis a serem adicionados a matriz, antes sempre adicionava-se todos os níveis hierárquicos de cada vértice)

- v1.0 (29/04/2019)
  - As rotinas estão separadas em 3 programas, um para cada modelo de rede estudado: 
    - waxman.py (para modelos baseados em distância euclideana);
    - erdos_renyi.py (para modelos de grafos aleatórios);
    - watts_strogatz.py (para o estudo de redes no modelo de pequeno mundo de Watts-Strogatz e da interpolação entre um reticulado e um grafo totalmente aleatório).
  - Foi adicionada a integração parcial dos códigos com o pacote python-igraph para facilitar a obtenção de resultados rápidos como a distribuição de graus da rede, ou o diâmetro, etc.

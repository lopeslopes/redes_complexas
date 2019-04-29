# redes_complexas

Esse repositorio foi criado para organizar meus programas da iniciacao cientifica, onde estudo redes complexas.
Pretendo dividir melhor meu codigo ao longo do tempo.


- v0.1
  - Medidas basicas e medida de posicao relativa funcionando para redes com arestas sem peso e sem direcao.

- v0.2
  - Adicionado programa para trabalhar com redes espaciais de Waxman. 
  - Medida de posicao relativa alterada (agora pode-se selecionar o numero de niveis a serem adicionados a matriz, antes sempre adicionava-se todos os niveis hierarquicos de cada vertice)

- v1.0 (29/04/2019)
  - As rotinas estao separadas em 3 programas, um para cada modelo de rede estudado: 
    - waxman.py (para modelos baseados em distancia euclideana);
    - erdos_renyi.py (para modelos de grafos aleatorios);
    - ws.py (para o estudo de redes no modelo de pequeno mundo de Watts-Strogatz e da interpolacao entre um reticulado e um grafo totalmente aleatorio).
  - Foi adicionada a integracao parcial dos codigos com o pacote python-igraph para facilitar a obtencao de resultados rapidos como a distribuicao de graus da ede, ou o diametro, etc.

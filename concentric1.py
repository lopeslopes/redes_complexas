import random
import numpy as np
import igraph as ig
import math
import scipy
from multiprocessing.dummy import Pool as ThreadPool 

random.seed()

class node:
    def __init__(self, index):
        self.index = index
        self.neighbors = []

    def __repr__(self):
        return repr(self.index)

def ergraph(n, p):        ''' cria uma rede no modelo de Erdos-Renyi '''
    vertices = [node(i) for i in range(n)]
    edges = [(i, j) for i in range(n) for j in range(i) if random.random() < p]

    for (i, j) in edges:
        vertices[i].neighbors.append(vertices[j])
        vertices[j].neighbors.append(vertices[i])

    return vertices

def maior_componente(vertices):   ''' retorna o maior componente conectado da rede '''
    visitados = set()
    pilha = []
    grafos = []
    i = -1

    for nodo in vertices:
        if nodo not in visitados:
            i = i+1
            grafos.append([])
            pilha = [nodo]

        while len(pilha)>0:
            v = pilha.pop()
            if v not in visitados:
                grafos[i].append(v)
                visitados.add(v)
                for vizinho in v.neighbors:
                    if vizinho not in visitados:
                        pilha.append(vizinho)

    aux = max(grafos, key=lambda x: len(x))
    for i in range(len(aux)):
        aux[i].index = i
    
    return (aux)

def matriz_adj(vertices, n):  ''' retorna a matriz de adj. de uma rede criada pelo metodo acima '''
    matriz = []
    for i in range(n):
        matriz.append([])
        for j in range(n):
            matriz[i].append(False)

    for v in vertices:
        for vizinho in v.neighbors:
            matriz[v.index][vizinho.index] = True
            matriz[vizinho.index][v.index] = True

    return matriz

def bola_d(no, d):
    bola = set()
    anel = set()
    bola.add(no)
    for m in range(d):
        for v in bola:
            for vizinho in v.neighbors:
                if (vizinho not in bola):
                    anel.add(vizinho)
        
        bola.update(anel)
    
    return bola

def anel_d(no, d):
    anel = set()
    if (d==0):
        anel.add(no)
    else:
        bola = bola_d(no, d-1)
        for v in bola:
            for vizinho in v.neighbors:
                if (vizinho not in bola):
                    anel.add(vizinho)

    return anel

def k_d(no, d):
    bola = bola_d(no, d)
    anel = anel_d(no, d+1)
    grau = 0

    for b in bola:
        for a in anel:
            if (a in b.neighbors):
                grau += 1

    return grau

def distribuicao(vertices):
    graumax = 0
    n = len(vertices)
    for k in range(n):
        if graumax < len(vertices[k].neighbors):
            graumax = len(vertices[k].neighbors)

    dist_grau = [0 for k in range(graumax)]

    for k in range(graumax):
        for l in range(n):
            if k == len(vertices[l].neighbors):
                dist_grau[k]+=1

    f = open('dados', 'w')
    for k in range(graumax):
        f.write(str(dist_grau[k])+"\n")

    f.close()

def cc_d(no, d):
    if (d==0):
        coef = 0
    else:
        anel = anel_d(no, d)
        
        num_arestas = 0
        for a in anel:
            for v in a.neighbors:
                if v in anel:
                    num_arestas = num_arestas + 1
        
        num_vert = len(anel)
        coef = (2*num_arestas)/(num_vert*(num_vert-1))

    return coef

def convergencia(no, d):
    anel = anel_d(no, d+1)
    num_vert = len(anel)
    kd = k_d(no, d)
    taxa = kd/num_vert
    return taxa

def intraring(no, d):
    anel = anel_d(no, d)
    media = 0
    for a in anel:
        for vizinho in a.neighbors:
            if (vizinho in anel):
                media = media + 1
    if (len(anel) != 0):
        media = media/len(anel)
    else:
        media = 0
    return media

def interring(no, d):
    anel = anel_d(no, d)
    anel1 = anel_d(no, d+1)
    media = 0
    for a in anel:
        for vizinho in a.neighbors:
            if (vizinho in anel1):
                media = media + 1
    if (len(anel) != 0):
        media = media/len(anel)
    else:
        media = 0
    return media

def h_common(no, d): 
    anel = anel_d(no, d)
    media = 0
    for a in anel:
        media = media + len(a.neighbors)
    media = media/len(anel)
    return media

def pos_relativa(no, n):  ''' medida de posicao relativa para redes sem peso e sem direcao '''
    d = 0
    pos_ind = []

    for i in range (n):
        pos_ind.append(-1)
    
    anel = anel_d(no, d)
    
    while(len(anel) != 0):
        for a in anel:
            pos_ind[a.index] = d
        d = d + 1
        anel.clear()
        anel = anel_d(no, d)

    return pos_ind

def pca(matriz):          ''' faz o pca de uma matriz '''
    m_np = np.array(matriz)
    mt_np = np.transpose(m_np)

    ''' subtracao da media '''
    
    media = []

    for linha in mt_np:
        soma = 0
        for i in range(len(linha)):
            soma = soma + linha[i]
        soma = soma/i
        media.append(soma)

    xc = []
    i = 0

    for linha in mt_np:
        xc.append([])
        for valor in linha:
            xc[i].append(valor-media[i])
        i = i+1

    ''' divisao pelo desvio padrao '''

    sd = []

    for linha in xc:
        desvio = 0
        for valor in linha:
            desvio = desvio + valor*valor
        desvio = math.sqrt(desvio/(len(linha)-1))
        sd.append(desvio)

    for i in range(len(xc)):
        for j in range(len(xc[0])):
            xc[i][j] = xc[i][j]/sd[i]

    ''' fim '''
    
    xc_np = np.array(xc)
    xct_np = np.transpose(xc_np)

    k = (1/(len(matriz)-1))*(xc_np @ xct_np)

    valores, vetores = np.linalg.eig(k)
    
    indices = np.argsort(valores)
    ind = indices.tolist()
    ind.reverse()

    ordenados = []

    s = 0
    sc = 0

    for i in range(len(ind)):
        if (valores[ind[i]] > 0):
            ordenados.append(vetores[ind[i]])
            s = s + valores[ind[i]]
    
    for i in range(2):
        sc = sc + valores[ind[i]]
    conservacao = (sc/s)*100

    print("conservacao da variancia: " + str(conservacao))

    ord_np = np.array(ordenados)
    ord_np = np.transpose(ord_np)

    pca = ord_np @ mt_np

    ''' subtracao da media do pca '''

    media = []

    for linha in pca:
        soma = 0
        for valor in linha:
            soma = soma + valor
        soma = soma/len(linha)
        media.append(soma)

    for i in range(len(pca)):
        for j in range(len(pca[0])):
            pca[i][j] = (pca[i][j]-media[i])

    ''' divisao pelo desivo padrao do pca '''

    sd = []

    for linha in pca:
        desvio = 0
        for valor in linha:
            desvio = desvio + valor*valor
        desvio = math.sqrt(desvio/(len(linha)-1))
        sd.append(desvio)

    for i in range(len(pca)):
        for j in range(len(pca[0])):
            pca[i][j] = pca[i][j]/sd[i]

    ''' fim '''

    return pca

def pos_relativa_w(nos):  ''' medida de posicao relativa p/ redes com peso'''
    cont = 0
    pos_rel = []
    for vertice in nos.vs:
        pos_rel.append([])
        for i in range(len(nos.vs)):
            distancia = 0
            aux = nos.get_shortest_paths(vertice, to=i, weights=nos.es["weight"], mode="ALL", output="epath")
            for e in aux[0]:
                aux2 = nos.es[e]["weight"]
                distancia = distancia + aux2
            pos_rel[cont].append(distancia)
        cont = cont + 1

    return pos_rel 
    
def iris():               ''' faz o pca da base dados iris e divide em 3 arquivos para plotar com cores diferentes '''
    f = open('iris.dat', 'r')

    data = f.readline()
    data = f.readline()
    setosa = []
    versicolor = []
    virginica = []

    while(data != ''):
        aux2 = []
        aux = data.split(',')
        for i in range(4):
            aux2.append(float(aux[i]))

        if(aux[4]=='setosa\n'):
            setosa.append(aux2)
        elif(aux[4]=='versicolor\n'):
            versicolor.append(aux2)
        elif(aux[4]=='virginica\n'):
            virginica.append(aux2)
        data = f.readline()

    f.close()

    total = []
    total.extend(setosa)
    total.extend(versicolor)
    total.extend(virginica)

    total2 = []
    i = 0

    for t in total:
        total2.append([])
        for medida in t:
            total2[i].append(float(medida))
        i = i+1

    total_np = np.array(total2)

    pca2 = pca(total_np)
    teste1 = np.transpose(pca2)

    f = open('setosa.dat', 'w')
    g = open('versicolor.dat', 'w')
    h = open('virginica.dat', 'w')

    for i in range(len(teste1)):
        if(total2[i] in setosa):
            f.write(str(teste1[i][0]) + " " + str(teste1[i][1]) + "\n")
        elif(total2[i] in versicolor):
            g.write(str(teste1[i][0]) + " " + str(teste1[i][1]) + "\n")
        elif(total2[i] in virginica):
            h.write(str(teste1[i][0]) + " " + str(teste1[i][1]) + "\n")

    f.close()



''' main (geralmente pra escolher uma medida e gravar em arquivo para plotar)'''

g = open('sao_carlos.gml', 'r')
nos = ig.Graph.Read_GML(g)
g.close()

arestas_re = []
for i in range(1500):
    arestas_re.append(nos.es[i])

nos2 = nos.subgraph_edges(arestas_re, delete_vertices=True)

resultado = pca(pos_relativa_w(nos2))

f = open('pca.dat', 'w')

for linha in resultado:
    f.write(str(linha[0]) + " " + str(linha[1]) + "\n")


f.close()

''' fim main '''
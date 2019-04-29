import random
import numpy as np
import math

random.seed()

''' definicoes '''

class node:
    def __init__(self, index):
        self.index = index
        self.neighbors = []

    def __repr__(self):
        return repr(self.index)

def ergraph(n, p):
    vertices = [node(i) for i in range(n)]
    edges = [(i, j) for i in range(n) for j in range(i) if random.random() < p]

    for (i, j) in edges:
        vertices[i].neighbors.append(vertices[j])
        vertices[j].neighbors.append(vertices[i])

    return vertices

''' funcoes praticas '''

def rede_mathematica(grafo):
    r1 = open('matriz_adj.dat', 'w')

    m = matriz_adj(grafo, len(grafo))
    for linha in m:
        for l in linha:
            if(l):
                r1.write(str(1)+ ' ')
            else:
                r1.write(str(0)+ ' ')
        r1.write('\n')

    r1.close()

def arquivo(matriz, nome):
    arq = open(nome, "w")

    if(type(matriz)==type(matriz[0])):
        for i in range(len(matriz)):
            for j in range(len(matriz[0])-1):
                arq.write(str(matriz[i][j]) + " ")
            j =+ 1
            arq.write(str(matriz[i][j]) + "\n")
    else:
        for i in range(len(matriz)):
            arq.write(str(matriz[i]) + "\n")

    arq.close()

''' funcoes de calculo '''

def maior_componente(vertices):
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
    
    return aux

def matriz_adj(vertices, n):
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

def pca(matriz, flag):
    m_np = np.array(matriz)
    mt_np = np.transpose(m_np)

    ''' standardizacao dos dados '''
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
    
    for i in range(3):
        sc = sc + valores[ind[i]]
    conservacao = (sc/s)*100

    if (len(ordenados) < len(ind)):
        dif = len(ind) - len(ordenados)
        for i in range(dif):
            ordenados.append([0 for j in range(len(ind))])
    print("conservacao da variancia: " + str(conservacao))
    
    ord_np = np.array(ordenados)
    ord_np = np.transpose(ord_np)

    pca = ord_np @ mt_np

    if (flag):
        ''' standardizacao do resultado do pca '''
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

        sd = []

        for linha in pca:
            desvio = 0
            for valor in linha:
                desvio = desvio + valor*valor
            desvio2 = (desvio/(len(linha)-1))**(.5)
            sd.append(desvio2)

        for i in range(len(pca)):
            for j in range(len(pca[0])):
                pca[i][j] = pca[i][j]/sd[i]

    return pca

def m_distancias(grafo, num):
    matriz = []
    aux = set()

    for v1 in grafo:
        aux.clear()
        d = 0
        aux.add(v1)
        matriz.append([])
        matriz[v1.index] = [0 for i in range(len(grafo))]
        while(d <= num):
            d = d+1
            anel = anel_d(v1, d)
            for a in anel:
                matriz[v1.index][a.index] = d
            aux.update(anel)

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

def medida(grafo, num):
    mdist = m_distancias(grafo, num)
    m_vizinhos = []
    n_viz = []

    for i in range(len(mdist)):
        m_vizinhos.append([])
        num = 0
        for j in range(len(mdist[0])):
            if (mdist[i][j] != 0):
                m_vizinhos[i].append(j)
                num += 1
        n_viz.append(num)

    maximo = max(n_viz)

    for i in range(len(m_vizinhos)):
        if(len(m_vizinhos[i]) < maximo):
            for j in range(len(m_vizinhos[i]), maximo):
                m_vizinhos[i].append(0)

    m_vizinhos = np.array(m_vizinhos)

    return m_vizinhos

''' main '''

grafo0 = ergraph(1024, 0.01)
grafo = maior_componente(grafo0)

rede_mathematica(grafo)

m_viz = medida(grafo, 5)

dados = pca(m_viz, True)
np.transpose(dados)
dados2d = [[dados[i][0], dados[i][1]] for i in range(len(dados))]


arquivo(dados2d, "pca_er.dat")

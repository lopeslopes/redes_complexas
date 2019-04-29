import random
import numpy as np
import math
import igraph as ig
random.seed()


''' estrutura da rede '''

class node:
    def __init__(self, index, xpos, ypos):
        self.index = index
        self.neighbors = []
        self.xpos = xpos
        self.ypos = ypos

    def __repr__(self):
        return repr(self.index)

def waxman_graph(n, densidade_l, dist_tipica):
    vertices = []
    arestas = []
    for i in range(n):
        xpos = random.random()
        ypos = random.random()
        vertices.append(node(i, xpos, ypos))
    
    for i in range(n):
        for j in range(n):
            d = distancia(vertices[i], vertices[j])
            if (random.random() < (densidade_l * math.exp(-(d/dist_tipica)))) and (i!=j):
                arestas.append((i, j, d))
    
    for (i, j, d) in arestas:
        if vertices[j] not in vertices[i].neighbors:
            vertices[i].neighbors.append(vertices[j])
        if vertices[i] not in vertices[j].neighbors:
            vertices[j].neighbors.append(vertices[i])

    grafo = [vertices, arestas]
    return grafo

def maior_componente(grafo):
    visitados = set()
    grafos = []
    k = -1
    vertices = grafo[0]

    for i in range(len(vertices)):
        if vertices[i] not in visitados:
            k = k+1
            grafos.append([])
            pilha = [vertices[i]]

        while len(pilha)>0:
            v = pilha.pop()
            if v not in visitados:
                grafos[k].append(v)
                visitados.add(v)
                for vizinho in v.neighbors:
                    if vizinho not in visitados:
                        pilha.append(vizinho)

    aux_v = max(grafos, key=lambda x: len(x))

    for i in range(len(aux_v)):
        aux_v[i].index = i

    aux_as = set()
    aux_a = []
    
    for v1 in aux_v:
        for v2 in aux_v:
            if v2 in v1.neighbors:
                d = distancia(v1, v2)
                if ((v2.index, v1.index, d) not in aux_as):
                    aux_as.add((v1.index, v2.index, distancia(v1, v2)))
    
    for a in aux_as:
        aux_a.append(a)
    
    aux_g = [aux_v, aux_a]
    return aux_g

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

def matriz_adj(vertices, n):
    matriz = []
    for i in range(n):
        matriz.append([])
        for j in range(n):
            matriz[i].append(0)

    for v in vertices:
        for vizinho in v.neighbors:
            matriz[v.index][vizinho.index] = 1
            matriz[vizinho.index][v.index] = 1

    return matriz

''' etapas da medida '''

def m_distancias(grafo, num):
    matriz = []
    aux = set()

    for v1 in grafo[0]:
        aux.clear()
        d = 0
        aux.add(v1)
        matriz.append([])
        matriz[v1.index] = [0 for i in range(len(grafo[0]))]
        while(d < num):
            d = d+1
            anel = anel_d(v1, d)
            for a in anel:
                matriz[v1.index][a.index] = d
            aux.update(anel)

    return matriz

def medida_v1(grafo, num):
    '''essa medida ignora o grau hierarquico dos nos, armazena
    na ordem dos indices '''
    mdist = m_distancias(grafo, num)
    m_vizinhos = []
    n_viz = []

    for i in range(len(mdist)):
        m_vizinhos.append([])
        aux = 0
        for j in range(len(mdist[0])):
            if (mdist[i][j] != 0):
                m_vizinhos[i].append(j)
                aux += 1
        n_viz.append(aux)

    maximo = max(n_viz)

    for i in range(len(m_vizinhos)):
        if(len(m_vizinhos[i]) < maximo):
            for j in range(len(m_vizinhos[i]), maximo):
                m_vizinhos[i].append(0)

    m_vizinhos = np.array(m_vizinhos)
    return m_vizinhos

def medida_v2(grafo, num):
    '''essa medida leva em conta o grau hierarquico dos nos, os vizinhos
    mais proximos aparecem primeiro '''
    mdist = m_distancias(grafo, num)
    m_vizinhos = []
    n_viz = []

    for i in range(len(mdist)):
        m_vizinhos.append([])
        aux = 0
        for k in range(1, num):
            for j in range(len(mdist[0])):
                if (mdist[i][j] == k):
                    m_vizinhos[i].append(j)
                    aux += 1
        n_viz.append(aux)

    maximo = max(n_viz)

    for i in range(len(m_vizinhos)):
        if(len(m_vizinhos[i]) < maximo):
            for j in range(len(m_vizinhos[i]), maximo):
                m_vizinhos[i].append(0)

    m_vizinhos = np.array(m_vizinhos)
    return m_vizinhos

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

''' entrada e saida de redes '''

def gerar_gml(grafo, nome):
    arq = open(str(nome)+'.gml', 'w')
    arq.write('graph\n[\n    directed 0\n')
    for v in grafo[0]:
        arq.write('  node\n  [\n    id '+str(v.index)+'\n    waschain 0\n    posx '+str(v.xpos)+'\n    posy '+str(v.ypos)+'\n  ]\n')
    for e in grafo[1]:
        arq.write('  edge\n  [\n    source '+str(e[0])+'\n    target '+str(e[1])+'\n    roadtype null\n    dist '+str(e[2])+'\n  ]\n') 
    arq.write(']')
    arq.close()

def ler_gml(arquivo):
    tudo = arquivo.readlines()
    vertices = []
    arestas = []

    for i in range(len(tudo)):
        if(tudo[i]=='  node\n'):
            i = i+2
            aux1 = tudo[i][7:-1]
            i = i+2
            aux2 = tudo[i][9:-1]
            i = i+1
            aux3 = tudo[i][9:-1]
            v = node(index=int(aux1), xpos=float(aux2), ypos=float(aux3))
            vertices.append(v)
        elif(tudo[i]=='  edge\n'):
            i = i+2
            aux1 = tudo[i][11:-1]
            i = i+1
            aux2 = tudo[i][11:-1]
            i = i+2
            aux3 = tudo[i][11:-1]
            arestas.append((int(aux1), int(aux2), float(aux3)))

    for (i, j, d) in arestas:
        if vertices[j] not in vertices[i].neighbors:
            vertices[i].neighbors.append(vertices[j])
        if vertices[i] not in vertices[j].neighbors:
            vertices[j].neighbors.append(vertices[i])
    
    grafo = [vertices, arestas]
    return grafo

def rep_spline(pca2d, c_spline, passo):
    representacao = []
    for i in range(len(pca2d)):
        dists = []
        for j in range(len(c_spline)):
            d = (c_spline[j][0]-pca2d[i][0])**2 + (c_spline[j][1]-pca2d[i][1])**2
            dists.append(d)
        aux = dists.index(min(dists))
        aux = aux*passo
        min_dist = min(dists)
        representacao.append([i, aux, min_dist])
    return representacao

''' utilidades '''

def distancia(v1, v2):
    dist = ((v1.xpos - v2.xpos)**2 + (v1.ypos - v2.ypos)**2)**0.5
    return dist

def arquivo(matriz, nome):
    arq = open(nome, "w")

    if(type(matriz)==type(matriz)):
        for i in range(len(matriz)):
            for j in range(len(matriz[0])-1):
                arq.write(str(matriz[i][j]) + " ")
            j =+ 1
            arq.write(str(matriz[i][j]) + "\n")
    else:
        for i in range(len(matriz)):
            arq.write(str(matriz[i]) + "\n")

    arq.close()

def rede_mathematica(grafo, nome):
    r1 = open('adj'+str(nome)+'.dat', 'w')

    m = matriz_adj(grafo[0], len(grafo[0]))
    for linha in m:
        for l in linha:
            r1.write(str(l)+ ' ')
        r1.write('\n')

    r1.close()

    r2 = open('coord'+str(nome)+'.dat', 'w')

    for i in range(len(m)):
        r2.write(str(grafo[0][i].xpos)+' '+str(grafo[0][i].ypos)+'\n')

    r2.close()

''' main '''


''' abrir um grafo '''
f = open('0,04.gml', 'r')
grafo = ler_gml(f)
f.close()
'''
'''

aux = matriz_adj(grafo[0], len(grafo[0]))
g = ig.GraphBase.Adjacency(aux, mode="undirected")
print(len(g.get_diameter()))


''' tirar a medida e aplicar o pca '''
'''
m_vizinhos = medida_v1(grafo, 20)
dados = pca(m_vizinhos, True)
dados = np.transpose(dados)


dados3d = [[dados[i][j].real for j in range(3)] for i in range(len(dados))]
arquivo(dados3d, "pca_3d.dat")
'''



''' abrir um pca ja existente aqui no programa '''
'''
f = open('pca_saocarlos-10.dat', 'r')
tudo = f.readlines()
f.close()
pca2d = []
for ponto in tudo:
    ponto = ponto[0:-2]
    aux = ponto.split(' ')
    pca2d.append([float(aux[0]), float(aux[1])])
'''


''' selecionar pontos de controle e obter curva do spline no mathematica '''
''' depois usar curva do spline para representar o pca em uma dimensao '''
'''
f = open('spline.dat', 'r')
tudo = f.readlines()
f.close()
curva = []
for ponto in tudo:
    aux = ponto.split(',')
    curva.append([float(aux[0]), float(aux[1])])

passo = 0.0001
rep = rep_spline(pca2d, curva, passo)
'''



''' fatiar a representacao do pca e depois colorir a rede '''
'''
fatias = [[] for i in range(3)]
for ponto in rep:
    if(ponto[1] > 0.00 and ponto[1] < 0.33):
        fatias[0].append(ponto[0])
    elif(ponto[1] > 0.33 and ponto[1] < 0.66):
        fatias[1].append(ponto[0])
    elif(ponto[1] > 0.66 and ponto[1] < 1.00):
        fatias[2].append(ponto[0])

rede_f = [[] for i in range(len(fatias))]
for vertice in grafo[0]:
    for i in range(len(fatias)):
        if(vertice.index in fatias[i]):
            rede_f[i].append([vertice.xpos, vertice.ypos])

for i in range(len(fatias)):
    arquivo(rede_f[i], 'rede_teste_parte'+str(i)+'.dat')

teste = [[],[],[]]
for ponto in rep:
    if(ponto[0] in fatias[0]):
        teste[0].append([ponto[1], ponto[2]])
    elif(ponto[0] in fatias[1]):
        teste[1].append([ponto[1], ponto[2]])
    elif(ponto[0] in fatias[2]):
        teste[2].append([ponto[1], ponto[2]])

for i in range(len(teste)):
    arquivo(teste[i], 'pca_teste'+str(i)+'.dat')
'''
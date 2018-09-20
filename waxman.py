import random
import numpy as np
import math

random.seed()


class node:
    def __init__(self, index, xpos, ypos):
        self.index = index
        self.neighbors = []
        self.xpos = xpos
        self.ypos = ypos

    def __repr__(self):
        return repr(self.index)

def distancia(v1, v2):
    dist = ((v1.xpos - v2.xpos)**2 + (v1.ypos - v2.ypos)**2)**0.5
    return dist

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

def pca(matriz, flag):
    m_np = np.array(matriz)
    mt_np = np.transpose(m_np)

    if(flag):
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
    else:
        xc = matriz
    
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

    if (len(ordenados) < len(ind)):
        dif = len(ind) - len(ordenados)
        for i in range(dif):
            ordenados.append([0 for j in range(len(ind))])
    print("conservacao da variancia: " + str(conservacao))
    
    ord_np = np.array(ordenados)
    ord_np = np.transpose(ord_np)

    pca = ord_np @ mt_np

    ''' subtracao da media do pca

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

     divisao pelo desivo padrao do pca 

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

     fim '''

    return pca

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

def m_distancias(grafo, num):
    matriz = []
    aux = set()

    for v1 in grafo[0]:
        aux.clear()
        d = 0
        aux.add(v1)
        matriz.append([])
        matriz[v1.index] = [0 for i in range(len(grafo[0]))]
        while(d <= num):
            d = d+1
            anel = anel_d(v1, d)
            for a in anel:
                matriz[v1.index][a.index] = d
            aux.update(anel)

    return matriz

def maior_distancia(grafo):
    aux = set()
    dmax = 0

    for v1 in grafo[0]:
        aux.clear()
        d = 0
        aux.add(v1)
        while(len(aux) != len(grafo[0])):
            d = d+1
            anel = anel_d(v1, d)
            aux.update(anel)
        if (d>=dmax):
            dmax = d

    return dmax

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

def gerar_gml(grafo):
    arq = open('waxman.gml', 'w')
    arq.write('graph\n[\n    directed 0\n')
    for v in grafo[0]:
        arq.write('  node\n  [\n    id '+str(v.index)+'\n    posx '+str(v.xpos)+'\n    posy '+str(v.xpos)+'\n  ]\n')
    for e in grafo[1]:
        arq.write('  edge\n  [\n    source '+str(e[0])+'\n    target '+str(e[1])+'\n    dist '+str(e[2])+'\n  ]\n') 
    arq.write(']')
    arq.close()

def rede_mathematica(grafo):
    r1 = open('rede1.dat', 'w')

    m = matriz_adj(grafo[0], len(grafo[0]))
    for linha in m:
        for l in linha:
            r1.write(str(l)+ ' ')
        r1.write('\n')

    r1.close()

    r2 = open('rede2.dat', 'w')

    for i in range(len(m)):
        r2.write(str(grafo[0][i].xpos)+' '+str(grafo[0][i].ypos)+'\n')

    r2.close()

def d_borda(pca1):
    
    n_viz = []
    i = 0
    for linha1 in pca1:
        n_viz.append(0)
        for linha2 in pca1:
            if ((linha2[0]<(linha1[0]+0.08)) and linha2[0]>(linha1[0]-0.08) and (linha2[1]<linha1[1]+0.08) and (linha2[1]>linha1[1]-0.08)):
                n_viz[i] = n_viz[i] + 1
        i = i+1

    borda = set()
    for i in range(len(n_viz)):
        if (n_viz[i]<=2):
            borda.add(i)
    
    '''borda = set()
    for i in range(len(pca1)):
        if((pca1[i][0]**2 + pca1[i][1]**2) > 9):
            borda.add(i)'''

    return borda

def d_borda_rede(grafo):
    borda_r = set()
    naoborda = set()
    for v in grafo[0]:
        if (v.xpos<0.1):
            borda_r.add(v.index)
        elif (v.xpos>0.9):
            borda_r.add(v.index)
        elif (v.ypos<0.1):
            borda_r.add(v.index)
        elif (v.ypos>0.9):
            borda_r.add(v.index)
        else:
            naoborda.add(v.index)
    
    return borda_r


grafo1 = waxman_graph(500, 2, 0.05)
grafo = maior_componente(grafo1)

borda = d_borda_rede(grafo)
b = open('borda.dat', 'w')
for i in range(len(borda)-1):
    b.write(str(borda.pop()) + '|')
b.write(str(borda.pop()))
b.close()


'''f = open('sao_carlos.gml', 'r')
grafo = ler_gml(f)'''


rede_mathematica(grafo)

print(maior_distancia(grafo))

teste = 1
while(teste == 1):
    u = input("digite o numero de niveis na medida: ")
    m = m_distancias(grafo, int(u))

    dados = pca(m, True)

    borda2 = d_borda_rede(grafo)

    f1 = open('pca1.dat', 'w')
    f2 = open('pca2.dat', 'w')

    for i in range(len(grafo[0])):
        if(i in borda2):
            f2.write(str(dados[i][0]) + ' ' + str(dados[i][1]) + '\n')
        else:    
            f1.write(str(dados[i][0]) + ' ' + str(dados[i][1]) + '\n')

    f1.close()
    f2.close()

    teste = int(input("0 para sair, 1 para outra medida: "))
import random
import numpy as np
import math
import igraph as ig 


def vizinhanca_d(indice, matriz_adj, alcance):
    aux1 = set()
    aux2 = set()
    aux3 = set()
    aux1.add(indice)
    aux3.add(indice)

    vizinhanca = [0 for i in range(len(matriz_adj))]

    for i in range(1,alcance):
        for v in aux1:
            for j in range(len(matriz_adj[v])):
                if(matriz_adj[v][j] == 1):
                    if(j not in aux3):
                        aux2.add(j)
                        aux3.add(j)
                        vizinhanca[j] = i
        aux1.clear()
        aux1 = aux2.copy()
        aux2.clear()

    return vizinhanca

def m_distancias(matriz_adj, alcance):
    matriz = []
    for i in range(len(matriz_adj)):
        matriz.append(vizinhanca_d(i, matriz_adj, alcance))

    return matriz

def medida_v1(m_dist, alcance):
    '''essa medida ignora o grau hierarquico dos nos, armazena
    na ordem dos indices '''
    m_vizinhos = []
    n_viz = []

    for i in range(len(m_dist)):
        m_vizinhos.append([])
        aux = 0
        for j in range(len(m_dist[0])):
            if (m_dist[i][j] != 0):
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

''' main '''

g = ig.Graph.Read_GML(open("lattice.gml", "r"))
aux = g.degree_distribution(bin_width = 1)
print(aux)

for k in range(4):
    for i in range(1,10):
        m_adj = g.get_adjacency()
        m_adj = m_adj.data
        dados_v1 = pca(np.array(medida_v1(m_distancias(m_adj, 5), 5)), True)
        dados = np.transpose(dados_v1)
        print(dados.shape)
        dados2d = [[dados[j][0].real, dados[j][1].real] for j in range(len(dados))]
        arquivo(dados2d, "pca_ws_"+str(i)+"_"+str(k)+".dat")
        g.rewire_edges((i*10**(k-4)), loops=False, multiple=False)

aux = g.degree_distribution(bin_width = 1)
print(aux)
g.write_gml(open("post_rewire.gml", "w"))
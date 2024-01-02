library(igraph)

#Leemos TFs
tfs.data = read.table("targets.txt", as.is=T)
tfs = tfs.data$V1
tfs
length(tfs)

#Leemos los genes que se regulan por cada FT

genes.in.net=c()

#Se quedarian fuera TFs que no se regulan, deberíamos añadirlo externamente, para tener nodo fuente, en este caso ya 
#lo tenemos

for(i in 1:length(tfs))
{
  cur.tf.file=paste0(tfs[i], ".txt", sep="")
  cur.tf.data=read.table(cur.tf.file, as.is=T)
  cur.tf.targets= cur.tf.data$V1
  genes.in.net= c(genes.in.net, cur.tf.targets)
}

length(genes.in.net)
length(unique(genes.in.net))
head(genes.in.net)
#Es un mero vector con nuestro universo sin nada mas, debemos contruir la red
#partiriamos de matriz de adyacencia, matriz cuadrada n x n, con n igual a numero
#de genes que tengamos, si existe relacion entre genes 1 sino 0, en coexpresion
#en este caso es matriz dirigida, no es lo mismo hay que establecer orden, rellenamos
#en filas si gen i regula a j ponemos 1 no necesariamente indicara que j regula a i

adj.matrix=matrix(0, nrow = length(unique(genes.in.net)), ncol = length(unique(genes.in.net)))
dim(adj.matrix)

rownames(adj.matrix)=unique(genes.in.net)
colnames(adj.matrix)=unique(genes.in.net)


#Cuando aparezca un gen en la lista del TF dado significa que es regulado por este, fila del tf col del gen
#Podemos fijar cada TF y rellenar en esa fila las columna dada, para fijar los TF como lo hacemos, poes con
#un for iterando sobre lista de tfs

sum (tfs %in% genes.in.net)
#por eso no es necesario añadirlos, lostenemos, comprobamos que nuestros tfs estan en las lisyas de regulomas
#sino habria que añadirlos manualmente a genes.in.net

i=1
for(i in 1:length(tfs))
{
  cur.tf=tfs[i]
  cur.tf.file=read.table(paste0(cur.tf, ".txt"))
  cur.tf.vec=cur.tf.file$V1
  adj.matrix[cur.tf, cur.tf.vec]=1
  
}

sum(adj.matrix)
#2 1 o mas en una misma fila es que ese Tf regula esos genes, si es en una misma columna es que mismo gen se
#regula por distintos tfs.

#CREAR LA RED A PARTIR DE LA MATRIZ DE ADYACENCIA. 

library(igraph)

transcriptional.network=graph.adjacency(adj.matrix, mode= "directed")
transcriptional.network
class(transcriptional.network)

write.graph(transcriptional.network, "my_network.gml", format = "gml")


#COMENZAMOS A VER SI LA RED ES LIBRE DE ESCALA
#TRABAJAMOS CON GRADOS DE ALIDA, LOS QUE NOS INTERESAN, PUES INFO VA DE FUENTE A DIANA


net.degrees=degree(transcriptional.network, mode = "total")
head(net.degrees)
#pueden tener 0 porque tenemos modo out pero lo haremos con mode total

hist(net.degrees)

#hay nodos concentradores y el resto muy bajos, se cumple demasiado, muy radical, casi todos caen muy rapido

tfs.degree.distribution=degree_distribution(transcriptional.network)
power.law.fit(net.degrees)

#pvalor muy pequeño es que NO SE AJUSTA A RED LIBRE DE ESCALA JOJJOJOJOJOJOJOJOJOJJOJOJOJO
#En este caso es porque se desplaza en exceso a existencia de hubs
#lo hemos hecho con kolmogorov smirnov


#podiamos hacer regresion lineal
node.degree.freq=table(net.degrees)
node.degree.freq
#frecuencia con la que aparece un grado en nuestro datast
lm.res=lm(log(node.degree.freq) ~ as.numeric(names(node.degree.freq)))
summary(lm.res)
plot(as.numeric(names(node.degree.freq)), log(node.degree.freq))
#R

#SUBCONJUNTOS (vamos a ver motivos solo entre TFs)
transcriptional.network=read_graph("my_network.gml", format = "gml")

par(mar = c(0, 0, 0, 0))
plot.igraph(transcriptional.network, vertex.size = 2, vertex.label = NA, edge.width = 0.5, arrow.size = 0.5, edge.arrow.size = 0.1, layout = layout_with_fr, main = "Transcriptional Network", vertex.label.color = "black", vertex.frame.color = NA)

#MOTIVOS DE RED DE UN NODO, NO PODEMOS saber si PAR o NAR pues estamos sin ponderacion, sin datos
#de rnaseq

#diag vemos la diagonal de la matriz de adyacencia y sum suma los 1 de la matriz, con sumar los numeros (1 o 0 )
#tendremos numero de autorregulaciones, y para ver si es significativo, creamo redes aleatorias y vemos nuemro de veces
#de autorregulacion, contamos cuantas veces es superior ese numero de autoreglaciones, dividimos y vemos pvalor, pj 0.05
#menos de 5% de veces tiene mas veces e motivo que en nuestra red real y or tanto es  motivo. COMO CREAR. 

#numeor de nodos fuente de interacciones sea el mismo y que las conexiones vayan aleatorias, es decir que vayana mismo numero
#de genes los tf, misma topologia y miramos si regulaciones son iguales, + o - frecuentes en nuestra red respecto a aleatorias
#si es mayor en aleatorias es antimotivo de red


#El número de "1" es el número de autorregulaciones


tfs.adj = as.matrix(get.adjacency(transcriptional.network))
class(tfs.adj)
tfs.adj[1:6,1:6]
#Todos tienen autorregulación en los 6 primeros

sum(diag(tfs.adj))
#tenemos 14 autorregulaciones ahora veremos si es motivo de red
#creamos muchas redes aleatorias y vemos cuantas veces tienen esto. 

num.nodos.al=ncol(tfs.adj)
#deberiamos asignar TFs aleatoriamente sin que necesariamente coincidan con los de nuestra red

out.degree=degree(transcriptional.network, mode = "out")


size=sum(out.degree != 0)
#este caso es raro, normalmente no todos los nodos van a ser TFs, eso sería size
random.tfs=sample(x=1:num.nodos.al, size=sum(out.degree != 0), replace = F)
length(random.tfs)
#la longitus de random.tfs son el numero de TF en nuestra redes aleatorias que debe ser el mismo que en la original
#podemos elegir TFs en red, son los únicos con grado de salida, siendo fuente de la interaccion

random.adj = matrix(0, nrow=num.nodos.al, ncol=num.nodos.al)
dim(random.adj)
#dentro de cada nodo no todos tienen el mimso numero de conexiones de salida. 
out.tfs.true=out.degree[out.degree != 0]

for (i in 1:length(random.tfs))
{
  random.adj[i, sample(x=1:num.nodos.al,size = out.tfs.true[i], replace= F)] = 1
}
sum(random.adj)

sum(diag(random.adj))

#no poner pvalor 0
#__________________________________________________________________________________________________________________

#MOTIVOS DE DOS NODOS: CORREGULACIONES ES LO QUE VAMOS A VER AHORA
#LO QUE SE TIENE QUE CUMPLIR ES QUE EN (i,j)=(j,i)=1
#multiplica (i,j) * (j,i), solo los que den 1 son corregulacion
#si multiplicamos matriz por traspuesta, solo habra 1 donde haya corregulacion, aunque puede haber sitios donde haya
#1 y no lo sea, la diagonal, que son eutorregulaciones, esas no. 

true.cor =(sum(tfs.adj*t(tfs.adj)) - sum (diag(tfs.adj)))/2

##tenemos 28 corregulaicones, nos quedamos con los 1 y eliminamos aquellos que eran autorregulaciones de la diag,
#hay que dividir entre 2 para que no cuente por separado (i,j) y (j,i), es autorregulacion, se cuenta como una

#crearemos mayor numero de redes aleatorias para establecer p valor. 

#COMO GENERAR MUCHAS REDES ALEATORIAS

# Cálculo de autorregulaciones y correlaciones
true.num.auto = sum(diag(tfs.adj))
true.num.cor = (sum(tfs.adj * t(tfs.adj)) - sum(diag(tfs.adj))) / 2

random.auto = c()
random.correg = c()
random.3.nodes = matrix(0, nrow = 100, ncol = 16)  # Cambia 100 al número real de randomizaciones

number.randomizations = 100

for (j in 1:number.randomizations) {
  random.tfs = sample(x = 1:num.nodos.al, size = sum(out.degree != 0), replace = FALSE)
  random.adj = matrix(0, nrow = num.nodos.al, ncol = num.nodos.al)
  
  for (i in 1:length(random.tfs)) {
    random.adj[random.tfs[i], sample(x = 1:num.nodos.al, size = out.tfs.true[i], replace = FALSE)] = 1
  }
  
  current.num.auto = sum(diag(random.adj))
  current.num.cor = (sum(random.adj * t(random.adj)) - current.num.auto) / 2
  
  random.auto = c(random.auto, current.num.auto)
  random.correg = c(random.correg, current.num.cor)
  
  # Cálculo de motivos de 3 nodos
  random.graph = graph.adjacency(random.adj, mode = "directed")
  current.3.nodes = graph.motifs(graph = random.graph, size = 3)
  random.3.nodes[j,] = current.3.nodes
}

# Verificación de autorregulaciones y correlaciones
p_value_auto = sum(true.num.auto <= random.auto) / number.randomizations
p_value_cor = sum(true.num.cor <= random.correg) / number.randomizations

# Verificación de motivos de 3 nodos
true.3.nodes = graph.motifs(graph = transcriptional.network, size = 3)
signf.motif = apply(random.3.nodes, 2, function(col) sum(true.3.nodes <= col) / number.randomizations)

# Imprime los p-values y motivos significativos
cat("P-value para autorregulaciones:", p_value_auto, "\n")
cat("P-value para correlaciones:", p_value_cor, "\n")
cat("Motivos de 3 nodos significativos:", signf.motif, "\n")

# Visualización del subgrafo isomórfico
plot.igraph(graph.isocreate(size = 3, number = 14, directed = TRUE))





## Procesamiento de datos del archivo rst
## por Viri J. Quintero, email: viridiana.j.quintero@gmail.com

# Función para extraer el árbol con nodos
rst_get_tree <- function(file, output){
  dis_tree <- (grep("tree with node labels for Rod Page's TreeView", file))
  write(file[(dis_tree +1)], file = output)
}
# Función para obtener los valores referentes a un nodo
rst_get_node <- function(file, node){
  # Se obtienen cuatro archivos de esta función, el primero es un archivo fasta
  # con la secuencia ancestral reconstruida correspondiente al nodo de interés.
  # El segundo archivo es una tabla (TSV) que incluye el aminoácido de cada
  # posición reconstruido junto con la probabilidad correspondiente, seguido 
  # de la mejor probabilidad. El aminoácido para el peor caso plausible y la 
  # probabilidad de este sitio.
  # El tercer archivo es un archivo fasta con la secuencia correspondiente al
  # peor caso plausible.
  # El cuarto archivo muestra el error esperado o sitios erróneamente esperados
  # así como el número de sitios ambiguos reconstruidos.
  
  # Archivo
  rst_file <- file
  
  # Obtener posiciones y valores para trabajar
  dis_node <- (grep("Prob distribution at node", rst_file))
  pos_node <- (grep((paste("Prob distribution at node ", node, sep = "")), rst_file))
  num_int <- dis_node[grep(pos_node, dis_node)]+4
  num_sec <- dis_node[grep(pos_node, dis_node) +1]-2
  num_site <- num_sec - num_int
  
  # Variable para guardar datos y encabezado
  out_file <- character(length = num_site+2)
  out_prob <- 0
  amn_code <- c("A", "R", "N", "D", "L", "Q", "E", "G", "H", "I",
                "L", "K", "M", "F", "P", "S", "T", "W", "Y", "P")
  out_file[1] <- paste("Pos", "A", "R", "N", "D", "L", "Q", "E", "G", "H", 
                       "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "P", 
                       "Best_prob", "Worst_case", "Worst_prob", sep="  ")
  
  # Secuencia de aminoácidos del nodo de interés
  for_serc <- (paste("node #", node, sep = ""))
  pos_sec <- (grep(for_serc, rst_file))
  pos_fast <- regexpr(for_serc, rst_file[pos_sec])
  sec_ance <- substr(rst_file[pos_sec], attr(pos_fast, "match.length")+1, nchar(rst_file[pos_sec]))
  sec_ance <- gsub(" ", "", sec_ance)
  
  # Guardar archivo fasta
  nod_fast <- character(length = 2)
  nod_fast[1] <-paste(">nodo_", node, sep="")
  nod_fast[2] <- sec_ance
  bad_fast <- character(length = 2)
  bad_fast[1] <-paste(">nodo_", node, "_peor_caso_plausible", sep="")
  
  # Ciclo for para extraer la información
  num_cont <- 0
  amb_site <- 0
  for (i in (num_int):(num_sec)){
    num_cont <- num_cont + 1
    sec_site <- substr(sec_ance, num_cont, num_cont)
    sec_line <- rst_file[i]
    amino_A <- substr(sec_line, (unlist(regexpr("A\\(", sec_line))+2), (unlist(regexpr("A\\(", sec_line))+6))
    amino_R <- substr(sec_line, (unlist(regexpr("R\\(", sec_line))+2), (unlist(regexpr("R\\(", sec_line))+6))
    amino_N <- substr(sec_line, (unlist(regexpr("N\\(", sec_line))+2), (unlist(regexpr("N\\(", sec_line))+6))
    amino_D <- substr(sec_line, (unlist(regexpr("D\\(", sec_line))+2), (unlist(regexpr("D\\(", sec_line))+6))
    amino_L <- substr(sec_line, (unlist(regexpr("L\\(", sec_line))+2), (unlist(regexpr("L\\(", sec_line))+6))
    amino_Q <- substr(sec_line, (unlist(regexpr("Q\\(", sec_line))+2), (unlist(regexpr("Q\\(", sec_line))+6))
    amino_E <- substr(sec_line, (unlist(regexpr("E\\(", sec_line))+2), (unlist(regexpr("E\\(", sec_line))+6))
    amino_G <- substr(sec_line, (unlist(regexpr("G\\(", sec_line))+2), (unlist(regexpr("G\\(", sec_line))+6))
    amino_H <- substr(sec_line, (unlist(regexpr("H\\(", sec_line))+2), (unlist(regexpr("H\\(", sec_line))+6))
    amino_I <- substr(sec_line, (unlist(regexpr("I\\(", sec_line))+2), (unlist(regexpr("I\\(", sec_line))+6))
    amino_L <- substr(sec_line, (unlist(regexpr("L\\(", sec_line))+2), (unlist(regexpr("L\\(", sec_line))+6))
    amino_K <- substr(sec_line, (unlist(regexpr("K\\(", sec_line))+2), (unlist(regexpr("K\\(", sec_line))+6))
    amino_M <- substr(sec_line, (unlist(regexpr("M\\(", sec_line))+2), (unlist(regexpr("M\\(", sec_line))+6))
    amino_F <- substr(sec_line, (unlist(regexpr("F\\(", sec_line))+2), (unlist(regexpr("F\\(", sec_line))+6))
    amino_P <- substr(sec_line, (unlist(regexpr("P\\(", sec_line))+2), (unlist(regexpr("P\\(", sec_line))+6))
    amino_S <- substr(sec_line, (unlist(regexpr("S\\(", sec_line))+2), (unlist(regexpr("S\\(", sec_line))+6))
    amino_T <- substr(sec_line, (unlist(regexpr("T\\(", sec_line))+2), (unlist(regexpr("T\\(", sec_line))+6))
    amino_W <- substr(sec_line, (unlist(regexpr("W\\(", sec_line))+2), (unlist(regexpr("W\\(", sec_line))+6))
    amino_Y <- substr(sec_line, (unlist(regexpr("Y\\(", sec_line))+2), (unlist(regexpr("Y\\(", sec_line))+6))
    amino_V <- substr(sec_line, (unlist(regexpr("V\\(", sec_line))+2), (unlist(regexpr("V\\(", sec_line))+6))
    
    all_amino <- as.numeric(c(amino_A, amino_R, amino_N, amino_D, amino_L, amino_Q, 
                              amino_E, amino_G, amino_H, amino_I, amino_L, amino_K, amino_M, 
                              amino_F, amino_P, amino_S, amino_T, amino_W, amino_Y, amino_V)) 
    
    # Probabilidad mayor
    pos_1 <- which.max(all_amino)
    pro_1 <- all_amino[pos_1]
    all_amino[pos_1] <- 0
    # Segunda probabilidad
    pos_2 <- which.max(all_amino)
    pro_2 <- all_amino[pos_2]
    # Peor caso plausible
    if (pro_1 & pro_2 >= 0.2){
      amb_site <- amb_site +1
      bad_site <- amn_code[pos_2]
      bad_prob <- pro_2
    } else {
      bad_site <- amn_code[pos_1]
      bad_prob <- pro_1
    }
    bad_fast[2] <- paste(bad_fast[2], bad_site, sep = "")
    
    out_file[num_cont+1] <- paste(sec_site, amino_A, amino_R, amino_N, amino_D, amino_L, amino_Q, 
                                  amino_E, amino_G, amino_H, amino_I, amino_L, amino_K, amino_M, 
                                  amino_F, amino_P, amino_S, amino_T, amino_W, amino_Y, amino_V,
                                  pro_1, bad_site, bad_prob, sep="  ")
    out_prob <- out_prob +(1 - pro_1)
  }
  
  # Error esperado o sitios erróneamente esperados
  
  err_site <- out_prob / (length(num_site)+1)
  err_file <- paste("El error esperado o sitios erroneamente esperados es igual a: ", as.character(err_site), 
                    "\nEl número de sitios ambiguos es igual a: ", amb_site, sep = "")
  
  # Guardar archivo 
  out_nam1 <- paste("nodo_", node, ".tsv", sep="")
  out_nam2 <- paste("error_", node, ".txt", sep="")
  out_fas1 <- paste("nodo_", node, ".fas", sep="")
  out_fas2 <- paste("peor_caso_", node, ".fas", sep="")
  write(out_file[out_file != ""], file = out_nam1)
  write(err_file[err_file != ""], file = out_nam2)
  write(nod_fast[nod_fast != ""], file = out_fas1)
  write(bad_fast[bad_fast != ""], file = out_fas2)
}

## MODO DE USO
## Es necesario cargar a una variable el contenido del archivo rst
## posteriormente introducir esa variable a las distintas funciones

# Leer archivo rst
rst_file <- readLines("rst")

# Obtener árbol con nodos
# La función recibe el archivo rst y el nombre con el que se desea guardar
# rst_get_tree(file, output)
rst_get_tree(rst_file, "anc_node.tre")

# Obtener valores referentes al nodo de interés
# La función recibe el archivo rst y el nodo de interés.
# rst_get_node(file, #node)
rst_get_node(rst_file, 385)

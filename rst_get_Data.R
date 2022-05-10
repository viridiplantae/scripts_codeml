rst_get_node <- function(file, node){
  #
  rst_file <- readLines(file)
  
  dis_node <- (grep("Prob distribution at node", rst_file))
  pos_node <- (grep((paste("Prob distribution at node ", node, sep = "")), rst_file))
  num_int <- dis_node[grep(pos_node, dis_node)]+4
  num_sec <- dis_node[grep(pos_node, dis_node) +1]-2
  num_site <- num_sec - num_int
  
  #variable para guardar datos y encabezado
  out_file <- character(length = num_site+2)
  out_file[1] <- paste("Pos", "A", "R", "N", "D", "L", "Q", "E", "G", "H",
                       "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y",
                       "P", sep="  ")
  
  #Secuencia de aminoácidos del nodo de interés
  for_serc <- (paste("node #", node, sep = ""))
  pos_sec <- (grep(for_serc, rst_file))
  pos_fast <- regexpr(for_serc, rst_file[pos_sec])
    sec_ance <- substr(rst_file[pos_sec], attr(pos_fast, "match.length")+1, nchar(rst_file[pos_sec]))
  sec_ance <- gsub(" ", "", sec_ance)
  
  #Ciclo for para extraer la información
  num_cont <- 0
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
    out_file[num_cont+1] <- paste(sec_site, amino_A, amino_R, amino_N, amino_D, amino_L, amino_Q, 
                                  amino_E, amino_G, amino_H, amino_I, amino_L, amino_K, amino_M, 
                                  amino_F, amino_P, amino_S, amino_T, amino_W, amino_Y, amino_P, sep="  ")
  }
  
  #Guardar archivo 
  out_name <- paste("nodo_", node, ".tsv", sep="")
  write(out_file[out_file != ""], file = out_name)
}

#rst_get_node(file, #node)
rst_get_node("rst.txt", 385)

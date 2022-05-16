# scripts_codeml

Script en R con funciones para el procesamiento de datos del archivo rst obtenido de ejecutar codeml de paml.

## Modo de uso
Cargar a una variable el archivo rst y posteriormente introducir esa variable a las distintas funciones.

```R
  rst_file <- file
```

### Función para extraer el árbol con nodos
La función recibe el archivo rst previamente cargado y el nombre con el que se desea guardar el árbol.tre
```R
rst_get_tree(file, output)
```
### Función para obtener los valores referentes a un nodo
La función recibe el archivo rst y el nodo de interés.
```R
rst_get_tree(file, output)
```
Se obtienen cuatro archivos de esta función, el primero es un archivo fasta con la secuencia ancestral reconstruida correspondiente al nodo de interés. El segundo archivo es una tabla (TSV) que incluye el aminoácido de cada posición reconstruido junto con la probabilidad correspondiente, seguido de la mejor probabilidad. El aminoácido para el peor caso plausible y la probabilidad de este sitio. El tercer archivo es un archivo fasta con la secuencia correspondiente al peor caso plausible. El cuarto archivo muestra el error esperado o sitios erróneamente esperados así como el número de sitios ambiguos reconstruidos.
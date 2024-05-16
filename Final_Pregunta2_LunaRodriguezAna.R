
    #Creador: Ana Fernanda Luna Rodriguez
    # Instalar y cargar el paquete pacman
    #install.packages("pacman")
  library(pacman)
  
  # Cargar los paquetes necesarios
  p_load("readr", #para llamar las bases de datos
         "ggplot2", #para graficar
         "ggrepel", #para etiquetar datos en una gráfica
         "matrixTests", #Para realizar la prueba estadística
         "dplyr") #facilita el manejo de datos
  
  # Cargar los datos desde la URL
  datos <- read_csv(file="https://raw.githubusercontent.com/ManuelLaraMVZ/Transcript-mica/main/examen2")
  datos
  
  # Extracción de genes controles (referencia)
  controles <- datos %>% 
    filter(Condicion == "Control")
  head(controles)
  
  # Sacar promedios de los controles 
  promedio_controles <- controles %>% 
    summarise(Mean_C1 = mean(Cx1),
              Mean_C2 = mean(Cx2),
              Mean_C3 = mean(Cx3),
              Mean_T1 = mean(T1),
              Mean_T2 = mean(T2),
              Mean_T3 = mean(T3)) %>% 
    mutate(Gen = "promedio_controles") %>% 
    select(Gen, Mean_C1, Mean_C2, Mean_C3, Mean_T1, Mean_T2, Mean_T3)
  promedio_controles
  
  # Extraer los genes de la tabla de datos
  genes <- datos %>% 
    filter(Condicion == "Target") %>%
    select(-Condicion)
  head(genes)
  
  # Sacar el 2^-DCT
  DCT <- genes %>% 
    mutate(DCT_C1 = 2^-(Cx1 - promedio_controles$Mean_C1),
           DCT_C2 = 2^-(Cx2 - promedio_controles$Mean_C2),
           DCT_C3 = 2^-(Cx3 - promedio_controles$Mean_C3),
           DCT_T1 = 2^-(T1 - promedio_controles$Mean_T1),
           DCT_T2 = 2^-(T2 - promedio_controles$Mean_T2),
           DCT_T3 = 2^-(T3 - promedio_controles$Mean_T3)) %>%
    select(miRNA, DCT_C1, DCT_C2, DCT_C3, DCT_T1, DCT_T2, DCT_T3)
  DCT
  
  # Calcular promedios de 2^-DCT para genes
  promedio_genes <- DCT %>%
    mutate(Mean_DCT_Cx = (DCT_C1 + DCT_C2 + DCT_C3) / 3,
           Mean_DCT_Tx = (DCT_T1 + DCT_T2 + DCT_T3) / 3)%>%
    select (1,8,9)
  promedio_genes

  
  
  
 
  
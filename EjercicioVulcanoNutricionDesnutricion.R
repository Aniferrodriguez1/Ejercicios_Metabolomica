#Creador: Ana Fernanda Luna Rodriguez
#Grafica de Volcano 
install.packages("pacman")
library(pacman)
p_load("readr", #para lamar las bases de datos
       "ggplot2", #para graficar
       "ggrepel", #para etiquetar datos en una gráfica
       "matrixTests", #Para realizar la prueba estadística
       "dplyr") #facilita el manejo de datos

datos <- read_csv(file="https://raw.githubusercontent.com/ManuelLaraMVZ/Transcript-mica/main/DesnutridasvsEunutridas.csv")
datos

#extracción de genes controles (referencia)
controles <- datos %>% 
  filter(Condicion=="Control")
head(controles)

#sacar promedios de los controles 
promedio_controles <-  controles %>% 
  summarise(Mean_C1 = mean(Cx1),
            Mean_C2 = mean(Cx2),
            Mean_C3 = mean(Cx3),
            Mean_T1 = mean(T1),
            Mean_T2 = mean(T2),
            Mean_T3 = mean(T3))
mutate(Gen="promedio_controles") %>% 
  select(7,1,2,3,4,5,6)
promedio_controles

#######
#extraer los genes de la tabla de "datos
                  
                  genes <- datos %>% 
                    filter(Condicion=="Target") %>%
                    select(-2)
                  head(genes)
                  
                  ##############################
                  
                  #Sacar el 2^-DCT
                  
                  DCT <- genes %>% 
                    mutate(DCT_C1=2^-(Cx1-promedio_controles$Mean_C1),
                           DCT_C2=2^-(Cx2-promedio_controles$Mean_C2),
                           DCT_C3=2^-(Cx3-promedio_controles$Mean_C3),
                           DCT_T1=2^-(T1-promedio_controles$Mean_T1),
                           DCT_T2=2^-(T2-promedio_controles$Mean_T2),
                           DCT_T3=2^-(T3-promedio_controles$Mean_T3)) %>%
                    select(-2,-3,-4,-5,-6,-7)
                  DCT
                  
                  promedio_genes <- DCT %>%
                    mutate(Mean_DCT_Cx=(DCT_C1+DCT_C2+DCT_C3)/3,
                           Mean_DCT_Tx=(DCT_T1+DCT_T2+DCT_T3)/3)
                  
                  promedio_genes
                  
                  #######################################
                  
                  #Definir prueba estadística: 
                  
                  tvalue_gen <- row_t_welch(promedio_genes[,c("DCT_C1",
                                                              "DCT_C2",
                                                              "DCT_C3")],
                                            promedio_genes[,c("DCT_T1",
                                                              "DCT_T2",
                                                              "DCT_T3")])
                  View(tvalue_gen)
                  
                  FCyPV <- promedio_genes %>%
                    select(1,8,9) %>%
                    mutate(p_value=tvalue_gen$pvalue,
                           Fold_change = Mean_DCT_Tx/Mean_DCT_Cx) %>%
                    select(1,4,5)
                  
                  FCyPV
                  
                  Logs <- FCyPV %>%
                    mutate(LPV=-log10(p_value),
                           LFC=log2(Fold_change)) %>%
                    select(1,4,5)
                  Logs
                  
                  #########
                  
                  vulcano <- ggplot(Logs, aes(x = LFC, y = LPV)) +
                    geom_point(size = 2, color = "blue") +
                    theme_classic() +
                    labs(
                      title = "Análisis comparativo de miRNAs",
                      caption = "Creador: Ana Fernanda Luna",
                      x = "Log2(2-DDCT)",
                      y = "-Log10(valor de p)"
                    )
                  vulcano
                  
                  ###########################
                  
                  #Límites
                  
                  limite_p <- 0.05
                  
                  limite_FC <-1.5
                  
                  down_regulated <- Logs %>%
                    filter(LFC < -log2(limite_FC),
                           LPV > -log10(limite_p))
                  
                  down_regulated
                  
                  up_regulated <- Logs %>%
                    filter(LFC > log2(limite_FC),
                           LPV > -log10(limite_p))
                  
                  up_regulated
                  
                  top_down_regulated <- down_regulated %>%
                    arrange(desc(LPV)) %>%
                    head (5)
                  
                  top_down_regulated
                  
                  top_up_regulated <- up_regulated %>%
                    arrange(desc(LPV)) %>%
                    head (5)
                  
                  top_up_regulated
                  
                  #######################
                  #Mejorar la grafica 
                  
                  vulcano2 <- vulcano +
                    geom_hline(yintercept = -log10(limite_p), 
                               linetype = "dashed") +
                    geom_vline(xintercept = c(-log2(limite_FC), log2(limite_FC)), 
                               linetype = "dashed")
                  
                  
                  vulcano2
                  

                  vulcano3 <- vulcano2 +
                    geom_point(data = up_regulated,
                               aes(x = LFC, y = LPV),
                               alpha = 1,
                               size = 3,
                               color = "pink" ) +
                    geom_point(data = down_regulated,
                               aes(x = LFC, y = LPV),
                               alpha = 1,
                               size = 3,
                               color = "purple")
                  
                  vulcano3
                  
                  vulcano4 <- vulcano3 +
                    geom_point(data = up_regulated,
                               aes(x = LFC, y = LPV),
                               alpha = 0.75,  # Ajustar la transparencia si es necesario
                               size = 2,  # Reducir el tamaño de los puntos puede ayudar
                               color = "pink") +
                    geom_point(data = down_regulated,
                               aes(x = LFC, y = LPV),
                               alpha = 0.75,  # Ajustar la transparencia si es necesario
                               size = 2,  # Reducir el tamaño de los puntos puede ayudar
                               color = "purple") +
                    geom_label_repel(data = top_up_regulated,
                                     aes(x = LFC, y = LPV, label = Gen),
                                     max.overlaps = 100,
                                     box.padding = 0.35,  # Ajustar el espacio alrededor de las etiquetas
                                     point.padding = 0.5) +  # Ajustar el espacio alrededor de los puntos
                    geom_label_repel(data = top_down_regulated,
                                     aes(x = LFC, y = LPV, label = Gen),
                                     max.overlaps = 100,
                                     box.padding = 0.35,
                                     point.padding = 0.5)
                  
                  vulcano4
                  
                  ggsave("Vulcano4_NUTDESNUTR.jpeg",
                         plot = vulcano4,
                         height = 12, width = 10,
                         dpi = 300)
                  
#este codigo realiza un análisis comparativo de la expresión génica entre animales
#desnutridos y eunutridos, identifica los genes que muestra cambios significativos
#en su expresión y visualiza estos resultados en un gráfico de volcan. Las etiquetas 
#ayudan a identificar genes específicos que son de particular interés debido a sus
#cambios significativos en la expresión
                  
                  
                  
                  
                  
                  
                  
#Creadora: Ana Fernanda Luna
#Gráfica de puntos

install.packages("pacman")
library(pacman)
p_load("readr", #para lamar las bases de datos
       "ggplot2", #para graficar
       "dplyr") #facilira el manejo de datos
datos <- read_csv(file="https://raw.githubusercontent.com/ManuelLaraMVZ/Transcript-mica/main/datos_miRNAs.csv")
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

#Grafica 

grafica_puntos <- ggplot(promedio_genes,
                         mapping=aes(x=Mean_DCT_Cx,
                                     y=Mean_DCT_Tx))+
      geom_point(size=3,color="blue")+
      labs(title="Analisis de puntos",
           x="Condición Control (2^-DCT)",
           y="Tratamiento (2^-DCT)")+
  geom_smooth(method="lm",
              alpha=0.05,
              linewidth=1, span=1) + #para regresión
  theme_minimal()
grafica_puntos

guardar <- ggsave("grafica_dispersion.jpeg",
                  plot= grafica_puntos,
                  width=6, height=4,
                  dpi=300)


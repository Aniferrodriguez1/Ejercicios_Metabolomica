#Creador: Ana Fernanda Luna
# Llamada de bases de datos 

#Llamar la base de datos desde la computadora

install.packages("readr")
library("readr")

titulacion<-read_csv(file="datos_titulacion.csv")


#Llamar desde un repositorio (internet)

repositorio<-read_csv(file="https://raw.githubusercontent.com/ManuelLaraMVZ/titulacion_amino_acidos/main/datos_titulacion%20(2).csv")
head(repositorio) #para ver el encabezado muestra los primeros 6 datos
View(repositorio) #para ver la tabla

#gráfica
install.packages("ggplot2")
library(ggplot2)

grafica<- ggplot(repositorio,
                 aes(x=Volumen, #el aes lo que hace es la estetica de la grafica
                 y=pH)) + #poner los datos en paréntesis dentro del paréntesis y cuidar las comas
  geom_line()+
  labs(title="Titulación de cisteina",
       x="Volumen ácido (uL)",
       y="Valor de pH")+
  theme_dark() #para cambiar el fondo de color
  
grafica #para visualizarlo la gráfica en puntos geom_point #para visualizar la linea geom_line

ggsave("titulación_repertorio.jpeg", plot = grafica, width=6, height=4, dpi=500) #para guardarlo en la carpeta

                 
                 
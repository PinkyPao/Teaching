# =============================================================================
# Clase: Más allá de la esperanza de vida
# Medidas resumen de la tabla de mortalidad y dispersión de la vida
# Paola Vázquez-Castillo
# =============================================================================

library(tidyverse)
library(readxl)
library(patchwork)

# =============================================================================
# 1. Funciones
# =============================================================================

# -----------------------------------------------------------------------------
# 1.1 Construcción de la tabla de mortalidad
#     Recibe un data frame con columnas EDAD y mx (una sola tabla)
#     Devuelve la tabla de mortalidad completa
# -----------------------------------------------------------------------------

tabla.mort <- function(df, radix = 1) {
  df %>%
    arrange(EDAD)%>%
    mutate(
      # ax: Andreev & Kingkade (2015) para edad 0, 0.5 para el resto
      ax = case_when(
        EDAD > 0      ~ 0.5,
        mx >= 0.06891 ~ 0.31411,
        mx >= 0.01724 ~ 0.04667 + 3.88089 * mx,
        TRUE          ~ 0.14903 - 2.05527 * mx
      ),
      qx = mx / (1 + (1 - ax) * mx),
      qx = if_else(EDAD == max(EDAD) , 1, qx),   # cierre de la tabla
      px = 1 - qx,
      lx = radix * cumprod(lag(px, default = 1)),
      dx = lx * qx,
      Lx = if_else(EDAD == max(EDAD)  & mx > 0,
                   lx / mx,
                   lx - (1 - ax) * dx),
      Tx = rev(cumsum(rev(Lx))),
      ex = Tx / lx
    )
}

#################################################
######### Leer datos para México
#################################################

Nx<- read_excel("data/0_Pob_Mitad_1950_2070.xlsx")
Dx <- read_excel("data/1_Defunciones_1950_2070.xlsx")

# Comprobamos que se haya leido bien
# head(Nx)
# head(Dx)

# PASO 0. Creamos Mx que seran el input para nuestra funcion de la tabla de mortalidad
new_df <- Nx %>%
  left_join(Dx, by = c("AÑO", "ENTIDAD", "CVE_GEO", "EDAD", "SEXO")) %>%
  mutate(mx = if_else(POBLACION > 0, DEFUNCIONES / POBLACION, 0)) %>%
  select(AÑO, ENTIDAD, SEXO, EDAD,  mx) %>%
  arrange(AÑO, ENTIDAD, SEXO, EDAD) 

# Creamos las tablas de vida para cada año, entidad y sexo
lt <- new_df %>%
  group_by(AÑO, ENTIDAD, SEXO) %>%
  group_modify(~ tabla.mort(.x)) %>%
  ungroup()

# Comprobamos que esten bien 
# head(lt)

lt_clase<-lt %>% filter(AÑO==2019,
                        SEXO=="Mujeres",
                        ENTIDAD=="México")

###############################################################
######### Estimacion de medidas de posicion
###############################################################

#########
# Mediana
#########

lt_clase %>% 
  # Aqui usamos 0.5 porque radix es 1. En otros casos se puede 1) estandarizar lx y dx 
  # o en vez de restar 0.5 restar lx[1]/2
  mutate(Mediana_ent = EDAD[which.min(abs(lx - 0.5))],
         pos=which(EDAD == Mediana_ent),
         Mediana=Mediana_ent + ((0.5-lx[pos])/(lx[pos+1]-lx[pos]))) %>% 
  summarise(Mediana=unique(Mediana))

# Transformamos esto en una funcion
mediana<-function(lt){
  lt %>% 
    mutate(Mediana_ent = EDAD[ which.min(abs(lx - 0.5))],
           pos=which(EDAD ==Mediana_ent),
           Mediana=Mediana_ent + ((0.5-lx[pos])/(lx[pos+1]-lx[pos]))) %>% 
    summarise(Mediana=unique(Mediana))
}

# Probamos la funcion
mediana(lt_clase)

# Si quiero hacer para todas las tablas de mortalidad que tengo 
lt %>% 
  group_by(ENTIDAD, AÑO, SEXO) %>% 
  group_modify(~ mediana(.x)) %>%
  ungroup()

#########
# Moda
#########

lt_clase %>% 
  # Aqui usamos 0.5 porque radix es 1. En otros casos se puede 1) estandarizar lx y dx 
  # o en vez de restar 0.5 restar lx[1]/2
  mutate(Moda_ent = EDAD[which.max(dx)],
         pos=which(EDAD ==Moda_ent),
         Moda=Moda_ent + ((dx[pos]-dx[pos+1])/(2*dx[pos]-dx[pos+1]-dx[pos-1]))) %>% 
  summarise(Moda=unique(Moda))

# Transformar en una funcion

moda <- function(lt){
  lt %>% 
    # filtramos porque queremos la moda de mortalidad a la edad adulta
      filter(EDAD >= 10) %>% 
  mutate(Moda_ent = EDAD[which.max(dx)],
         pos=which(EDAD ==Moda_ent),
         Moda=Moda_ent + ((dx[pos]-dx[pos+1])/(2*dx[pos]-dx[pos+1]-dx[pos-1]))) %>% 
    summarise(Moda=unique(Moda))
}



moda(lt_clase)

# Si quiero hacer para todas las tablas de mortalidad que tengo 
lt %>% 
  group_by(ENTIDAD, AÑO, SEXO) %>% 
  group_modify(~ moda(.x)) %>%
  ungroup()


###############################################################
######### Estimacion de medidas de longevidad de la vida
###############################################################

#########
# Desviación estandar
#########

lt_clase %>% 
  mutate(SD=sqrt(sum(((EDAD-ex[1])^2)*dx))) %>% 
  summarise(SD=unique(SD))

ds<-function(lt){
  lt %>% 
    mutate(SD=sqrt(sum(((EDAD-ex[1])^2)*dx)))%>% 
    summarise(SD=unique(SD))
}

# Si quiero hacer para todas las tablas de mortalidad que tengo 
lt %>% 
  group_by(ENTIDAD, AÑO, SEXO) %>% 
  group_modify(~ ds(.x)) %>%
  ungroup()

#########
# e-dagger 
#########
lt_clase %>% 
  mutate(e_dagger=sum(ex*dx)) %>% 
  summarise(e_dagger=unique(e_dagger))

e_dagger<-function(lt){
  lt %>% 
    mutate(e_dagger=sum(ex*dx)) %>% 
    summarise(e_dagger=unique(e_dagger))
}

# Si quiero hacer para todas las tablas de mortalidad que tengo 
lt %>% 
  group_by(ENTIDAD, AÑO, SEXO) %>% 
  group_modify(~ e_dagger(.x)) %>%
  ungroup()

#########
# Gini 
#########
lt_clase %>% 
  mutate(Gini=1-(1/ex[1])*sum(lx^2)) %>% 
  summarise(Gini=unique(Gini))

gini<-function(lt){
  lt %>% 
    mutate(Gini=1-(1/ex[1])*sum(lx^2)) %>% 
    summarise(Gini=unique(Gini))
}

# Si quiero hacer para todas las tablas de mortalidad que tengo 
lt %>% 
  group_by(ENTIDAD, AÑO, SEXO) %>% 
  group_modify(~ gini(.x)) %>%
  ungroup()

#########
# Entropía 
#########
lt_clase %>% 
  mutate(H=(e_dagger(.))/ex[1]) %>% 
  summarise(H=unique(H))

lt_clase %>% 
  filter(EDAD>=1) %>% 
  mutate(H=-(sum(lx*log(lx)))/sum(lx)) %>% 
  summarise(H=unique(H))

entropia<-function(lt){
    lt %>% 
    mutate(H_lx=-(sum(lx*log(lx),na.rm = T))/sum(lx),
           H_e_dagger=(e_dagger(.))/ex[1]) %>% 
    summarise(H_lx=unique(H_lx),
              H_e_dagger=unique(H_e_dagger))
  }

# Si quiero hacer para todas las tablas de mortalidad que tengo 
lt %>% 
  group_by(ENTIDAD, AÑO, SEXO) %>% 
  group_modify(~ entropia(.x)) %>%
  ungroup() %>% 
  # Para ver la diferencia entre una approximación y otra 
  mutate(dif=H_lx-H_e_dagger$e_dagger) %>% 
  ggplot(aes(x=AÑO, y=dif, col=SEXO))+
  geom_point()






